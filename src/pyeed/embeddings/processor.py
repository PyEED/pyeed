"""
Main embedding processor for coordinating embedding operations.

Provides high-level interfaces for batch processing, single sequence processing,
and database operations with automatic device management and model loading.
"""

from typing import List, Union, Any, Literal, Optional
import torch
from torch.nn import DataParallel, Module
from loguru import logger
import numpy as np
from numpy.typing import NDArray
import time
from concurrent.futures import ThreadPoolExecutor
import os

from .factory import ModelFactory
from .base import BaseEmbeddingModel
from .models import ESM2EmbeddingModel, ESMCEmbeddingModel, ESM3EmbeddingModel, ProtT5EmbeddingModel
from .database import update_protein_embeddings_in_db
from .utils import free_memory
from pyeed.dbconnect import DatabaseConnector


class EmbeddingProcessor:
    """
    Main processor for handling protein embedding operations.
    
    Automatically manages device selection, model loading, and provides
    simplified interfaces for all embedding operations.
    """
    
    def __init__(self):
        self._models: dict[str, BaseEmbeddingModel] = {}
        self._devices: List[torch.device] = []
        self._initialize_devices()
    
    def _initialize_devices(self) -> None:
        """Initialize available devices for computation."""
        if torch.cuda.is_available():
            device_count = torch.cuda.device_count()
            self._devices = [torch.device(f"cuda:{i}") for i in range(device_count)]
            logger.info(f"Initialized {device_count} GPU device(s): {self._devices}")
        else:
            self._devices = [torch.device("cpu")]
            logger.warning("No GPU available, using CPU.")
    
    def get_available_devices(self) -> List[torch.device]:
        """Get list of available devices."""
        return self._devices.copy()
    
    def get_or_create_model(
        self, 
        model_name: str, 
        device: Optional[torch.device] = None
    ) -> BaseEmbeddingModel:
        """Get existing model or create new one on specified or best available device."""
        if device is None:
            device = self._devices[0]  # Use first available device
        
        key = f"{model_name}_{device}"
        if key not in self._models:
            self._models[key] = ModelFactory.create_model(model_name, device)
            logger.info(f"Loaded model {model_name} on {device}")
        return self._models[key]
    
    def calculate_batch_embeddings(
        self,
        data: List[tuple[str, str]],
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        batch_size: int = 16,
        num_gpus: Optional[int] = None,
        db: Optional[DatabaseConnector] = None,
        embedding_type: Literal["last_hidden_state", "all_layers", "first_layer", "final_embeddings"] = "last_hidden_state"
    ) -> Optional[List[NDArray[np.float64]]]:
        """
        Calculate embeddings for a batch of sequences with automatic device management.
        
        Args:
            data: List of (accession_id, sequence) tuples
            model_name: Name of the model to use
            batch_size: Batch size for processing
            num_gpus: Number of GPUs to use (None = use all available)
            db: Database connector for storing results (optional)
            embedding_type: Type of embedding to calculate:
                - "last_hidden_state": Use last hidden state (most common)
                - "all_layers": Average across all transformer layers  
                - "first_layer": Use first layer embedding
                - "final_embeddings": Robust option that works across all models (recommended for compatibility)
            
        Returns:
            List of embeddings if db is None, otherwise None (results stored in DB)
        """
        # Disable tqdm to prevent threading issues with multiple GPUs
        os.environ['DISABLE_TQDM'] = 'True'
        
        if not data:
            logger.info("No sequences to process.")
            return []
        
        # Determine number of GPUs to use
        available_gpus = len([d for d in self._devices if d.type == 'cuda'])
        if num_gpus is None:
            num_gpus = available_gpus
        else:
            num_gpus = min(num_gpus, available_gpus)
        
        if num_gpus == 0:
            devices_to_use = [torch.device("cpu")]
            num_gpus = 1
        else:
            devices_to_use = [torch.device(f"cuda:{i}") for i in range(num_gpus)]
        
        logger.info(f"Processing {len(data)} sequences using {num_gpus} device(s)")
        
        # Load models for each device
        models = []
        for device in devices_to_use:
            try:
                model = self.get_or_create_model(model_name, device)
                models.append(model)
            except Exception as e:
                if "tqdm" in str(e).lower() or "_lock" in str(e).lower():
                    logger.warning(f"Model loading failed on {device} due to threading issue. Reducing to single GPU mode.")
                    # Fall back to single GPU mode to avoid threading issues
                    devices_to_use = [devices_to_use[0]]
                    num_gpus = 1
                    models = [self.get_or_create_model(model_name, devices_to_use[0])]
                    break
                else:
                    raise e
        
        # Split data across devices
        gpu_batches = [
            data[i::num_gpus] for i in range(num_gpus)
        ]
        
        start_time = time.time()
        all_embeddings = []
        
        if num_gpus == 1:
            # Single device processing
            embeddings = self._process_batch_single_device(
                gpu_batches[0], models[0], batch_size, db, embedding_type
            )
            all_embeddings.extend(embeddings)
        else:
            # Multi-device parallel processing
            with ThreadPoolExecutor(max_workers=num_gpus) as executor:
                futures = []
                for i, gpu_data in enumerate(gpu_batches):
                    if not gpu_data:
                        continue
                    
                    futures.append(
                        executor.submit(
                            self._process_batch_single_device,
                            gpu_data,
                            models[i],
                            batch_size,
                            db,
                            embedding_type
                        )
                    )
                
                for future in futures:
                    embeddings = future.result()
                    all_embeddings.extend(embeddings)
        
        end_time = time.time()
        logger.info(f"Batch processing completed in {end_time - start_time:.2f} seconds")
        
        return all_embeddings if db is None else None
    
    def _process_batch_single_device(
        self,
        data: List[tuple[str, str]],
        model: BaseEmbeddingModel,
        batch_size: int,
        db: Optional[DatabaseConnector] = None,
        embedding_type: str = "last_hidden_state"
    ) -> List[NDArray[np.float64]]:
        """Process batch on a single device."""
        all_embeddings = []
        
        for batch_start in range(0, len(data), batch_size):
            batch_end = min(batch_start + batch_size, len(data))
            batch = data[batch_start:batch_end]
            
            accessions, sequences = zip(*batch)
            current_batch_size = len(sequences)
            
            while current_batch_size > 0:
                try:
                    # Calculate embeddings based on type
                    if embedding_type == "last_hidden_state":
                        # no batching for last hidden state
                        embeddings_batch = [
                            model.get_single_embedding_last_hidden_state(seq)
                            for seq in sequences[:current_batch_size]
                        ]
                    elif embedding_type == "all_layers":
                        embeddings_batch = [
                            model.get_single_embedding_all_layers(seq)
                            for seq in sequences[:current_batch_size]
                        ]
                    elif embedding_type == "first_layer":
                        embeddings_batch = [
                            model.get_single_embedding_first_layer(seq)
                            for seq in sequences[:current_batch_size]
                        ]
                    elif embedding_type == "final_embeddings":
                        embeddings_batch = [
                            model.get_final_embeddings(seq)
                            for seq in sequences[:current_batch_size]
                        ]
                    else:
                        raise ValueError(f"Unknown embedding_type: {embedding_type}")
                    
                    # Store in database if provided
                    if db is not None:
                        update_protein_embeddings_in_db(
                            db, list(accessions[:current_batch_size]), embeddings_batch
                        )
                    
                    all_embeddings.extend(embeddings_batch)
                    break  # Successful execution
                
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
                    current_batch_size = max(1, current_batch_size // 2)
                    logger.warning(f"Reduced batch size to {current_batch_size} due to OOM error.")
        
        return all_embeddings
    
    def calculate_single_embedding(
        self,
        sequence: str,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        embedding_type: Literal["last_hidden_state", "all_layers", "first_layer", "final_embeddings"] = "last_hidden_state",
        device: Optional[torch.device] = None
    ) -> NDArray[np.float64]:
        """
        Calculate embedding for a single sequence.
        
        Args:
            sequence: Protein sequence
            model_name: Name of the model to use
            embedding_type: Type of embedding to calculate
            device: Specific device to use (optional)
            
        Returns:
            Embedding as numpy array
        """
        model = self.get_or_create_model(model_name, device)
        
        if embedding_type == "last_hidden_state":
            return model.get_single_embedding_last_hidden_state(sequence)
        elif embedding_type == "all_layers":
            return model.get_single_embedding_all_layers(sequence)
        elif embedding_type == "first_layer":
            return model.get_single_embedding_first_layer(sequence)
        elif embedding_type == "final_embeddings":
            return model.get_final_embeddings(sequence)
        else:
            raise ValueError(f"Unknown embedding_type: {embedding_type}")
    
    def calculate_database_embeddings(
        self,
        db: DatabaseConnector,
        batch_size: int = 16,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        num_gpus: Optional[int] = None,
        embedding_type: Literal["last_hidden_state", "all_layers", "first_layer", "final_embeddings"] = "last_hidden_state"
    ) -> None:
        """
        Calculate embeddings for all sequences in database that don't have embeddings.
        
        Args:
            db: Database connector
            batch_size: Batch size for processing
            model_name: Name of the model to use
            num_gpus: Number of GPUs to use (None = use all available)
            embedding_type: Type of embedding to calculate
        """
        # Retrieve sequences without embeddings
        query = """
        MATCH (p:Protein)
        WHERE p.embedding IS NULL AND p.sequence IS NOT NULL
        RETURN p.accession_id AS accession, p.sequence AS sequence
        """
        results = db.execute_read(query)
        data = [(result["accession"], result["sequence"]) for result in results]
        
        if not data:
            logger.info("No sequences to process.")
            return
        
        logger.info(f"Found {len(data)} sequences without embeddings")
        
        # Process using batch embedding method
        self.calculate_batch_embeddings(
            data=data,
            model_name=model_name,
            batch_size=batch_size,
            num_gpus=num_gpus,
            db=db,
            embedding_type=embedding_type
        )
    
    # Legacy compatibility methods (for backward compatibility with existing processor.py)
    def process_batches_on_gpu(
        self,
        data: List[tuple[str, str]],
        batch_size: int,
        model: Union[Any, DataParallel[Module]],
        tokenizer: Union[Any, None],
        db: DatabaseConnector,
        device: torch.device,
    ) -> None:
        """Legacy method for backward compatibility."""
        logger.warning("Using legacy process_batches_on_gpu method. Consider using calculate_batch_embeddings instead.")
        
        # Convert to new interface
        accessions, sequences = zip(*data)
        embedding_data = list(zip(accessions, sequences))
        
        # Use new method
        self.calculate_batch_embeddings(
            data=embedding_data,
            batch_size=batch_size,
            db=db
        )
    
    def get_batch_embeddings_unified(
        self,
        batch_sequences: List[str],
        model: Union[Any, DataParallel[Module]],
        tokenizer: Union[Any, None],
        device: torch.device = torch.device("cuda:0"),
        pool_embeddings: bool = True,
    ) -> List[NDArray[np.float64]]:
        """Legacy method for backward compatibility."""
        logger.warning("Using legacy get_batch_embeddings_unified method.")
        
        # Determine model type from the actual model instance
        base_model = model.module if isinstance(model, torch.nn.DataParallel) else model
        model_type = type(base_model).__name__
        
        # Map model class names to our model types
        if "ESMC" in model_type:
            embedding_model = ESMCEmbeddingModel("", device)
            embedding_model.model = base_model
            return embedding_model.get_batch_embeddings(batch_sequences, pool_embeddings)
        elif "ESM3" in model_type:
            embedding_model = ESM3EmbeddingModel("", device)
            embedding_model.model = base_model
            return embedding_model.get_batch_embeddings(batch_sequences, pool_embeddings)
        elif "T5Model" in model_type:
            embedding_model = ProtT5EmbeddingModel("", device)
            embedding_model.model = base_model
            embedding_model.tokenizer = tokenizer
            return embedding_model.get_batch_embeddings(batch_sequences, pool_embeddings)
        else:  # ESM-2 and other ESM models
            embedding_model = ESM2EmbeddingModel("", device)
            embedding_model.model = base_model
            embedding_model.tokenizer = tokenizer
            return embedding_model.get_batch_embeddings(batch_sequences, pool_embeddings)
    
    def calculate_single_sequence_embedding_last_hidden_state(
        self,
        sequence: str,
        device: torch.device = torch.device("cuda:0"),
        model_name: str = "facebook/esm2_t33_650M_UR50D",
    ) -> NDArray[np.float64]:
        """Legacy method for backward compatibility."""
        return self.calculate_single_embedding(sequence, model_name, "last_hidden_state", device)
    
    def calculate_single_sequence_embedding_all_layers(
        self,
        sequence: str,
        device: torch.device,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
    ) -> NDArray[np.float64]:
        """Legacy method for backward compatibility."""
        return self.calculate_single_embedding(sequence, model_name, "all_layers", device)
    
    def calculate_single_sequence_embedding_first_layer(
        self,
        sequence: str,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        device: torch.device = torch.device("cuda:0"),
    ) -> NDArray[np.float64]:
        """Legacy method for backward compatibility."""
        return self.calculate_single_embedding(sequence, model_name, "first_layer", device)
    
    def get_single_embedding_last_hidden_state(
        self, 
        sequence: str, 
        model: Any, 
        tokenizer: Any, 
        device: torch.device
    ) -> NDArray[np.float64]:
        """Legacy method for backward compatibility."""
        logger.warning("Using legacy get_single_embedding_last_hidden_state method.")
        return self._get_single_embedding_legacy(sequence, model, tokenizer, device, "last_hidden_state")
    
    def get_single_embedding_all_layers(
        self, 
        sequence: str, 
        model: Any, 
        tokenizer: Any, 
        device: torch.device
    ) -> NDArray[np.float64]:
        """Legacy method for backward compatibility."""
        logger.warning("Using legacy get_single_embedding_all_layers method.")
        return self._get_single_embedding_legacy(sequence, model, tokenizer, device, "all_layers")
    
    def get_single_embedding_first_layer(
        self, 
        sequence: str, 
        model: Any, 
        tokenizer: Any, 
        device: torch.device
    ) -> NDArray[np.float64]:
        """Legacy method for backward compatibility."""
        logger.warning("Using legacy get_single_embedding_first_layer method.")
        return self._get_single_embedding_legacy(sequence, model, tokenizer, device, "first_layer")
    
    def _get_single_embedding_legacy(
        self, 
        sequence: str, 
        model: Any, 
        tokenizer: Any, 
        device: torch.device,
        embedding_type: str
    ) -> NDArray[np.float64]:
        """Helper method for legacy single embedding methods."""
        # Determine model type and create appropriate embedding model
        base_model = model.module if isinstance(model, torch.nn.DataParallel) else model
        model_type = type(base_model).__name__
        
        if "ESMC" in model_type:
            embedding_model = ESMCEmbeddingModel("", device)
            embedding_model.model = base_model
        elif "ESM3" in model_type:
            embedding_model = ESM3EmbeddingModel("", device)
            embedding_model.model = base_model
        elif "T5Model" in model_type:
            embedding_model = ProtT5EmbeddingModel("", device)
            embedding_model.model = base_model
            embedding_model.tokenizer = tokenizer
        else:  # ESM-2 and other ESM models
            embedding_model = ESM2EmbeddingModel("", device)
            embedding_model.model = base_model
            embedding_model.tokenizer = tokenizer
        
        if embedding_type == "last_hidden_state":
            return embedding_model.get_single_embedding_last_hidden_state(sequence)
        elif embedding_type == "all_layers":
            return embedding_model.get_single_embedding_all_layers(sequence)
        elif embedding_type == "first_layer":
            return embedding_model.get_single_embedding_first_layer(sequence)
        else:
            raise ValueError(f"Unknown embedding_type: {embedding_type}")
    
    def cleanup(self) -> None:
        """Clean up all models and free memory."""
        for model in self._models.values():
            model.cleanup()
        self._models.clear()
        free_memory()


# Global processor instance
_processor = EmbeddingProcessor()


def get_processor() -> EmbeddingProcessor:
    """Get the global embedding processor instance."""
    return _processor 