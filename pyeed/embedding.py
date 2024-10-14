import gc

import torch
from transformers import EsmModel, EsmTokenizer


def get_batch_embeddings(sequences: list[str], batch_size: int = 16):
    # Load the ESM2 model and tokenizer
    model_name = "facebook/esm2_t33_650M_UR50D"
    model = EsmModel.from_pretrained(model_name)
    tokenizer = EsmTokenizer.from_pretrained(model_name)

    # Check if MPS (Metal Performance Shaders) is available and use it
    device = (
        torch.device("mps") if torch.backends.mps.is_built() else torch.device("cpu")
    )
    model = model.to(device)

    embedding_list = []
    model.eval()

    with torch.no_grad():
        # Process sequences in batches
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i : i + batch_size]

            # Tokenize the input sequences (must be a list of strings)
            inputs = tokenizer(
                batch, padding=True, truncation=True, return_tensors="pt"
            ).to(device)

            # Get model outputs
            outputs = model(**inputs)
            embeddings = outputs.last_hidden_state

            # Process each sequence in the batch
            for j in range(len(batch)):
                valid_token_mask = inputs["attention_mask"][j].bool()
                seq_embeddings = embeddings[j][valid_token_mask].mean(dim=0).cpu()
                embedding_list.append(seq_embeddings)

    return embedding_list


def free_memory():
    gc.collect()  # Python garbage collection
    if torch.backends.mps.is_built():
        torch.mps.empty_cache()
    elif torch.cuda.is_available():
        torch.cuda.empty_cache()


if __name__ == "__main__":
    free_memory()
