"""Test script showing different return dtypes."""

import asyncio

from pyeed.embeddings.esm2_async import embed_proteins_async
from pyeed.embeddings.pooling import mean_pooling


async def test_return_dtypes():
    """Demonstrate different return dtype options."""

    sequences = ["MKTVRQERLK", "ARNDCEQGH"]
    accessions = ["P00001", "P00002"]

    print("=" * 70)
    print("Return Dtype Options")
    print("=" * 70)

    # 1. Float32 (default, most compatible)
    print("\n1. Float32 (default)")
    result = await embed_proteins_async(
        sequences,
        accessions,
        pooling_methods=mean_pooling,
        return_dtype="float32",
    )
    emb = result["mean_pooling"][0]
    print(f"   Dtype: {emb.dtype}, Shape: {emb.shape}")
    print(f"   Range: [{emb.min():.3f}, {emb.max():.3f}]")
    print(f"   ✓ Full precision, most compatible")

    # 2. Float16 (recommended for Milvus)
    print("\n2. Float16 (recommended for Milvus ANN)")
    result = await embed_proteins_async(
        sequences,
        accessions,
        pooling_methods=mean_pooling,
        model_dtype="bfloat16",
        return_dtype="float16",
    )
    emb = result["mean_pooling"][0]
    print(f"   Dtype: {emb.dtype}, Shape: {emb.shape}")
    print(f"   Range: [{emb.min():.3f}, {emb.max():.3f}]")
    print(f"   ✓ 2x memory savings, Milvus FLOAT16_VECTOR ready")

    # 3. BFloat16 (Google's format)
    print("\n3. BFloat16 (Google format)")
    result = await embed_proteins_async(
        sequences,
        accessions,
        pooling_methods=mean_pooling,
        model_dtype="bfloat16",
        return_dtype="bfloat16",
    )
    emb = result["mean_pooling"][0]
    print(f"   Dtype: {emb.dtype}, Shape: {emb.shape}")
    print(f"   Range: [{emb.min():.3f}, {emb.max():.3f}]")
    print(f"   ✓ Better dynamic range than float16")

    # 4. Int8 (compact storage, hardcoded quantization)
    print("\n4. Int8 (compact storage with hardcoded scale)")
    result = await embed_proteins_async(
        sequences,
        accessions,
        pooling_methods=mean_pooling,
        return_dtype="int8",
    )
    emb = result["mean_pooling"][0]
    print(f"   Dtype: {emb.dtype}, Shape: {emb.shape}")
    print(f"   Range: [{emb.min()}, {emb.max()}]")
    print(f"   Sample: {emb[:5]}")
    print(f"   ✓ 4x memory savings, hardcoded scale=127")
    print(f"   ⚠ Not for Milvus ANN (use float16 instead)")

    # 5. Int16 (more precision than int8)
    print("\n5. Int16 (more precision, still compact)")
    result = await embed_proteins_async(
        sequences,
        accessions,
        pooling_methods=mean_pooling,
        return_dtype="int16",
    )
    emb = result["mean_pooling"][0]
    print(f"   Dtype: {emb.dtype}, Shape: {emb.shape}")
    print(f"   Range: [{emb.min()}, {emb.max()}]")
    print(f"   Sample: {emb[:5]}")
    print(f"   ✓ 2x memory savings, hardcoded scale=32767")

    print("\n" + "=" * 70)
    print("Recommendations:")
    print("  • Milvus/Vector DB: float16 or bfloat16")
    print("  • Compact storage: int8 or int16")
    print("  • Custom quantization: get float32 and quantize yourself")
    print("=" * 70)


if __name__ == "__main__":
    asyncio.run(test_return_dtypes())






