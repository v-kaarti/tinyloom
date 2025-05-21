# this file contains our compression scheme with *zstd enforced* as a set of helpers

import argparse, math, sys, gzip
import loompy as lp, numpy as np, scipy.sparse as sp
from bitarray import bitarray
import zstandard as zstd


# get min bits (note that we use signed bits, but in hindsight this was not necessary, but the code works and i dont want to touch it
def _min_bits(x: int) -> int:
    if x == 0:
        return 1
    for b in range(1, 33):
        if -(1 << (b - 1)) <= x < (1 << (b - 1)):
            return b
    return 32


# encode array with our succint delta scheme
def dvld_encode(arr: np.ndarray) -> bytes:
    arr = arr.astype(np.int32, copy=False)
    n = len(arr)
    deltas, widths, prev, max_b = [], [], 0, 1
    for i, v in enumerate(arr):
        d = int(v) if i == 0 else int(v) - prev
        b = _min_bits(d)
        deltas.append(d)
        widths.append(b)
        if b > max_b:
            max_b = b
        prev = int(v)
    hdr_bits = 0 if max_b == 1 else math.ceil(math.log2(max_b))
    bits = bitarray(endian="big")
    bits.frombytes(n.to_bytes(4, "big"))
    bits.frombytes(bytes([hdr_bits]))
    mask_cache, fmt_cache = {}, {}
    for d, b in zip(deltas, widths):
        if hdr_bits:
            bits.extend(f"{b-1:0{hdr_bits}b}")
        if b not in mask_cache:
            mask_cache[b] = (1 << b) - 1
            fmt_cache[b] = f"0{b}b"
        bits.extend(format(d & mask_cache[b], fmt_cache[b]))
    pad = (8 - len(bits) % 8) % 8
    bits.extend("0" * pad)
    return bits.tobytes()


# decode array with our succint delta scheme
def dvld_decode(buf: bytes) -> np.ndarray:
    bits, pos = bitarray(endian="big"), 0
    bits.frombytes(buf)
    n = int.from_bytes(bits[pos : pos + 32].tobytes(), "big")
    pos += 32
    hdr_bits = bits[pos : pos + 8].tobytes()[0]
    pos += 8
    out, prev = np.empty(n, dtype=np.int32), 0
    for i in range(n):
        b = int(bits[pos : pos + hdr_bits].to01(), 2) + 1 if hdr_bits else 1
        pos += hdr_bits if hdr_bits else 0
        u = int(bits[pos : pos + b].to01(), 2)
        pos += b
        d = u - (1 << b) if u >= 1 << (b - 1) else u
        val = d if i == 0 else prev + d
        out[i] = val
        prev = val
    return out


# csr -> complete scheme
def csr_to_delta_blob(csr: sp.csr_matrix) -> bytes:
    indptr_b, indices_b, data_b = map(dvld_encode, (csr.indptr, csr.indices, csr.data))
    return (
        csr.shape[0].to_bytes(4, "big")
        + csr.shape[1].to_bytes(4, "big")
        + len(indptr_b).to_bytes(4, "big")
        + indptr_b
        + len(indices_b).to_bytes(4, "big")
        + indices_b
        + len(data_b).to_bytes(4, "big")
        + data_b
    )


# coplete scheme -> csr
def delta_blob_to_csr(blob: bytes) -> sp.csr_matrix:
    p = 0
    rows = int.from_bytes(blob[p : p + 4], "big")
    p += 4
    cols = int.from_bytes(blob[p : p + 4], "big")
    p += 4

    def nxt():
        nonlocal p
        ln = int.from_bytes(blob[p : p + 4], "big")
        p += 4
        out = blob[p : p + ln]
        p += ln
        return out

    indptr, indices, data = map(dvld_decode, (nxt(), nxt(), nxt()))
    return sp.csr_matrix((data, indices, indptr), shape=(rows, cols))


# give it a spin !!!
def main():
    pa = argparse.ArgumentParser(
        description="enter a loom file to compare our compression scheme vs "
    )
    pa.add_argument("loom_file", help="where is .loom :interrobang:")
    args = pa.parse_args()

    ds = lp.connect(args.loom_file)
    dense = ds[:, :]
    ds.close()
    dense = np.rint(dense).astype(np.int32, copy=False)

    csr = sp.csr_matrix(dense)
    blob = csr_to_delta_blob(csr)

    z_blob = zstd.ZstdCompressor(level=19).compress(blob)
    de_blob = zstd.ZstdDecompressor().decompress(z_blob)
    if (csr != delta_blob_to_csr(de_blob)).nnz:
        print("reconstruction kaboomed")
        sys.exit(1)

    comp_size = len(z_blob)
    with open(args.loom_file, "rb") as f:
        orig = f.read()
    gzip_size = len(gzip.compress(orig))

    print(f"Custom compression size (delta+zstd): {comp_size} bytes")
    print(f"Gzip + loom size: {gzip_size} bytes")
    if comp_size < gzip_size:
        print("Our algorithm is smaller.")
    elif comp_size > gzip_size:
        print("Gzip + loom is smaller.")
    else:
        print("Sizes are equal.")


# we imported into other stuff for testing, validation, etc.
if __name__ == "__main__":
    main()
