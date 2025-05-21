# Tinyloom: A Succinct Compression Scheme for scRNA Matrices
 ### Authors
- Vijaykaarti Sundarapandiyan
- Faraz Ghahremani
### Tinyloom

Tinyloom is a compression scheme for scRNA matrices. Tinyloom uses the facts that a) these matrices are highly sparse b) these matrices have entries that require very few bits to store, with exception and c) existing compression techniques address these features for compression efficiency independently. Our method Tinyloom creates a succinct delta encoding for CSR compression, nested inside of Zstd compression. It achieves significantly smaller final size as compared to Loom format + Gzip. Note the "loom" part of the name is purely for fun.

#### References
 - Huge thanks to the Linarsson Lab for Loompy and sample scRNA matrices
	- https://linnarssonlab.org/loompy/
	- http://loom.linnarssonlab.org/
- https://pypi.org/project/bitarray/ (for efficient bitpacking)
