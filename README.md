# Parser

Parser contains the original implementation of the LZ parser and the random-access RLZ index.

# Cached RLZ Parser

This repository also contains an experimental cached RLZ parser wrapper:

```text
parser/cached_rlz_parser.hpp
```

The cached parser wraps the RLZ factorisation process and adds a cache for repeated suffix-array interval refinements. When the same refinement is needed again, the parser can reuse the cached interval instead of repeating the binary searches.

The comparison program is:

```text
parser/compare_rlz_cache.cpp
```

This program runs both the original RLZ parser and the cached RLZ parser, checks that they produce the same factorisation output, and reports timing and cache statistics such as hits, misses, number of entries, and hit rate.

# RLZ_Varki

Here the FM-based parser can be found, as well as its parallelized form.

# Powered 

This repository contains the RLZ-parser that uses powered backward search in FM-index of the reference, with a modified rank structure.

# PT16mer_RLZ

Contains implementation of the SA-based parser enhanced with the prefix table.

In this repository, the version combining the prefix table and powered backward search can be found.


# Refgen

Refgen contains a program that generates a reference string of arbitrary integer type for an input of that type.

# SA

SA contains a simple implementation of the prefix-doubling algorithm to compute suffix arrays of arbitrary integer type for inputs of arbitrary integer types.

# Benchmark dataset

In our experminets we used fles from

https://github.com/koeppl/phoni (1000 sequences of human Chromosome 19)
https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=4932 (1764 assemblies of S. cerevisiae)
https://zenodo.org/records/6577997 (3682 E.coli genomes)
