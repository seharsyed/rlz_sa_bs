# Fast RLZ-parser that uses powered backward search in FM-index

## TL/DR of what this is

This parser compresses sequences against the reference using a greedy strategy - by searching in a loop for the largest matching suffix of a prefix.
The largest parse - suffix of a prefix of the sequences - is computed using backward search in the FM-index of the reference. Moreover, the backward search is powered by joining several rank computations into one.

## Building indexes

When starting the project, run

```bash
make
```

Given a plain text BWT file `/path/to/bwt.txt` and a file containing its suffix array `/path/to/sa.txt`, to make the FM-index `bwt.bwt` with associated data `bwt_data.bwt` run the following in the repository root

```bash
$ ./make_bwt --rle -i /path/to/bwt.txt -sa /path/to/sa.txt -o bwt.bwt
```

This will create an index with block size of 256 symbols by splitting runs as necessary to store pairs (symbol, run length) in two bytes per run. Run `./make_bwt` for more information on how to generate different versions of the indexes.

For example,

```bash
$ ./make_bwt -i /path/to/bwt.txt -sa /path/to/sa.txt -o bwt.bwt
```

will create an index with no run length compression, where the ranks are computed using the extended [QuadRank](https://arxiv.org/abs/2602.04103) approach. In order for this to work correctly, the bwt text file must contain at most 4 various characters - we advise to create this file using [eBWT] (https://github.com/davidecenzato/PFP-eBWT) strategy, having the flag --GCA set for the SA computation.

For different block sizes, change the SMALL_BLOCK_SIZE in the Makefile.

Please run the command again if the number of dense blocks shows 0 for bigger file and the flag -r is on.

## Transforming BWT to meta characters

To transform the RLBWT file `/path/to/bwt.bwt` (the input must be run-length encoded) to a BWT file with metacharacters containing 4 letters, run

```bash
$ ./transform -i /path/to/bwt.rlbwt -o bwt_powered.txt
```
Run `./transform` for more information. To build the FM-index over this transformed file, simply follow the previous instructions, with `powered_bwt.txt` file as input. The SA file remains the same. Or run

```bash
$ ./transform -b -i /path/to/bwt.bwt -sa /path/to/sa.txt -o powered_bwt.txt
```

and the powered_bwt.bwt will contain the bwt structure over th meta-alphabet.

## RLZ-parser

Given the powered index of a reference file `powered_bwt.bwt` and a file containing the file names of the sequences to parse `sequences.txt` run:

```bash
$ ./rlz_parser -r powered_bwt.bwt -s sequences.txt -o ../parses/
```

This will create a binary file for each sequence containing the parses as pairs of 64-bit words (length, location). The binary files will be stored in the `../parses/` directory.
Run `./rlz_parser` for information.


## Decompression

Given the reference file `reference.txt` and a file containing the file names of the compressed sequences `compressed_sequences.txt` run:

```bash
$ ./decompress -r reference.txt -s compressed_sequences.txt -o ../decompressed/
```

This will recreate the sequences and save them in files in the suggested output directory.

## Example

To see an example of the files and commands they were created with see the `example` directory.

## Requirements

Compilation and execution has been tested on modern x86-64 systems and GCC supporting `-std=c++2a`. Code should be compatible with other compilers, but is not expected to compile correctly on compilers not supporting C++ 20. The project should also work on new ARM based apple systems.

