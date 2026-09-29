# Example usage


## TL/DR of what this is

This directory contains example files for all possible usage. In this file, the complete recipe for obtaining the parsed documents is described.


## Creating the bwt and sa files for the reference

In order to obtain the FM-index of the file `example_reference.txt` which contains the reference sequence, this command was ran in the directory containing the implementation  [eBWT] (https://github.com/davidecenzato/PFP-eBWT)

```bash
./pfpebwt -w 4 -p 10 --GCA path/to/RLZ-parser/example/example_reference.txt 
```

Note that the flags `-w` and `-p` can be omitted for larger sequences. The flag `--GCA` inables the computation of GCA (SA) file, which is crucial for next steps. Mind, however, that the input sequence must be in format as described in the example file - the file contains two sequence, but the second is empty.

This command created the files `example_reference.txt.da`, `example_reference.txt.ebwt`, `example_reference.txt.gca`, `example_reference.txt.I`, `example_reference.txt.log`, `example_reference.txt.info`.


## Building the FM-index

Using the files  `example_reference.txt.ebwt` and `example_reference.txt.gca` the FM-index is build with these commands


```bash
./make_bwt --rle -i ./example/example_reference.txt.ebwt -sa ./example/example_reference.txt.gca -o ./example/example_reference.bwt
```

The first command creates file meant for transforming the bwt into bwt with metacharacters. The second command creates file ready to be parsed against.

Sometimes, this command is required to be ran more than once in order to create the correct index.


## Creating FM-indeces with metacharacters

To create an index with 4 columns of the bwt matrix instead of one, run this


```bash
./transform -i ./example/example_reference.bwt -sa ./example/example_reference.txt.ebwt -o ./example/example_reference_four.txt
```

## Parsing the sequences

The example sequence to parse `example_sequence1.txt` and `example_sequence2.fa` are listed in the file `example_sequences.txt`. To create the parses for each sequence file, run

```bash
./rlz_parser -r ./example/example_reference_four.bwt -s ./example/sequences.txt -o ./example/parses/
```


## Decompressing the parsed files

To recreate the original sequences of the compressed files listed in `compressed_sequences.txt`, run

```bash
./decompress -r ./example/example_reference.txt -s ./example/compressed_sequences.txt -o ./example/decompressed/
```
