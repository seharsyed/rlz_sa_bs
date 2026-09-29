# Relative Lempel-Ziv (RLZ)

```
ooooooooo.   ooooo         oooooooooooo 
`888   `Y88. `888'        d'""""""d888' 
 888   .d88'  888               .888P   
 888ooo88P'   888              d888'    
 888`88b.     888            .888P      
 888  `88b.   888       o   d888'    .P 
o888o  o888o o888ooooood8 .8888888888P  
                                        v1.1.0
```

## Description

This software computes the Relative Lempel Ziv (RLZ) parse of the target sequence file using a reference file. By default, the software does character-level encoding.
However, character-level encoding can fail if the reference file does not contain all the unique characters in the sequence file. We have also provided an option to 
do the encoding at the bit-level. Doing the encoding this way prevents the issue that occurs when the sequence contains a character that is not present in the reference.
However, the parse might be larger than the character-level parse due to encoding the bit representation of the sequence file. The decompression expects that the files were originally ASCII (8 bit) encoded.

The software performs pattern matching by streaming the sequence file char by char (or bit by bit) against the FM-index of the reversed reference. This approach enables the software to simulate forward matching of the sequence file against the reference. The correct reference positions of the matches are obtained by applying an involution to the suffix array positions retrieved from the FM-index (which is built on the reversed reference text), a constant-time operation. 

Additionally, an executable has been added that supports sorting RLZ sequences directly in their compressed form. It takes an RLZ parse as input and produces the suffix array of the original sequence. Note that the suffix array is 0-based and does not include a terminal sentinel symbol, unlike conventional constructions.   

## RLZ Algorithm Workflow

To compress the target sequence file relative to a reference file, the software follows these steps:

Bit-level encoding:

0. Convert both the reference and sequence files to their "binary" representation.

Common steps (for all encoding types):

1. Reverse the reference sequence.

2. Build an FM-index from the reversed reference (in bit or character form).

3. Perform "backwards" (forward) match:
    Match each character or bit of the sequence against the reversed reference using the FM-index's backward matching capabilities (to simulate forward matching).
   
    3a. If a match is found, check if the next character or bit also matches.
   
    3b. If a match is found and it's the end of the sequence, apply involution to the suffix array (SA) position, subtract the length of the match from the position, and record the (pos, len) pair.
   
    3c. If a mismatch occurs, apply involution to the SA position, subtract the length of the match from the position, record the (prev pos, len - 1) pair, and restart the search from the mismatched character or bit.
   
4. Write all (pos, len) pairs sequentially to a file. This constitutes the RLZ parse.

> [!NOTE]
> The RLZ parse is computed with respect to the reference file.

## Prerequisites

- [CMake](https://cmake.org/) 3.15 or higher.
- GCC 9+
- C++17-compatible compiler.
- OpenMP

## Getting Started

### Building the Project

```
git clone https://github.com/rvarki/RLZ.git
cd RLZ
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make -j
```

### Running the project

After building the project, an executable named `rlz` will be created in the build directory. Run it with:
```
./rlz compress -r [reference file] -s [file to compress] [options] 
```

Additionally, an executable named `sort` will be created in the build directory. Run it with:

```
./sort rlz -r [reference file] -p [RLZ parse file] [options]
```

### Default Compression Example

In this section, we will go through a small example using the default character-level compression. In the data/dna directory, we have provided an example reference and target sequence file that were derived from DNA FASTA files.

1. To compress the sequence file with character compression, run the following command from the build directory

```
./rlz compress -r ../data/dna/dna_ref.txt -s ../data/dna/dna_seq.txt
```
This command will produce the following file in the data/dna directory: `dna_seq.txt.rlz`. The .rlz file contains the RLZ parse.

> [!NOTE]
> Multithreading is supported in the compression step with the -t [num. of threads] option which can significantly make the compression step faster. However, the RLZ parse is slightly different from what you would get if you run with a single thread. The reason is we cannot identify phrases that span where the file was split. Potentially might add an additional thread number of parse entries that would not exist if you ran with a single thread.

2. To decompress the file, run the following command

```
./rlz decompress -r ../data/dna/dna_ref.txt -p ../data/dna/dna_seq.txt.rlz 
```
This command should produce a file called `dna_seq.txt.out` in the data/dna directory. This is the decompressed sequence file.

3. Check to see if the file decompressed correctly
```
diff ../data/dna/dna_seq.txt ../data/dna/dna_seq.txt.out
```

There should be no output from this command if compressed and decompressed correctly. 


### Bit Compression Example

In this section, we will go through a small example using the bit-level compression option. In the data/english directory, we have provided an example reference and target sequence file that were derived from the English text in the [Pizza&Chili Corpus](https://pizzachili.dcc.uchile.cl/texts/nlang/).

1. To compress the sequence file, run the following command from the build directory

```
./rlz compress -r ../data/english/english_ref.txt -s ../data/english/english_seq.txt --bit
```

This command will produce the following file in the data/english directory: `english_seq.txt.rlz`. The .rlz file contains the RLZ parse.

2. To decompress the file, run the following command

```
./rlz decompress -r ../data/english/english_ref.txt -p ../data/english/english_seq.txt.rlz --bit
```
This command should produce a file called `english_seq.txt.out` in the data/english directory. This is the decompressed sequence file.

3. Check to see if the file decompressed correctly
```
diff ../data/english/english_seq.txt ../data/english/english_seq.txt.out
```

There should be no output from this command if compressed and decompressed correctly. 

### Sorting RLZ Parse Example

In this section, we show how to use the sort exectuable. We assume that you have already run the commands from the Default Compression example.

1. To sort the RLZ parse file, run the following command from the build directory

```
./sort rlz -r ../data/dna/dna_ref.txt -p ../data/dna/dna_seq.txt.rlz --naive
```

This command produces a file named `dna_seq.txt.rlz.sa` in the data/dna directory. The file contains the suffix array of the original text, computed by decomposing complete RLZ factors into their constituent complete and incomplete factors before sorting. Additional sorting options are available and can be viewed by passing the `-h` flag.

2. To sort the original sequence file, run the following command

```
./sort text -s ../data/dna/dna_seq.txt 
```

This command should produce a file called `dna_seq.txt.sa` in the data/dna directory. This contains the suffix array of the orginal text. This was found by directly sorting the sequence file itself.

3. Check if the files are the same

```
diff ../data/dna/dna_seq.txt.rlz.sa ../data/dna/dna_seq.txt.sa 
```

There should be no output from this command if the sorting was done correctly.

> [!NOTE]
> The suffix array is 0-based and does not assume the presence of a terminal sentinel symbol, unlike most constructions.

The previous sort command expands the RLZ representation into complete and incomplete factors to construct the full suffix array. Alternatively, it can construct a partial suffix array containing only suffixes represented by complete RLZ factors.

4. To sort only the complete RLZ factors, run the following command

```
./sort rlz -r ../data/dna/dna_ref.txt -p ../data/dna/dna_seq.txt.rlz --factors-only
```

This will overwrite the suffix array produced by the last sort command. 

5. To verify that the partial suffix array is a subset of the full suffix array, run the following command

```
awk 'NR==FNR{a[++n]=$1; next} $1==a[i+1]{i++} END{exit (i<n)}' ../data/dna/dna_seq.txt.rlz.sa ../data/dna/dna_seq.txt.sa  && echo "Ordered subset" || echo "Not an ordered subset"
```

The above command should output: "Ordered Subset"


> [!NOTE]
> To get more information from the tools. Run the commands with the -v option.

## License

This project is licensed under the GNU License - see the [LICENSE](https://github.com/rvarki/rlz/blob/main/LICENSE) file for details

## Dependencies
- [SDSL](https://github.com/simongog/sdsl-lite)
- [CLI11](https://github.com/CLIUtils/CLI11)
- [spdlog](https://github.com/gabime/spdlog)

## Acknowledgements

- [S. Kuruppu, S. J. Puglisi and J. Zobel, Relative Lempel-Ziv Compression of Genomes for Large-Scale Storage and Retrieval](http://dx.doi.org/10.1007/978-3-642-16321-0_20). Proc. 17th International Symposium on String Processing and Information Retrieval (SPIRE 2010) Lecture Notes in Computer Science, Volume 6393, (2010) pp. 201-206. 
