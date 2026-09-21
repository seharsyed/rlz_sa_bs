#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#include "rmq_tree.h"
#ifndef _OPENMP
#define _OPENMP ;
#endif
#include "libsais/src/libsais.h"

#define sequenceSeparator 2

// function that prints the instructions for using the tool
void print_help(char **argv) {
  std::cout << "Usage: " << argv[0]
            << " -r reference_file -c collection_file [options]" << std::endl;
  std::cout
      << "  Options: " << std::endl
      << "    -o <output file> : basename for the output files (default: "
         "reference_collection)"
      << std::endl
      << "    -b <buffer size> : size of the buffer for holding the temporary "
         "MS entries in MB (default: 100)"
      << std::endl
      << "    -t <number of threads> : number of threads to use (default: 8)"
      << std::endl
      << "    -n : do not write the .ms output file"
      << std::endl
      << "    -h : print this help" << std::endl;
  exit(-1);
}

struct Args {
  std::string reference;
  std::string collection;
  std::string outfile;
  int bufferText = 104857600;
  int numThreads = 8;
  bool noOutput = false;
};

// function for parsing the input arguments
void parseArgs(int argc, char **argv, Args &arg) {
  int c;
  // extern int optind;

  puts("==== Command line:");
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  puts("");

  std::string sarg;
  while ((c = getopt(argc, argv, "r:c:o:b:t:nh")) != -1) {
    switch (c) {
    case 'r':
      arg.reference.assign(optarg);
      break;
    case 'c':
      arg.collection.assign(optarg);
      break;
    case 'o':
      arg.outfile.assign(optarg);
      break;
    case 'b':
      arg.bufferText = 1048576 * std::stoi(optarg);
      break;
    case 't':
      arg.numThreads = std::stoi(optarg);
      break;
    case 'n':
      arg.noOutput = true;
      break;
    case 'h':
      print_help(argv);
      exit(-1);
      // fall through
    default:
      std::cout << "Unknown option \"" << (char)c << "\". Use -h for help."
                << std::endl;
      exit(-1);
    }
  }

  // the only input parameter is the file name
  if (argc < 4) {
    std::cout << "Invalid number of arguments" << std::endl;
    print_help(argv);
  }
  // set output files basename
  if (arg.outfile == "")
    arg.outfile = arg.reference + "_" + arg.collection;

  std::cout << "==== Parameters:" << std::endl;
  std::cout << "Reference: " << arg.reference << std::endl;
  std::cout << "Collection: " << arg.collection << std::endl;
  std::cout << "Output basename: " << arg.outfile << std::endl;
  std::cout << "Buffer size: " << arg.bufferText << std::endl;
  std::cout << "Number of threads: " << arg.numThreads << std::endl;
  std::cout << "====" << std::endl;
}

/* MS PROGRAM */
static std::string reference;
static int32_t *_SA;
static int32_t *_ISA;
static int32_t *_LCP;
static int32_t *_LRF;
static rmq_tree<int32_t> *_rmq;
static uint32_t _sizeReference;
static uint64_t _sn;

void constructISA(int32_t *sa, int32_t *isa, uint32_t _sizeReference) {
  fprintf(stdout, "\tComputing ISA...\n");
  for (uint32_t i = 0; i < _sizeReference; i++) {
    isa[sa[i]] = i;
  }
}

void constructISA_parallel(int32_t *sa, int32_t *isa, uint32_t _sizeReference) {
  fprintf(stdout, "\tComputing ISA...\n");
#pragma omp parallel for
  for (uint32_t i = 0; i < _sizeReference; i++) {
    isa[sa[i]] = i;
  }
}

std::pair<int, int> adjustInterval(int lo, int hi, int offset) {
  int psv = _rmq->psv(lo, offset);
  if (psv == -1) {
    psv = 0;
  } else {
    // psv++;
  }
  int nsv = _rmq->nsv(hi + 1, offset);
  if (nsv == -1) {
    nsv = _sizeReference - 1;
  } else {
    nsv--;
  }
  return {psv, nsv};
}

std::pair<int, int> contractLeft(int lo, int hi, int offset) {
  uint32_t suflo = _SA[lo];
  uint32_t sufhi = _SA[hi];
  if (suflo == _sizeReference - 1 ||
      sufhi == _sizeReference - 1) { // if true we must be at depth 1
    return std::make_pair(0, _sizeReference - 1); // root
  }
  uint32_t tmplo = _ISA[suflo + 1];
  uint32_t tmphi = _ISA[sufhi + 1];
  return adjustInterval(tmplo, tmphi, offset);
}

// Returns the leftmost occurrence of the element if it is present or (if not
// present) then returns -(x+1) where x is the index at which the key would be
// inserted into the array: i.e., the index of the first element greater than
// the key, or hi+1 if all elements in the array are less than the specified
// key.
inline int32_t binarySearchLB(int32_t lo, int32_t hi, uint32_t offset,
                              uint8_t c) {
  int32_t low = lo, high = hi;
  while (low <= high) {
    int32_t mid = (low + high) >> 1;
    uint8_t midVal = reference[_SA[mid] + offset];
    if (midVal < c) {
      low = mid + 1;
      __builtin_prefetch(&reference[_SA[(low + high) >> 1] + offset], 0, 0);
    } else if (midVal > c) {
      high = mid - 1;
      __builtin_prefetch(&reference[_SA[(low + high) >> 1] + offset], 0, 0);
    } else { // midVal == c
      if (mid == lo)
        return mid; // leftmost occ of key found
      uint8_t midValLeft = reference[_SA[mid - 1] + offset];
      if (midValLeft == midVal) {
        high = mid - 1; // discard mid and the ones to the right of mid
        __builtin_prefetch(&reference[_SA[(low + high) >> 1] + offset], 0, 0);
      } else {      // midValLeft must be less than midVal == c
        return mid; // leftmost occ of key found
      }
    }
  }
  return -(low + 1); // key not found.
}

inline int32_t binarySearchRB(int32_t lo, int32_t hi, uint32_t offset,
                              uint8_t c) {
  int32_t low = lo, high = hi;
  while (low <= high) {
    int32_t mid = (low + high) >> 1;
    uint8_t midVal = reference[_SA[mid] + offset];
    if (midVal < c) {
      low = mid + 1;
      __builtin_prefetch(&reference[_SA[(low + high) >> 1] + offset], 0, 0);
    } else if (midVal > c) {
      high = mid - 1;
      __builtin_prefetch(&reference[_SA[(low + high) >> 1] + offset], 0, 0);
    } else { // midVal == c
      if (mid == hi)
        return mid; // rightmost occ of key found
      uint8_t midValRight = reference[_SA[mid + 1] + offset];
      if (midValRight == midVal) {
        low = mid + 1; // discard mid and the ones to the left of mid
        __builtin_prefetch(&reference[_SA[(low + high) >> 1] + offset], 0, 0);
      } else {      // midValRight must be greater than midVal == c
        return mid; // rightmost occ of key found
      }
    }
  }
  return -(low + 1); // key not found.
}

void computeMatchingFactor(const std::string &collection, uint64_t i,
                           int32_t *pos, int32_t *len, int32_t &leftB,
                           int32_t &rightB) {
  uint32_t offset = *len;
  uint64_t j = i + offset;

  int32_t nlb = leftB, nrb = rightB, maxMatch;
  unsigned int match = _SA[nlb];
  while (j < collection.size()) { // scans the string from j onwards until a
                                  // maximal prefix is
                                  // found between reference and collection
    if (nlb == nrb) {
      if (reference[_SA[nlb] + offset] != collection[j]) {
        break;
      }
      leftB = nlb;
      rightB = nrb;
      maxMatch = nlb;
    } else { // refining the bucket in which the match is found, from left and
             // then from right
      nlb = binarySearchLB(nlb, nrb, offset, collection[j]);
      if (nlb < 0) {
        // auto tmp = true;
        maxMatch = -nlb - 1;
        if (maxMatch == nrb + 1) {
          maxMatch--;
        }
        match = _SA[maxMatch];
        break;
      }
      nrb = binarySearchRB(nlb, nrb, offset, collection[j]);

      leftB = nlb;
      rightB = nrb;
    }
    match = _SA[nlb];
    j++;
    offset++;
  }
  *pos = match;
  *len = offset;
}

void loadReferenceAndComputeDS(std::string refFileName) {
  errno = 0;
  FILE *infileRef = fopen(refFileName.c_str(), "r");
  if (!infileRef) {
    printf("Error opening file of base sequence %s, errno=%d\n",
           refFileName.c_str(), errno);
    exit(1);
  }
  printf("About to read ref\n");

  unsigned int n = 0;
  fseek(infileRef, 0, SEEK_END);
  n = ftell(infileRef) / sizeof(uint8_t);
  printf("Reference length: %u\n", n);
  fseek(infileRef, 0, SEEK_SET);
  if (n) {
    char *firstChar = new char[1];
    if (fread(firstChar, sizeof(char), 1, infileRef) != 1) {
      printf("Error reading first character of reference file\n");
      exit(1);
    }
    if (firstChar[0] == '>') {
      printf("Reading reference from fasta file\n");
      reference.reserve(n);
      std::ifstream streamInfileRef(refFileName);
      std::string line, content;
      while (std::getline(streamInfileRef, line).good()) {
        if (line.empty() || line[0] == '>') {
          reference += content;
          // reference += sequenceSeparator;
          std::string().swap(content);
        } else if (!line.empty()) {
          content += line;
        }
      }
      if (content.size())
        reference += content;
      std::string().swap(content);
      reference.resize(reference.size());
      streamInfileRef.close();
    } else {
      printf("Reading reference from binary file\n");
      fseek(infileRef, 0, SEEK_SET);
      reference.resize(n);
      printf("Resized string\n");
      if (n != fread(&reference[0], sizeof(uint8_t), n, infileRef)) {
        printf("Error reading %u bytes from file %s\n", n, refFileName.c_str());
        exit(1);
      }
      printf("Read everything from file\n");
    }
  } else {
    printf("Reference file is empty!\n");
    exit(1);
  }
  if (reference[reference.size() - 1] == '\n')
    reference.erase(reference.size() - 1, 1); // remove newline character

  reference += (char)0;
  printf("Appended 0 to string\n");

  _sizeReference = reference.size();
  fclose(infileRef);
  printf("Reference (size = %lu):\n\t", reference.size());

  auto t01 = std::chrono::high_resolution_clock::now();
  int32_t *sa = new int32_t[reference.size()];
  libsais(
      reinterpret_cast<unsigned char *>(const_cast<char *>(reference.c_str())),
      sa, reference.size(), 0, NULL);
  auto t02 = std::chrono::high_resolution_clock::now();
  printf(
      "Computing SA done in %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count());

  t01 = std::chrono::high_resolution_clock::now();
  _SA = sa;
  _ISA = new int32_t[reference.size()];
  _LRF = new int32_t[reference.size()];
  _LCP = new int32_t[reference.size()];
  t01 = std::chrono::high_resolution_clock::now();
  constructISA(_SA, _ISA, reference.size());
  t02 = std::chrono::high_resolution_clock::now();
  printf(
      "Computing ISA done in %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count());

  libsais_plcp(
      reinterpret_cast<unsigned char *>(const_cast<char *>(reference.c_str())),
      _SA, _LRF, reference.size());
  libsais_lcp(_LRF, _SA, _LCP, reference.size());

  for (uint32_t i = 0; i < reference.size(); i++) {
    _LRF[i] = std::max(_LCP[_ISA[i]], _LCP[_ISA[i] + 1]);
  }
  t02 = std::chrono::high_resolution_clock::now();
  printf(
      "Computing LCP and PLCP done in %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count());

  t01 = std::chrono::high_resolution_clock::now();
  fprintf(stderr, "\tComputing RMQ...\n");
  _rmq = new rmq_tree((int *)_LCP, reference.size(), 7);
  t02 = std::chrono::high_resolution_clock::now();
  printf(
      "Computing RMQ done in %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count());

  printf("Done!\n");
}

void loadReferenceAndComputeDS_parallel(std::string refFileName,
                                        uint32_t numThreads = 1) {
  omp_set_num_threads(numThreads);
  errno = 0;
  FILE *infileRef = fopen(refFileName.c_str(), "r");
  if (!infileRef) {
    printf("Error opening file of base sequence %s, errno=%d\n",
           refFileName.c_str(), errno);
    exit(1);
  }
  printf("About to read ref\n");

  unsigned int n = 0;
  fseek(infileRef, 0, SEEK_END);
  n = ftell(infileRef) / sizeof(uint8_t);
  printf("Reference length: %u\n", n);
  fseek(infileRef, 0, SEEK_SET);
  if (n) {
    char *firstChar = new char[1];
    if (fread(firstChar, sizeof(char), 1, infileRef) != 1) {
      printf("Error reading first character of reference file\n");
      exit(1);
    }
    if (firstChar[0] == '>') {
      printf("Reading reference from fasta file\n");
      reference.reserve(n);
      std::ifstream streamInfileRef(refFileName);
      std::string line, content;
      while (std::getline(streamInfileRef, line).good()) {
        if (line.empty() || line[0] == '>') {
          reference += content;
          // reference += sequenceSeparator;
          std::string().swap(content);
        } else if (!line.empty()) {
          content += line;
        }
      }
      if (content.size())
        reference += content;
      std::string().swap(content);
      reference.resize(reference.size());
      streamInfileRef.close();
    } else {
      printf("Reading reference from binary file\n");
      fseek(infileRef, 0, SEEK_SET);
      reference.resize(n);
      printf("Resized string\n");
      if (n != fread(&reference[0], sizeof(uint8_t), n, infileRef)) {
        printf("Error reading %u bytes from file %s\n", n, refFileName.c_str());
        exit(1);
      }
      printf("Read everything from file\n");
    }
  } else {
    printf("Reference file is empty!\n");
    exit(1);
  }
  if (reference[reference.size() - 1] == '\n')
    reference.erase(reference.size() - 1, 1); // remove newline character

  reference += (char)0;
  printf("Appended 0 to string\n");

  _sizeReference = reference.size();
  fclose(infileRef);
  printf("Reference (size = %lu):\n\t", reference.size());

  auto t01 = std::chrono::high_resolution_clock::now();
  int32_t *sa = new int32_t[reference.size()];
  libsais_omp(
      reinterpret_cast<unsigned char *>(const_cast<char *>(reference.c_str())),
      sa, reference.size(), 0, NULL, numThreads);
  auto t02 = std::chrono::high_resolution_clock::now();
  printf(
      "Computing SA done in %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count());

  t01 = std::chrono::high_resolution_clock::now();
  _SA = sa;
  _ISA = new int32_t[reference.size()];
  _LRF = new int32_t[reference.size()];
  _LCP = new int32_t[reference.size()];
  t01 = std::chrono::high_resolution_clock::now();
  constructISA_parallel(_SA, _ISA, reference.size());
  t02 = std::chrono::high_resolution_clock::now();
  printf(
      "Computing ISA done in %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count());

  libsais_plcp_omp(
      reinterpret_cast<unsigned char *>(const_cast<char *>(reference.c_str())),
      _SA, _LRF, reference.size(), numThreads);
  libsais_lcp_omp(_LRF, _SA, _LCP, reference.size(), numThreads);

#pragma omp parallel for
  for (uint32_t i = 0; i < reference.size(); i++) {
    _LRF[i] = std::max(_LCP[_ISA[i]], _LCP[_ISA[i] + 1]);
  }
  t02 = std::chrono::high_resolution_clock::now();
  printf(
      "Computing LCP and PLCP done in %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count());

  t01 = std::chrono::high_resolution_clock::now();
  fprintf(stderr, "\tComputing RMQ...\n");
  _rmq = new rmq_tree((int *)_LCP, reference.size(), 7);
  t02 = std::chrono::high_resolution_clock::now();
  printf(
      "Computing RMQ done in %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count());

  printf("Done!\n");
}

void computeMatchingStatistics(std::string collFileName, uint64_t bufferText,
                               std::string outputFileName = "") {

  std::ofstream streamOutfile(outputFileName, std::ios::out | std::ios::binary);
  uint64_t time = 0;
  auto t1 = std::chrono::high_resolution_clock::now();
  FILE *infile = fopen(collFileName.c_str(), "r");
  if (!infile) {
    fprintf(stderr, "Error opening file of sequence (%s)\n",
            collFileName.c_str());
    exit(1);
  }
  uint64_t _sizeCollection = 0;
  fseek(infile, 0L, SEEK_END);
  _sizeCollection = ftello(infile) / sizeof(uint8_t);
  // uint64_t prefix = std::min(prefix, _sizeCollection);
  _sn = _sizeCollection;
  printf("About to read sequence of size %lu with prefix %lu\n",
         _sizeCollection, _sn);
  fclose(infile);

  uint64_t capBufferWrite =
      (bufferText - _sizeReference) / sizeof(std::pair<uint32_t, uint32_t>) + 1;

  std::vector<std::pair<uint32_t, uint32_t>> matchingStatistics;
  matchingStatistics.reserve(2 * capBufferWrite);

  std::ifstream streamInfile(collFileName, std::ios::in);
  std::string line, content;
  content.reserve(_sizeReference * 2);
  uint64_t charactersRead = 0;
  uint64_t D = 1;
  while (std::getline(streamInfile, line)) {
    if (line.empty() || line[0] == '>') {
      content += sequenceSeparator;
      charactersRead++;
      int64_t i = 0;
      int32_t leftB = 0;
      int32_t rightB = _sizeReference - 1;
      int32_t pos = _sizeReference - 1, len = 0;
      D++;
      while (i < (static_cast<int64_t>(content.size())) -
                     1) { // stop before the separator
        computeMatchingFactor(content, i, &pos, &len, leftB, rightB);
        matchingStatistics.push_back(std::make_pair(pos, len));
        if (matchingStatistics.size() >= capBufferWrite) {
          streamOutfile.write(
              reinterpret_cast<char *>(matchingStatistics.data()),
              matchingStatistics.size() *
                  sizeof(std::pair<uint32_t, uint32_t>));
          matchingStatistics.clear();
          matchingStatistics.reserve(capBufferWrite);
        }
        len--;

        if (leftB == rightB) {
          while (len > _LRF[pos + 1]) {
            time++;
            matchingStatistics.push_back(std::make_pair(pos + 1, len));
            i++;
            len--;
            pos++;
          }
          if (matchingStatistics.size() >= capBufferWrite) {
            streamOutfile.write(
                reinterpret_cast<char *>(matchingStatistics.data()),
                matchingStatistics.size() *
                    sizeof(std::pair<uint32_t, uint32_t>));
            matchingStatistics.clear();
            matchingStatistics.reserve(capBufferWrite);
          }
          std::pair<int, int> interval =
              adjustInterval(_ISA[pos + 1], _ISA[pos + 1], len);

          leftB = interval.first;
          rightB = interval.second;
        } else {
          std::pair<int, int> interval = contractLeft(leftB, rightB, len);
          leftB = interval.first;
          rightB = interval.second;
        }
        i++;
      }
      len = 0;
      pos = _sizeReference - 1;
      matchingStatistics.push_back(std::make_pair(pos, len));
      content.erase(0, content.size());
    } else if (!line.empty()) {
      charactersRead += line.size();
      content += line;
    }
  }
  streamInfile.close();
  if (content.size() != 0) {
    // printf("Last sequence\n");
    content += sequenceSeparator;
    charactersRead++;
    if (charactersRead < _sizeCollection - 1) {
      _sizeCollection = charactersRead;
    }
    D++;
    int64_t i = 0;
    int32_t leftB = 0;
    int32_t rightB = _sizeReference - 1;
    int32_t pos = _sizeReference - 1, len = 0;
    while (i < (static_cast<int64_t>(content.size()))) {
      if (content[i] == sequenceSeparator) {
        leftB = 0;
        rightB = _sizeReference - 1;
        len = 0;
        pos = _sizeReference - 1;
        matchingStatistics.push_back(std::make_pair(pos, len));
      } else {
        computeMatchingFactor(content, i, &pos, &len, leftB, rightB);
        matchingStatistics.push_back(std::make_pair(pos, len));
        if (matchingStatistics.size() >= capBufferWrite) {
          streamOutfile.write(
              reinterpret_cast<char *>(matchingStatistics.data()),
              matchingStatistics.size() *
                  sizeof(std::pair<uint32_t, uint32_t>));
          matchingStatistics.clear();
          matchingStatistics.reserve(capBufferWrite);
        }
        // }
        len--;

        if (leftB == rightB) {
          while (len > _LRF[pos + 1]) {
            time++;
            matchingStatistics.push_back(std::make_pair(pos + 1, len));
            i++;
            len--;
            pos++;
          }
          if (matchingStatistics.size() >= capBufferWrite) {
            streamOutfile.write(
                reinterpret_cast<char *>(matchingStatistics.data()),
                matchingStatistics.size() *
                    sizeof(std::pair<uint32_t, uint32_t>));
            matchingStatistics.clear();
            matchingStatistics.reserve(capBufferWrite);
          }
          std::pair<int, int> interval =
              adjustInterval(_ISA[pos + 1], _ISA[pos + 1], len);
          leftB = interval.first;
          rightB = interval.second;
        } else {
          std::pair<int, int> interval = contractLeft(leftB, rightB, len);
          leftB = interval.first;
          rightB = interval.second;
        }
      }
      i++;
    }
    std::string().swap(content);
  }
  delete[] _SA;
  delete[] _ISA;
  delete[] _LCP;
  delete[] _LRF;
  delete _rmq;
  if (matchingStatistics.size() != 0) {
    streamOutfile.write(reinterpret_cast<char *>(matchingStatistics.data()),
                        matchingStatistics.size() *
                            sizeof(std::pair<uint32_t, uint32_t>));
    matchingStatistics.clear();
  }
  streamOutfile.close();

  auto t2 = std::chrono::high_resolution_clock::now();
  printf(
      "Time to compute matching statistics: %ld ms\n",
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
  printf("Times used euristics %ld\n", time);
}

void computeMatchingStatisticsParallel(std::string collFileName,
                                       uint64_t bufferText,
                                       std::string outputFileName = "",
                                       uint32_t numThreads = 16) {

  std::cout << "Number of threads: " << numThreads << '\n';
  omp_set_num_threads(numThreads);
  uint64_t time = 0;
  auto t1 = std::chrono::high_resolution_clock::now();
  FILE *infile = fopen(collFileName.c_str(), "r");
  if (!infile) {
    fprintf(stderr, "Error opening file of sequence (%s)\n",
            collFileName.c_str());
    exit(1);
  }
  uint64_t _sizeCollection = 0;
  fseek(infile, 0L, SEEK_END);
  _sizeCollection = ftello(infile) / sizeof(uint8_t);
  _sn = _sizeCollection;
  printf("About to read sequence of size %lu with prefix %lu\n",
         _sizeCollection, _sn);
  fclose(infile);

  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> matchingStatistics(
      numThreads);

  uint64_t numberSeqToRead = 0, capBufferWrite = 0;
  std::vector<uint64_t> parts;
  parts.reserve(numThreads + 1);
  if (bufferText / 2 < 2 * _sizeReference) {
    numberSeqToRead = 1;
    capBufferWrite =
        (bufferText - _sizeReference) / sizeof(std::pair<uint32_t, uint32_t>) +
        1;
    parts.push_back(0);
    for (uint32_t i = 1; i < numThreads; i++) {
      parts.push_back(std::round((uint64_t)_sizeReference / numThreads * i));
      // std::cout << "parts[" << i << "]: " << parts[i] << '\n';
    }
    parts.push_back(_sizeReference);
  } else {
    numberSeqToRead = bufferText / 2 / _sizeReference;
    capBufferWrite =
        (bufferText / 2) / sizeof(std::pair<uint32_t, uint32_t>) + 1;
  }
  // std::cout << "Number of sequences to read: " << numberSeqToRead << '\n';
  // std::cout << "Cap buffer write: " << capBufferWrite << '\n';

  for (uint32_t i = 0; i < numThreads; i++) {
    matchingStatistics[i].reserve(capBufferWrite / numThreads);
  }

  uint64_t charactersRead = 0;

  std::ofstream streamOutfile(outputFileName, std::ios::out | std::ios::binary);
  std::ifstream streamInfile(collFileName, std::ios::in);
  std::string line;
  while (charactersRead < _sn) {
    // open output files
    std::vector<std::ofstream> streamOutfiles(numThreads);
    for (uint32_t i = 0; i < numThreads; i++) {
      streamOutfiles[i].open(outputFileName + std::to_string(i),
                             std::ios::out | std::ios::binary);
      if (!streamOutfiles[i]) {
        fprintf(stderr, "Error opening file of sequence (%s)\n",
                (outputFileName + std::to_string(i)).c_str());
        exit(1);
      }
    }
    // uint64_t charactersToRead = std::min(bufferText, _sn - charactersRead);
    uint64_t charactersReadNow = 0;
    std::string content;
    content.reserve(numberSeqToRead * _sizeReference);
    uint32_t seqRead = 0;
    std::vector<uint64_t> posForSequenceSeparator;
    while (seqRead < numberSeqToRead &
           charactersRead + charactersReadNow < _sn) {
      std::getline(streamInfile, line);
      if (line[0] == '>') {
        if (content.size() == 0)
          continue;
        content += sequenceSeparator;
        posForSequenceSeparator.push_back(content.size());
        seqRead++;
      } else {
        content += line;
      }
      charactersReadNow += line.size() + 1;
    }
    // std::cout << "content.size(): " << content.size() << '\n';
    charactersRead += charactersReadNow;
    // std::cout << "charactersRead: " << charactersRead << '\n';

    if (numberSeqToRead != 1) {
      parts.clear();
      parts.push_back(0);
      if (numThreads < numberSeqToRead) {
        // std::cout << "Num threads < numberSeqToRead\n";
        // std::cout << "posForSequenceSeparator.size(): " <<
        // posForSequenceSeparator.size() << '\n';
        if (posForSequenceSeparator.size() > numThreads) {
          for (uint32_t x = 0; x < numThreads; x++) {
            parts.push_back(posForSequenceSeparator[std::round(
                posForSequenceSeparator.size() / numThreads * x)]);
            // std::cout << "parts[" << x << "]: " << parts[x] << '\n';
          }
          parts.push_back(
              posForSequenceSeparator[posForSequenceSeparator.size() - 1]);
        } else {
          // std::cout << "Num threads > numberSeqToRead\n";
          for (uint32_t x = 0; x < numThreads; x++) {
            parts.push_back(std::round(content.size() / numThreads * x));
            // std::cout << "parts[" << x << "]: " << parts[x] << '\n';
          }
          parts.push_back(content.size());
        }
      } else {
        // std::cout << "Num threads > numberSeqToRead\n";
        for (uint32_t x = 0; x < numThreads; x++) {
          parts.push_back(std::round(charactersReadNow / numThreads * x));
          // std::cout << "parts[" << x << "]: " << parts[x] << '\n';
        }
        parts.push_back(charactersReadNow);
      }
    } else {
      parts[numThreads] = charactersReadNow;
    }

#pragma omp parallel for
    for (uint32_t thread = 0; thread < numThreads; thread++) {
      // std::cout << omp_get_thread_num() << '\n';
      uint64_t i = parts[thread];
      int32_t leftB = 0;
      int32_t rightB = _sizeReference - 1;
      int32_t pos = _sizeReference - 1, len = 0;
      while (i < parts[thread + 1]) {
        if (content[i] == sequenceSeparator) {
          len = 0;
          pos = _sizeReference - 1;
          leftB = 0;
          rightB = _sizeReference - 1;
          matchingStatistics[omp_get_thread_num()].push_back(
              std::make_pair(pos, len));
          i++;
          continue;
        }
        computeMatchingFactor(content, i, &pos, &len, leftB, rightB);
        matchingStatistics[omp_get_thread_num()].push_back(
            std::make_pair(pos, len));
        if (matchingStatistics[omp_get_thread_num()].size() >=
            capBufferWrite / numThreads) {
          streamOutfiles[omp_get_thread_num()].write(
              reinterpret_cast<char *>(
                  &matchingStatistics[omp_get_thread_num()][0]),
              matchingStatistics[omp_get_thread_num()].size() *
                  sizeof(std::pair<uint32_t, uint32_t>));
          matchingStatistics[omp_get_thread_num()].clear();
          matchingStatistics[omp_get_thread_num()].reserve(capBufferWrite /
                                                           numThreads);
        }
        len--;

        if (leftB == rightB) {
          while ((len > _LRF[pos + 1]) & (i < parts[thread + 1])) {
            time++;
            matchingStatistics[omp_get_thread_num()].push_back(
                std::make_pair(pos + 1, len));
            i++;
            len--;
            pos++;
          }
          if (matchingStatistics[omp_get_thread_num()].size() >=
              capBufferWrite / numThreads) {
            streamOutfiles[omp_get_thread_num()].write(
                reinterpret_cast<char *>(
                    &matchingStatistics[omp_get_thread_num()][0]),
                matchingStatistics[omp_get_thread_num()].size() *
                    sizeof(std::pair<uint32_t, uint32_t>));
            matchingStatistics[omp_get_thread_num()].clear();
            matchingStatistics[omp_get_thread_num()].reserve(capBufferWrite /
                                                             numThreads);
          }
          std::pair<int, int> interval =
              adjustInterval(_ISA[pos + 1], _ISA[pos + 1], len);
          leftB = interval.first;
          rightB = interval.second;
        } else {
          std::pair<int, int> interval = contractLeft(leftB, rightB, len);
          leftB = interval.first;
          rightB = interval.second;
        }
        i++;
      }
      if (matchingStatistics[omp_get_thread_num()].size() != 0) {
        streamOutfiles[omp_get_thread_num()].write(
            reinterpret_cast<char *>(
                &matchingStatistics[omp_get_thread_num()][0]),
            matchingStatistics[omp_get_thread_num()].size() *
                sizeof(std::pair<uint32_t, uint32_t>));
        streamOutfiles[omp_get_thread_num()].close();
        matchingStatistics[omp_get_thread_num()].clear();
      }
    }
    content.clear();

    // unite the files
    for (uint32_t thread = 0; thread < numThreads; thread++) {
      std::ifstream streamInfileThread(outputFileName + std::to_string(thread),
                                       std::ios::in | std::ios::binary);
      if (!streamInfileThread) {
        fprintf(stderr, "Error opening file of sequence (%s)\n",
                (outputFileName + std::to_string(thread)).c_str());
        exit(1);
      }
      // std::cout << "file " << outputFileName + std::to_string(thread) << "
      // size: " << std::filesystem::file_size(outputFileName +
      // std::to_string(thread)) << '\n';
      if (std::filesystem::file_size(outputFileName + std::to_string(thread)) ==
          0) {
        std::cout << "file " << outputFileName + std::to_string(thread)
                  << " is empty\n";
        std::remove((outputFileName + std::to_string(thread)).c_str());
        continue;
      }
      streamOutfile << streamInfileThread.rdbuf();
      streamInfileThread.close();
    }
  }
  streamOutfile.close();
  for (uint32_t i = 0; i < numThreads; i++) {
    std::remove((outputFileName + std::to_string(i)).c_str());
  }

  std::cout << "Time to compute matching statistics: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::high_resolution_clock::now() - t1)
                   .count()
            << " ms\n";
}

void computeMatchingStatisticsParallelCMS(std::string collFileName,
                                          uint64_t bufferText,
                                          std::string outputFileName = "",
                                          uint32_t numThreads = 16) {

  std::cout << "Number of threads: " << numThreads << '\n';
  omp_set_num_threads(numThreads);
  uint64_t time = 0;
  auto t1 = std::chrono::high_resolution_clock::now();
  FILE *infile = fopen(collFileName.c_str(), "r");
  if (!infile) {
    fprintf(stderr, "Error opening file of sequence (%s)\n",
            collFileName.c_str());
    exit(1);
  }
  uint64_t _sizeCollection = 0;
  fseek(infile, 0L, SEEK_END);
  _sizeCollection = ftello(infile) / sizeof(uint8_t);
  _sn = _sizeCollection;
  printf("About to read sequence of size %lu with prefix %lu\n",
         _sizeCollection, _sn);
  fclose(infile);

  std::vector<std::vector<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>>>
      CMS(numThreads);

  uint64_t sizeElementCMS =
      sizeof(std::pair<uint32_t, std::pair<uint32_t, uint32_t>>);

  uint64_t numberSeqToRead = 0, capBufferWrite = 0;
  std::vector<uint64_t> parts;
  parts.reserve(numThreads + 1);
  if (bufferText / 2 < 2 * _sizeReference) {
    numberSeqToRead = 1;
    capBufferWrite = (bufferText - _sizeReference) / sizeElementCMS + 1;
    parts.push_back(0);
    for (uint32_t i = 1; i < numThreads; i++) {
      parts.push_back(std::round((uint64_t)_sizeReference / numThreads * i));
      // std::cout << "parts[" << i << "]: " << parts[i] << '\n';
    }
    parts.push_back(_sizeReference);
  } else {
    numberSeqToRead = bufferText / 2 / _sizeReference + 1;
    capBufferWrite = (bufferText / 2) / sizeElementCMS + 1;
  }
  // std::cout << "Number of sequences to read: " << numberSeqToRead << '\n';
  // std::cout << "Cap buffer write: " << capBufferWrite << '\n';

  for (uint32_t i = 0; i < numThreads; i++) {
    CMS[i].reserve(capBufferWrite / numThreads);
  }

  uint64_t charactersRead = 0;

  std::ofstream streamOutfile(outputFileName, std::ios::out | std::ios::binary);
  std::ifstream streamInfile(collFileName, std::ios::in);
  std::string line;
  uint64_t nCMSInBuffer = 0;
  std::vector<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>>
      CMStotalBuffer;
  CMStotalBuffer.reserve(capBufferWrite);
  while (charactersRead < _sn) {
    uint64_t charactersReadNow = 0;
    std::string content;
    content.reserve(numberSeqToRead * _sizeReference);
    uint32_t seqRead = 0;
    std::vector<uint64_t> posForSequenceSeparator;
    while (seqRead < numberSeqToRead &
           charactersRead + charactersReadNow < _sn) {
      std::getline(streamInfile, line);
      if (line[0] == '>') {
        if (content.size() == 0)
          continue;
        content += sequenceSeparator;
        posForSequenceSeparator.push_back(content.size());
        seqRead++;
      } else {
        content += line;
      }
      charactersReadNow += line.size() + 1;
    }
    // std::cout << "content.size(): " << content.size() << '\n';
    charactersRead += charactersReadNow;
    // std::cout << "charactersRead: " << charactersRead << '\n';

    if (numberSeqToRead != 1) {
      parts.clear();
      parts.push_back(0);
      if (numThreads < numberSeqToRead) {
        // std::cout << "Num threads < numberSeqToRead\n";
        // std::cout << "posForSequenceSeparator.size(): " <<
        // posForSequenceSeparator.size() << '\n';
        if (posForSequenceSeparator.size() > numThreads) {
          for (uint32_t x = 1; x < numThreads; x++) {
            parts.push_back(posForSequenceSeparator[std::round(
                posForSequenceSeparator.size() / numThreads * x)]);
            // std::cout << "parts[" << x << "]: " << parts[x] << '\n';
          }
          parts.push_back(
              posForSequenceSeparator[posForSequenceSeparator.size() - 1]);
        } else {
          // std::cout << "Num threads > numberSeqToRead\n";
          for (uint32_t x = 1; x < numThreads; x++) {
            parts.push_back(std::round(content.size() / numThreads * x));
            // std::cout << "parts[" << x << "]: " << parts[x] << '\n';
          }
          parts.push_back(content.size());
        }
      } else {
        // std::cout << "Num threads > numberSeqToRead\n";
        for (uint32_t x = 1; x < numThreads; x++) {
          parts.push_back(std::round(charactersReadNow / numThreads * x));
          // std::cout << "parts[" << x << "]: " << parts[x] << '\n';
        }
        parts.push_back(charactersReadNow);
      }
    } else {
      parts[numThreads] = charactersReadNow;
    }

#pragma omp parallel for reduction(+ : nCMSInBuffer)
    for (uint32_t thread = 0; thread < numThreads; thread++) {
      // std::cout << omp_get_thread_num() << '\n';
      uint64_t i = parts[thread];
      int32_t leftB = 0;
      int32_t rightB = _sizeReference - 1;
      int32_t pos = _sizeReference - 1, len = 0, prevPos = -2;
      while (i < parts[thread + 1]) {
        if (content[i] == sequenceSeparator) {
          len = 0;
          pos = _sizeReference - 1;
          leftB = 0;
          rightB = _sizeReference - 1;
          CMS[omp_get_thread_num()].push_back(
              std::make_pair(i, std::make_pair(pos, len)));
          i++;
          continue;
        }
        computeMatchingFactor(content, i, &pos, &len, leftB, rightB);
        if (prevPos + 1 != pos) {
          CMS[omp_get_thread_num()].push_back(
              std::make_pair(i, std::make_pair(pos, len)));
          nCMSInBuffer++;
        }
        len--;
        if (len < 0)
          len = 0;

        if (leftB == rightB) {
          while ((len > _LRF[pos + 1]) & (i < parts[thread + 1])) {
            time++;
            i++;
            len--;
            pos++;
          }
          std::pair<int, int> interval =
              adjustInterval(_ISA[pos + 1], _ISA[pos + 1], len);
          leftB = interval.first;
          rightB = interval.second;
        } else {
          std::pair<int, int> interval = contractLeft(leftB, rightB, len);
          leftB = interval.first;
          rightB = interval.second;
        }
        i++;
        prevPos = pos;
      }
    }
    content.clear();

    // unite the CMS into CMS total buffer
    for (uint32_t thread = 0; thread < numThreads; thread++) {
      CMStotalBuffer.insert(CMStotalBuffer.end(), CMS[thread].begin(),
                            CMS[thread].end());
      CMS[thread].clear();
    }

    if (CMStotalBuffer.size() >= capBufferWrite) {
      // write the CMS to file processing couples of entries at a time
      std::pair<uint32_t, std::pair<uint32_t, uint32_t>> firstHead, secondHead;
      std::vector<std::pair<uint32_t, uint32_t>> msBuffer;
      for (uint64_t i = 1; i < CMStotalBuffer.size(); i++) {
        firstHead = CMStotalBuffer[i - 1];
        secondHead = CMStotalBuffer[i];
        if (firstHead.second.second == 0) {
          msBuffer.push_back(firstHead.second);
          continue;
        }
        uint64_t offset = 0;
        for (uint64_t start = firstHead.first; start < secondHead.first;
             start++) {
          auto elementToWrite =
              std::make_pair(firstHead.second.first + offset,
                             firstHead.second.second - offset);
          offset++;
          msBuffer.push_back(elementToWrite);
        }
        streamOutfile.write(reinterpret_cast<char *>(msBuffer.data()),
                            msBuffer.size() *
                                sizeof(std::pair<uint32_t, uint32_t>));
        msBuffer.clear();
      }
      streamOutfile.write(reinterpret_cast<char *>(&secondHead.second),
                          sizeof(std::pair<uint32_t, uint32_t>));
      CMStotalBuffer.clear();
    }
  }
  std::cout << "Characters read: " << charactersRead << '\n';
  if (CMStotalBuffer.size() > 0) {
    // write the CMS to file processing couples of entries at a time
    std::pair<uint32_t, std::pair<uint32_t, uint32_t>> firstHead, secondHead;
    std::vector<std::pair<uint32_t, uint32_t>> msBuffer;
    for (uint64_t i = 1; i < CMStotalBuffer.size(); i++) {
      firstHead = CMStotalBuffer[i - 1];
      secondHead = CMStotalBuffer[i];
      if (firstHead.second.second == 0) {
        msBuffer.push_back(firstHead.second);
        continue;
      }
      uint64_t offset = 0;
      for (uint32_t start = firstHead.first; start < secondHead.first;
           start++) {
        auto elementToWrite = std::make_pair(firstHead.second.first + offset,
                                             firstHead.second.second - offset);
        offset++;
        msBuffer.push_back(elementToWrite);
      }
      streamOutfile.write(reinterpret_cast<char *>(msBuffer.data()),
                          msBuffer.size() *
                              sizeof(std::pair<uint32_t, uint32_t>));
      msBuffer.clear();
    }
    streamOutfile.write(reinterpret_cast<char *>(&secondHead.second),
                        sizeof(std::pair<uint32_t, uint32_t>));
    CMStotalBuffer.clear();
  }
  streamOutfile.close();

  std::cout << "Time to compute matching statistics: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::high_resolution_clock::now() - t1)
                   .count()
            << " ms\n";
}

int main(int argc, char **argv) {
  Args arg;
  parseArgs(argc, argv, arg);
  std::string referenceFileName = arg.reference;
  std::string collFileName = arg.collection;

  std::string outputFileName =
  arg.noOutput ? "/dev/null" : arg.outfile + ".ms";

  // std::cout << omp_get_max_threads() << '\n';
  if (arg.numThreads > 1) {
    std::cout << "Using " << arg.numThreads << " threads\n";
    // omp_set_num_threads(arg.numThreads);
    loadReferenceAndComputeDS_parallel(referenceFileName, arg.numThreads);
    computeMatchingStatisticsParallelCMS(collFileName, arg.bufferText,
                                         outputFileName, arg.numThreads);
  } else {
    loadReferenceAndComputeDS(referenceFileName);
    computeMatchingStatistics(collFileName, arg.bufferText,
                              outputFileName);
  }
  return 0;
}
