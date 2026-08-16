#!/usr/bin/env python3
"""
pt16_hitrate_check.py

Audits whether the PT16 miss rate is a data property or a code bug, and tells
you which k would give the best hit rate -- WITHOUT touching your C++.

Core idea
---------
For greedy RLZ, the k-mer at a phrase start is present in the reference iff that
phrase's longest-match length L satisfies L >= k. So:

        predicted PT16 hit rate for parameter k  =  P(L >= k)

This script reproduces the exact greedy longest-match parse using your suffix
array (the same computation as computeLZFactorAt in parser.hpp, incl. the
len<=1 -> literal rule), collects the phrase-length distribution, and reports
P(L >= k) for every candidate k at once.

How to read the result
-----------------------
* P(L >= 16) should land near the 8.4% your pt16 CSV reports.
    - If it matches  -> the low hit rate is real (E. coli divergence), NOT a
      lookup bug. The C++ path is exonerated.
    - If C++ reports FAR less than this -> the C++ lookup is dropping hits it
      should be finding (bucket index / L_ scan / packing). Chase that.
* Compare P(L >= 4), P(L >= 8), P(L >= 12), P(L >= 16) to pick k. Bigger hit
  rate is good, but remember smaller k resolves fewer characters per hit.

Inputs must be the SAME raw byte files your C++ reads (read_file<unsigned char>
for reference/inputs, read_file<unsigned int> for the SA -- native-endian
uint32). The "cleaned" .fna inputs are pure ACGT byte streams (no headers).

Usage
-----
  python3 pt16_hitrate_check.py \
      --reference   /path/to/reference.raw \
      --suffix-array /path/to/reference.sa \
      --input       /path/to/one_input.fna

  # or run over the first few entries of your C++ filenames list:
  python3 pt16_hitrate_check.py \
      --reference R --suffix-array SA --filenames files.txt --limit-files 3

  # start small if a full file is slow:
  python3 pt16_hitrate_check.py ... --input X.fna --max-bytes 1000000

No third-party dependencies (pure standard library).
"""

import argparse
import sys
from array import array


def load_reference(path):
    with open(path, "rb") as fh:
        ref = fh.read()          # bytes: O(1) integer indexing
    return ref


def load_suffix_array(path):
    sa = array("I")              # unsigned int == uint32, native endian (matches C++)
    if sa.itemsize != 4:
        sys.exit("error: array('I') is %d bytes here, expected 4 (uint32). "
                 "Adjust the SA loader for your platform." % sa.itemsize)
    with open(path, "rb") as fh:
        raw = fh.read()
    if len(raw) % 4 != 0:
        sys.exit("error: suffix-array file size is not a multiple of 4 bytes.")
    sa.frombytes(raw)
    return sa


def longest_match_len(ref, n_ref, sa, m, pattern, start, plen):
    """Length of the longest prefix of pattern[start:plen] that occurs in ref,
    via binary search over the suffix array (max LCP along the descent)."""
    lo = 0
    hi = m                      # exclusive
    best = 0
    max_l = plen - start
    while lo < hi:
        mid = (lo + hi) >> 1
        rp = sa[mid]
        pp = start
        l = 0
        # extend the common prefix between ref[sa[mid]:] and pattern[start:]
        while l < max_l and rp < n_ref and ref[rp] == pattern[pp]:
            l += 1
            rp += 1
            pp += 1
        if l > best:
            best = l
        if l == max_l:
            break               # whole remaining pattern matched
        if rp >= n_ref:
            lo = mid + 1         # suffix ran out -> it sorts before the pattern
        elif ref[rp] < pattern[pp]:
            lo = mid + 1
        else:
            hi = mid
    return best


def parse_file(ref, sa, input_path, ks, max_bytes):
    with open(input_path, "rb") as fh:
        pattern = fh.read()
    if max_bytes and len(pattern) > max_bytes:
        pattern = pattern[:max_bytes]

    n_ref = len(ref)
    m = len(sa)
    plen = len(pattern)

    total_phrases = 0
    total_covered = 0
    # per-k tallies: lookups_k = phrases with remaining >= k; hits_k = phrases with L >= k
    lookups = {k: 0 for k in ks}
    hits = {k: 0 for k in ks}

    i = 0
    while i < plen:
        remaining = plen - i
        L = longest_match_len(ref, n_ref, sa, m, pattern, i, plen)
        if L <= 1:
            L = 1               # literal rule (matches lzFactorize)
        for k in ks:
            if remaining >= k:
                lookups[k] += 1
                if L >= k:
                    hits[k] += 1
        total_phrases += 1
        total_covered += L
        i += L

    return {
        "path": input_path,
        "bytes": plen,
        "phrases": total_phrases,
        "avg_len": (total_covered / total_phrases) if total_phrases else 0.0,
        "lookups": lookups,
        "hits": hits,
    }


def read_filenames(path):
    files = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                files.append(line)
    return files


def print_report(results, ks):
    print()
    for r in results:
        print("file: %s" % r["path"])
        print("  input_bytes   : %d" % r["bytes"])
        print("  phrases       : %d" % r["phrases"])
        print("  avg phrase len: %.2f" % r["avg_len"])
        print("  predicted PT16 hit rate  = P(L >= k) = hits / lookups")
        for k in ks:
            lk = r["lookups"][k]
            hk = r["hits"][k]
            rate = (hk / lk) if lk else 0.0
            print("    k=%-3d  hit_rate = %7.4f%%   (%d / %d)"
                  % (k, 100.0 * rate, hk, lk))
        print()

    if len(results) > 1:
        print("=== combined over %d files ===" % len(results))
        tot_bytes = sum(r["bytes"] for r in results)
        tot_phr = sum(r["phrases"] for r in results)
        tot_cov = sum(r["avg_len"] * r["phrases"] for r in results)
        print("  input_bytes   : %d" % tot_bytes)
        print("  phrases       : %d" % tot_phr)
        print("  avg phrase len: %.2f" % (tot_cov / tot_phr if tot_phr else 0.0))
        for k in ks:
            lk = sum(r["lookups"][k] for r in results)
            hk = sum(r["hits"][k] for r in results)
            rate = (hk / lk) if lk else 0.0
            print("    k=%-3d  hit_rate = %7.4f%%   (%d / %d)"
                  % (k, 100.0 * rate, hk, lk))
        print()

    print("Interpretation:")
    print("  * Compare k=16 above to the hit rate in your pt16 CSV (~8.4%).")
    print("    Matches  -> low hit rate is the data, not a lookup bug.")
    print("    C++ much lower -> real bug in the C++ lookup path.")
    print("  * Use the k=4/8/12/16 rates to choose k (higher rate = more hits,")
    print("    but smaller k resolves fewer characters and narrows less).")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, help="reference raw byte file")
    ap.add_argument("--suffix-array", required=True, help="suffix array (uint32) file")
    ap.add_argument("--input", action="append", default=[],
                    help="an input file to parse (repeatable)")
    ap.add_argument("--filenames", help="text file listing input paths (like the C++ --filenames)")
    ap.add_argument("--limit-files", type=int, default=1,
                    help="when using --filenames, process only the first N (default 1)")
    ap.add_argument("--max-bytes", type=int, default=0,
                    help="cap each input to the first N bytes (0 = whole file)")
    ap.add_argument("--k", type=int, nargs="+", default=[4, 8, 12, 16],
                    help="candidate k values (default: 4 8 12 16)")
    args = ap.parse_args()

    ks = sorted(set(args.k))

    inputs = list(args.input)
    if args.filenames:
        inputs.extend(read_filenames(args.filenames)[:args.limit_files])
    if not inputs:
        ap.error("provide at least one --input or a --filenames list")

    sys.stderr.write("loading reference ...\n")
    ref = load_reference(args.reference)
    sys.stderr.write("loading suffix array ...\n")
    sa = load_suffix_array(args.suffix_array)
    if len(sa) != len(ref):
        sys.stderr.write("warning: |SA|=%d differs from |ref|=%d "
                         "(sentinel/off-by-one?)\n" % (len(sa), len(ref)))

    results = []
    for path in inputs:
        sys.stderr.write("parsing %s ...\n" % path)
        results.append(parse_file(ref, sa, path, ks, args.max_bytes))

    print_report(results, ks)


if __name__ == "__main__":
    main()