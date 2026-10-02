#!/usr/bin/env python3

"""
Estimate PT16 hit rates from greedy RLZ phrase lengths.

For a candidate k, a PT16-style lookup at a phrase start is a hit when
the longest RLZ match length L satisfies L >= k.

Reports predicted hit rates for one or more k values.
"""

import argparse
import sys
from array import array


def load_reference(path):
    with open(path, "rb") as fh:
        return fh.read()


def load_suffix_array(path):
    sa = array("I")

    if sa.itemsize != 4:
        sys.exit(
            "error: array('I') is %d bytes; expected 4 bytes."
            % sa.itemsize
        )

    with open(path, "rb") as fh:
        raw = fh.read()

    if len(raw) % 4 != 0:
        sys.exit(
            "error: suffix-array file size is not a multiple of 4 bytes."
        )

    sa.frombytes(raw)
    return sa


def longest_match_len(ref, n_ref, sa, m, pattern, start, plen):
    lo = 0
    hi = m
    best = 0
    max_l = plen - start

    while lo < hi:
        mid = (lo + hi) >> 1

        rp = sa[mid]
        pp = start
        length = 0

        while (
            length < max_l
            and rp < n_ref
            and ref[rp] == pattern[pp]
        ):
            length += 1
            rp += 1
            pp += 1

        if length > best:
            best = length

        if length == max_l:
            break

        if rp >= n_ref:
            lo = mid + 1
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

    lookups = {k: 0 for k in ks}
    hits = {k: 0 for k in ks}

    i = 0

    while i < plen:
        remaining = plen - i

        length = longest_match_len(
            ref,
            n_ref,
            sa,
            m,
            pattern,
            i,
            plen,
        )

        if length <= 1:
            length = 1

        for k in ks:
            if remaining >= k:
                lookups[k] += 1

                if length >= k:
                    hits[k] += 1

        total_phrases += 1
        total_covered += length
        i += length

    return {
        "path": input_path,
        "bytes": plen,
        "phrases": total_phrases,
        "avg_len": (
            total_covered / total_phrases
            if total_phrases
            else 0.0
        ),
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

    for result in results:
        print("file:", result["path"])
        print("  input_bytes   :", result["bytes"])
        print("  phrases       :", result["phrases"])
        print("  avg phrase len: %.2f" % result["avg_len"])

        for k in ks:
            lookups = result["lookups"][k]
            hits = result["hits"][k]

            hit_rate = hits / lookups if lookups else 0.0

            print(
                "  k=%-3d  hit_rate=%7.4f%%  hits=%d  lookups=%d"
                % (
                    k,
                    100.0 * hit_rate,
                    hits,
                    lookups,
                )
            )

        print()

    if len(results) > 1:
        print("combined:")

        total_bytes = sum(
            result["bytes"]
            for result in results
        )

        total_phrases = sum(
            result["phrases"]
            for result in results
        )

        total_covered = sum(
            result["avg_len"] * result["phrases"]
            for result in results
        )

        avg_len = (
            total_covered / total_phrases
            if total_phrases
            else 0.0
        )

        print("  input_bytes   :", total_bytes)
        print("  phrases       :", total_phrases)
        print("  avg phrase len: %.2f" % avg_len)

        for k in ks:
            lookups = sum(
                result["lookups"][k]
                for result in results
            )

            hits = sum(
                result["hits"][k]
                for result in results
            )

            hit_rate = hits / lookups if lookups else 0.0

            print(
                "  k=%-3d  hit_rate=%7.4f%%  hits=%d  lookups=%d"
                % (
                    k,
                    100.0 * hit_rate,
                    hits,
                    lookups,
                )
            )

        print()


def main():
    parser = argparse.ArgumentParser(
        description="Estimate PT16 hit rates from RLZ phrase lengths."
    )

    parser.add_argument(
        "--reference",
        required=True,
    )

    parser.add_argument(
        "--suffix-array",
        required=True,
    )

    parser.add_argument(
        "--input",
        action="append",
        default=[],
    )

    parser.add_argument(
        "--filenames",
    )

    parser.add_argument(
        "--limit-files",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--max-bytes",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--k",
        type=int,
        nargs="+",
        default=[4, 8, 12, 16],
    )

    args = parser.parse_args()

    ks = sorted(set(args.k))

    inputs = list(args.input)

    if args.filenames:
        inputs.extend(
            read_filenames(args.filenames)[
                :args.limit_files
            ]
        )

    if not inputs:
        parser.error(
            "provide --input or --filenames"
        )

    print("Loading reference...", file=sys.stderr)
    ref = load_reference(args.reference)

    print("Loading suffix array...", file=sys.stderr)
    sa = load_suffix_array(args.suffix_array)

    if len(sa) != len(ref):
        print(
            "warning: |SA|=%d, |reference|=%d"
            % (len(sa), len(ref)),
            file=sys.stderr,
        )

    results = []

    for path in inputs:
        print(
            "Parsing:",
            path,
            file=sys.stderr,
        )

        results.append(
            parse_file(
                ref,
                sa,
                path,
                ks,
                args.max_bytes,
            )
        )

    print_report(results, ks)


if __name__ == "__main__":
    main()