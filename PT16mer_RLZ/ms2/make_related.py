#!/usr/bin/env python3
"""A DNA base and a set of files closely related to it (Python 3.6+).

  make_related.py K P            generate a random base of K characters
  make_related.py FILE P         use FILE (any existing sequence) as the base
  make_related.py --base FILE P  the same (a K given as well is ignored: the
                                 base is as long as the file)

  options: [--count 10] [--out-dir DIR] [--seed S]
           [--substitution W] [--insertion W] [--deletion W] [--max-indel N]

Writes to DIR (default: related_K_P/, or related_<FILE name>_P/):

  base.txt           the generated base: K uniformly random A/C/G/T
                     characters (few repeats); not written with --base
  related_0.txt ...  --count copies of the base, each position mutated
                     independently with probability P
  list.txt           the related files' absolute paths (for --filenames)

Every base position is mutated with the same probability P, independently.
This is done exactly, without a coin flip per character: the number of
unmutated positions before the next mutated one is drawn from the
geometric distribution with parameter P, which is the same thing.

A mutated position gets one of (weights --substitution / --insertion /
--deletion, default 0.8 / 0.1 / 0.1):

  substitution   the character becomes a different one of A/C/G/T;
  insertion      1 to --max-indel random characters are inserted before it
                 (the character itself stays);
  deletion       it and up to --max-indel - 1 following characters are
                 removed (a mutation drawn inside a deletion is dropped:
                 that character is gone).

Indel lengths are geometric: each extra character half as likely. Files
are plain ACGT with no newline. Each related file has its own seed derived
from --seed, so a run is reproducible.
"""

import argparse
import math
import os
import random
import sys

ACGT = b"ACGT"
TO_ACGT = bytes(ACGT[b % 4] for b in range(256))  # byte value -> A/C/G/T


def random_dna(rng, length):
    if length == 0:
        return b""
    return rng.getrandbits(8 * length).to_bytes(length, "little").translate(TO_ACGT)


def mutated_positions(rng, n, p):
    """Every position in 0..n-1 independently with probability p, in order."""
    if p <= 0:
        return
    if p >= 1:
        yield from range(n)
        return

    log_q = math.log1p(-p)
    position = -1
    while True:
        # Unmutated positions before the next mutated one: Geometric(p).
        u = 1.0 - rng.random()  # in (0, 1]
        position += 1 + int(math.log(u) / log_q)
        if position >= n:
            return
        yield position


def indel_length(rng, max_indel):
    length = 1
    while length < max_indel and rng.random() < 0.5:
        length += 1
    return length


def mutate(base, rng, p, weights, max_indel):
    """A copy of `base` with every position mutated with probability p."""
    n = len(base)
    counts = {"substitution": 0, "insertion": 0, "deletion": 0,
              "inserted": 0, "deleted": 0, "dropped": 0}
    kinds = ["substitution", "insertion", "deletion"]

    pieces = []
    at = 0  # next base position not yet copied

    for position in mutated_positions(rng, n, p):
        if position < at:  # inside an earlier deletion
            counts["dropped"] += 1
            continue

        pieces.append(base[at:position])
        kind = rng.choices(kinds, weights=weights)[0]
        counts[kind] += 1

        if kind == "substitution":
            old = base[position]
            pieces.append(bytes([rng.choice([c for c in ACGT if c != old])]))
            at = position + 1
        elif kind == "insertion":
            length = indel_length(rng, max_indel)
            pieces.append(random_dna(rng, length))
            counts["inserted"] += length
            at = position  # the character itself is kept
        else:
            length = min(indel_length(rng, max_indel), n - position)
            counts["deleted"] += length
            at = position + length

    pieces.append(base[at:])
    return b"".join(pieces), counts


def write(path, data):
    with open(path, "wb") as f:
        f.write(data)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("values", nargs="*", metavar="K|FILE P",
                        help="K, the length of a generated base, or FILE, an "
                             "existing base (or give it with --base); then "
                             "P, the probability that a position is mutated "
                             "(e.g. 0.0001 = 0.01%%)")
    parser.add_argument("--base", default=None,
                        help="use this file as the base instead of "
                             "generating one")
    parser.add_argument("--count", type=int, default=10,
                        help="number of related files (default 10)")
    parser.add_argument("--out-dir", default=None,
                        help="output directory (default related_K_P/)")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--substitution", type=float, default=0.8)
    parser.add_argument("--insertion", type=float, default=0.1)
    parser.add_argument("--deletion", type=float, default=0.1)
    parser.add_argument("--max-indel", type=int, default=10)
    # Positionals may also come after the options (argparse on Python 3.6
    # does not mix them, so the leftovers are collected here).
    args, extra = parser.parse_known_args()
    unknown = [x for x in extra if x.startswith("-")]
    if unknown:
        parser.error("unrecognized arguments: " + " ".join(unknown))
    values = args.values + extra

    forms = ("give K P (generate a base), FILE P or --base FILE P "
             "(use an existing base)")

    # A value naming an existing file is the base.
    files = [v for v in values if os.path.isfile(v)]
    if files:
        if args.base is not None or len(files) > 1:
            parser.error("more than one base file given; " + forms)
        args.base = files[0]
        values = [v for v in values if v != files[0]]

    try:
        if args.base is not None:
            # K P with a base file: K is not needed (the base is the file),
            # so it is ignored rather than refused.
            if len(values) == 2:
                ignored_k = int(values[0])
                values = values[1:]
            else:
                ignored_k = None
            if len(values) != 1:
                parser.error(forms)
            k, p = None, float(values[0])
        else:
            if len(values) != 2:
                parser.error(forms)
            k, p = int(values[0]), float(values[1])
    except ValueError:
        parser.error("K must be an integer and P a number "
                     "(or a base file that does not exist?); " + forms)
    args.p = p

    if k is not None and k < 0:
        parser.error("K must not be negative")
    if not 0.0 <= p <= 1.0:
        parser.error("P must be between 0 and 1")

    if args.base is not None:
        name = os.path.splitext(os.path.basename(args.base))[0]
        out_dir = args.out_dir or "related_{}_{}".format(name, p)
        os.makedirs(out_dir, exist_ok=True)

        with open(args.base, "rb") as f:
            base = f.read()
        print("base: {} ({} characters, read)".format(args.base, len(base)),
              file=sys.stderr)
        if ignored_k is not None and ignored_k != len(base):
            print("note: K = {} ignored; the base is the file ({} characters)"
                  .format(ignored_k, len(base)), file=sys.stderr)
    else:
        out_dir = args.out_dir or "related_{}_{}".format(k, p)
        os.makedirs(out_dir, exist_ok=True)

        base = random_dna(random.Random(args.seed), k)
        base_path = os.path.join(out_dir, "base.txt")
        write(base_path, base)
        print("base: {} ({} characters, seed {})".format(base_path, len(base),
                                                         args.seed),
              file=sys.stderr)

    weights = [args.substitution, args.insertion, args.deletion]
    paths = []

    for k in range(args.count):
        rng = random.Random(args.seed * 1000003 + k + 1)
        data, c = mutate(base, rng, args.p, weights, args.max_indel)

        path = os.path.join(out_dir, "related_{}.txt".format(k))
        write(path, data)
        paths.append(os.path.abspath(path))

        mutations = c["substitution"] + c["insertion"] + c["deletion"]
        print("related: {} ({} characters, {} mutations = {:.4f}% of positions: "
              "{} subst, {} ins (+{}), {} del (-{}); {} fell in a deletion)".format(
                  path, len(data), mutations + c["dropped"],
                  100.0 * (mutations + c["dropped"]) / max(1, len(base)),
                  c["substitution"], c["insertion"], c["inserted"],
                  c["deletion"], c["deleted"], c["dropped"]), file=sys.stderr)

    list_path = os.path.join(out_dir, "list.txt")
    write(list_path, ("\n".join(paths) + "\n").encode())
    print("list: {} ({} files)".format(list_path, len(paths)), file=sys.stderr)


if __name__ == "__main__":
    main()
