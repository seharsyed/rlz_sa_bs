#!/usr/bin/env python3
"""Synthetic test data (Python 3.6+): one random DNA sequence and files closely related to it.

The base is k uniformly random A/C/G/T characters, so it has almost no
internal long repeats (low LCP within the file). Each related file is a copy of
the base with a few random edits -- substitutions, insertions and
deletions -- so the files are highly similar to the base and to each other.

  # the base only
  gen_related.py base --length 50000000 --out base.txt

  # related files from an existing base
  gen_related.py related --base base.txt --count 10 --rate 0.0001 \\
      --out-dir related/

  # both at once (base + related + a list file for --filenames)
  gen_related.py all --length 50000000 --count 10 --rate 0.0001 --out-dir data/

--rate is the chance that an edit starts at a given position (0.0001 =
0.01%); an edit is a substitution, an insertion or a deletion (chosen by
--substitution/--insertion/--deletion weights), and an indel is 1 to
--max-indel characters long (shorter ones more likely). So a file differs
from the base in about rate * length places, and its length changes a
little. Files are written as plain ACGT with no newline. Each related file
gets a different seed derived from --seed, so a run is reproducible.
"""

import argparse
import math
import os
import random
import sys

ACGT = b"ACGT"

# Maps every byte value to one of A, C, G, T (value mod 4).
TO_ACGT = bytes(ACGT[b % 4] for b in range(256))


def random_dna(rng: random.Random, length: int) -> bytes:
    # getrandbits rather than randbytes (Python 3.9+), for older Pythons.
    if length == 0:
        return b""
    raw = rng.getrandbits(8 * length).to_bytes(length, "little")
    return raw.translate(TO_ACGT)


def event_count(rng: random.Random, n: int, rate: float) -> int:
    """How many edits: Binomial(n, rate), drawn cheaply."""
    if n == 0 or rate <= 0:
        return 0
    if hasattr(rng, "binomialvariate"):  # Python 3.12+
        return rng.binomialvariate(n, rate)
    mean = n * rate
    if mean < 30:  # Poisson (Knuth), close to the binomial for small rates
        limit, k, p = math.exp(-mean), 0, 1.0
        while True:
            p *= rng.random()
            if p <= limit:
                return k
            k += 1
    sd = math.sqrt(mean * (1 - rate))
    return max(0, round(rng.gauss(mean, sd)))


def indel_length(rng: random.Random, max_indel: int) -> int:
    """1..max_indel, geometric-ish: each extra character half as likely."""
    length = 1
    while length < max_indel and rng.random() < 0.5:
        length += 1
    return length


def mutate(base: bytes, rng: random.Random, rate: float, weights, max_indel: int):
    """A copy of `base` with random edits; returns (bytes, counts)."""
    n = len(base)
    m = event_count(rng, n, rate)
    positions = sorted(rng.sample(range(n), m)) if m else []

    counts = {"substitution": 0, "insertion": 0, "deletion": 0,
              "inserted": 0, "deleted": 0}
    kinds = ["substitution", "insertion", "deletion"]

    pieces = []
    at = 0  # next base position not yet copied

    for position in positions:
        if position < at:  # inside an earlier deletion
            continue

        pieces.append(base[at:position])
        kind = rng.choices(kinds, weights=weights)[0]
        counts[kind] += 1

        if kind == "substitution":
            old = base[position]
            new = rng.choice([c for c in ACGT if c != old])
            pieces.append(bytes([new]))
            at = position + 1
        elif kind == "insertion":
            length = indel_length(rng, max_indel)
            pieces.append(random_dna(rng, length))
            counts["inserted"] += length
            at = position  # the base character itself is kept
        else:  # deletion
            length = min(indel_length(rng, max_indel), n - position)
            counts["deleted"] += length
            at = position + length

    pieces.append(base[at:])
    return b"".join(pieces), counts


def write(path: str, data: bytes) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "wb") as f:
        f.write(data)


def make_base(length: int, seed: int, out: str) -> bytes:
    data = random_dna(random.Random(seed), length)
    write(out, data)
    print(f"base: {out} ({len(data)} characters, seed {seed})", file=sys.stderr)
    return data


def make_related(base: bytes, args, out_dir: str) -> list:
    weights = [args.substitution, args.insertion, args.deletion]
    paths = []

    for k in range(args.count):
        rng = random.Random(args.seed * 1_000_003 + k + 1)
        data, counts = mutate(base, rng, args.rate, weights, args.max_indel)

        path = os.path.join(out_dir, f"{args.prefix}{k}.txt")
        write(path, data)
        paths.append(os.path.abspath(path))

        edits = counts["substitution"] + counts["insertion"] + counts["deletion"]
        print(f"related: {path} ({len(data)} characters, {edits} edits = "
              f"{100.0 * edits / max(1, len(base)):.4f}%: "
              f"{counts['substitution']} subst, {counts['insertion']} ins "
              f"(+{counts['inserted']}), {counts['deletion']} del "
              f"(-{counts['deleted']}))", file=sys.stderr)

    return paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # (required= for subparsers needs Python 3.7+; checked below instead)
    sub = parser.add_subparsers(dest="command")

    def related_options(p):
        p.add_argument("--count", type=int, required=True,
                       help="number of related files")
        p.add_argument("--rate", type=float, default=0.0001,
                       help="chance an edit starts at a position (default "
                            "0.0001 = 0.01%%)")
        p.add_argument("--substitution", type=float, default=0.8,
                       help="weight of substitutions (default 0.8)")
        p.add_argument("--insertion", type=float, default=0.1,
                       help="weight of insertions (default 0.1)")
        p.add_argument("--deletion", type=float, default=0.1,
                       help="weight of deletions (default 0.1)")
        p.add_argument("--max-indel", type=int, default=10,
                       help="longest insertion/deletion (default 10)")
        p.add_argument("--prefix", default="related_",
                       help="file name prefix (default related_)")
        p.add_argument("--out-dir", required=True)

    p_base = sub.add_parser("base", help="write the random base sequence")
    p_base.add_argument("--length", type=int, required=True)
    p_base.add_argument("--out", required=True)
    p_base.add_argument("--seed", type=int, default=1)

    p_related = sub.add_parser("related", help="related files from a base")
    p_related.add_argument("--base", required=True)
    p_related.add_argument("--seed", type=int, default=1)
    related_options(p_related)

    p_all = sub.add_parser("all", help="base + related files + list file")
    p_all.add_argument("--length", type=int, required=True)
    p_all.add_argument("--seed", type=int, default=1)
    related_options(p_all)

    args = parser.parse_args()

    if args.command is None:
        parser.error("choose a command: base, related or all")

    if args.command == "base":
        make_base(args.length, args.seed, args.out)
        return

    if args.command == "related":
        with open(args.base, "rb") as f:
            base = f.read()
        make_related(base, args, args.out_dir)
        return

    # all
    base_path = os.path.join(args.out_dir, "base.txt")
    base = make_base(args.length, args.seed, base_path)
    paths = make_related(base, args, args.out_dir)

    list_path = os.path.join(args.out_dir, "list.txt")
    write(list_path, ("\n".join(paths) + "\n").encode())
    print(f"list: {list_path} ({len(paths)} files)", file=sys.stderr)


if __name__ == "__main__":
    main()
