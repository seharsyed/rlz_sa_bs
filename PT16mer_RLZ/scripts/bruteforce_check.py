#!/usr/bin/env python3

import argparse
import struct


def load_sequence(path):
    with open(path, "rb") as f:
        return f.read().replace(b"\n", b"").replace(b"\r", b"")


def load_factors(path):
    factors = []

    with open(path, "rb") as f:
        count = struct.unpack("<Q", f.read(8))[0]

        for _ in range(count):
            input_pos = struct.unpack("<Q", f.read(8))[0]
            ref_pos = struct.unpack("<Q", f.read(8))[0]
            length = struct.unpack("<Q", f.read(8))[0]

            factors.append((input_pos, ref_pos, length))

    return factors


def brute_force_longest_match(reference, input_seq, input_pos):
    best_ref_pos = 0
    best_length = 0

    for ref_pos in range(len(reference)):
        length = 0

        while (
            input_pos + length < len(input_seq)
            and ref_pos + length < len(reference)
            and input_seq[input_pos + length] == reference[ref_pos + length]
        ):
            length += 1

        if length > best_length:
            best_length = length
            best_ref_pos = ref_pos

    return best_ref_pos, best_length


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--reference", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--factors", required=True)
    parser.add_argument("--max-factors", type=int, default=50)

    args = parser.parse_args()

    reference = load_sequence(args.reference)
    input_seq = load_sequence(args.input)
    factors = load_factors(args.factors)

    number_to_check = min(args.max_factors, len(factors))

    print("reference_length =", len(reference))
    print("input_length     =", len(input_seq))
    print("factor_count     =", len(factors))
    print("checking         =", number_to_check)

    for factor_index in range(number_to_check):
        input_pos, pt16_ref_pos, pt16_length = factors[factor_index]

        # Literals are not reference copies.
        if pt16_length == 1:
            continue

        oracle_ref_pos, oracle_length = brute_force_longest_match(
            reference,
            input_seq,
            input_pos
        )

        # PT16 factor must actually match.
        pt16_match = reference[
            pt16_ref_pos:pt16_ref_pos + pt16_length
        ]

        input_match = input_seq[
            input_pos:input_pos + pt16_length
        ]

        if pt16_match != input_match:
            print()
            print("FAIL: invalid PT16 reference occurrence")
            print("factor_index =", factor_index)
            print("input_pos    =", input_pos)
            print("pt16_ref_pos =", pt16_ref_pos)
            print("pt16_length  =", pt16_length)
            return

        # PT16 must have the true longest-match length.
        if pt16_length != oracle_length:
            print()
            print("FAIL: factor is not maximal")
            print("factor_index      =", factor_index)
            print("input_pos         =", input_pos)
            print("pt16_ref_pos      =", pt16_ref_pos)
            print("pt16_length       =", pt16_length)
            print("oracle_ref_pos    =", oracle_ref_pos)
            print("oracle_length     =", oracle_length)
            return

    print()
    print("PASS")
    print("All checked PT16 factors have the true longest-match length.")
    print("All checked PT16 reference positions are valid.")


if __name__ == "__main__":
    main()