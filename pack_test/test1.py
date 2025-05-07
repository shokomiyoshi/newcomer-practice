#!/usr/bin/env python3
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from pack_test.hello import hello_test1


def main():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Test script for uv packaging.",
    )
    group = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "-i",
        "--input1",
        help="argument for input1",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-j",
        "--input2",
        help="argument for input2",
        type=str,
    )
    group.add_argument(
        "-s",
        "--switch1",
        help="Turn on switch1. Switch1 and Switch2 are mutually exclusive.",
        action="store_true",
    )
    group.add_argument(
        "-t",
        "--switch2",
        help="Turn on switch2. Switch1 and Switch2 are mutually exclusive.",
        action="store_true",
    )
    args = parser.parse_args()
    hello_test1(args.input1, args.input2, args.switch1, args.switch2)


if __name__ == "__main__":
    main()
