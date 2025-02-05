from __future__ import annotations

import argparse


def optics_cli() -> None:
    parser = argparse.ArgumentParser(
        prog="legend-pygeom-optics",
        description="%(prog)s command line interface",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    g4gps_parser = subparsers.add_parser(
        "g4gps", help="build evt file from remage hit file"
    )
    g4gps_parser.add_argument(
        "--type",
        choices=("arb_file", "macro"),
        help="file that contains a list of detector ids that are part of the input file. default: %(default)e",
        default="macro",
    )
    g4gps_parser.add_argument(
        "spectrum", choices=("lar_emission", "pen_emission", "fiber_emission")
    )
    g4gps_parser.add_argument("output", help="output file")

    args = parser.parse_args()

    if args.command == "g4gps":
        if args.spectrum == "lar_emission":
            from legendoptics.lar import g4gps_lar_emissions_spectrum as g4gps_spec
        elif args.spectrum == "pen_emission":
            from legendoptics.pen import g4gps_pen_emissions_spectrum as g4gps_spec
        elif args.spectrum == "fiber_emission":
            from legendoptics.fibers import g4gps_fiber_emissions_spectrum as g4gps_spec

        g4gps_spec(args.output, args.type == "macro")
