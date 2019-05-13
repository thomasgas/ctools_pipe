import sys
import argparse

from background_sim import simulate_background


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate some CTA Science.')

    parser.add_argument('-b', '--background',
                        dest='background',
                        action='store_const',
                        const=sum,
                        help='Do the simulation of the background. Need yaml input file.')

    parser.add_argument('config_file',
                        metavar='background.yaml',
                        type=str,
                        help='yaml configuration file for background simulation')

    # add jobs scheduler: common to all scripts

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        print(" ------------------------- ")
        print("|      Need more help?    |")
        print("| gasparethomas@gmail.com |")
        print(" ------------------------- ")
        print(" \ ")
        print("  \ ")
        print("    /¯¯\ ")
        print("    @  @")
        print("    || |/")
        print("    |\_/|")
        print("    \___/")
        print()
        sys.exit(0)

    if args.background:
        simulate_background(args.config_file)

