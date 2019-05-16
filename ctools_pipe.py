import sys
import argparse
import yaml
import subprocess

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
        infile = args.config_file
        config_in = yaml.safe_load(open(infile))
        realizations = config_in['sim']['realizations']
        execution = config_in['exe']

        # launches job using local python
        if execution['mode'] == "local":
            for counter in range(realizations):
                print(f"process {counter} started")
                p = subprocess.Popen(
                    ['python', 'background_sim.py', infile, str(counter + 1)],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                # if everything goes well, the output is None
                # check this just for the first job
                if counter == 0:
                    (result, error) = p.communicate()
                    print(result, error)
        elif execution['mode'] == "bsub":
            details = execution['details']
            exec_string = "bsub "
            if details['queue']['name']  != "N/A":
                exec_string+=f"-q {details['queue']['name']} "
            if details['queue']['flags'] != "N/A":
                exec_string+=f"{details['queue']['flags']} "
            if details['mail'] != "N/A":
                exec_string+=f"-u {details['mail']} "
            if execution['others'] != "N/A":
                exec_string+=execution['others']
            
            print(exec_string)
            for counter in range(realizations):
                p = subprocess.Popen([*exec_string.split(" "),
                                      "python", "background_sim.py", infile, str(counter + 1)],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE
                                     )
                
                if counter == 0 and execution['debug'] == "yes":
                    (result, error) = p.communicate()
                    print(result, error)


