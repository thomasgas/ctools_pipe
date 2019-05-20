import sys
import argparse
import yaml
import subprocess
from utils import create_path
from environs import Env

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="python ctools_pipe.py",
        description='Simulate some CTA Science.',
    )

    parser.add_argument("-j",
                        "--jobs",
                        metavar='jobs.yaml',
                        type=str,
                        required=True,
                        help="MANDATORY: input a jobs.yaml configuration file. ")

    parser.add_argument("-b",
                        "--background",
                        metavar='background*.yaml',
                        type=str,
                        required=False,
                        help="Needs a background*.yaml configuration file.")

    parser.add_argument("-m",
                        "--models",
                        type=str,
                        metavar='model_input.yaml',
                        required=False,
                        help="Needs a model_input.yaml configuration file.")

    try:
        args = parser.parse_args()
    except:
        # parser.print_help()
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
        infile = args.background
        config_in = yaml.safe_load(open(infile))
        realizations = config_in['sim']['realizations']

        in_jobs = args.jobs
        jobs_exe = yaml.safe_load(open(in_jobs))
        execution = jobs_exe['exe']

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

            # create string for jobs submission
            exec_string = f"{execution['mode']} "

            if details['queue']['name'] != "N/A":
                exec_string += f"-q {details['queue']['name']} "
            if details['queue']['flags'] != "N/A":
                exec_string += f"{details['queue']['flags']} "
            if details['mail'] != "N/A":
                exec_string += f"-u {details['mail']} "
            if execution['others'] != "N/A":
                exec_string += f"{execution['others']} "
            
            print(exec_string)
            for counter in range(realizations):
                p = subprocess.Popen([*exec_string.split(" "),
                                      "python", "background_sim.py", infile, str(counter + 1)],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE
                                     )
                # print stderr ad stdout for first job
                if counter == 0 and execution['debug'] == "yes":
                    (result, error) = p.communicate()
                    print(result, error)

        elif execution['mode'] == "qsub":
            details = execution['details']
            # create scripts to lunch on a folder
            folder_launch = create_path(f"{execution['path']}/tmp/scripts")

            conda_path = create_path(execution['conda']['conda_path'])
            env = Env()
            for counter in range(realizations):
                script_name = f"launch_{str(counter).zfill(3)}.sh"
                out_file_name = folder_launch + "/" + script_name
                file_out = open(out_file_name, "w")

                # need to export also the env variables, if they are used
                if execution['path'].startswith("$"):
                    env_folder_name = execution['path'][1:].split('/', 1)[0]
                    evaluate_folder = env(env_folder_name)
                    file_out.write(f'export {env_folder_name}="{evaluate_folder}"\n')

                ctools_pipe_path=create_path(execution['software_path'])
                env_name = execution['conda']['env_name']
                caldb_path = create_path(execution['caldb'])

                file_out.write(f'export CALDB="{caldb_path}"\n')
                file_out.write(f'export PATH="{conda_path}/bin:$PATH"\n')
                file_out.write(f'export PATH="{conda_path}/lib:$PATH"\n')
                file_out.write(f'export PYTHON_EGG_CACHE="/lapp_data/cta/gasparetto/python_cache"\n')
                file_out.write(f'source activate {env_name}\n')
                file_out.write(f'python {ctools_pipe_path}/background_sim.py {ctools_pipe_path}/{infile} {str(counter + 1)} \n')
                file_out.write('source deactivate\n')
                file_out.close()

                exec_string = f"{execution['mode']} -V -j oe "
                if details['output'] != "N/A":
                    exec_string += f"-o {details['output']} "
                if details['queue'] != "N/A":
                    exec_string += f"-q {details['queue']} "
                if details['mail'] != "N/A":
                    exec_string += f"-M {details['mail']} "
                if details['flags'] != "N/A":
                    exec_string += f"-m {details['flags']} "
                if execution['others'] != "N/A":
                    exec_string += f"{execution['others']} "

                exec_string += f"{out_file_name}"
                print(exec_string)
                p = subprocess.Popen(exec_string.split(" "),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE
                                     )
                # print stderr ad stdout for first job
                if counter == 0 and execution['debug'] == "yes":
                    (result, error) = p.communicate()
                    print(result, error)

    if args.models:
        infile = args.models
        config_in = yaml.safe_load(open(infile))

        in_jobs = args.jobs
        jobs_exe = yaml.safe_load(open(in_jobs))
        execution = jobs_exe['exe']
        if execution['mode'] == "local":
            p = subprocess.Popen(
                ['python', 'model_creation.py', infile]
            )
            # if everything goes well, the output is None
            # check this just for the first job
            (result, error) = p.communicate()
            print(result, error)
