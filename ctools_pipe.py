import sys
import argparse
import yaml
import subprocess
from utils import create_path
from environs import Env
import glob
from irf_handler import IRFPicker
import pandas as pd


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

    parser.add_argument("-s",
                        "--simulation",
                        type=str,
                        metavar='simulation.yaml',
                        required=False,
                        help="Needs a simulation.yaml configuration file.")


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

    # input for jobs submission
    in_jobs = args.jobs
    jobs_exe = yaml.safe_load(open(in_jobs))
    execution = jobs_exe['exe']

    if args.background:
        infile = args.background
        config_in = yaml.safe_load(open(infile))
        realizations = config_in['sim']['realizations']

        # launches job using local python
        if execution['mode'] == "local":
            for counter in range(realizations):
                print(f"process {counter} started")
                p = subprocess.Popen(
                    ['python', 'background_sim.py', infile, in_jobs, str(counter + 1)],
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

            env_path = create_path(execution['env_path'])

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
                if counter == 0 and execution['debug'] == "yes":
                    exec_string+="-e error.log -o output.log "

                p = subprocess.Popen([*exec_string.split(" "),
                                      f"{env_path}/bin/python", "background_sim.py", infile, in_jobs, str(counter + 1)],
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
            folder_launch = create_path(f"{execution['path']}/tmp/back/scripts")

            conda_path = create_path(execution['conda']['conda_path'])
            env = Env()
            list_script_bash = []

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
                python_cache = create_path(execution['python_cache'])

                file_out.write(f'export CALDB="{caldb_path}"\n')
                file_out.write(f'export PATH="{conda_path}/bin:$PATH"\n')
                file_out.write(f'export PATH="{conda_path}/lib:$PATH"\n')
                file_out.write(f'export PYTHON_EGG_CACHE="{python_cache}"\n')
                file_out.write(f'source activate {env_name}\n')
                file_out.write(f'python {ctools_pipe_path}/background_sim.py {ctools_pipe_path}/{infile} {ctools_pipe_path}/{in_jobs} {str(counter + 1)} \n')
                file_out.write('source deactivate\n')
                file_out.close()

                list_script_bash.append(out_file_name)

            for file_bash in list_script_bash:
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

                exec_string += f"{file_bash}"
                print(exec_string)

                #p = subprocess.Popen(exec_string.split(" "),
                #                     stdout=subprocess.PIPE,
                #                     stderr=subprocess.PIPE
                #                     )
                # print stderr ad stdout for first job
                #if counter == 0 and execution['debug'] == "yes":
                #    (result, error) = p.communicate()
                #    print(result, error)

    if args.models:
        # model creation is done only locally since it's not a heavy computation
        infile = args.models
        config_in = yaml.safe_load(open(infile))

        if execution['mode'] == "local":
            p = subprocess.Popen(
                ['python', 'model_creation.py', infile, in_jobs]
            )
            # if everything goes well, the output is None
            # check this just for the first job
            (result, error) = p.communicate()
            print(result, error)

    if args.simulation:
        # input for simulation
        in_simu = args.simulation
        config_in = yaml.safe_load(open(in_simu))
        ctobssim_input = config_in['ctobssim']
        realizations = ctobssim_input['realizations']

        # load models
        models_path = create_path(ctobssim_input['models_in']['xml_path'])
        models_list = glob.glob(f"{models_path}/*")
        if len(models_list) == 0:
            print("No input model")
            sys.exit()

        max_models = config_in['source']['max_sources']

        ctools_pipe_path = create_path(jobs_exe['exe']['software_path'])

        # select IRF to choose the proper background folder
        irf = IRFPicker(config_in, ctools_pipe_path)
        name_irf = irf.irf_pick()

        # load backgrounds
        backgrounds_path = create_path(ctobssim_input['background_path'])
        fits_background_list = glob.glob(f"{backgrounds_path}/{name_irf}/background*.fits")
        if len(fits_background_list) == 0:
            print(f"No background for IRF {name_irf}")
            sys.exit()

        if config_in['source']['type'] == "GW":
            point_path = create_path(config_in['source']['pointings_path'])
            event_list = pd.read_csv(glob.glob(f"{point_path}/BNS*.txt")[0], sep=" ")
            events = event_list[:max_models]
            runs = events['run']
            mergers = events['MergerID']
            IRF_zenith = events['Mean Altitude']
            observatory = events['Observatory']
            opt_point_path = f"{point_path}/optimized_pointings"

        if execution['mode'] == "local":
            # loop over the models
            for counter, model in enumerate(models_list[:max_models]):
                # loop over the realizations for each model
                for sim in range(realizations):
                    path_background_to_use = fits_background_list[sim]
                    p = subprocess.Popen(
                        ['python',
                         'simulation_analysis.py',
                         in_simu,
                         in_jobs,
                         glob.glob(f"{model}/*xml")[0],
                         str(sim + 1),
                         path_background_to_use
                         ],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    # if everything goes well, the output is None
                    # check this just for the first job
                    #if counter == 0 and sim == 0:
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
            # loop over the models
            for counter, model in enumerate(models_list[:max_models]):
                # loop over the realizations for each model
                for sim in range(realizations):
                    path_background_to_use = fits_background_list[sim]
                    p = subprocess.Popen(
                        [*exec_string.split(" "),
                         "python",
                         'simulation_analysis.py',
                         in_simu,
                         in_jobs,
                         glob.glob(f"{model}/*xml")[0],
                         str(sim + 1),
                         path_background_to_use
                         ],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    # print stderr ad stdout for first job
                    if execution['debug'] == "yes" and counter == 0 and sim == 0:
                        (result, error) = p.communicate()
                        print(result, error)

        elif execution['mode'] == "qsub":
            details = execution['details']
            # create scripts to lunch on a folder
            folder_launch = create_path(f"{execution['path']}/tmp/sim/scripts")

            conda_path = create_path(execution['conda']['conda_path'])
            env = Env()
            list_script_bash = []
            # loop over all the models
            for model in models_list[:max_models]:
                # loop over the realizations for each model
                for sim in range(realizations):
                    model_to_out = glob.glob(f"{model}/*xml")[0]

                    path_background_to_use = fits_background_list[sim]

                    script_name = f"sim_launch_{model}_{str(sim).zfill(3)}.sh"
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
                    python_cache = create_path(execution['python_cache'])

                    file_out.write(f'export CALDB="{caldb_path}"\n')
                    file_out.write(f'export PATH="{conda_path}/bin:$PATH"\n')
                    file_out.write(f'export PATH="{conda_path}/lib:$PATH"\n')
                    file_out.write(f'export PYTHON_EGG_CACHE="{python_cache}"\n')
                    file_out.write(f'source activate {env_name}\n')
                    file_out.write(f'python {ctools_pipe_path}/simulation_analysis.py {ctools_pipe_path}/{in_simu} {ctools_pipe_path}/{in_jobs} {model_to_out}  {str(sim + 1)} {path_background_to_use} \n')
                    file_out.write('source deactivate\n')
                    file_out.close()

                    list_script_bash.append(out_file_name)

            for file_bash in list_script_bash:
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

                exec_string += f"{file_bash}"
                print(exec_string)
