# ctools_pipe
Pipelines for GRB/GW simulations task

### package needed (on top of ctools environment) 
#### (this will soon be changed with a setup.py file)
- yaml: `pip install PyYAML`
- astropy: `pip install astropy`
- environs: `pip install environs`
- `wget https://raw.githubusercontent.com/HESOFTS/sexten_2017/master/model_creator/scriptModel_variable.py`

## Description

The scripts are developed and tested inside a ctools Anaconda environment in order ot use gammalib/ctools functionalities. I personally advise to use the Anaconda package. It can be also easily used in a job submitting system (such as LSF or others), provided that the right path are given. This set of scripts is built in order to handle this correctly.

**NOTE**: remember to load the IRFs in the ctools environment. Just download them from the webpage https://www.cta-observatory.org/science/cta-performance/ and **unpack the FITS version in the parent folder of the $CALDB folder** (visible from inside the environment). Other supported IRFs are those that can be downloaded from redmine.

### what's inside
**1. simulation of background** to be saved and reused afterwards. The configuration file has to be prepared according to the job submitting system.
```
python ctools_pipe.py --background background.yaml --jobs jobs*.yaml
```
or
```
python ctools_pipe.py -b background.yaml -j jobs*.yaml
```

**2. creation of XML models** (ctools compliant) for sources with variable spectra over time. The idea is to put all the sources in the same point in space, but each of them have also a lightcurve attached, which is shaped as a squared wave so that it's zero everywhere (the source is OFF), expect for a time window in which the source is ON.
The creation of the model have to be done **AFTER** having sourced the ctools environment because the scrit is using ctools and gammalib, that are installed inside that environment.
```
python ctools_pipe.py --models model_input.yaml --jobs jobs_local.yaml
```
or
```
python ctools_pipe.py -m model_input.yaml -j jobs_local.yaml
```


**3. the simulation of the source.** Every source realization is attached to a background which was previously simulated: this saves time when thousands of simulations have to be done.

#### HELP:
You can always type: 
```
python ctools_pipe.py --help
```
to get the help message on how to use the script.

### LAPP (MUST) usage: READ CAREFULLY!!! 
- change the `batch_example_lapp.sh` to your `batch_LAPP.sh` after having set the `LAPP_APP_SHARED` folder properly.
- `source batch_LAPP.sh` in order to se the proper environmental variables
- change the `variables_example.sh` properly and create your `variable.sh` (see a `variable_LAPP.sh` example file).
- `source variable.sh`

`batch_LAPP.sh` and `variable.sh` have to be created once and then they just have to be sourced at the beginning of each new session.

**How to launch a job on MUST:**
- activate the environment: `source activate ctools_p3`
- launch command, i.e. `python ctools_pipe.py --background background.yaml --jobs jobs_qsub.yaml`. This will **NOT** launch a job but it will print on screen the string to be used for submitting the job (try with few files...you'll see). you can pipe it to a `launch_multiple.sh`, for example, like this:
```
python ctools_pipe.py --background background.yaml --jobs jobs_qsub.yaml > launch_multiple.sh
```
- EXIT from the environment: `source deactivate`
- submit the jobs, i.e. `source launch_multiple.sh`

...et voilà!!!
