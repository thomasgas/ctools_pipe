# ctools_pipe
Pipelines for GRB/GW simulations task

### package needed (on top of ctools environment) 
#### (this will soon be changed with a setup.py file)
- yaml: `pip install PyYAML`
- astropy: `pip install astropy`
- environs: `pip install environs`

## Description

The scripts are developed and tested inside a ctools Anaconda environment in order ot use gammalib/ctools functionalities. I personally advise to use the Anaconda package. It can be also easily used in a job submitting system (such as LSF or others), provided that the right path are given. This set of scripts is built in order to handle this correctly.

### what's inside
1. simulation of background to be saved and reused afterwards. The configuration file has to be prepared according to the job submitting system.
```
python ctools_pipe.py --background background.yaml --jobs jobs*.yaml
```
or
```
python ctools_pipe.py -b background.yaml -j jobs*.yaml
```

2. creation of XML models (ctools compliant) for sources with variable spectra over time. The idea is to put all the sources in the same point in space, but each of them have also a lightcurve attached, which is shaped as a squared wave so that it's zero everywhere (the source is OFF), expect for a time window in which the source is ON.
3. the simulation chain of the source. Every source realization is attached to a background which was previously simulated: this saves time when thousands of simulations have to be done.

#### HELP:
You can always type: 
```
python ctools_pipe.py --help
```
to get the help message on how to use the script.

### LAPP (MUST) usage
- `source batch_example_lapp.sh` after having set the `LAPP_APP_SHARED` folder properly.
- change the `variables_example.sh` properly and create your `variable.sh`.
- `source variable.sh`

`batch_example_lapp.sh` and `variable.sh` have to be created once and then they just have to be sourced at the beginning of each job submission.

Then all the variables are set and one can proceed to lunch the jobs. The environmental variables are needed because when the job is submitted, all environmental variables are lost (if not passes) and the relative paths are also lost.
