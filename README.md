# ctools_pipe
Pipelines for GRB/GW simulations task

The scripts are developed and tested inside a ctools Anaconda environment in order ot use gammalib/ctools functionalities. I personally advise to use the Anaconda package. It can be also easily used in a job submitting system (such as LSF or others), provided that the right path are given. This set of scripts is built in order to handle this correctly.

### what's inside
1. simulation of background to be saved and reused afterwards
2. creation of XML models (ctools compliant) for sources with variable spectra over time. The idea is to put all the sources in the same point in space, but each of them have also a lightcurve attached, which is shaped as a squared wave so that it's zero everywhere (the source is OFF), expect for a time window in which the source is ON.
3. the simulation chain of the source. Every source realization is attached to a background which was previously simulated: this saves time when thousands of simulations have to be done.
