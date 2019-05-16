#To be sourced ONLY FOR LAPP USER WHICH SEND JOB FROM LAPPSL

# Change to your shared home
LAPP_APP_SHARED="/lapp_data/cta/[YOUR_LAPP_NAME]/whatever/path/you/want"

# Define the PBS_HOME , by default equal to the LAPP_APP_SHARED
export PBS_O_HOME=$LAPP_APP_SHARED

# overload the grid variable because it'is false with a local acces
export MPI_SHARED_HOME=$LAPP_APP_SHARED

# Define the initial path equal to the HOME is standard situation but
# as we don't mount the HOME directory on the jobmanager client we have to define it
# by defaut is equal to the LAPP_APP_SHARED
# If we use some other PBS global variable, PBS_O_INITDIR will be the root path

export PBS_O_INITDIR=$LAPP_APP_SHARED

# Inportant variable to optimize the cluster usage
# if this variable is define, the variable TMPDIR will be define on the jobmanager client
# with a local path to the local disque and can be use as local path on the jobmanager client
# the benefit to use this variable when the sharing file systeme is not requiered is :
# - not overloading the sharing filesysteme
# - using the local disk ( better performencies)
# - this area is correctly remove when the job is finnish

# Be careful the variable to define on the submitter host (lappsl) is
#  TMP_DIR but the variable to use on the jobmanager client ( on your job)  is TMPDIR

export TMP_DIR=/var/spool/pbs/tmpdir

# force the usage of the -V options to export to jobmanager server the previous variables
alias qsub='qsub -V'
