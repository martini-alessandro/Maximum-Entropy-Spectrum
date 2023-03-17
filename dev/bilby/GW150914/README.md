## Things you may want to change

- The `mesa_m` parameter, set in the prior
- The number of CPUs available, set in `request-cpus`
- The path to the analysis executable bilby_pipe_mesa.py, needs to be updated

## To run on a cluster under HTCondor

$ bilby_pipe dynesty_config_mesa.ini
$ condor_submit_dag outdir_mesa_A/submit/dag_GW150914_mesa.submit

## To run directly
$ bilby_pipe dynesty_config_mesa.ini
$ bash outdir_mesa_A/submit/bash_GW150914_mesa.sh

