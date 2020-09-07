# Condor Scripts
When running code on Bristol's HPC grid, they don't like if I run code directly on the submit node

Instead submit jobs to the worker nodes with the scripts in this dir

#### Usage- Pull Study:
First build target pull-study and move it to this directory :^)

Then run

`condor_submit pull.job`

from a machine with access to the right condor config (i.e. sc01)


NOTE: this is all obviously a massive hack, if you look in pull.job I copy shared libraries around which probably isn't good
Any given job might just fail to run because of Reasons... if you run 10 jobs one might work?

