Code for Tianlai Project

If you need help for setting up Github on your computer, please follow this link

https://help.github.com/articles/generating-ssh-keys#platform-linux

To use the python modules, you need some python libs:

numpy

scipy

mpi4py

matplotlib

and add the project path to PYTHONPATH

export PYTHONPATH="PATH/TO/project_TL/:"$PYTHONPATH

please find the pipeline in input/

to run the pipeline, please use the manager.py

$python manager.py yours.pipe
