# compaction-frost
This repository contains a compaction model for frost heave / frozen fringe. The model solves a time-dependent system of PDEs with one spatial dimension (vertical). The model is described in a manuscript that will be linked to as soon as possible.

# Dependencies
The code requires [FEniCSx](https://fenicsproject.org). The jupyter notebooks can be run 
through a [Docker](https://www.docker.com) container with the command:

`docker run --init -ti -p 8888:8888 -v $(pwd):/home/fenics/shared -w /home/fenics/shared dolfinx/lab:stable`

Some of the notebooks rely on [ipyparallel](https://ipyparallel.readthedocs.io/en/latest/), which is installed via pip in the relevant notebooks.

# Running the code
The jupyter notebooks in the *notebooks* directory show how to run the code and reproduce figures from a manuscript.

# Contents
The *notebooks* directory contains the jupyter notebooks.

The *source* directory contains the source code, which is arranged as follows:

1. **solvers.py** contains the finite element methods (e.g., weak forms, discretization, etc.) and time-stepping loops for the solving the PDEs. 

2. **constitutive.py** sets constitutive relations for the material behaviour such as the sediment consolidation law, temperature profile, saturation, and permeability. These functions are relied on the PDEs/weak forms.

3. **params.py** sets most of the physical and nondimensional parameters in the problem. The parameters that are not set here are varied between simulations (passed to the time_stepping function) as shown in the jupyter notebook.

4. **post_process.py** contains a interpolation function that helps with saving and plotting the results (converts from dolfinx functions to numpy arrays).

5. **wrapper.py** contains functions that are used in the ipyparallel codes for varying a parameter over a range of values in parallel.