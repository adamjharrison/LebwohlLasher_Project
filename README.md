# LebwohlLasher_Project

This repository contains code for running the Lebwohl-Lasher Model of liquid crystals. Different versions of this code implement various strategies to speedup the serial python code

# Installation Instruction

### Prerequisites
This code requires different dependencies:
- Python 3.12+
- Numba
- Numpy
- Cython
- Distutils
- mpi4py
- MPI
- OpenMP

### Clone the repository
Clone this repository using the code below:

```
git clone https://github.com/adamjharrison/LebwohlLasher_Project
cd LebwohlLasher_Project
```

# Running the code
For Serial,Numba and Numpy implementations:

```
python LebowhlLasher<>.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG> <SAVEFILE>
```

For Cython implementations(current setup files are for use on Apple Silicon):

```
python setup_LebowhlLasher_cython<>.py build_ext -fi
python run_LebowhlLasher_cython<>.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG> <SAVEFILE>
```

For Cython-MPI(current setup files are for use on Apple Silicon):

```
python setup_LebowhlLasher_cython_mpi.py build_ext -fi
mpiexec -np <> python run_LebowhlLasher_cython_mpi.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG> <SAVEFILE>
```

For MPI:

```
mpiexec -np <> python LebowhlLasher_mpi4py.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG> <SAVEFILE>
```
