# X-MASS

## Introduction

A Python tool (**X-MASS**) that allows the massive set of ABSCO tables to be calculated using the HAPI software packagee with complete utilization of the parameters’ accuracy in HITRAN, including sophisticated line shapes

## Features

- Utilizes HITRAN's advanced line shapes and parameters for accuracy.
- Generates ABSCO tables in HDF5 format.
- Provides customizable input grids for pressure, temperature, and volume mixing ratios.
- Efficient and effectively compressed output files for rapid retrieval.

## Requirements

- Python 3.x
- [HAPI](https://github.com/hitranonline/hapi) (HITRAN Application Programming Interface)
- NumPy
- matplotlib
- [h5py](https://h5py.org)
- multiprocessing/asyncio/[nest_asyncio](https://pypi.org/project/nest-asyncio/)

## Installation

[to be updated]

Clone the repository and install the necessary dependencies:

```
bash
git clone https://github.com/VladimirMakhnev/X-MASS.git
cd x-mass 
```

## How to use 

1. Download the package into the working folder.

	*N.B.: Beware of lack of free space on the drive!*
	
	*N.B.: All input files should be in the working directory!*

2. Setup the **filenames.inp** which contains info about all other files and method of calculation: 

| **File**                              | **Description**                    |
|---------------------------------------|------------------------------------|
| **params.inp**                        | a file with parameters             |
| **pres_pRT.inp**                      | a file with pressure profile       |
|                                       | (in atm)                           |
| **temps.inp**                         | a file with temperature profile    |
|                                       | (in K)                             |
| **vms.inp**                           | a file with volume mixing ratio    |
| **wn.inp**                            | a file with wavenumber grid length |
| **01.H2O.SDV.HITRAN2020.25wing.hdf5** | a filename of resulting HDF5 file  |
| **MULTITHREADING**                    | a method of parallelization        |

### Options of calculation method

**PLAIN** -- the sequential calculation of cross-sections one at the time;

**PC** -- the ```asyncio``` calculation;

**MULTITHREADING** -- the ```multiprocessing``` calculation.


3. Setup the **params.inp** which contains info about details of calculations. 

The first 16 lines of **params.inp** file are consistent with ABSCO table:
 [manual](https://docserver.gesdisc.eosdis.nasa.gov/public/project/OCO/ABSCO_UsersGuide_20170724_corr2_v5.0.pdf). 
The next lines are:

| **Line**                              | **Description**                    |
|---------------------------------------|------------------------------------|
| **Wavenumber**                        | a number of grid points            |
| **Number_cores**                      | amount of computer cores allocated |
|                                       | for calculations (MULTITHREADING)  |
| **Profile**                           | a line profile used in calculations|
| **Profile_group**                     | a line profile parameter group     |

4. Prepare **pres_pRT.inp** pressure input file -- a column of N_p pressures.

	*N.B.: all pressures are in atm!*

5. Prepare **temps.inp** file -- a 2-D grid of N_p x N_t temperature values. 

    *N.B.: N_t values is needed for each pressure value! *

6. Prepare **vms.inp** file -- a column of N_vms volume mixing ratio. 

7. Setup number of grid points in **wn.inp**.

    *N.B.: the value should be the same as in **params.inp**!*

0. Launch the main.py script:

`python main.py `

## Credits

If you use X-MASS in your research or software development, please cite it using the following reference:
V.Yu. Makhnev, I.E. Gordon, L.S. Rothman, R.J. Hargreaves (in prep.)


















