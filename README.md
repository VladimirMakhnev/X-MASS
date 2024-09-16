# X-MASS

## Introduction

A Python tool (X-MASS) that allows the massive set of ABSCO tables to be calculated using the HAPI software packagee with complete utilization of the parameters’ accuracy in HITRAN, including sophisticated line shapes

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
- multiprocessing/asyncio/nest_asyncio

## Installation

[to be updated]

Clone the repository and install the necessary dependencies:

```
bash
git clone https://github.com/VladimirMakhnev/X-MASS.git
cd x-mass 
```

## How to use 

Prepare input files at the directory with scripts.
Launch the main.py script:

`python main.py `

where **filenames.inp** contains info about all other files and method of calculation: 

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


## Params file content

The configuration of **params.inp** file is consistent with ABSCO table [manual](https://docserver.gesdisc.eosdis.nasa.gov/public/project/OCO/ABSCO_UsersGuide_20170724_corr2_v5.0.pdf). 
First 11 lines correspond to the attributes from the manual, while next line describes: 
















