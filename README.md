# X-MASS

A Python tool (X-MASS) that allows the massive set of ABSCO tables to be calculated using the HAPI software packagee with complete utilization of the parameters’ accuracy in HITRAN, including sophisticated line shapes

## Features

- Utilizes HITRAN's advanced line shapes and parameters for accuracy.
- Generates ABSCO tables in HDF5 format.
- Provides customizable input grids for pressure, temperature, and volume mixing ratios.
- Efficient and effectively compressed output files for rapid retrieval.

## Requirements

- Python 3.x
- [HAPI](github.com/hitranonline/hapi) (HITRAN Application Programming Interface)
- NumPy
- matplotlib
- [h5py](h5py.org)
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

Launch the main.py script in folder with all input files:

`pythoon main.py params.inp pres_pRT.inp temps.inp vms.inp wn.inp 01.H2O.SDV.HITRAN2020.25wing.hdf5 MULTITHREADING `

where: 

| **File**                              | **Description**                    |
|---------------------------------------|------------------------------------|
| **params.inp**                        | a file with parameters             |
| **pres_pRT.inp**                      | a file with pressure profile       |
| **temps.inp**                         | a file with temperature profile    |
| **vms.inp**                           | a file with volume mixing ratio    |
| **wn.inp**                            | a file with wavenumber grid length |
| **01.H2O.SDV.HITRAN2020.25wing.hdf5** | a filename of resulting HDF5 file  |
| **MULTITHREADING**                    | a method of parallelization        |

## Options of calculation method

**PLAIN** -- the sequential calculation of cross-sections one at the time;

**PC** -- the ```asyncio``` calculation;

**MULTITHREADING** -- the ```multiprocessing``` calculation.










### Basic Markdown Format Rules:
- **Headers:** Use `#`, `##`, `###` to create header levels.
- **Emphasis:** Use `*text*` or `_text_` for italics, and `**text**` for bold.
- **Lists:** Use `-` or `*` for unordered lists, and `1.` for ordered lists.
- **Links:** `[Link text](URL)` to create links.
- **Images:** `![Alt text](Image URL)` for images.
- **Code blocks:** Use triple backticks (```) for code blocks, and single backticks (`) for inline code.```
- **Tables:** Create tables using pipes (`|`) and hyphens for column separation.






