# X-MASS

A Python tool (**X-MASS**) for the batch calculation of absorption cross-section (ABSCO) tables from the HITRAN database using the [HAPI](https://github.com/hitranonline/hapi) library. X-MASS computes cross-sections over a pressure × temperature × volume-mixing-ratio grid and stores them in the NASA ABSCO HDF5 table format, making full use of the line-shape parameters provided in HITRAN — including the speed-dependent Voigt and Hartmann-Tran profiles and first-order (Rosenkranz) line mixing — wherever these parameters are available in the database (the Voigt profile is used for the remaining lines).

## Features

- HITRAN2020 line-by-line data fetched directly from [HITRANonline](https://hitran.org).
- Advanced line profiles via HAPI: Voigt, speed-dependent Voigt, Hartmann-Tran, Lorentz; optional first-order (Rosenkranz) line mixing for Voigt/SDVoigt.
- 25 cm⁻¹ absolute wing cutoff (the recommended community standard and the MT_CKD "local line" convention), with an optional Lorentz-pedestal ("plinth") removal for use together with the MT_CKD continuum without double counting.
- NASA-ABSCO-format 4-D HDF5 output (`P × T × VMR × ν`), gzip+shuffle compressed, written incrementally — memory use stays at one spectrum regardless of grid size.
- Multiprocessing parallelism with per-worker line-list caching; a failed grid point is reported and left as a NaN row instead of killing the run.
- Command-line interface with input validation (`--validate-only`) and full run logging to `output.log`.

## Requirements

- Python ≥ 3.8 (Windows and Linux are both supported)
- `numpy`, `h5py`, `matplotlib`, `nest_asyncio`, `hitran-api` (HAPI; tested with 1.3.0.0) — see `requirements.txt`
- Internet access to hitran.org at run time (the line list is fetched on every run)

## Installation

```bash
git clone https://github.com/VladimirMakhnev/X-MASS.git
cd X-MASS
python -m venv venv
source venv/bin/activate          # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Quickstart

A complete small example (CO fundamental band, 2000–2200 cm⁻¹, 2×2×2 grid, ~15 s on 4 cores) lives in `examples/CO_quickstart`:

```bash
cd examples/CO_quickstart
python ../../main.py                 # runs the calculation -> CO_quickstart.hdf5
python check_output.py               # verifies the result against reference values
python ../plot_table.py CO_quickstart.hdf5   # prints the structure, saves a plot
```

## Usage

```
python main.py [CONFIG] [options]

positional:
  CONFIG            master input file (default: filenames.inp)

options:
  --method {PLAIN,PC,MULTITHREADING}   override the calculation method
  --cores N         override Number_cores from the params file
  --keep-dat        also write the legacy per-point .dat text files
  --validate-only   check all inputs, print a run summary, exit
  --quiet           no console output (everything goes to output.log)
  --version         show version
```

All input files are read from the current working directory (run X-MASS from the directory that holds your input set).

### The master input file (`filenames.inp`)

One entry per line, in this order:

| Line | Content                                    |
|------|--------------------------------------------|
| 1    | ignored (historical script name)           |
| 2    | parameters file (`params.inp`)             |
| 3    | pressure grid file (atm, one per line)     |
| 4    | temperature grid file (K, Np rows × Nt columns) |
| 5    | broadener VMR file (one value per line)    |
| 6    | wavenumber file (number of grid points)    |
| 7    | output HDF5 file name                      |
| 8    | method: `PLAIN`, `PC`, or `MULTITHREADING` |

### The parameters file (`params.inp`)

`key:value` pairs, one per line. The first 16 lines follow the ABSCO table
convention (see the [ABSCO User Guide](https://docserver.gesdisc.eosdis.nasa.gov/public/project/OCO/ABSCO_UsersGuide_20170724_corr2_v5.0.pdf));
lines 1–16 are positional, the switches below them are recognized by name:

| Key                | Meaning                                                        |
|--------------------|----------------------------------------------------------------|
| `version`, `addl_ident`, `gas_name`, `comment`, `Gas_Q_abs_Absorbtion` | metadata copied into the HDF5 attributes |
| `wn_begin`, `wn_end` | spectral range (cm⁻¹)                                        |
| `Gas_index`        | HITRAN molecule number (e.g. `05` = CO)                        |
| `Pressure`, `Temperature`, `Broadener_Q_brd_VMR` | grid dimensions (informational) |
| `broadener_name`   | stored in the broadener dataset attribute (e.g. `air`)         |
| `Broadener_Index`  | broadener index for the dataset names                          |
| `Wavenumber`       | number of wavenumber grid points (must match `wn.inp`)         |
| `Number_cores`     | worker count for `MULTITHREADING`                              |
| `Profile`          | `Voigt`, `SDVoigt`, `HT`, or `Lorentz`                         |
| `Profile_group`    | HAPI parameter group to fetch (e.g. `160-char`)                |
| `Line_mixing`      | `ON`/`OFF`: first-order (Rosenkranz) line mixing for Voigt/SDVoigt, where HITRAN provides the parameters |
| `Remove_pedestal`  | `ON`/`OFF`: subtract the 25 cm⁻¹ Lorentz pedestal (MT_CKD-consistent local-line tables) |
| `Keep_dat`         | `ON`/`OFF`: also write the legacy per-point `.dat` text files  |
| `Legacy_VMS_names` | `ON`/`OFF`: name the broadener dataset `Broadener_XX_VMS` (pre-0.9 convention) instead of the ABSCO-conformant `..._VMR` |

### Calculation methods

- **MULTITHREADING** (recommended) — `multiprocessing` pool; each worker loads the line list once, results stream back and are written to HDF5 as they finish.
- **PLAIN** — sequential, same code path in a single process.
- **PC** — `asyncio` compatibility mode (GIL-bound; no speedup for CPU-bound work).

## Output format

The HDF5 file follows the ABSCO table layout:

| Dataset               | Shape               | Units          |
|-----------------------|---------------------|----------------|
| `Gas_XX_Absorption`   | (Np, Nt, Nvmr, Nwn) | cm²/molecule   |
| `Pressure`            | (Np,)               | atm            |
| `Temperature`         | (Np, Nt)            | K              |
| `Broadener_YY_VMR`    | (Nvmr,)             | dimensionless  |
| `Wavenumber`          | (Nwn,)              | cm⁻¹           |
| `Gas_Index`, `Broadener_Index` | scalars    |                |

Failed grid points (e.g. NaN inputs) are reported at the end of the run, left as NaN rows (`missing_value` attribute), and make the process exit with a nonzero code. The absorption dataset is chunked per grid point and compressed with shuffle+gzip.

## Line-shape conventions

Cross-sections use an absolute wing cutoff of 25 cm⁻¹, the de facto Earth remote-sensing standard, matching both the recently recommended community wing-cutoff standard (Gharib-Nezhad et al. 2024, RASTI 3, 44) and the "local line shape" definition of the MT_CKD water-vapor continuum (Mlawer et al. 2023, JQSRT 306, 108645). Lines centered up to 25 cm⁻¹ outside the requested range are included, so edge points are complete. With `Remove_pedestal:ON` the Lorentz pedestal value at ±25 cm⁻¹ is subtracted across each line's window, producing MT_CKD-consistent local-line tables.

Note that beyond-Voigt parameters are available in HITRAN for a subset of molecules and transitions; where they are absent, the calculation falls back to the Voigt parametrization (this is a property of the database, not of the code).

## Troubleshooting

- **`KeyError` in workers with older HAPI**: do not call `cache2storage()` after `fetch_by_ids()` — in HAPI ≤ 1.3.0.0 it corrupts the stored header of tables fetched with extra parameter groups (X-MASS avoids this internally).
- **HAPI banner prints even with `--quiet`**: HAPI prints its banner at import time; it is harmless.
- **No internet**: the line list is fetched from hitran.org on every run; offline runs are currently not supported.
- **`Number_cores` larger than the machine**: the run aborts; use `--cores` to override without editing files.

## Credits

If you use X-MASS in your research or software development, please cite:
V.Yu. Makhnev, I.E. Gordon, L.S. Rothman, R.J. Hargreaves (in prep.)
