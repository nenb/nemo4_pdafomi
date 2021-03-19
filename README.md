# NEMO4 PDAFOMI

NEMO4 PDAFOMI contains the [PDAFOMI](http://pdaf.awi.de/trac/wiki) attachments to the [NEMO4 ocean model](https://www.nemo-ocean.eu/).

These attachments are currently a **work in progress** and should be used with caution.

<hr>

## Requirements

In addition to the requirements for running PDAFOMI and NEMO, this project also requires [pre-commit](https://pre-commit.com/) and [FORD](https://github.com/Fortran-FOSS-Programmers/ford) for full functionality. These requirements can be installed into your local Python environment with  `pip`:
``` python
python -m pip install -r requirements.txt
```

Once `pre-commit` is installed, please install the required git hooks with `pre-commit install`.

## Project Layout

Attachments are currently located at `src/nemo_r4.0.4/ext/PDAF`. There are two subdirectories: `pdaf_bindings` contains the PDAF interface and callback routines, and `nemo_src` contains modifications to the original NEMO4 source code.

## Documentation

The documentation for this project is built with [FORD](https://github.com/Fortran-FOSS-Programmers/ford). Documentation is available [here](https://nenb.github.io/nemo4_pdafomi/).
