# turbulence-hyporheic-clogging

[![DOI](https://zenodo.org/badge/639960484.svg)](https://zenodo.org/badge/latestdoi/639960484)

## Installation

### Requirements:

    - OpenFOAM 7
    
### Building:

    - Go to `src/libs` and run `wmake`
    - Go to `src/solvers/cloggingFoam` and run `wmake`

## Running:

    - surface-water simulations run with `pimpleFoam`
    - subsurface simulations run with `cloggingFoam`


