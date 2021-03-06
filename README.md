# PowerModelsDistribution.jl

<img src="https://lanl-ansi.github.io/PowerModelsDistribution.jl/dev/assets/logo.svg" align="left" width="200" alt="PowerModelsDistribution logo">

Release:
[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://lanl-ansi.github.io/PowerModelsDistribution.jl/stable/)

Dev:
[![Build Status](https://travis-ci.org/lanl-ansi/PowerModelsDistribution.jl.svg?branch=master)](https://travis-ci.org/lanl-ansi/PowerModelsDistribution.jl)
[![codecov](https://codecov.io/gh/lanl-ansi/PowerModelsDistribution.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/PowerModelsDistribution.jl)
[![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://lanl-ansi.github.io/PowerModelsDistribution.jl/latest/)

PowerModelsDistribution.jl is an extention package of PowerModels.jl for Steady-State Power Distribution Network Optimization.  It is designed to enable computational evaluation of emerging power network formulations and algorithms in a common platform.  The code is engineered to decouple problem specifications (e.g. Power Flow, Optimal Power Flow, ...) from the power network formulations (e.g. AC, linear-approximation, SOC-relaxation, ...).
This enables the definition of a wide variety of power network formulations and their comparison on common problem specifications.

## Core Problem Specifications

- Power Flow (pf)
- Optimal Power Flow (opf), for the Bus Injection Model (BIM) as well as the Branch Flow Model (BFM)
- Continuous load shed, minimum load delta (mld), for the Branch Flow Model (LPLinUBFPowerModel), AC Polar (ACPPowerModel), and Network Flow Approximation (NFAPowerModel)

## Core Network Formulations

- AC (polar and rectangular coordinates)
- SDP BFM relaxation
- SOC BFM and BIM relaxation (W-space)
- Linear approximation (LinDist3Flow and simplified unbalanced DistFlow)

## Network Data Formats

- Matlab ".m" files (extended for three-phase)
- OpenDSS ".dss" files

**Warning:** This package is under active development and may change drastically without warning.

## Development

Community-driven development and enhancement of PowerModelsDistribution are welcome and encouraged. Please fork this repository and share your contributions to the master with pull requests.

## Acknowledgments

This code has been developed as part of the Advanced Network Science Initiative at Los Alamos National Laboratory.  The primary developers are David Fobes(@pseudocubic) and Carleton Coffrin(@ccoffrin) with support from the following contributors:

- Sander Claeys (@sanderclaeys) KU Leuven, transformer models and ACR formulation
- Frederik Geth (@frederikgeth) CSIRO, Distribution modeling advise

## License

This code is provided under a BSD license as part of the Multi-Infrastructure Control and Optimization Toolkit (MICOT) project, LA-CC-13-108.
