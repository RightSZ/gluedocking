# gluedocking

## Overview

gluedocking is an R package designed to streamline molecular docking workflows. It provides a comprehensive set of functions for downloading protein and ligand structures, preparing docking files, running AutoDock Vina docking calculations, and analyzing docking results. The package wraps common molecular docking tools, enabling researchers to efficiently perform virtual screening and drug discovery studies within the R environment.

## Features

-   Download protein structures from the RCSB PDB database
-   Download ligand structures from the PubChem database
-   Prepare receptor and ligand files using AutoDock Tools
-   Convert molecular file formats using OpenBabel
-   Automatically calculate docking box parameters
-   Generate AutoDock Vina configuration files
-   Run molecular docking calculations
-   Parse and analyze docking results

Find out more at <https://github.com/RightSZ/gluedocking>

## Installation

### Prerequisites

Before using gluedocking, you need to install several external tools:

1.  **Python** - Required for running AutoDock Tools scripts
2.  **MGLTools** - Contains prepare_receptor4.py and prepare_ligand4.py scripts
3.  **OpenBabel** - For molecular file format conversion
4.  **AutoDock Vina** - For molecular docking calculations

You can install the development version of gluedocking from [GitHub](https://github.com/RightSZ/gluedocking) with:

``` r
# Install the development version from GitHub
if(!require(devtools)) install.packages("devtools")
devtools::install_github("RightSZ/gluedocking")
```

## License

This package is licensed under the MIT License. See the LICENSE file for details.
