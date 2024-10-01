# RETMAP

> [!WARNING]
> This software is still in alpha development, so certain features may not work as
> intended and the documentation may not be up to date. If you run into any problems, 
> please contact the maintainer or open a new issue. 

RETMAP (REfined Thermodynamical Model for Adsorption Prediction) is a Python 
application used for calculating adsorption properties, such as isobars and isotherms, 
using experimental or simulated data. For this, it makes use of Polanyi's adsorption 
potential theory, which is an effective method for predicting adsorption properties 
in the neighborhood of recorded data. 

The full documentation for the software can be found 
[here](https://stflavian.github.io/RETMAP/).

## Quick setup

> [!IMPORTANT]
> Most testing has been done using Python 3.12 on Windows and Linux. If you run into
> any issues when installing the application, please let us know by opening a new 
> issue.

RETMAP is an application built using Python 3.12, but should be compatible with all
supported versions of Python. In terms of dependencies, RETMAP is built using 
Numpy, Scipy, and Matplotlib, which are automatically installed alongside the 
application. RETMAP can be installed directly from PyPI using:

> pip install retmap

The application is platform independent, and should behave identically on all 
operating systems (Windows, MacOS, Linux, etc.). More details about the installation 
can be found [here](https://stflavian.github.io/RETMAP/installation/).

### Running RETMAP

In the working folder containing the configuration file, run:

> retmap config.in

## Terms of use

RETMAP is an open-source project licensed under the MIT licence. Usage of this software
in a scientific work should be cited accordingly.

## Authorship

This project started development during Flavian's Bachelor's End Project at Eindhoven 
University of Technology.

**Authors**: F. Stavarache [1], A. Luna-Triguero [2, 3], S. Calero [1, 3], and J. M. Vicent-Luna [1, 3]

[1] Materials Simulation & Modelling, Department of Applied Physics, Eindhoven 
University of Technology

[2] Energy Technology, Department of Mechanical Engineering, Eindhoven University of 
Technology

[3] Eindhoven Institute for Renewable Energy Systems (EIRES), Eindhoven University of 
Technology
