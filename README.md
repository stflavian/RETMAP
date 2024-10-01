# RETMAP

RETMAP (REfined Thermodynamical Model for Adsorption Prediction) is a Python application 
for calculating adsorption properties, such as isobars and isotherms, using experimental 
or simulated data. For this, it makes use of Polanyi's adsorption potential model, which
is an effective method for predicting adsorption properties in the neighbourhood of
known data. 

The documentation for the software can be found at https://stflavian.github.io/adsorpyon/.

## Quick setup

RETMAP is an application built using Python 3.7, meaning that it should be
compatible with most Python versions in use at the moment. In terms of libraries,
RETMAP is built using Numpy, Scipy, and Matplotlib, all of which should be installed
on the machine for the application to work. 

To decrease the chances of errors, it is recommended that the application is run inside 
a Python [virtual environment](https://docs.python.org/3/library/venv.html) containing the exact versions of the libraries 
listed in the `requirements.txt` file. Also, since the application was created in Linux, 
little to no testing has been done on Windows. For this reason, we encourage you to 
install RETMAP inside [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) or to contact us in case 
any unexpected bugs are encountered.

### Setting it up on Linux (Ubuntu, Debian, Fedora, Arch, etc.)

> git clone https://github.com/stflavian/adsorpyon.git \
> cd adsorpyon \
> pip install -r requirements.txt \
> vim /examples/run/run.py

Change the path in line 4 from `"/home/user/path/to/adsorpyon"` to the path in which 
you downloaded RETMAP, e.g. `"/home/john/Desktop/adsorpyon"`.

### Running RETMAP

In the working folder run:

> ./run.py

## Terms of use

RETMAP is an open-source project licenced under the MIT licence. Usage of this software
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
