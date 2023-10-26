# Installation

Adsorpyon is an application built using Python 3.7, meaning that it should be
compatible with most Python versions in use at the moment. In terms of libraries,
adsorpyon is built using Numpy, Scipy, and Matplotlib, all of which should be installed
on the machine for the application to work. 

To decrease the chances of errors, it is recommended that the application is run inside 
a Python [virtual environment](https://docs.python.org/3/library/venv.html) containing 
the exact versions of the libraries listed in the `requirements.txt` file. Also, since 
the application was created in Linux, little to no testing has been done on Windows. 
For this reason, we encourage you to install adsorpyon inside 
[Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) 
or to [contact us](contact.md) in case any unexpected bugs are encountered.

### Install Python and download adsorpyon

- Download and install Python 3.7 or a more recent version from the 
[official website](https://www.python.org/downloads/).
- Verify the installation by opening a terminal and typing `python --version`. The 
output should display the version of Python you just installed.
- Use `git clone https://github.com/stflavian/adsorpyon.git` or download the **.zip**
file containing the code to your desired location and unzip the archive. 

### Create a virtual environment

This step is not necessary in case you already have Numpy, Scipy, and Matplotlib 
installed globally on your Python version. Still, we **highly** recommend that you 
install these libraries inside a virtual environment to minimize the chance of errors.

- Using the terminal, create a virtual environment for your project by running 
`python -m venv env` in your project directory. This will create a new directory called 
env that contains a copy of the Python interpreter and pip package manager.
- Activate the virtual environment by running `source env/bin/activate` on Linux/MacOS or 
`env\Scripts\activate.bat` on Windows.

>**Note**
>
>Make sure that the Python virtual environment is active during the next step of the 
installation and when using the code. When the environment is active, its name is 
displayed before the terminal prompt inside parenthesis, e.g. `(env) John@PC: ~`.


### Install requirements and set up the run.py file
- Install dependencies by running `pip install -r requirements.txt` in your project 
directory. This will install all the packages listed in the `requirements.txt` file.
- Open the `run.py` file found in the `/examples/run` folder.
- Change the path in line 4 from `"/home/user/path/to/adsorpyon"` to the path in which 
you downloaded adsorpyon, e.g. `"/home/john/Desktop/adsorpyon"` on Linux or 
`C:\Users\John\Desktop\adsorpyon` on Windows.

>**Note**
>
> Make sure to use the updated version of the `run.py` when creating new working 
folders. 

