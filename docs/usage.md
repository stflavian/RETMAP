# How to use

Adsorpyon is a command line application, meaning that it is intended to be used solely
from a terminal. 

- Copy the updated `run.py` file to your working folder.
- Copy the adsorption data that needs to be processed to the working folder.
- Create the properties file for the molecule used to obtain the adsorption data.  
- Create the `config.in` file in the same folder and add the necessary information needed
for the application to run.
- Activate the virtual environment containing the needed libraries for the 
application to run.[^1]
- Run the application using `./run.py` on Linux/MacOS or `python run.py` on any 
platform.

>**NOTE**
> 
> Depending on the settings you use in `config.in` the application will create separate
folders for plots and written output. To keep the working folder clean it is
recommended that you also place the adsorption data and the properties file in a folder
named `Input`.

There is no notification or update tracker implemented in the application at the 
moment. The application will stop once it completes all the tasks mentioned in 
`config.in`, or once it encounters an error in the program. The application creates
files and plots in sequence, meaning that files will be created as the code is running.

Since the configuration file is written in plain text, it is easy to integrate the 
application inside other shell scripts to automate the processing procedure. 

[^1]: This step is not necessary in case you already have Numpy, Scipy, and Matplotlib 
installed globally on your Python version. Still, we **highly** recommend that you 
install these libraries inside a virtual environment to minimize the chance of errors.