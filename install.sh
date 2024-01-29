#!/bin/bash

ADSORPYON_PATH=$(pwd)

python=$(python --version | awk '{print $1}')
version=$(python --version | awk '{print $2}')

if [[ $python != "Python" ]] ; then
  echo "Error: Python was not found on this system!"
  exit 1
else
  echo "Python was found on the system ..."
  major=$(echo "$version" | cut -d. -f1)
  minor=$(echo "$version" | cut -d. -f2)

  if [[ $major != 3 || $minor -lt 7 ]] ; then
    echo "Error: Python version found is too old; needed at least 3.7, found $major.$minor"
    exit 1
  else
    touch install.log
    echo "Creating virtual environment using Python version found ($major.$minor) ..."
    python -m venv env >> install.log
    source env/bin/activate >> install.log
    echo "Upgrading pip ..."
    python -m pip install --upgrade pip >> install.log
    echo "Installing Python dependencies ..."
    pip install numpy scipy matplotlib >> install.log
    echo "Creating run file ..."
    touch run.sh >> install.log
    {
    echo $"#!/bin/bash"
    echo $"source $ADSORPYON_PATH/env/bin/activate"
    echo $"python $ADSORPYON_PATH/adsorpyon/main.py"
    echo $"deactivate"
    } >> run.sh
    chmod +x run.sh
    echo "Installation completed successfully!"
    exit 0
  fi
fi

