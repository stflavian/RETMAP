$ADSORPYON_PATH=$pwd

$python_version = python --version
$python = $python_version.Split(' ')[0]
$version = $python_version.Split(' ')[1]

Write-Output $python $version

if ($python -ne "Python") {
  Write-Output "Error: Python was not found on this system!"
  Exit
} else {
  Write-Output "Python was found on the system ..."
  $major = [int]$version.Split('.')[0]
  $minor = [int]$version.Split('.')[1]
  if ($major -ne 3 -or $minor -lt 7) {
    Write-Output "Error: Python version found is too old; needed at least 3.7, found $major.$minor"
    Exit
  } else {
      Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
      New-Item -Name "install.log" -ItemType File
      Write-Output "Creating virtual environment using Python version found ($major.$minor) ..."
      python -m venv env | Out-File -Append "install.log"
      .\env\Scripts\Activate.ps1 | Out-File -Append "install.log"
      Write-Output "Upgrading pip ..."
      python -m pip install --upgrade pip | Out-File -Append "install.log"
      Write-Output "Installing Python dependencies ..."
      pip install -r requirements.txt | Out-File -Append "install.log"
      Write-Output "Creating run file ..."
      New-Item -Name "run.ps1" -ItemType File | Out-File -Append "install.log"

      Write-Output ". $ADSORPYON_PATH\env\Scripts\Activate.ps1" | Out-File -Append "run.ps1"
      Write-Output "python $ADSORPYON_PATH\adsorpyon\main.py" | Out-File -Append "run.ps1"
      Write-Output "deactivate" | Out-File -Append "run.ps1"

      Write-Output "Installation completed successfully!"
      deactivate
  }
}
