# Introduction
Aroma is a utility package for evaluating aromatic properties via NICS
methods, preferably by NICS-pi,zz. It is designed as a “Plug-In” utility for the
computational chemistry packages; Gaussian 09, Gaussian 16, Orca 5 and
Orca 6. NBO 6 and NBO 7 are also supported for the purpose of CMO-NICS
calculations (using the NCS procedure). 


This package built in collaboration and active guidance of [Prof Dr. Amnon Stanger](https://chemistry.technion.ac.il/team/%d7%90%d7%9e%d7%a0%d7%95%d7%9f-%d7%a9%d7%98%d7%a0%d7%92%d7%a8/) at Technion and Scientific Interactive team. 

(c) Technion 

(c) Scientific Interactive


# Building 

## Requirments 
Python >= 3.9 

## Create a venv
python -m venv venv 
source venv/bin/activate

## Install Deps
pip install -r requirements.txt 

## Create Windows Installer
(use git bash)
bash build_win.sh

## Create Linux Installer
bash build_linux.sh

## Create MacOS Insatller
bash build_macos.sh

