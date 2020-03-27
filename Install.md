# auto-qchem installation instructions

### Prerequisite for Windows users only
Install openbabel binaries (comes with a GUI)
https://github.com/openbabel/openbabel/releases/download/openbabel-2-4-1/OpenBabel-2.4.1.exe

### Installation with Miniconda (Windows, MacOS)

Install the miniconda with python=3.7 suitable for your OS: [miniconda installation webpage](https://docs.conda.io/en/latest/miniconda.html)

In anaconda prompt (Windows) or terminal (MacOS) create a conda environment. You can use a different environment name, in this instructions we will use ```autoqchem```.
For Windows 'Anaconda Powershell Prompt' has more functionality than 'Anaconda Prompt', e.g. tab-completion, but either
one works fine.
```bash
conda create --name autoqchem python=3.7
```
Activate the environment
```bash
conda activate autoqchem
```

Install mainstream python packages
```bash
conda install jupyter pandas scipy matplotlib pymongo pyyaml fabric xlrd appdirs openpyxl
```
Answer "y" to the questions. If you'd like to automatically say yes to all question you can
add ```-y``` at the end of the command, e.g. ```conda install jupyter -y```

Install openbabel v2.4.1 (v3.0.0 is available, but crashes on many structures) and rdkit
```bash
conda install -c conda-forge openbabel=2.4.1
conda install -c rdkit rdkit
```

Install imolecule (allows 3d rotateable renderings within jupyter), a specific version is required here
```bash
python -m pip install imolecule==0.1.13
```

Clone autoqchem from github. If you don't have ```git```, please install ```git``` for your OS [git installation webpage](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git). Note that
 the last argument to the command below is a local directory to which auto-qchem will be installed, please update it to 
 a suitable one. It can be anywhere, but it needs to either be empty or not yet existing. For the purpose of this instructions I will use ```some_directory```.
  A good example would be ```~/software/github/auto-qchem```.
```bash
git clone https://github.com/PrincetonUniversity/auto-qchem.git some_directory
```
Navigate to the directory where auto-qchem github repository is installed (a file ```setup.py``` should be present in that directory) and run the setup script
```bash
cd some_directory
python setup.py install
```

That's it!

### Verify the installation
Verify the installation by running a special template notebook ```framework_functionality_test.ipynb```. It doesn't do a whole lot, but 
runs through some key functions.

Navigate to notebooks directory
```bash
cd some_directory/notebooks
```
Fire up jupyter-notebook
```bash
jupyter-notebook framework_functionality_test.ipynb
```
Run the notebook!


