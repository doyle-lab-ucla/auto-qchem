# auto-qchem installation instructions

### For Windows users only
install openbabel binary package, the python openbabel package for windows only includes bindings to binaries
https://github.com/openbabel/openbabel/releases/download/openbabel-2-4-1/OpenBabel-2.4.1.exe

### Installation with miniconda (Windows, MacOS)

Install the miniconda with python=3.7 suitable for your OS: [miniconda installation webpage](https://docs.conda.io/en/latest/miniconda.html)

In anaconda prompt (Windows) or terminal (MacOS) create a conda environment executing
```bash
conda create --name autoqchem python=3.7
```
activate the environment
```bash
conda activate autoqchem
```

install mainstream python packages
```bash
conda install jupyter pandas scipy matplotlib pymongo pyyaml fabric xlrd appdirs
```

install openbabel v2.4.1 (v3.0.0 is available, but crashes on many structures)
```bash
conda install -c conda-forge openbabel=2.4.1
```

install imolecule (allows 3d rotateable renderings within jupyter), a specific version is required here
```bash
python -m pip install imolecule==0.1.13
```

clone autoqchem from github, if you don't have git, please install git for your OS [git installation webpage](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git). Note that the last argument to this command is a local directory to which auto-qchem will be installed, please update it to a suitable one. For the purpose of this instructions I will use ```some_directory```
```bash
git clone https://github.com/PrincetonUniversity/auto-qchem.git some_directory
```
navigate to the directory where auto-qchem github repository is installed (a file ```setup.py``` should be present in that directory) and run the setup script
```bash
cd some_directory
python setup.py install
```

et voila!

