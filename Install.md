# auto-qchem installation instructions

### Prerequisite for Windows users only
Install openbabel binaries (comes with a GUI)
https://github.com/openbabel/openbabel/releases/download/openbabel-2-4-1/OpenBabel-2.4.1.exe

### Installation with Miniconda (Windows, MacOS)

Install the miniconda with python=3.7 suitable for your OS: [miniconda installation webpage](https://docs.conda.io/en/latest/miniconda.html)

In anaconda prompt (Windows) or terminal (MacOS) create a conda environment executing. You can use a different environment name, in this instructions we will use ```autoqchem```
```bash
conda create --name autoqchem python=3.7
```
Activate the environment
```bash
conda activate autoqchem
```

Install mainstream python packages
```bash
conda install jupyter pandas scipy matplotlib pymongo pyyaml fabric xlrd appdirs
```
Answer "y" to the questions. If you'd like to automatically say yes to all question you can
add ```-y``` at the end of the command, e.g. ```conda install jupyter -y```

Install openbabel v2.4.1 (v3.0.0 is available, but crashes on many structures)
```bash
conda install -c conda-forge openbabel=2.4.1
```

Install imolecule (allows 3d rotateable renderings within jupyter), a specific version is required here
```bash
python -m pip install imolecule==0.1.13
```

Clone autoqchem from github, if you don't have git, please install git for your OS [git installation webpage](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git). Note that the last argument to this command is a local directory to which auto-qchem will be installed, please update it to a suitable one. For the purpose of this instructions I will use ```some_directory```
```bash
git clone https://github.com/PrincetonUniversity/auto-qchem.git some_directory
```
Navigate to the directory where auto-qchem github repository is installed (a file ```setup.py``` should be present in that directory) and run the setup script
```bash
cd some_directory/auto-qchem
python setup.py install
```

et voila!



