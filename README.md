# Welcome to Auto-QChem!

### Quick links

Installation instructions https://github.com/PrincetonUniversity/auto-qchem/blob/master/Install.md

DB interface user guide https://github.com/PrincetonUniversity/auto-qchem/blob/master/DB.md

Functional documentation https://princetonuniversity.github.io/auto-qchem

### Update your version

Open your terminal (bash or Anaconda prompt) and activate your python environment 

```conda activate autoqchem```

Navigate to the source code of your repository (top level directory where auto-qchem is installed). This 
directory shall contain a ```setup.py``` file. To check if it's there execute ```ls```.

```
cd your_auto_qchem_directory
```

Update your code from github and re-install the package

```
git pull
python setup.py install
```

### Run notebooks

Template notebooks are stored in the auto-qchem repository under in ```notebooks``` directory

```
cd your_auto_qchem_directory
cd notebooks
```

Start a jupyter notebook 

```
jupyter-notebook
```

A new tab will open in your system web browser from which you can run the notebooks. An excellent documentation 
on jupyter notebooks and how to run them exists
 [here](https://jupyter-notebook.readthedocs.io/en/stable/examples/Notebook/Running%20Code.html).
