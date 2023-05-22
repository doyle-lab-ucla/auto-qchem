
<img src="https://github.com/beef-broccoli/misc-files/blob/9332ec68f7f798a3c2819dad9a0d2280769985ee/autoqchem.png" alt="logo" width="600" align="center"/>

### Link to publication

[Auto-QChem is published!](https://pubs.rsc.org/en/content/articlelanding/2022/re/d2re00030j#!divCitation). Free pdf is accessible [here](https://drive.google.com/file/d/1M8Ydqlk5Kbc_8WoR5dAm_JIbf2IBJTlU/view?usp=share_link)

### Note to external user
Auto-QChem currently supports Slurm scheduler at Princeton via slurm_manager.py; and SGE/UGE-type scheduler at UCLA via sge_manager.py. 

If you are an external user and your computational cluster uses either slurm or sge scheduler, you just need to make some minor changes to either .py files to make sure you can log in with credentials at your institution.

If you are an external user and your computational cluster uses other schedulers, you can still modify either .py files to adapt to your cluster, but significant changes might be required and we unfortunatey won't be able to help without access to your cluster. 

### Quick links

[Installation instructions](https://github.com/doyle-lab-ucla/auto-qchem/blob/master/Install.md)

[DB interface user guide](https://github.com/doyle-lab-ucla/auto-qchem/blob/master/DB.md)

[Descriptor database](https://autoqchem.org)

[Functional documentation](https://doyle-lab-ucla.github.io/auto-qchem)

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
