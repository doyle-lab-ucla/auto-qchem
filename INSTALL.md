# auto-qchem installation instructions

### New installation instructions for auto-qchem

1. Install conda: https://www.anaconda.com. I find it easier to just download the installer and install like how you would with any software.
2. Conda allows you to create separate software development environments for different purposes. We will create an environment just for auto-qchem. Open terminal and run this, and answer "y" to the prompt:
```
conda create --name autoqchem python=3.8 ipython
```
3. With conda environments, you need to activate them before you use. Run this in terminal: 
 ```
 conda activate autoqchem 
 ```
4. We have some separate packages to install first. (OpenBabel does not work on pip). Run the following commands in terminal one by one, and answer "y": 
 ```
conda install -c conda-forge openbabel
 ```
5. Finally, we can install auto-qchem directly through pip 
 ```
pip install auto-qchem
 ```
6. That's it, you are ready to use auto-qchem. 

---

### If you want to start over...
Sometimes things don't work, and you might want to delete the autoqchem environment and start over with a clean install. To remove the conda environment and all packages associated with it, run:

**(Optional)** if you are in the autoqchem environment, deactivate first:
 ```
conda deactivate
 ```
then run this to remove autoqchem environment:
 ```
conda remove -n autoqchem --all
 ```
And you are ready to start over at step 2 in the new installation instructions. 

---

### If you haven't used auto-qchem for a while...
Chances are there have been updates for the software.

**(Optional)** if you are NOT in the autoqchem environment, activate first:
 ```
conda activate autoqchem
 ```

**(Optional)** you can check the current version of auto-qchem:
 ```
conda list auto-qchem
 ```

Run this to update to the newest version:
 ```
pip install auto-qchem --upgrade
 ```

---

---

---


### Prerequisite for Windows users only
I have not tested for Windows. 

Install openbabel binaries (comes with a GUI)
https://github.com/openbabel/openbabel/releases/download/openbabel-2-4-1/OpenBabel-2.4.1.exe


