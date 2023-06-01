# auto-qchem installation instructions

### New installation instructions for auto-qchem v1.2.6
*v1.2.6 runs in python 3.7 environment, which seems to be the more stable version.*
*You should use this version, unless your mac use Apple M chips. In that case, please install auto-qchem v1.3.0 (see below)* 
1. Install conda: https://www.anaconda.com. I find it easier to just download the installer and install like how you would with any software.
2. Conda allows you to create separate software development environments for different purposes. We will create an environment just for auto-qchem. Open terminal and run this, and answer "y" to the prompt:
```
conda create --name autoqchem python=3.7 ipython
```
3. With conda environments, you need to activate them before you use. Run this in terminal: 
 ```
 conda activate autoqchem 
 ```
4. We have some separate packages to install first.
(These packages are a little finicky, to ensure expected behavior it's best to install them as specified here). Run the following commands in terminal one by one, and answer "y": 
 ```
conda install -c conda-forge openbabel=2.4.1
conda install -c conda-forge py3dmol
pip install rdkit-pypi
 ```
5. Finally, we can install auto-qchem directly through pip 
(Previously, we had to build from sources, but auto-qchem is now a proper python package. We also add a dash for the package name, just to make it clearer)
 ```
pip install auto-qchem==1.2.6
 ```
6. That's it, you are ready to use auto-qchem. 

---

### New installation instructions for auto-qchem v1.3.0
*This version runs on python 3.8, which is the oldest python version Apple M chips support.*
*This version also avoids vulnerability issue in numpy<=1.21. However, conformer generation might break more easily since we have to run openbabel 3.1*
*This version has not been extensively tested.*

1. Install conda: https://www.anaconda.com. I find it easier to just download the installer and install like how you would with any software.
2. Conda allows you to create separate software development environments for different purposes. We will create an environment just for auto-qchem. Open terminal and run this, and answer "y" to the prompt:
```
conda create --name autoqchem python=3.8 ipython
```
3. With conda environments, you need to activate them before you use. Run this in terminal: 
 ```
 conda activate autoqchem 
 ```
4. We have some separate packages to install first.
(These packages are a little finicky, to ensure expected behavior it's best to install them as specified here). Run the following commands in terminal one by one, and answer "y": 
 ```
conda install -c conda-forge openbabel
 ```
5. Finally, we can install auto-qchem directly through pip 
(Previously, we had to build from sources, but auto-qchem is now a proper python package. We also add a dash for the package name to differ it from environment name)
 ```
pip install auto-qchem==1.3.0
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

**Current versions: v1.2.6, v1.3.0**

**(Optional)** if you are NOT in the autoqchem environment, activate first:
 ```
conda activate autoqchem
 ```

**(Optional)** you can check the current version of auto-qchem:
 ```
conda list auto-qchem
 ```

If you are running auto-qchem v1.2.x and want to keep running the more stable v1.2.x, run this to update:
 ```
pip install auto-qchem==1.2.6 --upgrade
 ```

If you are running auto-qchem v1.3.x, or want to run the newer v1.3.x,  run this to update:
 ```
pip install auto-qchem==1.3.0 --upgrade
 ```

---

---

---


### Prerequisite for Windows users only
Install openbabel binaries (comes with a GUI)
https://github.com/openbabel/openbabel/releases/download/openbabel-2-4-1/OpenBabel-2.4.1.exe


