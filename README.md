<p align="center">
  <img src="https://github.com/beef-broccoli/misc-files/blob/9332ec68f7f798a3c2819dad9a0d2280769985ee/autoqchem.png" alt="logo" width="600">
</p>

## Quick links

[Installation instructions](https://github.com/doyle-lab-ucla/auto-qchem/blob/master/INSTALL.md)

[Database user instructions](https://github.com/doyle-lab-ucla/auto-qchem/blob/master/DB.md)

[Code base documentation](https://doyle-lab-ucla.github.io/auto-qchem)

---
[Repo with example jupyter notebooks](https://github.com/doyle-lab-ucla/auto-qchem-notebook-examples)

[Database link](https://autoqchem.org)

[Auto-QChem paper (publisher)](https://pubs.rsc.org/en/content/articlelanding/2022/re/d2re00030j#!divCitation), [free pdf](https://drive.google.com/file/d/1M8Ydqlk5Kbc_8WoR5dAm_JIbf2IBJTlU/view?usp=share_link)

## Note to external user
Auto-QChem currently supports Slurm scheduler at Princeton via slurm_manager.py; and SGE/UGE-type scheduler at UCLA via sge_manager.py. 

If you are an external user and your computational cluster uses either slurm or sge scheduler, you just need to make some minor changes to either .py files to make sure you can log in with credentials at your institution.

If you are an external user and your computational cluster uses other schedulers, you can still modify either .py files to adapt to your cluster, but significant changes might be required and we unfortunatey won't be able to help without access to your cluster.
