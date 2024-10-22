from setuptools import setup

setup(
    name='auto-qchem',
    version='1.3.15',
    packages=['autoqchem'],
    data_files=['config.yml'],
    url='https://github.com/doyle-lab-ucla/auto-qchem',
    license='GPL',
    author='Andrzej Zuranski, Benjamin Shields, Jason Wang, Winston Gee',
    description='auto-qchem',
    long_description='automated dft calculation management software',
    install_requires=['numpy>=1.22',
                      'pandas>=1.3',
                      'pyyaml>=6.0',
                      'scipy>=1.7',
                      'fabric>=2.6',
                      'paramiko>=3.1',
                      'pymongo>=3.10',
                      'appdirs>=1.4',
                      'ipywidgets>=8.0',
                      'py3Dmol>=1.8',
                      'jupyterlab>=3.4',
                      'notebook>=6.4',
                      'ipywidgets>=8.0',
                      'xlrd>=2.0',
                      'openpyxl>=3.0',
                      'rdkit',
                      'tqdm',
                      'matplotlib>=3.5'
                      ],
    python_requires='>=3.8'
)
