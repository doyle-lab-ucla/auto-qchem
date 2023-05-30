from setuptools import setup

setup(
    name='auto-qchem',
    version='1.2.0',
    packages=['autoqchem'],
    data_files=['config.yml'],
    url='https://github.com/PrincetonUniversity/auto-qchem',
    license='GPL',
    author='Andrzej Zuranski, Benjamin Shields',
    author_email='zuranski@princeton.edu, bjs4@princeton.edu',
    description='auto-qchem',
    install_requires=['numpy==1.21.6',
                      'pandas==1.3.5',
                      'pyyaml==6.0',
                      'scipy==1.7.3',
                      'fabric==2.6.0',
                      'paramiko==3.1.0',
                      'pymongo==4.3.2',
                      'appdirs==1.4.4',
                      # 'rdkit',
                      'ipywidgets==8.0.6',
                      'py3Dmol==1.8.1'
                      ]
)
