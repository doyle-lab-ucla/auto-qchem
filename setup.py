from setuptools import setup

setup(
    name='auto-qchem',
    version='0.1',
    packages=['autoqchem'],
    data_files=['config.yml'],
    url='https://github.com/PrincetonUniversity/auto-qchem',
    license='GPL',
    author='Andrzej Zuranski, Benjamin Shields',
    author_email='zuranski@princeton.edu, bjs4@princeton.edu',
    description='auto-qchem',
    install_requires=['numpy',
                      'pandas',
                      'pyyaml',
                      'scipy',
                      'fabric',
                      'paramiko',
                      'pymongo',
                      'appdirs']
)
