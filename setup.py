from setuptools import setup

setup(
    name='auto-qchem',
    version='0.1',
    packages=['autoqchem'],
    data_files=['config.yml'],
    url='https://github.com/PrincetonUniversity/auto-qchem',
    license='GPL',
    author='AndrzejZuranski',
    author_email='zuranski@princeton.edu',
    description='auto-qchem',
    install_requires=['numpy',
                      'pandas',
                      'scipy',
                      'yaml',
                      'fabric',
                      'paramiko',
                      'openbabel',
                      'IPython',
                      'pymongo']
)
