# Auto-QChem Database User Guide

### 1. Query available molecules

Navigate to landing page [link](http://128.112.127.147)

![landing page](images/landing_page.png)

#### 1.1 Query Form

##### a) Fields
Query form has 2 fields, both optional:
* __Select tags (multiple choice)__ - tags are used to mark a collection of molecules,
if you select multiple tags, molecules for all tags will be displayed. If left blank
all molecules will be displayed
* __SMARTS substructure__ - queries the structures using the SMARTS query (SMILES strings
are a subset of SMARTS), quick reference to the SMARTS query language can be looked up 
here: [https://www.daylight.com/dayhtml_tutorials/languages/smarts/index.html](https://www.daylight.com/dayhtml_tutorials/languages/smarts/index.html) 

##### b) Buttons
There are two buttons ```Query``` and ```Export```.
 * __Query__ - queries the DB and displays the table of queried molecules
 * __Export__ - shall be used after the ```Query```, it will download the displayed 
 table as an ```.xlsx``` file.
 
 Below a result of an example query on a datasets with 1166 molecules, with SMARTS query for anhydrides.
 ![query result](images/query_result.png)
 
 
 #### 1.2 Descriptors Lookup
 
 For each entry in the table a link to a descriptors lookup called ```descriptors``` is available in the rightmost 
 column. It will display the the qchem descriptors for that particule molecule. If the molecule contains multiple
 conformations, a "Boltzmann" average of all descriptors is shown.
 
 ![descriptors lookup](images/descriptors_lookup.png)
 
 ### 2. Descriptors extraction
 
 Once molecules have been queried, their descriptors can be extracted into an ```.xlsx``` file by toggling the
 ```Download descriptors``` bar and filling the form.
 
 ![download form](images/download_form.png)
 
 #### 2.1 Download Form

##### a) Fields
All fields are required

* __Descriptor Presets (multiple choice)__ - the following presets are available, choose as many as needed:
    * __Global__ - molecule level descriptors, e.g. homo energy, dipole moment, molecular weight, etc.
    * __MinMaxAtomic__ - atomic level descriptors minimum and maximum over the atoms within the molecule, e.g. 
    buried volume, Mulliken charge, NMR shift, etc.
    * __SubstructureAtomic__ - atomic level descriptors for each of the atoms in the SMARTS substrcture used for
     filtering
    * __SubstructureCore__ - atomic level descriptors for the common core of atoms within the dataset, the common
    core is determined using the MCS procedure from rdkit 
    [reference](https://www.rdkit.org/docs/source/rdkit.Chem.MCS.html)
    If substructure has been used for filtering, the common core will include the substructure
    * __Transitions__ - top 10 excited state transitions ordered by their oscillation strength
* __Conformer option (single choice)__ - choose from the following options:
    * __Boltzmann__ - Boltzmann average
    * __Max__ - lowest energy conformer (maximum weight conformer)
    * __Min__ - highest energy conformer (minimum weight conformer)
    * __Mean__ - arithmetic average
    * __Std__ - standard deviation over the conformers
    * __Any__ - randomily chosen conformer

##### b) Buttons
    
* __Download__ - download the descriptors to an ```.xlsx``` file. Note: when extracting descriptors for many molecules
this operation can take up to few minutes, depending on the server load.
     
