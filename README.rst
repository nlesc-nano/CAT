################################################################################
Compound Attachment/Analysis Tool
################################################################################

CAT is a collection of tools designed for the construction, and subsequent analysis, of various chemical compounds.

Installation
============

- Download miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_. 

- Create a new virtual environment using the following commands:

  - ``conda create -n CAT`` 

- Activate the new virtual environment
  
  - ``source activate CAT``

To exit the virtual environment type  ``source deactivate``.
    
    
.. _dependecies:

Dependencies installation
-------------------------

- Type in your terminal:

  ``conda activate qmflows``  

Using the conda environment the following packages should be installed:    


- install rdkit_ using conda:

  - ``conda install -y -q --name qmflows -c rdkit rdkit``

- install HDF5_ using conda:

  - ``conda install -y -q --name qmflows -c anaconda h5py``
    
    
.. _installation:

Package installation
--------------------
Finally install the package:

- Update PLAMS_ using pip:
  - ``pip install git+https://github.com/SCM-NV/PLAMS@master#egg=plams-1.2``
    
- Install **CAT** using pip:
  - ``pip install git+https://github.com/SCM-NV/qmflows@master#egg=qmflows-0.3.0``

Now you are ready to use **CAT**.  
