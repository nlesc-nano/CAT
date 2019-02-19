
.. image:: https://travis-ci.org/BvB93/CAT.svg?branch=master
   :target: https://travis-ci.org/BvB93/CAT
.. image:: https://readthedocs.org/projects/cat/badge/?version=latest
   :target: https://cat.readthedocs.io/en/latest
.. image:: https://img.shields.io/badge/python-3.6-blue.svg
   :target: https://www.python.org

#################################
Compound Attachment/Analysis Tool
#################################

**CAT** is a collection of tools designed for the construction, and subsequent analysis, of various chemical compounds.
Further information is provided in the documentation_.

Installation
============

- Download miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_. 

- Create a new virtual environment using the following commands:

  - ``conda create -n CAT`` 

- The virtual environment can be enabled and disabled by, respectively, typing:

  - Enable: ``conda activate CAT`` 
  
  - Disable: ``conda deactivate``
    
    
.. _dependecies:

Dependencies installation
-------------------------

Using the conda environment the following packages should be installed:    

- rdkit_: ``conda install -y --name CAT -c rdkit rdkit``

- HDF5_: ``conda install -y --name CAT -c anaconda h5py``
    
    
.. _installation:

Package installation
--------------------
Finally, install **CAT** using pip:
   
- **CAT**: ``pip install git+https://github.com/BvB93/CAT@master#egg=CAT-0.2.0``

Now you are ready to use **CAT**. 

Input files
============

An example input file, input.py, is located in the example_ directory.

.. _documentation: https://cat.readthedocs.io/en/latest/
.. _miniconda: http://conda.pydata.org/miniconda.html
.. _anaconda: https://www.continuum.io/downloads
.. _installConda: https://docs.anaconda.com/anaconda/install/
.. _HDF5: http://www.h5py.org/ 
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _QMFlows: https://github.com/SCM-NV/qmflows
.. _example: https://github.com/BvB93/CAT/tree/master/CAT/examples
