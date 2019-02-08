
.. image:: https://travis-ci.org/SCM-NV/qmflows.svg?branch=master
   :target: https://travis-ci.org/SCM-NV/qmflows 
.. image:: https://img.shields.io/badge/python-3.6-blue.svg

#################################
Compound Attachment/Analysis Tool
#################################

CAT is a collection of tools designed for the construction, and subsequent analysis, of various chemical compounds.
Further information is provided in the documentation_.

Installation
============

- Download miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_. 

- Create a new virtual environment using the following commands:

  - ``conda create -n CAT`` 

- The virtual environment can be enabled and disabled by typing:

  - Enable: ``conda activate CAT`` 
  
  - Disable: ``conda deactivate``
    
    
.. _dependecies:

Dependencies installation
-------------------------

Using the conda environment the following packages should be installed:    

- install rdkit_ using conda:

  - ``conda install -y --name CAT -c rdkit rdkit``
  

- install HDF5_ using conda:

  - ``conda install -y --name CAT -c anaconda h5py``
    
    
.. _installation:

Package installation
--------------------
Finally install the package:

- Update PLAMS_ using pip:
  - ``pip install git+https://github.com/SCM-NV/PLAMS@master#egg=plams-1.2``
    
- Install QMFlows_ using pip:
  - ``pip install git+https://github.com/SCM-NV/qmflows@master#egg=qmflows-0.3.0``
 
- Install **CAT** using pip:
  - ``pip install git+https://github.com/BvB93/CAT@master#egg=CAT-0.1.0``

Now you are ready to use **CAT**. 

Input files
============

An example input file, input.py, is located in the example_ directory.

.. _documentation: https://cat.readthedocs.io/en/latest/
.. _miniconda: http://conda.pydata.org/miniconda.html
.. _anaconda: https://www.continuum.io/downloads
.. _installConda: http://conda.pydata.org/docs/install/quick.html
.. _Noodles: http://nlesc.github.io/noodles/
.. _HDF5: http://www.h5py.org/ 
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _jupyter-notebook: http://jupyter.org/
.. _tutorial-qmflows: https://github.com/SCM-NV/qmflows/tree/master/jupyterNotebooks
.. _examples: https://github.com/SCM-NV/qmflows/tree/master/src/qmflows/examples
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _QMFlows: https://github.com/SCM-NV/qmflows
.. _example: https://github.com/BvB93/CAT/tree/master/CAT/examples
