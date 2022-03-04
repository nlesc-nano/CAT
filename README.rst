.. image:: https://github.com/nlesc-nano/CAT/workflows/Build%20with%20Conda/badge.svg
   :target: https://github.com/nlesc-nano/CAT/actions?query=workflow%3A%22Build+with+Conda%22
.. image:: https://readthedocs.org/projects/cat/badge/?version=latest
   :target: https://cat.readthedocs.io/en/latest/
.. image:: https://zenodo.org/badge/169708827.svg
   :target: https://zenodo.org/badge/latestdoi/169708827
.. image:: https://badge.fury.io/py/nlesc-CAT.svg
   :target: https://badge.fury.io/py/nlesc-CAT

|

.. image:: https://img.shields.io/badge/python-3.6-blue.svg
   :target: https://docs.python.org/3.6/
.. image:: https://img.shields.io/badge/python-3.7-blue.svg
   :target: https://docs.python.org/3.7/
.. image:: https://img.shields.io/badge/python-3.8-blue.svg
   :target: https://docs.python.org/3.8/
.. image:: https://img.shields.io/badge/python-3.9-blue.svg
   :target: https://docs.python.org/3.9/
.. image:: https://img.shields.io/badge/python-3.10-blue.svg
   :target: https://docs.python.org/3.10/

###############################
Compound Attachment Tool 0.11.0
###############################

**CAT** is a collection of tools designed for the construction of various chemical compounds.
Further information is provided in the documentation_.

Package installation
--------------------
**CAT** can be installed via pip as following:

- **CAT**: ``pip install nlesc-CAT --upgrade``

Note that, while not strictly necessary, it is recommended to first create a conda environment:

- Download and install miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Create a new virtual environment:  ``conda create --name CAT python``

- Activate the environment:: ``conda activate CAT``

Input files
============

Running **CAT** and can be done with the following command:
``init_cat my_settings.yaml``. The user merely has to provide a yaml_ file
with the job settings, settings which can be tweaked and altered to suit ones
purposes (see example1_). Alternatively, **CAT** can be run like a regular
python script, bypassing the command-line interface
(*i.e.* ``python input.py``, see example2_).

An extensive description of the various available settings is available in
the documentation_.


.. _yaml: https://yaml.org/
.. _documentation: https://cat.readthedocs.io/en/latest/
.. _miniconda: http://conda.pydata.org/miniconda.html
.. _anaconda: https://www.continuum.io/downloads
.. _installConda: https://docs.anaconda.com/anaconda/install/
.. _HDF5: http://www.h5py.org/
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _QMFlows: https://github.com/SCM-NV/qmflows
.. _example1: https://github.com/BvB93/CAT/blob/master/examples/input_settings.yaml
.. _example2: https://github.com/BvB93/CAT/blob/master/examples/input.py
