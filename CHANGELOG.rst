###########
Change Log
###########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

0.4.5
*****

* All raw input scripts are now stored in the structures.hdf5 file.


0.4.4
*****

* Split CAT_database.py into database.py and database_functions.py.
* Unoptimized starting structures are now exported to the database.
* Added the sphinx autosummary extension.


0.4.3
*****

* Improved interaction between the database and BDE module.
* Cleaned up BDE module.
* HDF5 indices are now always sorted when itneraction with the database.


0.4.2
*****

* Numerous bug fixes.
* A couple of code-style changes.


0.4.1
*****

* COSMO-RS calculations now allow for COSMO-surface construction
  at the DFT level.


0.4.0
*****

* Introduction of the CAT.Database class.
* Central object of CAT has been changed into a dataframe of
  molecules rather than lists molecules.
* Updated a number of tests.


0.3.3
*****

* Changed qmflows template import syntax (see: https://github.com/SCM-NV/qmflows/pull/132).
* Changed yaml loader.


0.3.2
*****

* Further (minor) updates and bug fixes to the database interaction.
* Overhaul of the bond dissociation energy (BDE) module.
* Job settings are now stored in the database.


0.3.0
*****

* Massive overhaul of the CAT database interaction.
* Moved functions related to functiona group recognizition to
  CAT.attachment.ligand_anchoring.py.
* Multiple minor bug fixes.


[Unreleased]
************

Added
-----

* Empty Python project directory structure.
