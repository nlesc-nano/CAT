##########
Change Log
##########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

0.5.1
*****

* Documentation update.


0.5.0
*****

* CAT_ has been split into 3 seperate packages (see https://github.com/nlesc-nano/CAT/issues/39):

  * CAT_: A collection of tools designed for the automatic construction of composite chemical compounds.
  * nano-CAT_: A collection of tools for the analysis of nanocrystals.
  * data-CAT_: A databasing framework for the Compound Attachment Tools package (CAT_).

* Docstrings have been changed into NumPy style.
* Added typehints.
* Added the CAT.SettingsDataFrame and CAT.SettingsSeries classes.
* Added more tests.
* Cleaned up all input-parsing related modules.
* Custom function groups (*i.e.* SMILES_ strings) can now be specified in the input
  under the optional.ligand.functional_groups key (see https://github.com/nlesc-nano/CAT/issues/13).

.. _CAT: https://github.com/nlesc-nano/CAT/
.. _nano-CAT: https://github.com/nlesc-nano/nano-CAT/
.. _data-CAT: https://github.com/nlesc-nano/data-CAT/
.. _SMILES: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system


0.4.6
*****

* Added an interface between MongoDB_ and the CAT.Database_ class (see https://github.com/nlesc-nano/CAT/issues/11).

.. _MongoDB: https://www.mongodb.com/
.. _CAT.Database: https://cat.readthedocs.io/en/latest/7_database.html


0.4.5
*****

* All raw input scripts are now stored in the structures.hdf5 file
  (see: https://github.com/nlesc-nano/CAT/issues/36).


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
