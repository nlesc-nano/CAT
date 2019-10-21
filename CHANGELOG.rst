##########
Change Log
##########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

0.6.5
*****
* Updated Nano-CAT to 0.2.4: https://github.com/nlesc-nano/nano-CAT/pull/20.
* Updated Data-CAT to 0.1.5: https://github.com/nlesc-nano/data-CAT/pull/17.
* Import assertions from AssertionLib_ rather than CAT_.
* Simplified to ``AsArray()`` context manager.
* Added the ``["keep_files"]`` option for quantum dot optimizations.
* Removed ``CRSJob()`` and ``CRSResults()``; import them from PLAMS_ instead.
* WiP: Introduction of a new CAT template system (``WorkFlow()``).


0.6.4
*****
* Moved the ligand bulkiness workflow from the `ligand` to the `qd` block in the CAT input. See `nano-CAT`_ 0.2.3.
* Updated the formula for the ligand bulkiness calculation. See `nano-CAT`_ 0.2.3.


0.6.3
*****
* Fixed a bug where hypervalent atoms where assigned incorrect atomic charges.


0.6.2
*****
* Added multiple improvements (and bug fixes) to the ligand conformation optimizer.
* Added a context manager for the `plams.Molecule.as_array()` method.
* Added an optimizer for the ligand vector.
* Updated the ligand bulkiness workflow in `nano-CAT`_ 0.2.2.


0.6.1
*****
* Added a workflow for calculating ligand bulkiness in `nano-CAT`_ 0.2.1.


0.6.0
*****
* Implemented an interface to MATCH_ (Multipurpose Atom-Typer for CHARMM) in Nano-CAT.
* Added a workflow for creating CP2K input files with the MATCH-assigned atom types & charges.
* Updated the handling of assertions, see ``CAT.assertions.assertion_manager``.


0.5.5
*****
* Lowered Python version requirement from >=3.7 to >=3.6.


0.5.4
*****
* Minor updates to the logger.
* Cleaned up CAT.jobs.py.
* ``check_sys_var()`` is now only called if an ADF-specific Job is requirest.
* Job hashes are now stored in (and retrieved from) $JN.hash files (plain text).
* Added a permanent Database_ instance to .optional.database.db.
* Parsing of functional group SMILES_ strings is now carried out during the Schema_ validation.
* Updated Data-CAT_ to 0.1.2; changed status from pre-alpha to alpha
  (see https://github.com/nlesc-nano/data-CAT/pull/13).



0.5.3
*****
* Moved Molecule to file exporting (*i.e.* .xyz and .pdb creation) from data-CAT_ to CAT_.
* Molecules can now be exported to .mol and .mol2 formats (in addition to .pdb and .xyz format).
* Increased the clarity of many exceptions (see https://github.com/nlesc-nano/CAT/issues/45).
* Updated the documentation.
* Introduced a proper logger (see https://github.com/nlesc-nano/CAT/issues/46).
* Updated data-CAT_ to 0.1.1 (https://github.com/nlesc-nano/data-CAT/pull/12) and
  nano_CAT_ to 0.1.2 (https://github.com/nlesc-nano/nano-CAT/pull/10).


0.5.2
*****
* Added more tests.
* Added a more explicit error message to ``_smiles_to_rdmol()``.


0.5.1
*****
* Documentation update.
* Updated to the ligand dissociation module in nano-CAT_ (see https://github.com/nlesc-nano/nano-CAT/issues/1).
* Added the ``keep_files`` keyword to the cosmo-rs and ligand dissociation workflows.
  Default value: ``True``.
* See https://github.com/nlesc-nano/nano-CAT/pull/9.


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


0.4.6
*****
* Added an interface between MongoDB_ and the CAT.Database_ class (see https://github.com/nlesc-nano/CAT/issues/11).


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
* Empty Python project directory structure.


.. _AssertionLib: https://github.com/nlesc-nano/AssertionLib
.. _CAT: https://github.com/nlesc-nano/CAT
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _MATCH: http://brooks.chem.lsa.umich.edu/index.php?page=match&subdir=articles/resources/software
.. _CP2K: https://www.cp2k.org/
.. _Database: https://cat.readthedocs.io/en/latest/7_database.html#class-api
.. _Schema: https://github.com/keleshev/schema
.. _nano-CAT: https://github.com/nlesc-nano/nano-CAT/
.. _data-CAT: https://github.com/nlesc-nano/data-CAT/
.. _SMILES: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
.. _MongoDB: https://www.mongodb.com/
.. _CAT.Database: https://cat.readthedocs.io/en/latest/7_database.html
