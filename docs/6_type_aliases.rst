.. _Type Aliases:

Type Aliases
============

Aliases are available for a large number of job types,
allowing one to pass a :class:`str` instead of a :class:`type` object, thus simplifying
the input settings for **CAT**. Aliases are insensitive towards capitalization
(or lack thereof).

A comprehensive list of :class:`plams.Job<scm.plams.core.basejob.Job>` subclasses and their respective
aliases (*i.e.* :class:`str`) is presented below.

Aliases
~~~~~~~

-   |ADFJob| = ``"adf"`` = ``"adfjob"``

-   |AMSJob| = ``"ams"`` = ``"amsjob"``

-   |UFFJob| = ``"uff"`` = ``"uffjob"``

-   |BANDJob| = ``"band"`` = ``"bandjob"``

-   |DFTBJob| = ``"dftb"`` = ``"dftbjob"``

-   |MOPACJob| = ``"mopac"`` = ``"mopacjob"``

-   |ReaxFFJob| = ``"reaxff"`` = ``"reaxffjob"``

-   |Cp2kJob| = ``"cp2k"`` = ``"cp2kjob"``

-   |ORCAJob| = ``"orca"`` = ``"orcajob"``

-   |DiracJob| = ``"dirac"`` = ``"diracjob"``

-   |GamessJob| = ``"gamess"`` = ``"gamessjob"``

-   |DFTBPlusJob| = ``"dftbplus"`` = ``"dftbplusjob"``

-   |CRSJob| = ``"crs"`` = ``"cosmo-rs"`` = ``"crsjob"``


.. |ADFJob| replace:: :class:`ADFJob<scm.plams.interfaces.adfsuite.adf.ADFJob>`
.. |AMSJob| replace:: :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`
.. |UFFJob| replace:: :class:`UFFJob<scm.plams.interfaces.adfsuite.uff.UFFJob>`
.. |BANDJob| replace:: :class:`BANDJob<scm.plams.interfaces.adfsuite.band.BANDJob>`
.. |DFTBJob| replace:: :class:`DFTBJob<scm.plams.interfaces.adfsuite.dftb.DFTBJob>`
.. |MOPACJob| replace:: :class:`MOPACJob<scm.plams.interfaces.adfsuite.mopac.MOPACJob>`
.. |ReaxFFJob| replace:: :class:`ReaxFFJob<scm.plams.interfaces.adfsuite.reaxff.ReaxFFJob>`
.. |Cp2kJob| replace:: :class:`Cp2kJob<scm.plams.interfaces.thirdparty.cp2k.Cp2kJob>`
.. |ORCAJob| replace:: :class:`ORCAJob<scm.plams.interfaces.thirdparty.orca.ORCAJob>`
.. |DiracJob| replace:: :class:`DiracJob<scm.plams.interfaces.thirdparty.dirac.DiracJob>`
.. |GamessJob| replace:: :class:`GamessJob<scm.plams.interfaces.thirdparty.gamess.GamessJob>`
.. |DFTBPlusJob| replace:: :class:`DFTBPlusJob<scm.plams.interfaces.thirdparty.dftbplus.DFTBPlusJob>`
.. |CRSJob| replace:: :class:`CRSJob<scm.plams.interfaces.adfsuite.crs.CRSJob>`
