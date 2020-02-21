"""
CAT.recipes
===========

A number of recipes constructed using the CAT and Nano-CAT packages.

Examples
--------
.. code:: python

    >>> from CAT.recipes import bulk_workflow
    >>> from CAT.recipes import replace_surface
    >>> from CAT.recipes import dissociate_surface, row_accumulator
    ...

"""

try:
    from nanoCAT import recipes as _recipes
    from nanoCAT.recipes import *

    __all__ = _recipes.__all__
    del _recipes

except ImportError as ex:
    __all__ = []
    raise ImportError("Usage of the 'CAT.recipes' module requires the Nano-CAT package: "
                      "'https://github.com/nlesc-nano/nano-CAT'") from ex
