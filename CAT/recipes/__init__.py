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
    import nanoCAT

except ImportError:
    import warnings
    warnings.simplefilter("always", category=ImportWarning)
    warnings.warn("Usage of the 'CAT.recipes' module requires the Nano-CAT package: "
                  "'https://github.com/nlesc-nano/nano-CAT'", ImportWarning)
    del warnings

    __all__ = []

else:
    del nanoCAT

    from .bulk import bulk_workflow
    from .mark_surface import replace_surface
    from .surface_dissociation import dissociate_surface, row_accumulator

    __all__ = ['bulk_workflow', 'replace_surface', 'dissociate_surface', 'row_accumulator']
