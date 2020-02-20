"""
CAT.recipes
===========

A number of recipes using the CAT and Nano-CAT packages.

Examples
--------
.. code:: python

    >>> from CAT.recipes import bulk_workflow
    >>> from CAT.recipes import replace_surface
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

    __all__ = ['bulk_workflow', 'replace_surface']
