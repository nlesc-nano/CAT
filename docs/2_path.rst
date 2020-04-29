.. _Path:

path
====

Default Settings
~~~~~~~~~~~~~~~~

.. code::

    path: null

|

Arguments
~~~~~~~~~

.. attribute:: path

    :Parameter:     * **Type** - :class:`str` or :class:`NoneType`
                    * **Default value** – ``None``

    The path were all working directories are/will be stored.
    To use the current working directory, use one of the following values:
    ``None``, ``"."``, ``""`` or ``"path_to_workdir"``.

    .. note::
        The yaml format uses ``null`` rather than ``None`` as in Python.
