#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# qmflows documentation build configuration file, created by
# sphinx-quickstart on Wed Nov  8 12:07:40 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
import datetime

sys.path.insert(0, os.path.abspath('..'))

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '2.4'


# Add any Sphinx extension module names here, as strings.
# They can be extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.duration'
]


# This value controls how to represents typehints. The setting takes the following values:
#     'signature' – Show typehints as its signature (default)
#     'none' – Do not show typehints
# New in version 2.1.
autodoc_typehints = 'none'


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']


# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string: source_suffix = ['.rst', '.md']
source_suffix = '.rst'


# The master toctree document.
master_doc = 'index'


# General information about the project.
project = 'CAT'
_year = str(datetime.datetime.now().year)
author = 'B. F. van Beek & J. Belic'
copyright = f'{_year}, {author}'


# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the built documents.
release = '0.10.2'  # The full version, including alpha/beta/rc tags.
version = release.rsplit('.', maxsplit=1)[0]


# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = [
    '_build',
    'Thumbs.db',
    '.DS_Store'
]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {'includehidden': False}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory.
# They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
# This is required for the alabaster theme
# refs: http://alabaster.readthedocs.io/en/latest/installation.html#sidebars
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
        'donate.html',
    ]
}


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'catdoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}


# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc,
     'CAT.tex',
     'CAT Documentation',
     author,
     'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc,
     'CAT',
     'CAT Documentation',
     [author],
     1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc,
     'CAT',
     'CAT Documentation',
     author,
     'A collection of tools designed for the construction, and subsequent analysis, of various chemical compounds.',
     'Miscellaneous'),
]


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'h5py': ('https://docs.h5py.org/en/latest/', None),
    'pandas': ('https://pandas.pydata.org/pandas-docs/stable/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'plams': ('https://www.scm.com/doc/plams/', None),
    'qmflows': ('https://qmflows.readthedocs.io/en/latest/', None),
    'FOX': ('https://auto-fox.readthedocs.io/en/latest/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
    'pymongo': ('https://pymongo.readthedocs.io/en/stable/', None),
    'nanoutils': ('https://nano-utils.readthedocs.io/en/latest/', None),
    'rdkit': ('https://www.rdkit.org/docs/', None),
}


# This value selects if automatically documented members are sorted alphabetical (value 'alphabetical'),
# by member type (value 'groupwise') or by source order (value 'bysource').
autodoc_member_order = 'bysource'


# True to parse NumPy style docstrings.
# False to disable support for NumPy style docstrings.
# Defaults to True.
napoleon_numpy_docstring = True

# True to use the :ivar: role for instance variables.
# False to use the .. attribute:: directive instead.
# Defaults to False.
napoleon_use_ivar = False


# True to parse NumPy style docstrings.
# False to disable support for NumPy style docstrings.
# Defaults to True.
napoleon_google_docstring = False


# True to use the .. admonition:: directive for the Example and Examples sections.
# False to use the .. rubric:: directive instead. One may look better than the other depending on what HTML theme is used.
# Defaults to False.
napoleon_use_admonition_for_examples = True


# True to use the .. admonition:: directive for Notes sections.
# False to use the .. rubric:: directive instead.
#  Defaults to False.
napoleon_use_admonition_for_notes = True


# True to use the .. admonition:: directive for References sections.
# False to use the .. rubric:: directive instead.
# Defaults to False.
napoleon_use_admonition_for_references = True

# This value contains a list of modules to be mocked up.
# This is useful when some external dependencies are not met at build time and break the building
# process. You may only specify the root package of the dependencies themselves and
# omit the sub-modules:
autodoc_mock_imports = ["rdkit"]

# A string of reStructuredText that will be included at the end of every source file that is read.
# This is a possible place to add substitutions that should be available in every file (another being rst_prolog).
rst_epilog = """
.. |plams.Molecule| replace:: :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>`
.. |plams.Atom| replace:: :class:`plams.Atom<scm.plams.mol.atom.Atom>`
.. |plams.Bond| replace:: :class:`plams.Bond<scm.plams.mol.bond.Bond>`
.. |Atom.properties| replace:: :class:`plams.Atom.properties<scm.plams.mol.atom.Atom>`
.. |Atom.symbol| replace:: :class:`plams.Atom.symbol<scm.plams.mol.atom.Atom>`
.. |Atom.mass| replace:: :class:`plams.Atom.mass<scm.plams.mol.atom.Atom>`
.. |Sequence| replace:: :class:`Sequence<collections.abc.Sequence>`
.. |Iterable| replace:: :class:`Iterable<collections.abc.Iterable>`
.. |Any| replace:: :class:`Any<object>`
.. |OrderedDict| replace:: :class:`OrderedDict<collections.OrderedDict>`
"""
