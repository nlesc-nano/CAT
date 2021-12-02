#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import importlib

from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit CAT/__version__.py
version = {}
with open(os.path.join(here, 'CAT', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst', encoding='utf-8') as readme_file:
    readme = readme_file.read()

docs_require = [
    'sphinx>=2.4',
    'sphinx_rtd_theme',
    'data-CAT>=0.7.0',
    'nano-CAT>=0.7.0',
    'auto-FOX>=0.10.0',
]

tests_require = [
    'pytest>=5.4.0',
    'pytest-cov',
    'pytest-mock',
]
if sys.version_info >= (3, 7):
    tests_require += docs_require

# Check if rdkit is manually installed (as it is not available via pypi)
try:
    importlib.import_module("rdkit")
except ModuleNotFoundError:
    print(
        "`CAT` requires the `rdkit` package: https://anaconda.org/conda-forge/rdkit",
        file=sys.stderr,
    )

setup(
    name='nlesc-CAT',
    version=version['__version__'],
    description=('A collection of tools designed for the automatic '
                 'construction of chemical compounds.'),
    long_description=f'{readme}\n\n',
    long_description_content_type='text/x-rst',
    author=['Bas van Beek', 'Jelena Belic'],
    author_email='b.f.van.beek@vu.nl',
    url='https://github.com/nlesc-nano/CAT',
    packages=[
        'CAT',
        'CAT.attachment',
        'CAT.data',
        'CAT.dye',
        'CAT.data.coskf',
        'CAT.data.templates',
        'CAT.data.CORES',
        'CAT.data.LIGANDS',
        'CAT.data_handling',
        'CAT.workflows'
    ],
    package_dir={'CAT': 'CAT'},
    package_data={
        'CAT': [
            'data/templates/*yaml',
            'data/coskf/*coskf',
            'data/coskf/misc/*coskf',
            'data/CORES/*xyz',
            'data/LIGANDS/*xyz',
            'workflows/workflow_yaml.yaml',
            'workflows/*.pyi',
            'py.typed'
        ]
    },
    entry_points={
        'console_scripts': ['init_cat=CAT.data_handling.entry_points:main']
    },
    include_package_data=True,
    license='GNU Lesser General Public License v3 or later',
    zip_safe=False,
    keywords=[
        'quantum-mechanics',
        'molecular-mechanics',
        'science',
        'chemistry',
        'python-3',
        'python-3-6',
        'python-3-7',
        'python-3-8',
        'python-3-9',
        'automation',
        'scientific-workflows'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Typing :: Typed',
    ],
    test_suite='tests',
    python_requires='>=3.6',
    install_requires=[
        'Nano-Utils>=0.4.3',
        'h5py',
        'numpy',
        'scipy',
        'pandas',
        'pyyaml>=5.1',
        'schema!=0.7.5',
        'AssertionLib>=2.2.3',
        'plams>=1.5.1',
        'contextlib2>=0.6.0; python_version=="3.6"',
        'typing-extensions>=3.7.4.2',
        'qmflows>=0.11.0',
    ],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=tests_require,
    extras_require={
        'test': tests_require,
        'doc': docs_require,
    }
)
