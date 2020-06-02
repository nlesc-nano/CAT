#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit CAT/__version__.py
version = {}
with open(os.path.join(here, 'CAT', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst', encoding='utf-8') as readme_file:
    readme = readme_file.read()

setup(
    name='CAT',
    version=version['__version__'],
    description=('A collection of tools designed for the automatic '
                 'construction of chemical compounds.'),
    long_description=readme + '\n\n',
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
            'data/CORES/*xyz',
            'data/LIGANDS/*xyz',
            'workflows/workflow_yaml.yaml',
            'py.typed',
            '*.pyi'
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
        'automation',
        'scientific-workflows'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'
    ],
    test_suite='tests',
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'scipy',
        'pandas<1.0.0',
        'pyyaml>=5.1',
        'schema',
        'AssertionLib>=2.2.3',
        'plams@git+https://github.com/SCM-NV/PLAMS@a5696ce62c09153a9fa67b2b03a750913e1d0924',
        'qmflows@git+https://github.com/SCM-NV/qmflows@master',
    ],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=[
        'pytest',
        'pytest-cov',
        'pytest-mock',
        'pycodestyle',
        'data-CAT@git+https://github.com/nlesc-nano/data-CAT@devel',
        'nano-CAT@git+https://github.com/nlesc-nano/nano-CAT@devel'
    ],
    extras_require={
        'test': ['pytest',
                 'pytest-cov',
                 'pytest-mock',
                 'pycodestyle',
                 'data-CAT@git+https://github.com/nlesc-nano/data-CAT@devel',
                 'nano-CAT@git+https://github.com/nlesc-nano/nano-CAT@devel'],
        'doc': ['sphinx>=2.0', 'sphinx_rtd_theme', 'sphinx-autodoc-typehints']
    }
)
