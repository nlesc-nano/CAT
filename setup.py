#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit CAT/__version__.py
version = {}
with open(os.path.join(here, 'CAT', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='CAT',
    version=version['__version__'],
    description='A collection of tools designed for the automatic construction, and subsequent analysis, of chemical compounds.',
    long_description=readme + '\n\n',
    author=['Bas van Beek'],
    author_email='b.f.van.beek@vu.nl',
    url='https://github.com/BvB93/CAT',
    packages=[
        'CAT',
        'CAT.analysis',
        'CAT.attachment',
        'CAT.data',
        'CAT.data.coskf',
        'CAT.data.templates',
        'CAT.data_handling'
    ],
    package_dir={'CAT': 'CAT'},
    package_data={
        'CAT': [
            'data/templates/*json',
            'data/templates/*yaml',
            'data/coskf/*coskf'
        ]
    },
    entry_points={
          'console_scripts': ['init_cat=CAT.data_handling.input_parser:main']
    },
    include_package_data=True,
    license='GNU Lesser General Public License v3 or later',
    zip_safe=False,
    keywords=[
        'quantum-mechanics',
        'science',
        'chemistry',
        'python-3',
        'python-3-6',
        'automation',
        'scientific-workflows'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
        'License :: OSI Approved :: GNU Lesser General Public License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
    ],
    test_suite='tests',
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'pyyaml',
        'schema',
        'plams@git+https://github.com/SCM-NV/PLAMS@master#egg=plams-1.2',
        'qmflows@git+https://github.com/SCM-NV/qmflows@master'
    ],
    setup_requires=[
        'pytest-runner',
        'sphinx',
        'sphinx_rtd_theme',
        'recommonmark'
    ],
    tests_require=[
        'pytest',
        'pytest-cov',
        'pycodestyle',
    ],
    extras_require={
        'test': ['pytest', 'pytest-cov', 'pytest-mock', 'nbsphinx', 'pycodestyle'],
        'doc': ['sphinx', 'sphinx_rtd_theme', 'nbsphinx']
    }
)
