"""
CAT.logger
==========

Index
-----
.. currentmodule:: CAT.logger
.. autosummary::
    logger

API
---
.. autodata:: logger
    :annotation: = logging.RootLogger

"""

import logging

__all__ = ['logger']

#: A logger for CAT
logger = logging.getLogger()
logging.basicConfig(level=logging.DEBUG,
                    format='[%(asctime)s] %(levelname)s: %(message)s',
                    datefmt='%H:%M:%S')
