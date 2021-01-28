"""Placeholder."""

from .attachment import distribution
from .attachment.distribution import (
    distribute_idx as distribute_idx,
    uniform_idx as uniform_idx,
)

# flake8: noqa: E402

__all__ = distribution.__all__
__doc__ = distribution.__doc__
del distribution
