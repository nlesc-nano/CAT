"""Tests for :mod:`CAT.settings_dataframe`."""

import numpy as np

from CAT.assertion.assertion_manager import assertion
from CAT.frozen_settings import FrozenSettings
from CAT.settings_dataframe import (SettingsDataFrame, SettingsSeries)

_DICT = {'a': True, 'b': False, 'c': [1, 2, 3, 4]}
DF = SettingsDataFrame(np.random.rand(10, 3), settings=_DICT)
SERIES = SettingsSeries(np.random.rand(10), settings=_DICT)


def test_df_and_series() -> None:
    """Tests for :class:`.SettingsDataFrame` and :class:`.SettingsSeries`."""
    settings = FrozenSettings(_DICT)

    assertion.eq(DF.settings, settings)
    assertion.eq(SERIES.settings, settings)
    assertion.eq(DF[0].settings, settings)
    assertion.eq(SERIES.to_frame().settings, settings)
    assertion.eq(DF.copy().settings, settings)
    assertion.eq(SERIES.copy().settings, settings)
