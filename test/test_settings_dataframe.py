"""Tests for :mod:`CAT.settings_dataframe`."""

import numpy as np

from CAT.assertion_functions import assert_value
from CAT.frozen_settings import FrozenSettings
from CAT.settings_dataframe import (SettingsDataFrame, SettingsSeries)

_DICT = {'a': True, 'b': False, 'c': [1, 2, 3, 4]}
DF = SettingsDataFrame(np.random.rand(10, 3), settings=_DICT)
SERIES = SettingsSeries(np.random.rand(10), settings=_DICT)


def test_df_and_series() -> None:
    """Tests for :class:`.SettingsDataFrame` and :class:`.SettingsSeries`."""
    settings = FrozenSettings(_DICT)

    assert_value(DF.settings, settings)
    assert_value(SERIES.settings, settings)
    assert_value(DF[0].settings, settings)
    assert_value(SERIES.to_frame().settings, settings)
    assert_value(DF.copy().settings, settings)
    assert_value(SERIES.copy().settings, settings)
