"""
CAT.workflows.workflow_dicts
============================

A module for loading :class:`WorkFlow` templates.

"""


import pathlib

import yaml
import numpy as np

import CAT

__all__ = ['finalize_templates']

OPT = ('opt', '')
MOL = ('mol', '')
HDF5_INDEX = ('hdf5 index', '')
JOB_SETTINGS_QD_OPT = ('job_settings_QD_opt', '')
ASA_INT = ('ASA', 'E_int')
ASA_STRAIN = ('ASA', 'E_strain')
ASA_E = ('ASA', 'E')
SETTINGS1 = ('settings', '1')
SETTINGS2 = ('settings', '2')
SETTINGS_ASA = ('settings', 'ASA 1')

FrozenSettings = CAT.frozen_settings.FrozenSettings

_TEMPLATE_UPDATE: dict = {
    'asa': {'import_columns': {ASA_INT: np.nan, ASA_STRAIN: np.nan, ASA_E: np.nan},
            'export_columns': (SETTINGS_ASA, ASA_INT, ASA_STRAIN, ASA_E)},
    'ligand_opt': {'import_columns': {HDF5_INDEX: -1, OPT: False},
                   'export_columns': (HDF5_INDEX, OPT, SETTINGS1, SETTINGS2)},
    'qd_attach': {'import_columns': {HDF5_INDEX: -1, OPT: False},
                  'export_columns': (HDF5_INDEX,)},
    'qd_opt': {'import_columns': {HDF5_INDEX: -1, OPT: False},
               'export_columns': (HDF5_INDEX, OPT, JOB_SETTINGS_QD_OPT, SETTINGS1, SETTINGS2)},
}


def load_templates() -> dict:
    """Load the templates from ``CAT/data/workflow_dicts/workflow.yaml``."""
    base = pathlib.Path(CAT.__path__[0])
    path = base / 'data' / 'workflow_dicts' / 'workflow.yaml'

    with open(path, 'r') as f:
        return yaml.load(f, Loader=yaml.FullLoader)


def finalize_templates() -> FrozenSettings:
    """Update the templates using :data:`._TEMPLATE_UPDATE`."""
    templates = load_templates()
    for k, v1 in templates.items():
        v2 = _TEMPLATE_UPDATE[k]
        v1.update(v2)

        # Convert all lists stored under the 'template' key into tuples
        for key, list_ in v1['template'].items():
            v1['template'][key] = tuple(list_)

    return FrozenSettings(templates)
