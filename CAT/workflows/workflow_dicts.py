import pathlib

import yaml
import numpy as np

import CAT

__all__ = ['finalize_templates']

ASA_INT = ('ASA', 'E_int')
ASA_STRAIN = ('ASA', 'E_strain')
ASA_E = ('ASA', 'E')
SETTINGS_ASA = ('settings', 'ASA 1')

FrozenSettings = CAT.frozen_settings.FrozenSettings


def load_templates():
    base = pathlib.Path(CAT.__path__[0])
    path = base / 'data' / 'workflow_dicts' / 'workflow.yaml'

    with open(path, 'r') as f:
        return yaml.load(f, Loader=yaml.FullLoader)


def finalize_templates():
    template_update = {
        'asa': {'import_columns': {ASA_INT: np.nan, ASA_STRAIN: np.nan, ASA_E: np.nan},
                'export_columns': (SETTINGS_ASA, ASA_INT, ASA_STRAIN, ASA_E)}
    }

    templates = load_templates()
    for k, v1 in templates.items():
        v2 = template_update[k]
        v1.update(v2)
    return FrozenSettings(templates)
