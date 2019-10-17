import pathlib
import yaml

import CAT

__all__ = ['ASA']


def load_templates():
    base = pathlib.Path(CAT.__path__[0])
    path = base / 'data' / 'workflow_dicts' / 'workflow.yaml'

    with open(path, 'r') as f:
        return yaml.load(f, Loader=yaml.FullLoader)


ASA, *_ = sorted(load_templates().values(), key=str)
