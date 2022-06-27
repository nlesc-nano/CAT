from typing import Mapping, Tuple, List, Type
from .workflow import WorkFlow as _WorkFlow

__all__: List[str] = ...

WorkFlow: Type[_WorkFlow] = ...

OPT: Tuple[str, str] = ...
MOL: Tuple[str, str] = ...
FORMULA: Tuple[str, str] = ...
HDF5_INDEX: Tuple[str, str] = ...
JOB_SETTINGS_QD_OPT: Tuple[str, str] = ...
JOB_SETTINGS_CRS: Tuple[str, str] = ...
JOB_SETTINGS_BDE: Tuple[str, str] = ...
JOB_SETTINGS_ASA: Tuple[str, str] = ...
JOB_SETTINGS_CDFT: Tuple[str, str] = ...
CDFT_MU: Tuple[str, str] = ...
CDFT_CHI: Tuple[str, str] = ...
CDFT_ETA: Tuple[str, str] = ...
CDFT_S: Tuple[str, str] = ...
CDFT_GAMMA: Tuple[str, str] = ...
CDFT_OMEGA: Tuple[str, str] = ...
CDFT_NUCLEOFUGE: Tuple[str, str] = ...
CDFT_ELECTROFUGE: Tuple[str, str] = ...
CDFT_W_MINUS: Tuple[str, str] = ...
CDFT_W_PLUS: Tuple[str, str] = ...
CDFT_ELECTROPHILICITY: Tuple[str, str] = ...
CDFT_DELTAF_MINUS: Tuple[str, str] = ...
CDFT_DELTAF_PLUS: Tuple[str, str] = ...
CDFT_MU_MINUS: Tuple[str, str] = ...
CDFT_MU_PLUS: Tuple[str, str] = ...
ASA_INT: Tuple[str, str] = ...
ASA_STRAIN: Tuple[str, str] = ...
ASA_E: Tuple[str, str] = ...
SETTINGS1: Tuple[str, str] = ...
SETTINGS2: Tuple[str, str] = ...
SETTINGS_SOLV1: Tuple[str, str] = ...
SETTINGS_SOLV2: Tuple[str, str] = ...
SETTINGS_ASA: Tuple[str, str] = ...
SETTINGS_BDE1: Tuple[str, str] = ...
SETTINGS_BDE2: Tuple[str, str] = ...
SETTINGS_CDFT: Tuple[str, str] = ...
V_BULK: Tuple[str, str] = ...
CONE_ANGLE: Tuple[str, str] = ...
BRANCH_DISTANCE: Tuple[str, str] = ...
BRANCH_SIZE: Tuple[str, str] = ...
