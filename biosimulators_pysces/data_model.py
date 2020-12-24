""" Data model for mapping KiSAO terms to PySCeS integrators and their settings

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-12-23
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_utils.data_model import ValueType

__all__ = [
    'KISAO_ALGORITHM_MAP',
]

KISAO_ALGORITHM_MAP = {
    'KISAO_0000019': {
        'id': 'CVODE',
        'settings': {
            "KISAO_0000211": {
                "id": "cvode_abstol",
                "name": "Absolute tolerance",
                "type": ValueType.float,
                "value": 1e-15,
            },
            "KISAO_0000571": {
                "id": "cvode_abstol_factor",
                "name": "Absolute tolerance adjustment factor",
                "type": ValueType.float,
                "value": 1e-6,
            },
            "KISAO_0000209": {
                "id": "cvode_reltol",
                "name": "Relative tolerance",
                "type": ValueType.float,
                "value": 1e-9,
            },
            "KISAO_0000570": {
                "id": "cvode_auto_tol_adjust",
                "name": "Auto reduce tolerances",
                "type": ValueType.boolean,
                "value": True,
            },
            "KISAO_0000415": {
                "id": "cvode_mxstep",
                "name": "Maximum number of steps",
                "type": ValueType.integer,
                "value": 1000,
            },
        },
    },
    'KISAO_0000088': {
        'id': 'LSODA',
        'settings': {
            'KISAO_0000211': {
                'id': 'lsoda_atol',
                "name": "Absolute tolerance",
                'type': ValueType.float,
                'default': 1.0e-12,
            },
            'KISAO_0000209': {
                'id': 'lsoda_rtol',
                "name": "Relative tolerance",
                'type': ValueType.float,
                'default': 1.0e-7,
            },
            'KISAO_0000415': {
                'id': 'lsoda_mxstep',
                "name": "Maximum number of steps",
                'type': ValueType.integer,
                'default': None,
            },
            'KISAO_0000559': {
                'id': 'lsoda_h0',
                "name": "Step size to be attempted for the first step",
                'type': ValueType.float,
                'default': None,
            },
            'KISAO_0000467': {
                'id': 'lsoda_hmax',
                "name": "Maximum absolute step size allowed",
                'type': ValueType.float,
                'default': None,
            },
            'KISAO_0000485': {
                'id': 'lsoda_hmin',
                "name": "Minimum absolute step size allowed",
                'type': ValueType.float,
                'default': None,
            },
            'KISAO_0000219': {
                'id': 'lsoda_mxordn',
                "name": "Maximum order to be allowed for the non-stiff (Adams) method",
                'type': ValueType.integer,
                'default': 12,
            },
            'KISAO_0000220': {
                'id': 'lsoda_mxords',
                "name": "Maximum order to be allowed for the stiff (BDF) method",
                'type': ValueType.integer,
                'default': 5,
            },
        },
    },
}
