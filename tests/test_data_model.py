from biosimulators_pysces import data_model
from biosimulators_utils.simulator.data_model import SoftwareInterface
import json
import os
import unittest


class DataModelTestCase(unittest.TestCase):

    def test_data_model_matches_specifications(self):
        with open(os.path.join(os.path.dirname(__file__), '..', 'biosimulators.json'), 'r') as file:
            specs = json.load(file)

        self.assertEqual(
            set(data_model.KISAO_ALGORITHM_MAP.keys()),
            set(alg_specs['kisaoId']['id'] for alg_specs in specs['algorithms']
                if SoftwareInterface.biosimulators_docker_image.value in alg_specs['availableSoftwareInterfaceTypes']))

        for alg_specs in specs['algorithms']:
            if SoftwareInterface.biosimulators_docker_image.value not in alg_specs['availableSoftwareInterfaceTypes']:
                continue

            alg_kisao_id = alg_specs['kisaoId']['id']
            alg_props = data_model.KISAO_ALGORITHM_MAP[alg_kisao_id]
            alg_param_props = alg_props['settings']

            self.assertEqual(
                set(alg_param_props.keys()),
                set(param_specs['kisaoId']['id'] for param_specs in alg_specs['parameters']
                    if SoftwareInterface.biosimulators_docker_image.value in param_specs['availableSoftwareInterfaceTypes']))

            for param_specs in alg_specs['parameters']:
                param_kisao_id = param_specs['kisaoId']['id']
                param_props = alg_param_props[param_kisao_id]
                self.assertEqual(param_props['type'].value, param_specs['type'])
