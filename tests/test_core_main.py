""" Tests of the command-line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_pysces import __main__
from biosimulators_pysces import core
from biosimulators_utils.combine import data_model as combine_data_model
from biosimulators_utils.combine.exceptions import CombineArchiveExecutionError
from biosimulators_utils.combine.io import CombineArchiveWriter
from biosimulators_utils.report import data_model as report_data_model
from biosimulators_utils.report.io import ReportReader
from biosimulators_utils.simulator.exceptions import AlgorithmDoesNotSupportModelFeatureException
from biosimulators_utils.simulator.exec import exec_sedml_docs_in_archive_with_containerized_simulator
from biosimulators_utils.simulator.specs import gen_algorithms_from_specs
from biosimulators_utils.sedml import data_model as sedml_data_model
from biosimulators_utils.sedml.io import SedmlSimulationWriter
from biosimulators_utils.sedml.utils import append_all_nested_children_to_doc
from biosimulators_utils.warnings import BioSimulatorsWarning
from kisao.exceptions import AlgorithmCannotBeSubstitutedException
from kisao.warnings import AlgorithmSubstitutedWarning
from unittest import mock
import datetime
import dateutil.tz
import numpy
import numpy.testing
import os
import pysces
import shutil
import tempfile
import unittest


class CliTestCase(unittest.TestCase):
    DOCKER_IMAGE = 'ghcr.io/biosimulators/biosimulators_pysces/pysces:latest'
    NAMESPACES = {
        'sbml': 'http://www.sbml.org/sbml/level2/version4',
    }

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_exec_sed_task_successfully(self):
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=os.path.join(os.path.dirname(__file__), 'fixtures', 'biomd0000000002.xml'),
                language=sedml_data_model.ModelLanguage.SBML.value,
                changes=[],
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000088',
                    changes=[
                        sedml_data_model.AlgorithmParameterChange(
                            kisao_id='KISAO_0000209',
                            new_value='1e-8',
                        ),
                    ],
                ),
                initial_time=5.,
                output_start_time=10.,
                output_end_time=20.,
                number_of_points=20,
            ),
        )

        variables = [
            sedml_data_model.Variable(id='time', symbol=sedml_data_model.Symbol.time, task=task),
            sedml_data_model.Variable(
                id='AL',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='AL']",
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
            sedml_data_model.Variable(
                id='BLL',
                target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="BLL"]',
                target_namespaces=self.NAMESPACES,
                task=task),
            sedml_data_model.Variable(
                id='IL',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='IL']",
                target_namespaces=self.NAMESPACES,
                task=task),
        ]

        variable_results, _ = core.exec_sed_task(task, variables)

        self.assertTrue(sorted(variable_results.keys()), sorted([var.id for var in variables]))
        self.assertEqual(variable_results[variables[0].id].shape, (task.simulation.number_of_points + 1,))
        numpy.testing.assert_almost_equal(
            variable_results['time'],
            numpy.linspace(task.simulation.output_start_time, task.simulation.output_end_time, task.simulation.number_of_points + 1),
        )

        for results in variable_results.values():
            self.assertFalse(numpy.any(numpy.isnan(results)))

    def test_exec_sed_task_error_handling(self):
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            task = sedml_data_model.Task(
                model=sedml_data_model.Model(
                    source=os.path.join(os.path.dirname(__file__), 'fixtures', 'biomd0000000002.xml'),
                    language=sedml_data_model.ModelLanguage.SBML.value,
                    changes=[],
                ),
                simulation=sedml_data_model.UniformTimeCourseSimulation(
                    algorithm=sedml_data_model.Algorithm(
                        kisao_id='KISAO_0000001',
                        changes=[
                            sedml_data_model.AlgorithmParameterChange(
                                kisao_id='KISAO_0000209',
                                new_value='2e-8',
                            ),
                        ],
                    ),
                    initial_time=5.,
                    output_start_time=10.,
                    output_end_time=20.,
                    number_of_points=20,
                ),
            )

            variables = [
                sedml_data_model.Variable(
                    id='time', symbol=sedml_data_model.Symbol.time, task=task),
                sedml_data_model.Variable(
                    id='AL',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='AL']",
                    target_namespaces=self.NAMESPACES,
                    task=task),
                sedml_data_model.Variable(
                    id='BLL',
                    target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="BLL"]',
                    target_namespaces=self.NAMESPACES,
                    task=task),
                sedml_data_model.Variable(
                    id='IL',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='IL']",
                    target_namespaces=self.NAMESPACES,
                    task=task),
            ]

            # Configure task
            task.model.source = os.path.join(self.dirname, 'bad-model.xml')
            with open(task.model.source, 'w') as file:
                file.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
                file.write('<sbml2 xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">')
                file.write('  <model id="model">')
                file.write('  </model>')
                file.write('</sbml2>')
            with self.assertRaisesRegex(ValueError, 'could not be imported'):
                core.exec_sed_task(task, [])
            task.model.source = os.path.join(os.path.dirname(__file__), 'fixtures', 'biomd0000000002.xml')

            task.simulation.algorithm.kisao_id = 'KISAO_0000448'
            with self.assertRaisesRegex(AlgorithmCannotBeSubstitutedException, 'No algorithm can be substituted'):
                core.exec_sed_task(task, variables)
            task.simulation.algorithm.kisao_id = 'KISAO_0000088'

            task.simulation.algorithm.changes[0].kisao_id = 'KISAO_0000001'
            with self.assertRaisesRegex(NotImplementedError, 'is not supported'):
                core.exec_sed_task(task, variables)
            task.simulation.algorithm.changes[0].kisao_id = 'KISAO_0000209'

            task.simulation.algorithm.changes[0].new_value = 'two e minus 8'
            with self.assertRaisesRegex(ValueError, 'is not a valid'):
                core.exec_sed_task(task, variables)
            task.simulation.algorithm.changes[0].new_value = '2e-8'

            task.simulation.output_end_time = 20.1
            with self.assertRaisesRegex(NotImplementedError, 'must specify an integer number of time points'):
                core.exec_sed_task(task, variables)
            task.simulation.output_end_time = 20.

            variables[0].symbol += '*'
            with self.assertRaisesRegex(NotImplementedError, 'symbols are not supported'):
                core.exec_sed_task(task, variables)
            variables[0].symbol = sedml_data_model.Symbol.time

            variables[1].target = "/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter[@id='kf_0']"
            with self.assertRaisesRegex(ValueError, 'targets could not be recorded'):
                core.exec_sed_task(task, variables)
            variables[1].target = "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='AL']"

        # algorithm substition
        task = sedml_data_model.Task(
            model=sedml_data_model.Model(
                source=os.path.join(os.path.dirname(__file__), 'fixtures', 'biomd0000000002.xml'),
                language=sedml_data_model.ModelLanguage.SBML.value,
                changes=[],
            ),
            simulation=sedml_data_model.UniformTimeCourseSimulation(
                algorithm=sedml_data_model.Algorithm(
                    kisao_id='KISAO_0000088',
                    changes=[
                        sedml_data_model.AlgorithmParameterChange(
                            kisao_id='KISAO_0000209',
                            new_value='1e-8',
                        ),
                    ],
                ),
                initial_time=5.,
                output_start_time=10.,
                output_end_time=20.,
                number_of_points=20,
            ),
        )

        variables = []

        task.simulation.algorithm.changes[0].new_value = 'not a number'
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaisesRegex(ValueError, 'is not a valid'):
                core.exec_sed_task(task, variables)

        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            with self.assertWarnsRegex(BioSimulatorsWarning, 'ignored because it is not a valid float'):
                core.exec_sed_task(task, variables)

        task.simulation.algorithm.changes[0].kisao_id = 'KISAO_9999999'
        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaisesRegex(NotImplementedError, 'is not supported'):
                core.exec_sed_task(task, variables)

        with mock.patch.dict('os.environ', {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            with self.assertWarnsRegex(BioSimulatorsWarning, 'ignored because it is not supported'):
                core.exec_sed_task(task, variables)

    def test_exec_sedml_docs_in_combine_archive_successfully(self):
        doc, archive_filename = self._build_combine_archive()

        out_dir = os.path.join(self.dirname, 'out')
        core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                                report_formats=[
                                                    report_data_model.ReportFormat.h5,
                                                    report_data_model.ReportFormat.csv,
                                                ],
                                                bundle_outputs=True,
                                                keep_individual_outputs=True)

        self._assert_combine_archive_outputs(doc, out_dir)

    def _build_combine_archive(self, algorithm=None):
        doc = self._build_sed_doc(algorithm=algorithm)

        archive_dirname = os.path.join(self.dirname, 'archive')
        if not os.path.isdir(archive_dirname):
            os.mkdir(archive_dirname)

        model_filename = os.path.join(archive_dirname, 'model_1.xml')
        shutil.copyfile(
            os.path.join(os.path.dirname(__file__), 'fixtures', 'biomd0000000002.xml'),
            model_filename)

        sim_filename = os.path.join(archive_dirname, 'sim_1.sedml')
        SedmlSimulationWriter().run(doc, sim_filename)

        updated = datetime.datetime(2020, 1, 2, 1, 2, 3, tzinfo=dateutil.tz.tzutc())
        archive = combine_data_model.CombineArchive(
            contents=[
                combine_data_model.CombineArchiveContent(
                    'model_1.xml', combine_data_model.CombineArchiveContentFormat.SBML.value, updated=updated),
                combine_data_model.CombineArchiveContent(
                    'sim_1.sedml', combine_data_model.CombineArchiveContentFormat.SED_ML.value, updated=updated),
            ],
            updated=updated,
        )
        archive_filename = os.path.join(self.dirname,
                                        'archive.omex' if algorithm is None else 'archive-{}.omex'.format(algorithm.kisao_id))
        CombineArchiveWriter().run(archive, archive_dirname, archive_filename)

        return (doc, archive_filename)

    def _build_sed_doc(self, algorithm=None):
        if algorithm is None:
            algorithm = sedml_data_model.Algorithm(
                kisao_id='KISAO_0000088',
                changes=[
                    sedml_data_model.AlgorithmParameterChange(
                        kisao_id='KISAO_0000209',
                        new_value='1e-8',
                    ),
                ],
            )

        doc = sedml_data_model.SedDocument()
        doc.models.append(sedml_data_model.Model(
            id='model_1',
            source='model_1.xml',
            language=sedml_data_model.ModelLanguage.SBML.value,
            changes=[],
        ))
        doc.simulations.append(sedml_data_model.UniformTimeCourseSimulation(
            id='sim_1_time_course',
            algorithm=algorithm,
            initial_time=0.,
            output_start_time=0.1,
            output_end_time=0.2,
            number_of_points=20,
        ))
        doc.tasks.append(sedml_data_model.Task(
            id='task_1',
            model=doc.models[0],
            simulation=doc.simulations[0],
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_time',
            variables=[
                sedml_data_model.Variable(
                    id='var_time',
                    symbol=sedml_data_model.Symbol.time,
                    task=doc.tasks[0],
                ),
            ],
            math='var_time',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_AL',
            variables=[
                sedml_data_model.Variable(
                    id='var_AL',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='AL']",
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_AL',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_BLL',
            variables=[
                sedml_data_model.Variable(
                    id='var_BLL',
                    target='/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id="BLL"]',
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_BLL',
        ))
        doc.data_generators.append(sedml_data_model.DataGenerator(
            id='data_gen_IL',
            variables=[
                sedml_data_model.Variable(
                    id='var_IL',
                    target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='IL']",
                    target_namespaces=self.NAMESPACES,
                    task=doc.tasks[0],
                ),
            ],
            math='var_IL',
        ))
        doc.outputs.append(sedml_data_model.Report(
            id='report_1',
            data_sets=[
                sedml_data_model.DataSet(id='data_set_time', label='Time', data_generator=doc.data_generators[0]),
                sedml_data_model.DataSet(id='data_set_AL', label='AL', data_generator=doc.data_generators[1]),
                sedml_data_model.DataSet(id='data_set_BLL', label='BLL', data_generator=doc.data_generators[2]),
                sedml_data_model.DataSet(id='data_set_IL', label='IL', data_generator=doc.data_generators[3]),
            ],
        ))

        append_all_nested_children_to_doc(doc)

        return doc

    def _assert_combine_archive_outputs(self, doc, out_dir):
        self.assertEqual(set(['reports.h5', 'reports.zip', 'sim_1.sedml']).difference(set(os.listdir(out_dir))), set())

        report = doc.outputs[0]

        # check HDF report
        report_results = ReportReader().run(report, out_dir, 'sim_1.sedml/report_1', format=report_data_model.ReportFormat.h5)

        self.assertEqual(sorted(report_results.keys()), sorted([d.id for d in doc.outputs[0].data_sets]))

        sim = doc.tasks[0].simulation
        self.assertEqual(len(report_results[report.data_sets[0].id]), sim.number_of_points + 1)
        numpy.testing.assert_almost_equal(
            report_results[report.data_sets[0].id],
            numpy.linspace(sim.output_start_time, sim.output_end_time, sim.number_of_points + 1),
        )

        for data_set_result in report_results.values():
            self.assertFalse(numpy.any(numpy.isnan(data_set_result)))

        # check CSV report
        report_results = ReportReader().run(report, out_dir, 'sim_1.sedml/report_1', format=report_data_model.ReportFormat.csv)

        self.assertEqual(sorted(report_results.keys()), sorted([d.id for d in doc.outputs[0].data_sets]))

        sim = doc.tasks[0].simulation
        self.assertEqual(len(report_results[report.data_sets[0].id]), sim.number_of_points + 1)
        numpy.testing.assert_almost_equal(
            report_results[report.data_sets[0].id],
            numpy.linspace(sim.output_start_time, sim.output_end_time, sim.number_of_points + 1),
        )

        for data_set_result in report_results.values():
            self.assertFalse(numpy.any(numpy.isnan(data_set_result)))

    def test_algorithm_substitution(self):
        doc, archive_filename = self._build_combine_archive()

        out_dir = self.dirname

        # NONE
        env = {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}

        original_pysces_model = pysces.model

        def pysces_model(filename):
            model = original_pysces_model(filename)
            model.__events__ = True
            return model

        with mock.patch.dict(os.environ, env):
            with mock.patch('pysces.model', side_effect=pysces_model):
                with self.assertRaises(CombineArchiveExecutionError):
                    core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                                            report_formats=[
                                                                report_data_model.ReportFormat.h5,
                                                                report_data_model.ReportFormat.csv,
                                                            ],
                                                            bundle_outputs=True,
                                                            keep_individual_outputs=True)

        # SAME FRAMEWORK
        env = {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}
        with mock.patch.dict(os.environ, env):
            with mock.patch('pysces.model', side_effect=pysces_model):
                with mock.patch.object(pysces.PyscesModel.PysMod, 'Simulate', side_effect=Exception('Stop')):
                    with self.assertRaisesRegex(CombineArchiveExecutionError, 'Stop$'):
                        with self.assertWarns(AlgorithmSubstitutedWarning):
                            core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                                                    report_formats=[
                                                                        report_data_model.ReportFormat.h5,
                                                                        report_data_model.ReportFormat.csv,
                                                                    ],
                                                                    bundle_outputs=True,
                                                                    keep_individual_outputs=True)

    def test_exec_sedml_docs_in_combine_archive_with_all_algorithms(self):
        for alg in gen_algorithms_from_specs(os.path.join(os.path.dirname(__file__), '..', 'biosimulators.json')).values():
            doc, archive_filename = self._build_combine_archive(algorithm=alg)

            out_dir = os.path.join(self.dirname, alg.kisao_id)
            core.exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                                    report_formats=[
                                                        report_data_model.ReportFormat.h5,
                                                        report_data_model.ReportFormat.csv,
                                                    ],
                                                    bundle_outputs=True,
                                                    keep_individual_outputs=True)
            self._assert_combine_archive_outputs(doc, out_dir)

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: ')

    def test_exec_sedml_docs_in_combine_archive_with_cli(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        env = self._get_combine_archive_exec_env()

        with mock.patch.dict(os.environ, env):
            with __main__.App(argv=['-i', archive_filename, '-o', out_dir]) as app:
                app.run()

        self._assert_combine_archive_outputs(doc, out_dir)

    def _get_combine_archive_exec_env(self):
        return {
            'REPORT_FORMATS': 'h5,csv',
            'BUNDLE_OUTPUTS': '1',
            'KEEP_INDIVIDUAL_OUTPUTS': '1',
        }

    def test_exec_sedml_docs_in_combine_archive_with_docker_image(self):
        doc, archive_filename = self._build_combine_archive()
        out_dir = os.path.join(self.dirname, 'out')
        docker_image = self.DOCKER_IMAGE
        env = self._get_combine_archive_exec_env()

        exec_sedml_docs_in_archive_with_containerized_simulator(
            archive_filename, out_dir, docker_image, environment=env, pull_docker_image=False)

        self._assert_combine_archive_outputs(doc, out_dir)

    def test_exec_sedml_docs_in_combine_archive_with_docker_image(self):
        archive_filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'Parmar-BMC-Syst-Biol-2017-iron-distribution.omex')
        out_dir = os.path.join(self.dirname, 'out')
        docker_image = self.DOCKER_IMAGE
        env = {
            'REPORT_FORMATS': 'h5'
        }

        exec_sedml_docs_in_archive_with_containerized_simulator(
            archive_filename, out_dir, docker_image, environment=env, pull_docker_image=False)

        report = sedml_data_model.Report(
            data_sets=[
                sedml_data_model.DataSet(id='data_set_time', label='time'),
                sedml_data_model.DataSet(id='data_set_FeDuo', label='FeDuo'),
            ]
        )

        report_results = ReportReader().run(report, out_dir, 'Parmar2017_Deficient_Rich_tracer.sedml/report_1',
                                            format=report_data_model.ReportFormat.h5)

        self.assertEqual(set(report_results.keys()), set(['data_set_time', 'data_set_FeDuo']))

        self.assertEqual(len(report_results['data_set_time']), 300 + 1)
        numpy.testing.assert_almost_equal(
            report_results['data_set_time'],
            numpy.linspace(0., 5100., 300 + 1),
        )

        for data_set_result in report_results.values():
            self.assertFalse(numpy.any(numpy.isnan(data_set_result)))
