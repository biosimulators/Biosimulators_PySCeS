""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-12-23
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import KISAO_ALGORITHM_MAP
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog  # noqa: F401
from biosimulators_utils.plot.data_model import PlotFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults  # noqa: F401
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, UniformTimeCourseSimulation,  # noqa: F401
                                                  Variable, Symbol)
from biosimulators_utils.utils.core import validate_str_value, parse_value
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.exec import exec_sed_doc
from biosimulators_utils.simulator.data_model import AlgorithmSubstitutionPolicy, ALGORITHM_SUBSTITUTION_POLICY_LEVELS
from biosimulators_utils.simulator.exceptions import AlgorithmDoesNotSupportModelFeatureException
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.simulator.warnings import AlgorithmSubstitutedWarning
from biosimulators_utils.utils.core import raise_errors_warnings
import functools
import os
cwd = os.getcwd()  # because PySCeS changes the working directory
import pysces  # noqa: E402
os.chdir(cwd)
import tempfile  # noqa: E402
import warnings  # noqa: E402


__all__ = ['exec_sedml_docs_in_combine_archive', 'exec_sed_task']


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                       report_formats=None, plot_formats=None,
                                       bundle_outputs=None, keep_individual_outputs=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        report_formats (:obj:`list` of :obj:`ReportFormat`, optional): report format (e.g., csv or h5)
        plot_formats (:obj:`list` of :obj:`PlotFormat`, optional): report format (e.g., pdf)
        bundle_outputs (:obj:`bool`, optional): if :obj:`True`, bundle outputs into archives for reports and plots
        keep_individual_outputs (:obj:`bool`, optional): if :obj:`True`, keep individual output files

    Returns:
        :obj:`CombineArchiveLog`: log
    """
    sed_doc_executer = functools.partial(exec_sed_doc, exec_sed_task)
    return exec_sedml_docs_in_archive(sed_doc_executer, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      report_formats=report_formats,
                                      plot_formats=plot_formats,
                                      bundle_outputs=bundle_outputs,
                                      keep_individual_outputs=keep_individual_outputs)


def exec_sed_task(task, variables, log=None):
    ''' Execute a task and save its results

    Args:
       task (:obj:`Task`): task
       variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
       log (:obj:`TaskLog`, optional): log for the task

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log
    '''
    log = log or TaskLog()

    model = task.model
    sim = task.simulation

    raise_errors_warnings(validation.validate_task(task),
                          error_summary='Task `{}` is invalid.'.format(task.id))
    raise_errors_warnings(validation.validate_model_language(task.model.language, ModelLanguage.SBML),
                          error_summary='Language for model `{}` is not supported.'.format(model.id))
    raise_errors_warnings(validation.validate_model_change_types(task.model.changes, ()),
                          error_summary='Changes for model `{}` are not supported.'.format(model.id))
    raise_errors_warnings(validation.validate_model_changes(task.model),
                          error_summary='Changes for model `{}` are invalid.'.format(model.id))
    raise_errors_warnings(validation.validate_simulation_type(task.simulation, (UniformTimeCourseSimulation, )),
                          error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
    raise_errors_warnings(validation.validate_simulation(task.simulation),
                          error_summary='Simulation `{}` is invalid.'.format(sim.id))
    raise_errors_warnings(validation.validate_data_generator_variables(variables),
                          error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))
    target_x_paths_to_sbml_ids = validation.validate_variable_xpaths(variables, task.model.source, attr='id')

    # Get the current working directory because PySCeS opaquely changes it
    cwd = os.getcwd()

    # Read the model
    interfaces = pysces.PyscesInterfaces.Core2interfaces()
    fid, model_filename = tempfile.mkstemp(suffix='.psc')
    os.close(fid)
    try:
        interfaces.convertSBML2PSC(sbmlfile=task.model.source, pscfile=model_filename)
    except Exception as exception:
        os.remove(model_filename)
        raise ValueError('Model at {} could not be imported:\n  {}'.format(
            task.model.source, str(exception).replace('\n', '\n  ')))
    model = pysces.model(model_filename)
    os.remove(model_filename)

    # Load the algorithm specified by `simulation.algorithm.kisao_id`
    sim = task.simulation
    integrator = KISAO_ALGORITHM_MAP.get(sim.algorithm.kisao_id, None)
    if integrator is None:
        raise NotImplementedError("".join([
            "Algorithm with KiSAO id '{}' is not supported. ".format(sim.algorithm.kisao_id),
            "Algorithm must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                '{}: {}'.format(kisao_id, algorithm['id'])
                for kisao_id, algorithm in KISAO_ALGORITHM_MAP.items())),
        ]))
    model.mode_integrator = integrator['id']

    # Apply the algorithm parameter changes specified by `task.simulation.algorithm.changes`
    for change in sim.algorithm.changes:
        setting = integrator['settings'].get(change.kisao_id, None)
        if setting is None:
            raise NotImplementedError("".join([
                "Algorithm parameter with KiSAO id '{}' is not supported. ".format(change.kisao_id),
                "Parameter must have one of the following KiSAO ids:\n  - {}".format('\n  - '.join(
                    '{}: {} ({})'.format(kisao_id, setting['id'], setting['name'])
                    for kisao_id, setting in integrator['settings'].items())),
            ]))

        value = change.new_value
        if not validate_str_value(value, setting['type']):
            raise ValueError("'{}' is not a valid {} value for parameter {}".format(
                value, setting['type'].name, change.kisao_id))

        parsed_value = parse_value(value, setting['type'])
        model.__settings__[setting['id']] = parsed_value

    # override algorithm choice if there are events
    if integrator['id'] == 'LSODA' and model.__events__:
        model.mode_integrator = 'CVODE'
        substitution_policy = get_algorithm_substitution_policy()
        if (
            ALGORITHM_SUBSTITUTION_POLICY_LEVELS[substitution_policy]
            >= ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.SIMILAR_VARIABLES]
        ):
            warnings.warn('CVODE (KISAO_0000019) will be used rather than LSODA (KISAO_0000088) because the model has events',
                          AlgorithmSubstitutedWarning)
        else:
            raise AlgorithmDoesNotSupportModelFeatureException('LSODA cannot execute the simulation because the model has events')

    if model.mode_integrator == 'CVODE':
        model.__settings__['cvode_return_event_timepoints'] = False

    # setup time course
    model.sim_start = sim.initial_time
    model.sim_end = sim.output_end_time
    model.sim_points = (
        sim.number_of_points
        * (sim.output_end_time - sim.initial_time)
        / (sim.output_end_time - sim.output_start_time)
        + 1
    )
    if model.sim_points != int(model.sim_points):
        raise NotImplementedError('Time course must specify an integer number of time points')

    # execute simulation
    model.Simulate()

    # extract results
    variable_results = VariableResults()
    unpredicted_symbols = []
    unpredicted_targets = []
    results, labels = model.data_sim.getAllSimData(lbls=True)
    labels = {label: i_label for i_label, label in enumerate(labels)}

    for variable in variables:
        if variable.symbol:
            if variable.symbol == Symbol.time:
                i_result = labels['Time']
                variable_results[variable.id] = results[:, i_result][-(sim.number_of_points + 1):]
            else:
                unpredicted_symbols.append(variable.symbol)

        else:
            sbml_id = target_x_paths_to_sbml_ids[variable.target]
            i_result = labels.get(sbml_id, None)
            if i_result is None:
                unpredicted_targets.append(variable.target)
            else:
                variable_results[variable.id] = results[:, i_result][-(sim.number_of_points + 1):]

    if unpredicted_symbols:
        raise NotImplementedError("".join([
            "The following variable symbols are not supported:\n  - {}\n\n".format(
                '\n  - '.join(sorted(unpredicted_symbols)),
            ),
            "Symbols must be one of the following:\n  - {}".format(Symbol.time),
        ]))

    if unpredicted_targets:
        raise ValueError(''.join([
            'The following variable targets could not be recorded:\n  - {}\n\n'.format(
                '\n  - '.join(sorted(unpredicted_targets)),
            ),
            'Targets must have one of the following ids:\n  - {}'.format(
                '\n  - '.join(sorted(set(labels.keys()).difference(set(['Time'])))),
            ),
        ]))

    # restore working directory
    os.chdir(cwd)

    # log action
    log.algorithm = 'KISAO_0000019' if model.mode_integrator == 'CVODE' else 'KISAO_0000088'

    arguments = {}
    for key, val in model.__settings__.items():
        if model.mode_integrator == 'CVODE':
            if key.startswith('cvode_'):
                arguments[key] = val
        else:
            if key.startswith('lsoda_'):
                arguments[key] = val

    log.simulator_details = {
        'method': model.Simulate.__module__ + '.' + model.Simulate.__name__,
        'arguments': arguments,
    }

    # return results and log
    return variable_results, log
