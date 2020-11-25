""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from Biosimulations_utils.simulation.data_model import TimecourseSimulation, SimulationResultsFormat
from Biosimulations_utils.simulator.utils import exec_simulations_in_archive
import os


__all__ = ['exec_combine_archive', 'exec_simulation']


def exec_combine_archive(archive_file, out_dir):
    """ Execute the SED tasks defined in a COMBINE archive and save the outputs

    Args:
        archive_file (:obj:`str`): path to COMBINE archive
        out_dir (:obj:`str`): directory to store the outputs of the tasks
    """
    exec_simulations_in_archive(archive_file, exec_simulation, out_dir, apply_model_changes=True)


def exec_simulation(model_filename, model_sed_urn, simulation, working_dir, out_filename, out_format):
    ''' Execute a simulation and save its results

    Args:
       model_filename (:obj:`str`): path to the model
       model_sed_urn (:obj:`str`): SED URN for the format of the model (e.g., `urn:sedml:language:sbml`)
       simulation (:obj:`TimecourseSimulation`): simulation
       working_dir (:obj:`str`): directory of the SED-ML file
       out_filename (:obj:`str`): path to save the results of the simulation
       out_format (:obj:`SimulationResultsFormat`): format to save the results of the simulation (e.g., `HDF5`)
    '''
    # check that model is encoded in SBML
    if model_sed_urn != "urn:sedml:language:sbml":
        raise NotImplementedError("Model language with URN '{}' is not supported".format(model_sed_urn))

    # check that simulation is a time course simulation
    if not isinstance(simulation, TimecourseSimulation):
        raise NotImplementedError('{} is not supported'.format(simulation.__class__.__name__))

    # check that model parameter changes have already been applied (because handled by :obj:`exec_simulations_in_archive`)
    if simulation.model_parameter_changes:
        raise NotImplementedError('Model parameter changes are not supported')

    # check that the desired output format is supported
    if out_format != SimulationResultsFormat.HDF5:
        raise NotImplementedError("Simulation results format '{}' is not supported".format(out_format))

    # Read the model located at `os.path.join(working_dir, model_filename)` in the format
    # with the SED URN `model_sed_urn`.

    # Apply the model parameter changes specified by `simulation.model_parameter_changes`

    # Load the algorithm specified by `simulation.algorithm`

    # Apply the algorithm parameter changes specified by `simulation.algorithm_parameter_changes`

    # Simulate the model from `simulation.start_time` to `simulation.end_time`

    # Save a report of the results of the simulation with `simulation.num_time_points` time points
    # beginning at `simulation.output_start_time` to `out_filename` in `out_format` format.
    # This should save all of the variables specified by `simulation.model.variables`.

    pass
