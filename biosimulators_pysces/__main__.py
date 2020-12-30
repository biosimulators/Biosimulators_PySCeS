""" BioSimulators-compliant command-line interface to the `PySCeS <http://pysces.sourceforge.net/>`_ simulation program.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-12-23
:Copyright: 2020, BioSimulators Team
:License: MIT
"""

from ._version import __version__
from .core import exec_sedml_docs_in_combine_archive
from biosimulators_utils.simulator.cli import build_cli
from biosimulators_utils.simulator.data_model import AlgorithmSubstitutionPolicy
from biosimulators_utils.simulator.environ import ENVIRONMENT_VARIABLES
import os
cwd = os.getcwd()  # because PySCeS changes the working directory
import pysces  # noqa: E402
os.chdir(cwd)

App = build_cli('pysces', __version__,
                'PySCeS', pysces.__version__, 'http://pysces.sourceforge.net/',
                exec_sedml_docs_in_combine_archive,
                environment_variables=[
                    ENVIRONMENT_VARIABLES[AlgorithmSubstitutionPolicy]
                ])


def main():
    with App() as app:
        app.run()
