""" BioSimulators-compliant command-line interface to the `PySCeS <http://pysces.sourceforge.net/>`_ simulation program.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-12-23
:Copyright: 2020, BioSimulators Team
:License: MIT
"""

from . import get_simulator_version
from ._version import __version__
from .core import exec_sedml_docs_in_combine_archive
from biosimulators_utils.simulator.cli import build_cli

App = build_cli('biosimulators-pysces', __version__,
                'PySCeS', get_simulator_version(), 'http://pysces.sourceforge.net/',
                exec_sedml_docs_in_combine_archive)


def main():
    with App() as app:
        app.run()
