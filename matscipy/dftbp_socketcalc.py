# ======================================================================
# matscipy - Python materials science tools
# https://github.com/libAtoms/matscipy
#
# Copyright (2020) Alexander Held
#                  Thomas Reichenbach
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

"""The dftbsocket module.

This python module provides an ASE calculator for DFTB+.
It operates DFTB+ through the i-PI interface protocol, which is implemented in DFTB+.
i-PI can be obtained from https://github.com/i-pi/i-pi/tree/master.
The i-PI dependency is not (yet) compatible with python 3,
so the DFTB+ wrapper will not work with python 3 until i-PI is fixed.

A usage example is found in the example directory.
Please mind that this implementation is NOT WELL TESTED!!!

Make sure to use the same socket specification in the ASE calculator and in the dftb+ input file!
Unix sockets and network sockets are supported (only unix sockets tested so far).
Please also read the Socket Driver section in the DFTB+ user manual.

Known issues:
-Stress tensor is not tested. Sign and/or unit could be wrong
-Not much testing was done, except to see if it runs without errors
-No parallel tests were done
-i-PI officialy only allows for triangular cell matrices, but it seems to not have a hard-coded limit.
 The full cell seems to be passed through the i-PI interface.
-There could be a transpose-error in the cell, I am not sure what the definition of the i-PI cell is
-Network sockets not tested


Author: Alexander Held
"""

from __future__ import (
    division,
    absolute_import,
    print_function,
    unicode_literals
)
from ase.calculators.calculator import Calculator, all_changes
from ipi.engine.forcefields import FFSocket
from ipi.interfaces.sockets import InterfaceSocket
from ipi.utils.mathtools import invert_ut3x3
from ase.units import Hartree, Bohr
from ase import Atoms
import subprocess
import time
import numpy as np
from io import open


class FakeIPIAtoms(object):
    def __init__(self, ase_atoms):
        self.ase_atoms = ase_atoms

    @property
    def q(self):
        return self.ase_atoms.positions.flatten() / Bohr


class FakeIPICell(object):
    def __init__(self, ase_atoms):
        self.ase_atoms = ase_atoms

    def array_pbc(self, ipi_pbcpos):
        ase_atoms = Atoms(
            positions=(ipi_pbcpos * Bohr).reshape(-1, 3),
            cell=self.ase_atoms.get_cell(),
            pbc=self.ase_atoms.get_pbc()
        )
        ase_atoms.wrap()
        ipi_pbcpos[:] = FakeIPIAtoms(ase_atoms).q

    @property
    def h(self):
        # I think we have to transpose, but I'm not quite sure
        h = self.ase_atoms.get_cell().T / Bohr
        # Triangular form seems not to be hard-coded to i-PI,
        # so we do not have to assert this here:
        # assert h[1, 0] == 0, 'triangular form expected'
        # assert h[2, 0] == 0, 'triangular form expected'
        # assert h[2, 1] == 0, 'triangular form expected'
        return h

    @property
    def ih(self):
        return invert_ut3x3(self.h)


class DFTBPlusSocketCalc(Calculator):
    implemented_properties = ['energy', 'forces', 'stress', 'extra']

    # taken from default of ipi.inputs.forcefields.InputFFSocket
    default_parameters = {
        'pars': {},
        'name': 'DFTBPlus',
        'latency': 1.0,  # decrease if you expect less than 1 seconds for a single step
        'dopbc': True,  # wrap atoms to cell before sending to client?
        'active': np.array([-1]),
        'address': 'localhost',
        'port': 65535,
        'slots': 4,
        'mode': 'inet',  # also 'unix' is possible for unix sockets
        'timeout': 0.0,
        'match_mode': 'auto',
        'command': ('dftb+', ),  # could also be s.th. like ('srun', 'dftb+')
        'stdout': 'dftb+.stdout.log',
        'stderr': 'dftb+.stderr.log',
        'logmode': 'a',  # append to log files
        'bufsize': -1,  # -1 is system default buffering
    }
    ffsocket = None
    dftbplus_process = None
    dftbplus_stdout = None
    dftbplus_stderr = None

    def calculate(
            self,
            atoms=None,
            properties=['energy'],
            system_changes=all_changes
    ):
        Calculator.calculate(self, atoms, properties, system_changes)
        if self.initialized and ({'numbers', 'pbc'} & set(system_changes)):
            self.shutdown()
        self._initialize_if_needed()
        request = self.ffsocket.queue(FakeIPIAtoms(atoms), FakeIPICell(atoms))
        while request["status"] != "Done":
            returncode = self.dftbplus_process.poll()
            if returncode is not None:
                self.ffsocket.release(request)
                self.shutdown()  # raises on nonzero return code
                assert returncode == 0, 'Returncode should be zero here. Bug?'
                raise EnvironmentError(
                    'DFTB+ terminated with exit status 0 before handing back result!'
                )
            time.sleep(self.ffsocket.latency)
        ipi_potential, ipi_force, ipi_virial, ipi_extra = request["result"]
        self.ffsocket.release(request)
        self.results = {
            'energy': ipi_potential * Hartree,
            'forces': ipi_force.reshape(-1, 3) * Hartree / Bohr,
            'stress': ipi_virial * Hartree / atoms.get_volume(),  # sign correct???
            'extra': ipi_extra,  # I have no idea, what DFTB+ sends here, maybe magmoms or something
        }

    @property
    def initialized(self):
        return self.ffsocket is not None

    def _initialize_if_needed(self):
        if self.initialized:
            return
        p = self.parameters
        self.ffsocket = FFSocket(
            pars=p['pars'],
            name=p['name'],
            latency=p['latency'],
            dopbc=p['dopbc'],
            active=p['active'],
            interface=InterfaceSocket(
                address=p['address'],
                port=p['port'],
                slots=p['slots'],
                mode=p['mode'],
                timeout=p['timeout'],
                match_mode=p['match_mode']
            ),
        )
        self.ffsocket.run()
        self.dftbplus_stdout = open(p['stdout'], p['logmode'], encoding='utf-8')
        self.dftbplus_stderr = open(p['stderr'], p['logmode'], encoding='utf-8')
        self.dftbplus_process = subprocess.Popen(
            p['command'],
            bufsize=p['bufsize'],
            stdout=self.dftbplus_stdout,
            stderr=self.dftbplus_stderr
        )

    def set(self, **kwargs):
        changed_parameters = Calculator.set(self, **kwargs)
        if bool(changed_parameters):
            self.shutdown()
            self.reset()

    def shutdown(self):
        if not self.initialized:
            return
        try:
            self.ffsocket.stop()
            returncode = self.dftbplus_process.wait()
            self.dftbplus_stdout.close()
            self.dftbplus_stderr.close()
            if returncode != 0:
                raise subprocess.CalledProcessError(returncode, self.parameters['command'][0])
        finally:
            self.ffsocket = None
            self.dftbplus_process = None
            self.dftbplus_stdout = None
            self.dftbplus_stderr = None

    def __del__(self):
        # just in case someone forgot to call shutdown
        self.shutdown()
