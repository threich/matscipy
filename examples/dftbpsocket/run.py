from __future__ import (
    division,
    absolute_import,
    print_function,
    unicode_literals
)
from ase.io import read
from ase.optimize import FIRE
from matscipy.dftbp_socketcalc import DFTBPlusSocketCalc

atoms = read('dftb_in.hsd')  # ASE can read atoms from DFTB+ input file
atoms.center()
# Socket (address for mode "unix" or address and port for mode "inet") has to match what is written in DFTB+ input file!
# See also Socket driver in DFTB+ manual.
atoms.calc = DFTBPlusSocketCalc(
    mode='unix',
    address='zundel',
    latency=1e-2,  # should be smaller than the expected time for a single geometry step
    timeout=6e2,  # value taken from example
    command=('/PATH/TO/DFTBPlus/bin/dftb+', ),  # command tuple as interpreted by python's subprocess.popen
)

try:
    opt = FIRE(atoms, trajectory='relax.traj')
    opt.run(fmax=0.01)
finally:
    # shutdown is also called through atoms.calc.__del__ of garbage collector, but we do it explicitely here
    atoms.calc.shutdown()
