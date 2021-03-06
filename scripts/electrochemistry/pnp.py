#!/usr/bin/env python3
# ======================================================================
# matscipy - Python materials science tools
# https://github.com/libAtoms/matscipy
#
# Copyright (2019) Johannes Hoermann, University of Freiburg
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
"""Solves 1D Poisson-Nernst-Planck system

Copyright 2019 IMTEK Simulation
University of Freiburg

Authors:
  Johannes Hoermann <johannes.hoermann@imtek-uni-freiburg.de>
"""
import datetime, logging, os, sys
import numpy as np

from matscipy.electrochemistry import PoissonNernstPlanckSystem

def main():
    """Solve Poisson-Nernst-Planck system and store distribution.

    Specify quantities in SI units at this command line interface."""

    import argparse

    # in order to have both:
    # * preformatted help text and ...
    # * automatic display of defaults
    class ArgumentDefaultsAndRawDescriptionHelpFormatter(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class = ArgumentDefaultsAndRawDescriptionHelpFormatter)

    parser.add_argument('outfile', metavar='OUT', default=None, nargs='?',
                        help='binary numpy .npz or plain text .txt output file')

    # physical system parameters
    #parser.add_argument('--box','-b', default=[50.0e-10,50.0e-10,100.0e-10],
    #                    nargs=3, required=False, type=float, dest="box",
    #                    metavar=('X','Y','Z'), help='Box dimensions (m)')

    parser.add_argument('--concentrations','-c',
                        default=[0.1,0.1], type=float, nargs='+',
                        metavar='c', required=False, dest="concentrations",
                        help='Ion species concentrations c (mol m^-3, mM)')

    parser.add_argument('--charges','-z',
                        default=[1,-1], type=float, nargs='+',
                        metavar='z', required=False, dest="charges",
                        help='Ion species number charges z')

    parser.add_argument('--potential','-u',
                        default=0.05, type=float,
                        metavar='U', required=False, dest="potential",
                        help='Potential drop from left to right dU (V)')

    parser.add_argument('--length','-l',
                        default=100.0e-9, type=float,
                        metavar='L', required=False, dest="length",
                        help='Domain length (m)')

    parser.add_argument('--temperature','-T',
                        default=298.15, type=float,
                        metavar='T', required=False, dest="temperature",
                        help='Temperature (K)')

    parser.add_argument('--relative-permittivity','--epsilon-r', '--eps',
                        default=79.0, type=float,
                        metavar='eps', required=False,
                        dest="relative_permittivity",
                        help='Relative permittivity')

    parser.add_argument('--compact-layer','--stern-layer', '--lambda-s',
                        default=0.0, type=float,
                        metavar='L', required=False,
                        dest="lambda_S",
                        help='Stern or compact layer thickness (for Robin BC)')

    parser.add_argument('--boundary-conditions','-bc',
                        default='cell', type=str,
                        metavar='BC', required=False,
                        dest="boundary_conditions",
                        choices=(
                            'interface', # open half-space
                            'cell', # 1D electorchemical cell with zero flux BC
                            'cell-stern', # 1D cell with linear compact layer regime
                            'cell-stern-explicit', # same as cell-stern
                            'cell-robin', # 1D cell with implict compact layer by Robin BC
                            'cell-stern-implicit', # same as cell-robin
                            ),
                        help='Boundary conditions')

    # technical settings
    parser.add_argument('--segments','-N',
                        default=200, type=int,
                        metavar='N', required=False,
                        dest="segments",
                        help='Number of discretization segments')

    parser.add_argument('--maximum-iterations','--maxit',
                        default=20, type=int,
                        metavar='N', required=False,
                        dest="maxit",
                        help='Maximum number of Newton iterations')

    parser.add_argument('--absolute-tolerance',
                        default=1e-8, type=float,
                        metavar='e', required=False,
                        dest="absolute_tolerance",
                        help='Absolute tolerance Newton solver convergence criterion')

    # output settings
    # parser.add_argument('--output-interval',
    #                     default=1, type=int,
    #                     metavar='N', required=False,
    #                     dest="outinterval",
    #                     help='Print log messages every Nth Newton step')

    parser.add_argument('--convergence-stats', default=False, required=False,
                        action='store_true', dest="convergence",
                        help='Record and store Newton solver convergence statistics')

    parser.add_argument('--debug', default=False, required=False,
                        action='store_true', dest="debug", help='debug flag')
    parser.add_argument('--verbose', default=False, required=False,
                        action='store_true', dest="verbose", help='verbose flag')
    parser.add_argument('--log', required=False, nargs='?', dest="log",
                        default=None, const='pnp.log', metavar='LOG',
                        help='Write log file pnp.log, optionally specify log file name')

    try:
        import argcomplete
        argcomplete.autocomplete(parser)
        # This supports bash autocompletion. To enable this, pip install
        # argcomplete, activate global completion, or add
        #      eval "$(register-python-argcomplete lpad)"
        # into your .bash_profile or .bashrc
    except ImportError:
        pass

    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.verbose:
        loglevel = logging.INFO
    else:
        loglevel = logging.WARNING

    # PoissonNernstPlanckSystem makes extensive use of Python's logging module

    # logformat  = ''.join(("%(asctime)s",
    #  "[ %(filename)s:%(lineno)s - %(funcName)s() ]: %(message)s"))
    logformat  = "[ %(filename)s:%(lineno)s - %(funcName)s() ]: %(message)s"

    logging.basicConfig(level=loglevel,
                        format=logformat)

    # explicitly modify the root logger (necessary?)
    logger = logging.getLogger()
    logger.setLevel(loglevel)

    # remove all handlers
    for h in logger.handlers: logger.removeHandler(h)

    # create and append custom handles
    ch = logging.StreamHandler()
    formatter = logging.Formatter(logformat)
    ch.setFormatter(formatter)
    ch.setLevel(loglevel)
    logger.addHandler(ch)

    if args.log:
        fh = logging.FileHandler(args.log)
        fh.setFormatter(formatter)
        fh.setLevel(loglevel)
        logger.addHandler(fh)

    # set up system
    pnp = PoissonNernstPlanckSystem(
        c =         np.array(args.concentrations, dtype=float),
        z =         np.array(args.charges, dtype=float),
        L =         float(args.length),
        T =         float(args.temperature),
        delta_u =   float(args.potential),
        lambda_S =  float(args.lambda_S),
        relative_permittivity = float(args.relative_permittivity) )

    # technical settings
    # pnp.output  = args.convergence # makes Newton solver display convergence plots
    pnp.N       = args.segments # uniformlya distanced grid points
    pnp.maxit   = args.maxit # maximum number of Newton iterations
    # pnp.outfreq = args.outinterval
    pnp.e       = args.absolute_tolerance # absolute tolerance

    if args.boundary_conditions == 'cell':
        pnp.useStandardCellBC()
    elif args.boundary_conditions in ('cell-stern','cell-stern-explicit'):
        pnp.useSternLayerCellBC(implicit=False)
    elif args.boundary_conditions in ('cell-robin','cell-stern-implicit'):
        pnp.useSternLayerCellBC(implicit=True)
    elif args.boundary_conditions == 'interface':
        pnp.useStandardInterfaceBC()
    else:
        raise ValueError("Boundary conditions '{}' not implemented!".format(
            args.boundary_conditions ))

    pnp.init()
    pnp.solve()

    extra_kwargs = {}
    if args.convergence:
        extra_kwargs.update( {
            'convergence_step_absolute': pnp.convergenceStepAbsolute,
            'convergence_step_relative': pnp.convergenceStepRelative,
            'convergence_residual_absolute': pnp.convergenceResidualAbsolute } )


    if not args.outfile:
        outfile = sys.stdout
        format  = 'txt'
    else:
        outfile = args.outfile
        _, format = os.path.splitext(outfile)

    if format == '.npz':
        np.savez(file=outfile,
            x=pnp.grid, u=pnp.potential, c=pnp.concentration, **extra_kwargs)
    else: # elif format == '.txt'
        comment = '\n'.join((
            'Poisson-Nernst-Planck system, generated on {:s}, {:s}'.format(
                os.uname().nodename,
                datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') ),
            'All quantities in SI units.',
            'grid x (m), potential u (V), concentration c (mM)'))
        header = comment + '\n' + '{:20s} {:22s} '.format('x','u') + ' '.join(
            [ '{:22s}'.format('c{:02d}'.format(k)) for k in range(pnp.M)])
        data = np.column_stack([pnp.grid, pnp.potential, pnp.concentration.T])
        np.savetxt(outfile, data, fmt='%22.15e', header=header)

    # write out final state as usual, but mark process failed if not converged
    if not pnp.converged:
        sys.exit(1)

if __name__ == '__main__':
    # Execute everything else
    main()
