#! /usr/bin/env python

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

from __future__ import (
    division,
    absolute_import,
    print_function,
    unicode_literals
)
import sys
import unittest
import matscipytest


class TestSlidingP(matscipytest.MatSciPyTestCase):
    @unittest.skipIf("ipi" not in sys.modules or sys.version_info.major >= 3,
                     'dftb+ socket requires python2 and i-pi')
    def test_instantiate(self):
        from matscipy.dftbp_socketcalc import DFTBPlusSocketCalc
        DFTBPlusSocketCalc(
            mode='unix',
            address='foo',
            latency=1e-2,
            timeout=6e2,
            command=('/does/not/exist______', ),
        )

if __name__ == '__main__':
    unittest.main()
