#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Test suit for the CP2K ASE calulator.

http://www.cp2k.org
Author: Ole Schütt <ole.schuett@mat.ethz.ch>
"""

from ase.test import NotAvailable
from ase.structure import molecule
from ase.calculators.calculator import get_calculator
from ase.units import GPa, Pascal, Hartree
from ase.optimize import BFGS

from cp2k import CP2K

#===============================================================================
def test_H2_LDA():
    calc = CP2K(label="test_H2_LDA")
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    energy = h2.get_potential_energy()/Hartree
    diff = abs((energy + 0.932722287302)/energy)
    assert(diff < 1e-10)
    print("passed test 'H2_LDA'")

#===============================================================================
def test_H2_PBE():
    calc = CP2K(xc='PBE', label="test_H2_PBE")
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    energy = h2.get_potential_energy()/Hartree
    diff = abs((energy + 0.961680073441)/energy)
    assert(diff < 1e-10)
    print("passed test 'H2_PBE'")

#===============================================================================
def test_H2_LS():
    inp = """&FORCE_EVAL
               &DFT
                 &QS
                   LS_SCF ON
                 &END QS
               &END DFT
             &END FORCE_EVAL"""
    calc = CP2K(label="test_H2_LS", inp=inp)
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    energy = h2.get_potential_energy()/Hartree
    diff = abs((energy + 0.932722212414)/energy)
    assert(diff < 1e-10)
    print("passed test 'H2_LS'")

#===============================================================================
def test_O2():
    calc = CP2K(label="test_O2")
    o2 = molecule('O2', calculator=calc)
    o2.center(vacuum=2.0)
    energy = o2.get_potential_energy()/Hartree
    diff = abs((energy + 29.669802288970264)/energy)
    assert(diff < 1e-10)
    print("passed test 'O2'")

#===============================================================================
def test_restart():
    calc = CP2K()
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    h2.get_potential_energy()
    calc.write("test_restart") # write a restart
    calc2 = CP2K(restart="test_restart") # load a restart
    assert not calc2.calculation_required(h2, ['energy'])
    print("passed test 'restart'")

#===============================================================================
def test_geopt():
    calc = CP2K(label="test_geopt")
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    dyn = BFGS(h2)
    dyn.run(fmax=0.05)
    dist = h2.get_distance(0, 1)
    diff = abs(dist - 1.36733746519)
    assert(diff < 1e-10)
    print("passed test 'geopt'")

#===============================================================================
def main():
    CP2K.command = "mpiexec -np 2 ./cp2k_shell.pdbg"

    test_H2_LDA()
    test_H2_PBE()
    test_H2_LS()
    test_O2()
    test_restart()
    test_geopt()

main()
#EOF
