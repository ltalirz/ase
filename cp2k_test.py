#!/usr/bin/python
# -*- coding: utf-8 -*-

from ase.test import NotAvailable
from ase.structure import molecule
from ase.calculators.calculator import get_calculator
from ase.units import GPa, Pascal, Hartree
from cp2k import CP2K

def main():
    CP2K.command = "mpiexec -np 4 /data/schuetto/svn/cp2k/cp2k/exe/Linux/cp2k_shell.pdbg"

    #calc = CP2K(label="somelabel", xc='LDA', max_scf=5, debug=False, txt="lala.out", basis_set='SZV-MOLOPT-SR-GTH')
    calc = CP2K(label="somelabel", xc='LDA', max_scf=5, debug=True, txt="lala.out", basis_set='SZV-MOLOPT-SR-GTH')
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    e2 = h2.get_potential_energy()
    print "Energy: [eV] ", e2
    print "Energy: [Ha] ", e2/Hartree
    print h2.get_stress(voigt=False)/GPa
    #return
    calc.set(xc='PBE')
    e2pbe = h2.get_potential_energy()
    h1 = h2.copy()
    del h1[1]
    h1.set_initial_magnetic_moments([1])
    h1.calc = calc
    e1pbe = h1.get_potential_energy()
    calc.set(xc='LDA')
    e1 = h1.get_potential_energy()
    print(2 * e1 - e2)
    print(2 * e1pbe - e2pbe)
    print e1, e2, e1pbe, e2pbe

    #calc.write() # write a restart
    #calc2 = CP2K(name) # load a restart
    #print calc.parameters, calc.results, calc.atoms
    #assert not calc.calculation_required(h1, ['energy'])
    #h1 = calc.get_atoms()
    #
    #print h1.get_potential_energy()
    #label = 'dir/' + name + '-h1'
    #calc = Calculator(label=label, atoms=h1, xc='LDA', **par)
    #
    #print h1.get_potential_energy()
    #print Calculator.read_atoms(label).get_potential_energy()

main()
#EOF
