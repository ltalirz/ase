#!/usr/bin/python
# -*- coding: utf-8 -*-

from ase.units import Hartree
from ase.calculators.cp2k import CP2K
from ase.lattice.cubic import FaceCenteredCubic
import os
import numpy as np

def run_cp2k(atoms, label):
    ntask = int(os.environ['PBS_NP'])
    print("ntasks = {}".format(ntask))
    CP2K.command = "mpiexec -np {} /share/apps/cp2k/R14666/exe/hypatia-libxc220/cp2k_shell.popt".format(ntask)
    
    calc = CP2K(xc='PBE', 
                txt=label+'.out',
                vdw='NONE', 
                print_level='LOW',
                label=label     )
    atoms.set_calculator(calc)

    energy = atoms.get_potential_energy() /Hartree
    print("a = {:.2f} Energy [Ha]: {:.6f}".format(a,energy/nk**3))
    f.write("a = {:.2f} Energy [Ha]: {:.6f}\n".format(a,energy/nk**3))
    

nk = 5

f = open('energies', 'a')
f.write('# unit cell,   energy/unit cell ({}x{}x{})\n'.format(nk,nk,nk))

for a in np.arange(4.18,4.23,0.01):
#for a in [4.18, 4.19, 4.20]:
    label = "a_{:.2f}".format(a)

    au = FaceCenteredCubic(
        size=(nk, nk, nk), 
        symbol='Au',
        latticeconstant=a,
    )
    run_cp2k(atoms=au, label=label)

f.close()
    

#    ##calc = CP2K(label="somelabel", xc='LDA', max_scf=5, debug=False, txt="lala.out", basis_set='SZV-MOLOPT-SR-GTH')
