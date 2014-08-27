# -*- coding: utf-8 -*-
"""This module defines an ASE interface to CP2K.

http://www.cp2k.org
Author: Ole Schütt <ole.schuett@mat.ethz.ch>
"""

import sys, os, re
import numpy as np
from os import path
from subprocess import Popen, PIPE
from tempfile import mkstemp, mktemp

import ase.io
from ase import Atoms
from ase.units import Rydberg, Hartree, Bohr
from ase.calculators.calculator import Calculator, all_changes, Parameters, ReadError


class CP2K(Calculator):
    """Class for doing CP2K calculations."""

    implemented_properties = ['energy', 'forces', 'stress']
    command = 'cp2k_shell.popt'

    default_parameters = dict(
        txt = "-",
        xc = 'LDA',
        basis_set = 'DZVP-MOLOPT-SR-GTH',
        pseudo_potential = 'auto',
        basis_set_file = "BASIS_MOLOPT",
        potential_file = "POTENTIAL",
        max_scf = 50,
        cutoff = 400 * Rydberg,
        charge = 0)

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='project001', atoms=None, debug=False, **kwargs):
        """Construct CP2K-calculator object."""

        self._debug  = debug
        self._force_env_id = None
        self._child  = None

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        # launch cp2k_shell child process
        self._pipe_fn = mktemp(prefix='cp2k_shell_', suffix=".fifo")
        os.mkfifo(self._pipe_fn) 
        cmd = self.command + " --redirect-comm="+self._pipe_fn
        if(self._debug): print cmd
        self._child = Popen(cmd.split(), stdin=PIPE, bufsize=1)
        self._pipe = open(self._pipe_fn, "r")
        assert(self._recv() == "* READY")


    def __del__(self):
        """Release force_env and terminate cp2k_shell child process"""
        self._release_force_env()
        if(self._child):
            self._send("EXIT")
            assert(self._child.wait() == 0) # child process exited properly?
            self._pipe.close()
            os.remove(self._pipe_fn)
            self._child = self._pipe_fn = self._pipe = None


    def reset(self):
        """Clear all information from old calculation."""
        Calculator.reset(self)
        self._release_force_env()


    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2, ...)."""
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()


    def write(self, label):
        """Write atoms, parameters and calculated properties into restart files."""
        self.atoms.write(label + '_restart.traj')
        self.parameters.write(label + '_params.ase')
        open(label+'_results.ase', "w").write(repr(self.results))


    def read(self, label):
        """Read atoms, parameters and calculated properties from restart files."""
        from numpy import array
        self.atoms = ase.io.read(label+'_restart.traj')
        self.parameters = Parameters.read(label + '_params.ase')
        self.results = eval(open(label+'_results.ase').read())


    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """Do the calculation."""

        Calculator.calculate(self, atoms, properties, system_changes)

        if("numbers" in system_changes or "initial_magmoms" in system_changes):
            self._release_force_env()

        if(self._force_env_id is None):
            self._create_force_env()

        n_atoms = len(self.atoms)
        if("cell" in system_changes):
            cell = self.atoms.get_cell()
            self._send("SET_CELL %d"%self._force_env_id)
            self._send(" ".join(["%.10f"%x for x in cell.flat]))
            assert(self._recv() == "* READY")

        if("positions" in system_changes):
            self._send("SET_POS %d"%self._force_env_id)
            self._send("%d"%(3*n_atoms))
            for pos in self.atoms.get_positions():
                self._send("%.10f   %.10f   %.10f"%(pos[0], pos[1], pos[2]))
            self._send("*END")
            assert(float(self._recv()) >= 0) # max change -> ignore
            assert(self._recv() == "* READY")

        self._send("EVAL_EF %d"%self._force_env_id)
        assert(self._recv() == "* READY")

        self._send("GET_E %d"%self._force_env_id)
        self.results['energy'] = float(self._recv()) * Hartree
        assert(self._recv() == "* READY")

        forces = np.zeros(shape=(n_atoms,3) )
        self._send("GET_F %d"%self._force_env_id)
        assert(int(self._recv()) == 3*n_atoms)
        for i in range(n_atoms):
            line = self._recv()
            forces[i,:] = [float(x) for x in line.split()]
        assert(self._recv() == "* END")
        assert(self._recv() == "* READY")
        self.results['forces'] = forces * Hartree / Bohr

        self._send("GET_STRESS %d"%self._force_env_id)
        line = self._recv()
        assert(self._recv() == "* READY")

        stress = np.array([float(x) for x in line.split()]).reshape(3,3)
        assert(np.all(stress == np.transpose(stress))) #should be symmetric
        # Convert 3x3 stress tensor to Voigt form as required by ASE
        stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                           stress[1, 2], stress[0, 2], stress[0, 1]])
        self.results['stress'] = stress * Hartree / Bohr**3

        self.write(self.label) #TODO: this should not be called so often


    def _create_force_env(self):
        assert(self._force_env_id is None)
        inp = self._generate_input()

        if(self._debug):
            print inp
        fd, inp_fn = mkstemp(suffix=".inp")
        f = os.fdopen(fd, "w")
        f.write(inp)
        f.close()

        label_dir = path.dirname(self.label)
        if(len(label_dir)>0 and not path.exists(label_dir)):
            print "Creating directory: "+label_dir
            os.makedirs(label_dir) # cp2k expects dirs to exist
        out_fn = self.parameters.txt
        if(out_fn=="-"): out_fn = "__STD_OUT__"
        self._send("LOAD %s %s"%(inp_fn, out_fn))
        self._force_env_id = int(self._recv())
        assert(self._force_env_id > 0)
        assert(self._recv() == "* READY")
        os.remove(inp_fn)


    def _release_force_env(self):
        if(self._force_env_id):
            self._send("DESTROY %d"%self._force_env_id)
            assert(self._recv() == "* READY")
            self._force_env_id = None


    def _generate_input(self):
        p = self.parameters

        output  = "!!!! Generated by ASE !!!!\n"
        output += "&GLOBAL\n"
        output += "PROJECT %s\n"%self.label
        output += "&END GLOBAL\n"

        output += "&FORCE_EVAL\n"
        output += "&PRINT\n"
        output += "  &STRESS_TENSOR ON\n"
        output += "  &END STRESS_TENSOR\n"
        output += "&END PRINT\n"

        output += "METHOD Quickstep\n"
        output += "STRESS_TENSOR ANALYTICAL\n"


        output += "&SUBSYS\n"

        # determine pseudo-potential
        potential = p.pseudo_potential
        if(p.pseudo_potential.lower() == "auto"):
            if(p.xc.upper() == "LDA"):
                potential = "GTH-PADE"
            elif(p.xc.upper() in ("PADE", "BP", "BLYP", "PBE",)):
                potential = "GTH-"+p.xc.upper()
            else:
                warn("No matching pseudo potential found, falling back to GTH-PBE", RuntimeWarning)
                potential = "GTH-PBE" # fall back

        # write atomic kinds
        valence_electrons = self._parse_basis_set()
        for elem in set(self.atoms.get_chemical_symbols()):
            output += "&KIND %s\n"%elem
            output += "  BASIS_SET %s\n"%p.basis_set
            q = valence_electrons[(elem, p.basis_set)]
            output += "  POTENTIAL %s-q%d\n"%(potential, q)
            output += "&END KIND\n"

        # write coords
        output += "&COORD\n"
        n_electrons = 0
        for elem, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions()):
            n_electrons += valence_electrons[(elem, p.basis_set)]
            output += "  %s  %.10f   %.10f   %.10f\n"%(elem, pos[0], pos[1], pos[2])
        output += "&END COORD\n"

        # write cell
        output += "&CELL\n"
        pbc = "".join([a for a,b in zip("XYZ",self.atoms.get_pbc()) if b])
        if(len(pbc)==0): pbc = "NONE"
        output += "  PERIODIC %s\n"%pbc
        cell = self.atoms.get_cell()
        output += "  A  %.10f   %.10f   %.10f\n"%(cell[0,0], cell[0,1], cell[0,2])
        output += "  B  %.10f   %.10f   %.10f\n"%(cell[1,0], cell[1,1], cell[1,2])
        output += "  C  %.10f   %.10f   %.10f\n"%(cell[2,0], cell[2,1], cell[2,2])
        output += "&END CELL\n"
        output += "&END SUBSYS\n"

        # write DFT-section
        output += "&DFT\n"
        output += "BASIS_SET_FILE_NAME %s\n"%p.basis_set_file
        output += "POTENTIAL_FILE_NAME %s\n"%p.potential_file
        #if(any(atoms.get_initial_magnetic_moments()!=0)):
        if(n_electrons%2 != 0): #TODO is this the proper way?
            output += "SPIN_POLARIZED .TRUE.\n"
        if(p.charge != 0):
            output += "CHARGE %d\n"%p.charge

        output += "&XC\n"
        output += "  &XC_FUNCTIONAL %s\n"%p.xc
        output += "  &END XC_FUNCTIONAL\n"
        output += "&END XC\n"

        output += "&MGRID\n"
        output += "  CUTOFF [eV] %.3f\n"%p.cutoff
        output += "&END MGRID\n"

        output += "&SCF\n"
        output += "  MAX_SCF %d\n"%p.max_scf
        output += "&END SCF\n"

        output += "&END DFT\n"
        output += "&END FORCE_EVAL\n"

        return(output)


    def _parse_basis_set(self):
        """Parses basis_sets and returns number of valence electrons"""
        valence_electrons = {}
        lines = open(self.parameters.basis_set_file).readlines()
        for line in lines:
            if(len(line.strip()) == 0): continue
            if(line.startswith("#")): continue
            if(line.strip()[0].isdigit()): continue

            basis_set_id = tuple(line.split()[:2])
            m = re.search("-q(\d+)",line)
            if(m):
                valence_electrons[basis_set_id] = int(m.group(1))

        return(valence_electrons)


    def _send(self, line):
        """Send a line to the cp2k_shell"""
        assert(self._child.poll() is None) # child process still alive?
        if(self._debug): print("Sending: "+line)
        self._child.stdin.write(line+"\n")


    def _recv(self):
        """Receive a line from the cp2k_shell"""
        assert(self._child.poll() is None) # child process still alive?
        line = self._pipe.readline().strip()
        if(self._debug): print("Received: "+line)
        return(line)

#EOF
