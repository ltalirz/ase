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


#===============================================================================
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
        charge = 0,
        print_level = 'MEDIUM',
        vdw = 'NONE',
        pair_potential = 'DFTD3',
        vdw_parameter_file = 'dftd3.dat',
    )


    #---------------------------------------------------------------------------
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='project001', atoms=None, debug=False, save_input=True, **kwargs):
        """Construct CP2K-calculator object."""

        self._debug  = debug
        self._save_input = save_input
        self._force_env_id = None
        self._child  = None

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        # launch cp2k_shell child process
        #self._pipe_fn = mktemp(prefix='cp2k_shell_', suffix=".fifo")
        self._pipe_fn = mktemp(prefix='cp2k_shell_', suffix=".fifo", dir=os.getcwd())
        os.mkfifo(self._pipe_fn)
        cmd = self.command + " --redirect-comm="+self._pipe_fn
        if(self._debug): print cmd
        self._child = Popen(cmd.split(), stdin=PIPE, bufsize=1)
        self._pipe = open(self._pipe_fn, "r")
        assert(self._recv() == "* READY")


    #---------------------------------------------------------------------------
    def __del__(self):
        """Release force_env and terminate cp2k_shell child process"""
        self._release_force_env()
        if(self._child):
            self._send("EXIT")
            assert(self._child.wait() == 0) # child process exited properly?
            self._pipe.close()
            os.remove(self._pipe_fn)
            self._child = self._pipe_fn = self._pipe = None


    #---------------------------------------------------------------------------
    def reset(self):
        """Clear all information from old calculation."""
        Calculator.reset(self)
        self._release_force_env()


    #---------------------------------------------------------------------------
    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2, ...)."""
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()


    #---------------------------------------------------------------------------
    def write(self, label):
        """Write atoms, parameters and calculated properties into restart files."""
        self.atoms.write(label + '_restart.traj')
        self.parameters.write(label + '_params.ase')
        open(label+'_results.ase', "w").write(repr(self.results))


    #---------------------------------------------------------------------------
    def read(self, label):
        """Read atoms, parameters and calculated properties from restart files."""
        from numpy import array
        self.atoms = ase.io.read(label+'_restart.traj')
        self.parameters = Parameters.read(label + '_params.ase')
        self.results = eval(open(label+'_results.ase').read())


    #---------------------------------------------------------------------------
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


    #---------------------------------------------------------------------------
    def _create_force_env(self):
        assert(self._force_env_id is None)
        inp = self._generate_input()

        if(self._debug):
            print inp
        if(self._save_input):
            fname = self.parameters.txt + '.inp'
            f = open(fname, "w")
            f.write(inp)
            f.close()
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


    #---------------------------------------------------------------------------
    def _release_force_env(self):
        if(self._force_env_id):
            self._send("DESTROY %d"%self._force_env_id)
            assert(self._recv() == "* READY")
            self._force_env_id = None


    #---------------------------------------------------------------------------
    def _generate_input(self):
        p = self.parameters

        root = parse_input("")
        root.add_keyword("GLOBAL", "PROJECT "+self.label)
        root.add_keyword("GLOBAL", "PRINT_LEVEL "+p.print_level)
        root.add_keyword("FORCE_EVAL", "METHOD Quickstep")
        root.add_keyword("FORCE_EVAL", "STRESS_TENSOR ANALYTICAL")
        root.add_keyword("FORCE_EVAL/PRINT/STRESS_TENSOR", "_SECTION_PARAMETERS_ ON")
        root.add_keyword("FORCE_EVAL/DFT", "BASIS_SET_FILE_NAME "+p.basis_set_file)
        root.add_keyword("FORCE_EVAL/DFT", "POTENTIAL_FILE_NAME "+p.potential_file)
        root.add_keyword("FORCE_EVAL/DFT/XC/XC_FUNCTIONAL", "_SECTION_PARAMETERS_ "+p.xc)
        root.add_keyword("FORCE_EVAL/DFT/MGRID", "CUTOFF [eV] %.3f"%p.cutoff)
        root.add_keyword("FORCE_EVAL/DFT/SCF",  "MAX_SCF %d"%p.max_scf)
        root.add_keyword("FORCE_EVAL/DFT/SCF/OT",  "PRECONDITIONER FULL_SINGLE_INVERSE")
        root.add_keyword("FORCE_EVAL/DFT/SCF/OT",  "MINIMIZER CG")

        if(any(self.atoms.get_initial_magnetic_moments()!=0)):
            raise(NotImplementedError("Magmetic moments not implemented"))

        if(p.charge != 0):
            root.add_keyword("FORCE_EVAL/DFT", "CHARGE %d"%p.charge)

        if(root.get_subsection("FORCE_EVAL/SUBSYS")):
            raise(Exception("Section SUBSYS exists already"))

        valence_electrons = self._parse_basis_set()

        root.add_keyword("FORCE_EVAL/DFT/XC/VDW_POTENTIAL", "POTENTIAL_TYPE "+p.vdw)
        root.add_keyword("FORCE_EVAL/DFT/XC/VDW_POTENTIAL/PAIR_POTENTIAL", "TYPE "+p.pair_potential)
        root.add_keyword("FORCE_EVAL/DFT/XC/VDW_POTENTIAL/PAIR_POTENTIAL", "PARAMETER_FILE_NAME "+p.vdw_parameter_file)
        root.add_keyword("FORCE_EVAL/DFT/XC/VDW_POTENTIAL/PAIR_POTENTIAL", "REFERENCE_FUNCTIONAL PBE")


        # write coords
        n_electrons = 0
        for elem, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions()):
            n_electrons += valence_electrons[(elem, p.basis_set)]
            line = "%s  %.10f   %.10f   %.10f"%(elem, pos[0], pos[1], pos[2])
            root.add_keyword("FORCE_EVAL/SUBSYS/COORD", line, unique=False)

        if(n_electrons%2 != 0): #TODO is this the proper way?
            root.add_keyword("FORCE_EVAL/DFT", "SPIN_POLARIZED .TRUE.")

        # write cell
        pbc = "".join([a for a,b in zip("XYZ",self.atoms.get_pbc()) if b])
        if(len(pbc)==0): pbc = "NONE"
        root.add_keyword("FORCE_EVAL/SUBSYS/CELL", "PERIODIC "+pbc)
        cell = self.atoms.get_cell()
        for i, a in enumerate("ABC"):
            line = "%s  %.10f   %.10f   %.10f"%(a, cell[i,0], cell[i,1], cell[i,2])
            root.add_keyword("FORCE_EVAL/SUBSYS/CELL", line)

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
        subsys = root.get_subsection("FORCE_EVAL/SUBSYS")
        for elem in set(self.atoms.get_chemical_symbols()):
            s = InputSection(name="KIND", params=elem)
            s.keywords.append("BASIS_SET "+p.basis_set)
            q = valence_electrons[(elem, p.basis_set)]
            s.keywords.append("POTENTIAL %s-q%d"%(potential, q))
            subsys.subsections.append(s)

        output_lines = ["!!! Generated by ASE !!!"] + root.write()
        return("\n".join(output_lines))


    #---------------------------------------------------------------------------
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


    #---------------------------------------------------------------------------
    def _send(self, line):
        """Send a line to the cp2k_shell"""
        assert(self._child.poll() is None) # child process still alive?
        if(self._debug): print("Sending: "+line)
        self._child.stdin.write(line+"\n")


    #---------------------------------------------------------------------------
    def _recv(self):
        """Receive a line from the cp2k_shell"""
        assert(self._child.poll() is None) # child process still alive?
        line = self._pipe.readline().strip()
        if(self._debug): print("Received: "+line)
        return(line)


#===============================================================================
class InputSection(object):
    def __init__(self, name, params=None):
        self.name = name.upper()
        self.params = params
        self.keywords = []
        self.subsections = []


    #---------------------------------------------------------------------------
    def write(self):
        output = []
        for k in self.keywords:
            output.append(k)
        for s in self.subsections:
            if(s.params):
                output.append("&%s %s"%(s.name, s.params))
            else:
                output.append("&%s"%s.name)
            for l in s.write():
                output.append("   %s"%l)
            output.append("&END %s"%s.name)
        return(output)


    #---------------------------------------------------------------------------
    def add_keyword(self, path, line, unique=True):
        parts = path.upper().split("/",1)
        candidates = [s for s in self.subsections if s.name==parts[0]]
        if(len(candidates) == 0):
            s = InputSection(name=parts[0])
            self.subsections.append(s)
            candidates = [s]
        elif(len(candidates) != 1):
            raise(Exception("Multiple %s sections found "%parts[0]))

        key = line.split()[0].upper()
        if(len(parts) > 1):
            candidates[0].add_keyword(parts[1], line, unique)
        elif(key == "_SECTION_PARAMETERS_"):
            if(candidates[0].params != None):
                raise(Exception("Section parameter of section %s already set"%(parts[0])))
            candidates[0].params = line.split(" ",1)[1].strip()
        else:
            existing_keys = [k.split()[0].upper() for k in candidates[0].keywords]
            if(unique and key in existing_keys):
                raise(Exception("Keyword %s already present in section %s"%(key, parts[0])))
            candidates[0].keywords.append(line)


    #---------------------------------------------------------------------------
    def get_subsection(self, path):
        parts = path.upper().split("/",1)
        candidates = [s for s in self.subsections if s.name==parts[0]]
        if(len(candidates) > 1):
            raise(Exception("Multiple %s sections found "%parts[0]))
        if(len(candidates) == 0):
            return(None)
        if(len(parts) == 1):
            return(candidates[0])
        return(candidates[0].get_subsection(parts[1]))


#===============================================================================
def parse_input(inp):
    root_section = InputSection("CP2K_INPUT")
    section_stack = [root_section]

    for line in inp.split("\n"):
        line = line.split("!", 1)[0].strip()
        if(len(line)==0): continue

        if(line.upper().startswith("&END")):
            s = section_stack.pop()
        elif(line[0] == "&"):
            parts = line.split(" ", 1)
            name = parts[0][1:]
            if(len(parts) > 1):
                s = InputSection(name=name, params=parts[1].strip())
            else:
                s = InputSection(name=name)
            section_stack[-1].subsections.append(s)
            section_stack.append(s)
        else:
            section_stack[-1].keywords.append(line)

    return(root_section)

#EOF
