#!/usr/bin/env python
# encoding: utf-8

"""
Object-oriented abstractions for GAUSSIAN input files
"""

from __future__ import print_function, division
# Python stdlib
from string import ascii_letters
import sys
import os
from collections import namedtuple
from datetime import datetime
# Chimera stuff
# Additional 3rd parties
# Own

if sys.version_info.major == 3:
    basestring = str

# QM_BASIS_SETS_EXT = ['d,f-Diffuse (+)', 'p,d,f-Diffuse (++)', 'd,f-Polarization (*)', 'p,d,f-Polarization (**)']
MM_FORCEFIELDS = {
    'General': ['Amber', 'Dreiding', 'UFF'],
    'Water': ['TIP3P']
}
MEM_UNITS = ('KB', 'MB', 'GB', 'TB', 'KW', 'MW', 'GW', 'TW')
JOB_TYPES = ('SP', 'Opt', 'IRC', 'IRCMax', 'Scan', 'Freq', 'Polar', 'ADMP',
             'Force', 'Stable', 'Volume')
JOB_OPTIONS = {
    'SP': (),
    'Opt': ('Min', 'TS'),
    'IRC': (),
    'IRCMax': (),
    'Scan': (),
    'Freq': ('Raman', 'NRaman', 'NNRaman', 'NoRaman', 'VCD', 'ROA'),
    'Polar': (),
    'ADMP': (),
    'Force': (),
    'Stable': (),
    'Volume': ()
}
QM_METHODS = ('AM1', 'PM3', 'PM3MM', 'PM6', 'PDDG', 'HF', 'DFT', 'CASSCF', 'MP2',
              'MP3', 'MP4(SDQ)', 'MP4(SDTQ)', 'MP5', 'QCISD', 'CCD', 'CCSD',
              'QCISD(T)', 'QCISD(TQ)', 'BD', 'EPT', 'CBS', 'W1', 'CIS', 'TD',
              'EOM', 'ZINDO', 'DFTB', 'CI', 'GVB', 'G1', 'G2', 'G2MP2', 'G3',
              'G3MP2', 'G3B3', 'G3MP2B3', 'G4', 'G4MP2')
QM_BASIS_SETS = ('STO-3G', '3-21G', '6-21G', '4-31G', '6-31G', "6-31G(d')",
                 "6-31G(d',p')", '6-311G', 'D95V', 'D95', 'SHC', 'CEP-4G',
                 'CEP-31G', 'CEP-121G', 'LanL2MB', 'LanL2DZ', 'SDD', 'SDDAll',
                 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z', 'cc-pV6Z')
QM_BASIS_SETS_EXT = ('', '+', '++', '*', '**')
QM_BASIS_SETS_WITH_ARGS = ('SDD', 'SHC')
QM_FUNCTIONALS = {
    'Pure':      ('VSXC', 'HCTH', 'HCTH93', 'HCTH147', 'HCTH407', 'tHCTH',
                  'M06L', 'B97D', 'B97D3', 'SOGGA11', 'M11L', 'N12', 'MN12L'),
    'Hybrid':    ('B3LYP', 'B3P86', 'B3PW91', 'B1B95', 'mPW1PW91', 'mPW1LYP',
                  'mPW1PBE', 'mPW3PBE', 'B98', 'B971', 'B972', 'PBE1PBE', 'B1LYP',
                  'O3LYP', 'BHandH', 'BHandHLYP', 'BMK', 'M06', 'M06HF', 'M062X',
                  'tHCTHhyb', 'APFD', 'APF', 'SOGGA11X', 'PBEh1PBE', 'TPSSh', 'X3LYP'),
    'RS hybrid': ('HSEH1PBE', 'OHSE2PBE', 'OHSE1PBE', 'wB97XD', 'wB97', 'wB97X',
                  'LC-wPBE', 'CAM-B3LYP', 'HISSbPBE', 'M11', 'N12SX', 'MN12SX')
}
QM_FUNCTIONALS_ALL = set(f for v in QM_FUNCTIONALS.values() for f in v)


class GaussianInputFile(object):

    """
    Object-oriented abstraction of a Gaussian input file.

    Represents an input file, section by section. Almost everything
    can be set directly from the class initialization, but some
    parameters will need manual assignment with special methods.

    Implemented sections are:

    - Link 0 (% commands)
    - Route (# lines)
    - Title
    - Molecule specification (see `pygaussian.GaussianAtom`)
    - Restraints (ModRedundant)
    - Extra basis sets

    It also offers some support for QM/MM jobs.

    More info about Gaussian input formats can be consulted
    in http://gaussian.com/input/.
    """

    def __init__(self, title='Untitled job', *args, **kwargs):
        self.title = title
        self._route = {}
        self._link = {}
        self._job = None
        self._charge = None
        self._multiplicity = None
        self._qm_method = None
        self._qm_functional = None
        self._qm_basis_set = None
        self._qm_basis_sets_extra = []
        self._mm_forcefield = None
        self._mm_forcefield_extra = None
        self._atoms = []
        self._restraints = []

        # Set and verify
        for k, v in kwargs.items():
            if hasattr(self, k) or isinstance(getattr(type(self), k, None), property):
                try:
                    setattr(self, k, v)
                except (TypeError, ValueError) as e:
                    print('! Could not set {} with value'
                          ' {} because {}'.format(k, v, e))
            else:
                print('! Keyword {} not recognized'.format(k))

    def __str__(self):
        return '\n'.join(self.build(joined=False)[1:])

    def build(self, joined=True):
        sections = [self.timestamp]
        # Link 0
        link = self.link
        if link:
            sections.append(link)
        # Route
        sections.append(self.route)
        sections.append('')
        # Title
        sections.append(self.title)
        sections.append('')
        # Charge, spin and atoms
        sections.append(self.system)
        sections.append('')
        # Variables and configuration
        restraints = self.restraints
        if restraints:
            sections.append(restraints)
            sections.append('')
        qm_basis_set_extra = self.qm_basis_set_extra
        if qm_basis_set_extra:
            sections.append(qm_basis_set_extra)
            sections.append('')
        sections.append('')
        sections.append('')
        if joined:
            return '\n'.join(sections)
        return sections
    
    @property
    def timestamp(self):
        now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        return '! Generated by pygaussian at ' + now
    
    # Link 0 section
    @property
    def link(self):
        s = []
        for keyword, options in self._link.items():
            if not options:
                s.append('%' + keyword)
            elif isinstance(options, (tuple, list)):
                s.append('%{}={}'.format(keyword, ','.join(options)))
            else:
                s.append('%{}={}'.format(keyword, options))
        return '\n'.join(s)

    def add_link_option(self, keyword, *options):
        self._link[keyword] = options

    @property
    def processors(self):
        return self._link.get('nprocshared')

    @processors.setter
    def processors(self, value):
        if value < 1:
            raise ValueError('processors must be greater than 0')
        self._link['nprocshared'] = int(value)

    @property
    def memory(self):
        return self._memory

    @memory.setter
    def memory(self, value):
        mem_units = 'GB'
        if isinstance(value, (list, tuple)):
            if len(value) == 1:
                value = value
            elif len(value) == 2:
                value, mem_units = value
            else:
                ValueError('Set memory with: <value>[, "units"]')
        if value <= 0:
            raise ValueError('Memory value must be greater than 0')
        if mem_units not in MEM_UNITS:
            raise ValueError('Memory unit not recognized.')
        self._memory = int(value), mem_units
        self._link['mem'] = '{}{}'.format(*self._memory)

    @property
    def checkpoint(self):
        return self._link.get('chk')

    @checkpoint.setter
    def checkpoint(self, value):
        # if not os.path.isfile(value):
        #     raise ValueError('File {} does not exist'.format(value))
        self._link['chk'] = value

    ########################################################
    @property
    def route(self):
        s = ['#p', self.modeling]
        for keyword, options in self._route.items():
            if not options:
                s.append(keyword)
            elif len(options) == 1:
                s.append('{}={}'.format(keyword, *options))
            else:
                s.append('{}=({})'.format(keyword, ','.join(options)))
        return ' '.join(s)

    def add_route_option(self, keyword, *options):
        self._route[keyword] = options

    @property
    def modeling(self):
        method = self.qm_method
        if method == 'DFT':
            method = self.qm_functional
        if not method:
            raise ValueError('Method must be set with `qm_method`. If '
                             'using DFT, `qm_functional` is also needed.')
        basis = 'gen' if self._qm_basis_sets_extra else self.qm_basis_set
        if not basis:
            raise ValueError('Basis set must be set with `qm_basis_set` '
                             'or `qm_basis_sets_extra`')
        forcefield = self.mm_forcefield
        if forcefield:
            if basis == 'gen':  # QMMM jobs usually specify ECP too
                basis = 'genecp'
            if self._mm_forcefield_extra:
                return 'oniom=({}/{}:{}/hardfirst)'.format(method, basis, forcefield)
            return'oniom=({}/{}:{})'.format(method, basis, forcefield)
        return '{}/{}'.format(method, basis)

    # Job type
    @property
    def job(self):
        job = self._job
        if job is None:
            raise ValueError('Job is not set!')
        return job

    @job.setter
    def job(self, value):
        if value not in JOB_TYPES:
            raise ValueError('Job must be either of {}'.format(JOB_TYPES))
        self._job = value
        self._route[value] = ()

    @property
    def job_options(self):
        return self._route.get(self._job, ())

    @job_options.setter
    def job_options(self, value):
        if not isinstance(value, (tuple, list)):
            value = [value]
        self._route[self._job] = value

    @property
    def freq(self):
        if 'freq' in self._route:
            freq_ = self._route.get('freq')
            if freq_:
                return freq_
            return True
        return False

    @freq.setter
    def freq(self, value):
        if self._job not in ('Opt', 'Polar'):
            raise ValueError('Freq can only be set if job is Opt or Polar.')
        if isinstance(value, basestring):
            self._route['freq'] = value
        elif value is True:
            self._route['freq'] = ()
        else:
            raise ValueError('freq must be str or bool')

    # QM model
    @property
    def qm_method(self):
        return self._qm_method

    @qm_method.setter
    def qm_method(self, value):
        if value not in QM_METHODS:
            raise ValueError('Method must be either of: '
                             '{}'.format(', '.join(QM_METHODS)))
        self._qm_method = value

    @property
    def qm_functional(self):
        return self._qm_functional

    @qm_functional.setter
    def qm_functional(self, value):
        if self.qm_method != 'DFT':
            raise ValueError('Funtionals can only be set if method == DFT')
        if value not in QM_FUNCTIONALS_ALL:
            raise ValueError('functional {} not recognized'.format(value))
        self._qm_functional = value

    @property
    def qm_basis_set(self):
        if self._qm_basis_sets_extra:
            return self._qm_basis_set, self._qm_basis_sets_extra
        return self._qm_basis_set

    @qm_basis_set.setter
    def qm_basis_set(self, value):
        if value.rstrip('+*') not in QM_BASIS_SETS:
            raise ValueError('Basis set {} not recognized. '
                             'Try with one of: {}'.format(value, ', '.join(QM_BASIS_SETS)))
        self._qm_basis_set = value

    @property
    def qm_basis_set_extra(self):
        return self._qm_basis_sets_extra

    def add_extra_basis_set(self, basis_set, elements, extra_args=None, position=None):
        basis_set = CustomBasisSet(basis_set, elements, extra_args=extra_args,
                                   position=position)
        self._qm_basis_sets_extra.append(basis_set)

    # MM model
    @property
    def mm_forcefield(self):
        return self._mm_forcefield

    @mm_forcefield.setter
    def mm_forcefield(self, value):
        if value not in MM_FORCEFIELDS:
            raise ValueError('Forcefield {} not recognized'.format(value))
        self._mm_forcefield = value


    def add_mm_forcefield(self, value):
        if value.endswith('.frcmod') and os.path.isfile(value):
            self._mm_forcefield_extra = import_from_frcmod(value)
        else:
            raise ValueError('Supply a .frcmod file to load new parameters')

    # System
    @property
    def system(self):
        return '\n'.join(['{} {}'.format(self.charge, self.spin)] + 
                         map(str, self.atoms))

    @property 
    def atoms(self):
        if not self._atoms:
            raise ValueError('System does not contain any atoms! Use .add_atom()!')
        if len(self._atoms) > 250000:
            raise ValueError('Max number of atoms is 250 000.')
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        if value:
            for a in value:
                self.add_atom(atom=a)

    def add_atom(self, atom=None, *args, **kwargs):
        if atom is None:
            atom = GaussianAtom(*args, **kwargs)
        if isinstance(atom, GaussianAtom):
            self._atoms.append(atom)
        else:
            raise TypeError('Provide either `atom` or options to construct '
                            'a GaussianAtom instance.')

    # ModRedundant restraints
    @property
    def restraints(self):
        return self._restraints

    def add_restraint(self, restraint, atoms, min_=None, max_=None, diag_elem=None,
                      nsteps=None, stepsize=None):
        if restraint not in ('A', 'F', 'B', 'K', 'R', 'D', 'H', 'S'):
            raise ValueError('Restraint {} not recognized'.format(restraint))
        if not (1 <= len(atoms) <= 4):
            raise ValueError('Supply between 1 and 4 atoms.')
        if not all(isinstance(i, int) or i == '*' for i in atoms):
            raise ValueError('Atoms must be int or *')

        prefixes = ['X', 'B', 'A', 'D']
        line = [prefixes[len(atoms) - 1]]
        line.extend(atoms)
        line.append(restraint)
        if restraint == 'H':
            if diag_elem is None:
                raise ValueError('diag_elem must be set if restraint == H')
            line.append(diag_elem)
        elif restraint == 'S':
            if None in (nsteps, stepsize):
                raise ValueError('nsteps and stepsize must be set if restraint == S')
            line.extend([nsteps, stepsize])
        if max_ is not None:
            if min_ is not None:
                line.append(min_)
            line.append(max_)
        self._restraints.append(' '.join(map(str, line)))

    @property
    def charge(self):
        if self._charge is None:
            raise ValueError('Please set charge or use .compute_charge()')
        charge = self.compute_charge()
        if self._charge is not None and self._charge != charge:
            print('! Registered charge ({}) does not match '
                  'computed charge from atoms ({})'.format(self._charge, charge))
        return self._charge

    @charge.setter
    def charge(self, value):
        self._charge = value

    def compute_charge(self):
        if self._atoms:
            return sum(a.charge for a in self.atoms if a.charge is not None)
    
    @property
    def multiplicity(self):
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, value):
        self._multiplicity = value


class GaussianAtom(object):

    """
    Object-oriented abstraction of a cartesian atom specification in
    Gaussian input files.

    It offers support for QM and MM atoms. More info is available at
    http://gaussian.com/molspec/.
    """

    def __init__(self, element, coordinates, atom_type=None, charge=None, freeze_code=None,
                 residue_number=None, residue_name=None, pdb_name=None, fragment=None,
                 iso=None, spin=None, zeff=None, qmom=None, nmagm=None, znuc=None,
                 oniom_layer=None, oniom_link=None, oniom_bonded=None, is_link=False):
        self._element = None
        self._coordinates = None
        self._atom_type = None
        self._charge = None
        self._freeze_code = None
        self._residue_number = None
        self._residue_name = None
        self._pdb_name = None
        self._fragment = None
        self._iso = None
        self._spin = None
        self._zeff = None
        self._qmom = None
        self._nmagm = None
        self._znuc = None
        self._oniom_layer = None
        self._oniom_link = None
        self._oniom_bonded = None
        self.is_link = bool(is_link)

        # Set and verify
        self.element = element
        self.coordinates = coordinates
        self.atom_type = atom_type
        self.charge = charge
        self.freeze_code = freeze_code
        self.residue_number = residue_number
        self.residue_name = residue_name
        self.pdb_name = pdb_name
        self.fragment = fragment
        self.iso = iso
        self.spin = spin
        self.zeff = zeff
        self.qmom = qmom
        self.nmagm = nmagm
        self.znuc = znuc
        self.oniom_layer = oniom_layer
        self.oniom_link = oniom_link

    @property
    def element(self):
        return self._element

    @element.setter
    def element(self, value):
        if value is None:
            self._element = None
            return
        if not value or value[0] not in ascii_letters:
            raise ValueError('Element cannot be empty and must start with a letter')
        self._element = value

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value):
        if value is None:
            if self.is_link:
                self._coordinates = None
                return
        try:
            if len(value) == 3:
                self._coordinates = tuple(float(v) for v in value)
                return
        except (ValueError, TypeError):
            pass
        raise ValueError('Coordinates must be 3-tuple of float (x, y, z)')

    @property
    def atom_type(self):
        return self._atom_type

    @atom_type.setter
    def atom_type(self, value):
        if value is None:
            self._atom_type = None
            return
        if not value or value[0] not in ascii_letters:
            raise ValueError('Atom_type cannot be empty and must start with a letter')
        self._atom_type = value

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, value):
        if value is None:
            self._charge = None
            return
        if self.atom_type is None:
            raise ValueError('Setting charges requires setting atom_type beforehand')
        try:
            self._charge = float(value)
        except TypeError:
            raise ValueError('Charge values must be float or float-like')

    @property
    def freeze_code(self):
        return self._freeze_code

    @freeze_code.setter
    def freeze_code(self, value):
        if value is None:
            self._freeze_code = None
            return
        if isinstance(value, int) and value <= 0:
            self._freeze_code = value
            return
        raise ValueError('Freeze_code must be int and <= 0')

    @property
    def residue_number(self):
        return self._residue_number

    @residue_number.setter
    def residue_number(self, value):
        if value is None:
            self._residue_number = None
            return
        if not isinstance(value, int):
            raise ValueError('residue_number must be int')
        self._residue_number = value

    @property
    def residue_name(self):
        return self._residue_name

    @residue_name.setter
    def residue_name(self, value):
        if value is None:
            self._residue_name = None
            return
        if isinstance(value, basestring) and value and value[0] not in ascii_letters:
            self._residue_name = value
            return
        raise ValueError('residue_name cannot be empty and must start with a letter')

    @property
    def pdb_name(self):
        return self._pdb_name

    @pdb_name.setter
    def pdb_name(self, value):
        if value is None:
            self._pdb_name = None
            return
        if isinstance(value, basestring) and value and value[0] not in ascii_letters:
            self._pdb_name = value
            return
        raise ValueError('pdb_name cannot be empty and must start with a letter')

    @property
    def fragment(self):
        return self._fragment

    @fragment.setter
    def fragment(self, value):
        if value is None:
            self._fragment = None
            return
        if not isinstance(value, int):
            raise ValueError('spin values must be int')
        self._fragment = value

    @property
    def iso(self):
        return self._iso

    @iso.setter
    def iso(self, value):
        if value is None:
            self._iso = None
            return
        try:
            self._iso = float(value)
        except ValueError:
            raise ValueError('iso must be float or float-like')

    @property
    def spin(self):
        return self._spin

    @spin.setter
    def spin(self, value):
        if value is None:
            self._spin = None
            return
        try:
            self._spin = float(value)
        except ValueError:
            raise ValueError('spin must be float or float-like')

    @property
    def zeff(self):
        return self._zeff

    @zeff.setter
    def zeff(self, value):
        if value is None:
            self._zeff = None
            return
        try:
            self._zeff = float(value)
        except ValueError:
            raise ValueError('zeff must be float or float-like')

    @property
    def qmom(self):
        return self._qmom

    @qmom.setter
    def qmom(self, value):
        if value is None:
            self._qmom = None
            return
        try:
            self._qmom = float(value)
        except ValueError:
            raise ValueError('qmom must be float or float-like')

    @property
    def nmagm(self):
        return self._nmagm

    @nmagm.setter
    def nmagm(self, value):
        if value is None:
            self._nmagm = None
            return
        try:
            self._nmagm = float(value)
        except ValueError:
            raise ValueError('nmagm must be float or float-like')

    @property
    def znuc(self):
        return self._znuc

    @znuc.setter
    def znuc(self, value):
        if value is None:
            self._znuc = None
            return
        try:
            self._znuc = float(value)
        except ValueError:
            raise ValueError('znuc must be float or float-like')

    @property
    def oniom_layer(self):
        return self._oniom_layer

    @oniom_layer.setter
    def oniom_layer(self, value):
        if value is None:
            self._oniom_layer = None
            return
        if value in ('H', 'h', 'M', 'm', 'L', 'l'):
            self._oniom_layer = value.upper()
            return
        raise ValueError('oniom_layer must be H, M, or L')

    @property
    def oniom_link(self):
        return self._oniom_link

    @oniom_link.setter
    def oniom_link(self, value):
        if value is None:
            self._oniom_link = None
            return
        if not isinstance(value, GaussianAtom):
            raise ValueError('oniom_link must be a GaussianAtom instance')

    @property
    def oniom_bonded(self):
        return self._oniom_bonded

    @oniom_bonded.setter
    def oniom_bonded(self, value):
        if value is None:
            self._oniom_bonded = None
            return
        try:
            self._oniom_bonded = int(value)
        except ValueError:
            raise ValueError('oniom_bonded must be int or int-like')

    @property
    def oniom_scale_factors(self):
        return self._oniom_scale_factors

    @oniom_scale_factors.setter
    def oniom_scale_factors(self, value):
        if value is None:
            self._oniom_scale_factors = value
            return
        try:
            if len(value) <= 3:
                self._oniom_scale_factors = tuple(float(v) for v in value)
                return
        except (ValueError, TypeError):
            pass
        raise ValueError('oniom_scale_factors must be tuple of float, three values max')

    #--------------------------------------
    # Helper methods
    #--------------------------------------
    @property
    def atom_spec(self):
        """
        Summarize atom information in a single string
        """
        line = [self.element]
        if self.atom_type:
            line.append('-{}'.format(self.atom_type))
        if self.charge is not None:
            line.append('-{}'.format(self.charge))

        return ''.join(line)

    @property
    def keywords(self):
        """
        Get a dict with all the keywords
        """
        keywords = ['residue_number', 'residue_name', 'pdb_name', 'fragment',
                    'iso', 'spin', 'zeff', 'qmom', 'nmagm', 'znuc']
        return {k: getattr(self, k) for k in keywords}

    @property
    def keywords_spec(self):
        """
        Summarize all keywords in a single string
        """
        keywords = []
        if self.residue_number is not None:
            keywords.append("{}={}".format("RESNum", self.residue_number))
        if self.residue_name is not None:
            keywords.append("{}={}".format("RESName", self.residue_name))
        if self.pdb_name is not None:
            keywords.append("{}={}".format("PDBName", self.pdb_name))
        if self.fragment is not None:
            keywords.append("{}={}".format("Fragment", self.fragment))
        if self.iso is not None:
            keywords.append("{}={}".format("Iso", self.iso))
        if self.spin is not None:
            keywords.append("{}={}".format("Spin", self.spin))
        if self.zeff is not None:
            keywords.append("{}={}".format("ZEff", self.zeff))
        if self.qmom is not None:
            keywords.append("{}={}".format("QMom", self.qmom))
        if self.nmagm is not None:
            keywords.append("{}={}".format("NMagM", self.nmagm))
        if self.znuc is not None:
            keywords.append("{}={}".format("ZNuc", self.znuc))

        if keywords:
            return '({})'.format(','.join(keywords))

    @property
    def coordinates_spec(self):
        if self.coordinates is not None:
            return ' '.join(map(str, self.coordinates))

    def __str__(self):
        # Atom element, name, charge
        line = [self.atom_spec]

        # Atom keywords
        keywords = self.keywords_spec
        if keywords:
            line.append(keywords)

        # Freeze code (rigidity)
        if self.freeze_code is not None:
            line.append(' {}'.format(self.freeze_code))

        # Coordinates
        coords = self.coordinates_spec
        if coords:
            line.append(' {}'.format(coords))

        # ONIOM config
        if self.oniom_layer:
            line.append(' {}'.format(self.oniom_layer))
        link = self.oniom_link
        if link:
            line.append(' {}'.format(link.atom_spec))
            if link.oniom_bonded:
                line.append(' {}'.format(link.oniom_bonded))
            if link.oniom_scale_factors:
                line.append(' {}'.format(' '.join(link.scale_factors)))

        return ''.join(map(str, line))


class CustomBasisSet(object):

    """
    Object-oriented abstraction of an extra basis set for Gaussian
    input files.

    It includes support for ebsel local databases as well as
    online databases from BSE and Cosmo.
    """

    def __init__(self, basis_set, elements, position=0, extra_args="",):
        if basis_set.rstrip('+*') not in QM_BASIS_SETS:
            raise ValueError('Basis set {} not recognized'.format(basis_set))
        if basis_set in QM_BASIS_SETS_WITH_ARGS and not extra_args:
            raise ValueError('{} requires extra_args'.format(basis_set))
        if position and len(elements) > 1:
            raise ValueError('Position can only be set if one element is suppled')
        self.basis_set = basis_set
        self.position = position
        self.extra_args = extra_args
        self.elements = elements

    def __str__(self):
        lines = [' '.join(self.elements), self.position]
        lines.append(self.basis_set)
        lines.append(self.extra_args)
        lines.append(' ****')
        return '\n'.join(lines)

    @classmethod
    def from_database(cls, basis_set, element, database='cosmologic-services'):
        """
        Import basis_set from online database.

        TODO:
            Support https://bse.pnl.gov/bse/portal/
            Support http://cosmologic-services.de/basis-sets/basissets.php
        """
        try:
            import requests
            from bs4 import BeautifulSoup
        except ImportError:
            raise ImportError('You need to install requests and beautifulsoup4 '
                              'to import basis sets from online databases.')
        api = {
            'cosmologic-services': {
                'header': '! Downloaded from http://cosmologic-services.de/basis-sets/basissets.php',
                'url': 'http://cosmologic-services.de/basis-sets/getbasis.php',
                'bs_parser': 'html.parser',
                'parser': cls._cosmo_parser,
                'data': {'basis': basis_set, element: element, 
                         'kind': 'Basis', 'format': 'Gaussian'}
            },
            'comp.chem.umn.edu': {
                'header': '! Downloaded from http://comp.chem.umn.edu/basissets/basis.cgi',
                'url': 'http://comp.chem.umn.edu/basissets/basis.cgi',
                'bs_parser': 'html.parser',
                'parser': cls._umn_parser,
                'data': {'basis_list': basis_set, 'element_list': element,
                         'format_list': 'Gaussian'}
            }
        }
        try:
            db = api[database]
        except KeyError:
            raise ValueError('Database {} not supported'.format(database))
        try:
            r = requests.post(db['url'], db['data'], timeout=10)
        except requests.exceptions.RequestException:
            raise
        if not r.ok:
            raise ValueError('Could not retrieve data from {}'.format(database))
        html = BeautifulSoup(r.content, db['bs_parser'])
        return '\n'.join([db['header'], db['parser'](html)])

    @staticmethod
    def _cosmo_parser(html):
        return html.find('pre').text

    @staticmethod
    def _umn_parser(html):
        text = str(html.find('body').contents[4])
        return text.replace('<br>', '\n').replace('</br>', '').replace(u'\xa0', u' ')

    @classmethod
    def from_bse(cls, basis_set, *elements):
        try:
            from ebsel.EMSL_local import EMSL_local
        except ImportError:
            raise ImportError('Access to Basis Set Exchange db requires ebsel package')

        db = EMSL_local(fmt="g94")
        try:
            basis = db.get_basis(basis_set, elements=elements)
        except UnboundLocalError:
            if basis_set not in [b for (b, d) in db.get_available_basis_sets()]:
                raise ValueError('Basis set {} not recognized'.format(basis_set))
            supported_elements = db.get_available_elements(basis_set)
            for element in elements:
                if element not in supported_elements:
                    raise ValueError('Basis set {} does not support element {}'.format(basis_set, element))
        else:
            # prepend a '-' to each element to prevent Gaussian errors if atom not present in system
            return '\n'.join([b.replace('****\n', '****\n-') for b in basis])


def import_from_frcmod(path):

    pass 

if __name__ == '__main__':
    atom = GaussianAtom(element='C', coordinates=(10, 10, 10), atom_type='CT', 
                        charge=1.0, residue_number=1, oniom_layer='H')
    print(atom)
