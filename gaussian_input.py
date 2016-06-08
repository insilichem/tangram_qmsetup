#!/usr/bin/env python
# encoding: utf-8

# Get used to importing this in your Py27 projects!
from __future__ import print_function, division
# Python stdlib
from string import ascii_letters
import sys
sys.tracebacklimit = 0
# Chimera stuff
# Additional 3rd parties
# Own

if sys.version_info.major == 3:
    basestring = str

"""
An object-oriented abstraction of a GAUSSIAN input file.
"""


class GaussianInputFile(object):

    def __init__(self, *args, **kwargs):
        pass


class GaussianAtom(object):

    def __init__(self, element, coordinates=None, atom_type=None, charge=None, freeze_code=None,
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

if __name__ == '__main__':
    atom = GaussianAtom(element='C', coordinates=None, atom_type='CT', charge=1.0,
                        residue_number=1, oniom_layer='H')
    print(atom)
