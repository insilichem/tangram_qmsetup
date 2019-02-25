#!/usr/bin/env python
# encoding: utf-8


from __future__ import print_function, division
# Python stdlib
from tkFileDialog import askopenfilename, asksaveasfilename
from copy import copy
from traceback import print_exc
import os
import json
# Chimera stuff
import chimera
from chimera.baseDialog import NotifyDialog
from WriteMol2 import chimera2sybyl
# Additional 3rd parties
# Own
from pygaussian import GaussianAtom, GaussianInputFile
from pygaussian import JOB_OPTIONS, QM_FUNCTIONALS
try:
    from bondorder.core import assign_bond_orders
except ImportError:
    print('! Tangram BondOrder not installed. Bond order perception disabled!')
    def assign_bond_orders(*args):
        pass


class Controller(object):

    def __init__(self, gui=None, model=None, *args, **kwargs):
        self.gui = gui
        self.gui.controller = self
        self.model = model

        # Flags
        self._basis_set_dialog = None
        self._layers_dialog = None
        self._mmtypes_dialog = None
        self._modredundant_dialog = None

        # Tie everything up
        self.set_mvc()

    def set_mvc(self):
        """
        Tie GUI to specific actions and events in the controller
        """
        # Main buttons
        for name, button in self.gui.buttonWidgets.items():
            command = getattr(self, '_cmd_' + name, None)
            if command:
                button.configure(command=command)

        # Button actions
        buttons = ('ui_layers', 'ui_solvent_btn', 'ui_qm_basis_per_atom',
                   'ui_redundant_btn', 'ui_checkpoint_btn','ui_mm_set_types_btn')
        for btn in buttons:
            button = getattr(self.gui, btn)
            command = getattr(self, '_cmd' + btn[2:], None)
            if command:
                button.configure(command=command)

        # Event callbacks & variable tracing
        with_callback = ('ui_job', 'ui_job_options', 'ui_calculation',
                         'ui_qm_methods', 'ui_qm_functional_type',
                         'ui_qm_functionals', 'ui_qm_basis_kind', 'ui_qm_basis_ext',
                         'ui_mm_water_forcefield', 'ui_memory_units')
        for name in with_callback:
            item = getattr(self.gui, name)
            command = getattr(self, '_cb' + name[2:], None)
            if command:
                item.configure(command=command)

        variables = ('var_molecule_replicas', 'var_job', 'var_job_options',
                     'var_calculation', 'var_qm_method',
                     'var_qm_functional', 'var_qm_functional_type', 'var_qm_basis_set',
                     'var_qm_basis_kind', 'var_qm_basis_ext',
                     'var_mm_forcefield', 'var_mm_water_forcefield',
                     'var_mm_residues', 'var_mm_external',
                     'var_charge_qm', 'var_charge_mm','var_multiplicity_qm', 
                     'var_multiplicity_mm', 'var_title', 'var_checkpoint', 
                     'var_checkpoint_path', 'var_nproc', 'var_memory', 
                     'var_memory_units', 'var_mm_from_mol2')
        for name in variables:
            var = getattr(self.gui, name)
            command = getattr(self, '_trc' + name[3:], None)
            if command:
                var.trace('w', command)
                command()

        # Custom widgets with other callback keywords
        self.gui.ui_molecules['selectioncommand'] = self._cmd_molecules

        # Chimera triggers
        chimera.triggers.addHandler('Molecule', self._trg_molecule_changed, None)
        self._cmd_molecules()

    # Button actions start with _cmd
    def _cmd_Preview(self, *args):
        contents = ''
        self.gui.ui_preview.configure(text_state='normal')
        self.gui.ui_preview.clear()
        try:
            current_file = self.model.build_model_from_current_state(with_replicas=False)[0]
            contents = current_file.build(timestamp=True)
        except Exception as e:
            self.gui.status('Could not preview file due to {}: {}'.format(type(e).__name__,
                            str(e)[:50]), color='red', blankAfter=5)
            print_exc()
        else:
            self.gui.ui_preview.setvalue('\n'.join(contents.splitlines()[1:]))
            return contents
        finally:
            self.gui.ui_preview.configure(text_state='disabled')

    def _cmd_Copy(self, *args):
        contents = self._cmd_Preview()
        if contents:
            self.gui.uiMaster().clipboard_clear()
            self.gui.uiMaster().clipboard_append(contents)
            self.gui.status('Copied to clipboard!', blankAfter=5)

    def _cmd_Export(self, *args):
        state = self.model.state
        gfiles = self.model.build_model_from_current_state()
        try:
            contents = gfiles[0].build(timestamp=True)
        except:
            contents = ''
        if not contents:
            self.gui.status('Export failed! Check current config.', color='red', blankAfter=5)
            return
        path = asksaveasfilename(title='Choose destination (.com)',
                                 filetypes=[('Gaussian input', '*.com'), ('All', '*')],
                                 defaultextension='.com')
        if not path:
            return
        self.gui.status('Generating input files...')
        fn, ext = os.path.splitext(path)
        # First, dump GUI state
        with open('{}_state.json'.format(fn), 'w') as f:
            print(state)
            json.dump(state, f, default=lambda a: None) #, skipkeys=True
        # Second, export files
        for i, gfile in enumerate(gfiles):
            outpath = '{fn}{i}{ext}'.format(fn=fn, i=i+1 if i else '', ext=ext)
            with open(outpath, 'w') as f:
                contents = gfile.build(timestamp=True)
                f.write(contents)
        self.gui.status('Saved {} file(s) to {}!'.format(len(gfiles), path), blankAfter=5)

    def _cmd_Import(self, *args):
        path = askopenfilename(title='Choose .json state file',
                               filetypes=[('JSON files', '*.json'), ('All', '*')])
        if not path:
            return
        with open(path) as f:
            state = json.load(f)
        self.gui.load_state(state)

    def _cmd_Close(self, *args):
        del self.model
        self.gui.Close()
        del self

    def _cmd_molecules(self, *args):
        m = self.gui.ui_molecules.getvalue()
        if m is None:
            self.gui.buttonWidgets['Export']['state'] = 'disabled'
            self.gui.buttonWidgets['Preview']['state'] = 'disabled'
            self.gui.buttonWidgets['Copy']['state'] = 'disabled'
            return

        self.gui.buttonWidgets['Export']['state'] = 'normal'
        self.gui.buttonWidgets['Preview']['state'] = 'normal'
        self.gui.buttonWidgets['Copy']['state'] = 'normal'

        n_coordsets = len(m.coordSets)
        if n_coordsets > 1:
            self.gui.ui_replicas_chk['state'] = 'normal'
            self.gui.ui_replicas_chk['text'] = '{} frames'.format(n_coordsets)
        else:
            self.gui.ui_replicas_chk['text'] = 'No frames'
            self.gui.ui_replicas_chk['state'] = 'disabled'
            self.gui.var_molecule_replicas.set(0)

    def _cmd_layers(self, *args):
        if self._layers_dialog is None:
            from gui import ONIOMLayersDialog
            self._layers_dialog = ONIOMLayersDialog(self.gui._layers,
                                                    master=self.gui.uiMaster())
        self._layers_dialog.enter()

    def _cmd_qm_basis_per_atom(self, *args):
        if self._basis_set_dialog is None:
            from gui import BasisSetDialog
            self._basis_set_dialog = BasisSetDialog(self.gui._qm_basis_extra,
                                                    master=self.gui.uiMaster())
        self._basis_set_dialog.enter()

    def _cmd_mm_set_types_btn(self, *args):
        if self._mmtypes_dialog is None:
            from gui import MMTypesDialog
            self._mmtypes_dialog = MMTypesDialog(self.gui._mmtypes, self.gui.var_mm_forcefield,
                                                self.gui._mm_frcmod, master=self.gui.uiMaster())
        self._mmtypes_dialog.enter()

    def _cmd_checkpoint_btn(self, *args):
        path = asksaveasfilename()
        if path:
            self.gui.var_checkpoint_path.set(path)

    def _cmd_redundant_btn(self, *args):
        if self._modredundant_dialog is None:
            from gui import ModRedundantDialog
            molecule = self.gui.ui_molecules.getvalue()
            if molecule is None:
                raise chimera.UserError('No molecule selected!')
            self._modredundant_dialog = ModRedundantDialog(self.gui._restraints,
                                                           molecule.atoms,
                                                           master=self.gui.uiMaster(),
                                                           callback=self._cb_after_modredundant)
        self._modredundant_dialog.enter()

    # Event callbacks start with _cb
    def _cb_after_modredundant(self):
        job_options = self.gui.var_job_options.get()
        if self.gui._restraints and 'modredundant' not in job_options:
            if job_options:
                job_options += ',modredundant'
            else:
                job_options = 'modredundant'
            self.gui.var_job_options.set(job_options)
        self._modredundant_dialog = None

    # Variables are traced with _trc methods
    def _trc_calculation(self, *args):
        value = self.gui.var_calculation.get()
        if value == 'ONIOM':
            self.gui.ui_layers['state'] = 'normal'
            self.gui.ui_mm_set_types_btn['state'] = 'normal'
            self.gui.ui_mm_water_forcefield['menubutton_state'] = 'normal'
            self.gui.ui_charges_mm['state'] = 'normal'
            self.gui.ui_multiplicity_mm['state'] = 'normal'
        else:  # == QM
            self.gui.ui_layers['state'] = 'disabled'
            self.gui.ui_mm_set_types_btn['state'] = 'disabled'
            self.gui.ui_mm_water_forcefield['menubutton_state'] = 'disabled'
            self.gui.ui_charges_mm['state'] = 'disabled'
            self.gui.ui_multiplicity_mm['state'] = 'disabled'

    def _trc_checkpoint(self, *args):
        if self.gui.var_checkpoint.get():
            self.gui.ui_checkpoint_fld.configure(state='normal')
            self.gui.ui_checkpoint_btn.configure(state='normal')
        else:
            self.gui.ui_checkpoint_fld.configure(state='disabled')
            self.gui.ui_checkpoint_btn.configure(state='disabled')

    def _trc_job(self, *args):
        value = self.gui.var_job.get()
        options = JOB_OPTIONS.get(value)
        self.gui.ui_job_options.clear()
        if options:
            size = len(options) * 20
            size = size if size <= 200 else 200
            self.gui.ui_job_options.component('scrolledlist').setlist(options)
            self.gui.ui_job_options.configure(scrolledlist_hull_height=len(options)*20)
        else:
            self.gui.ui_job_options.configure(scrolledlist_hull_height=0)

        if value == 'Opt':
            self.gui.ui_redundant_btn['state'] = 'normal'
        else:
            self.gui.ui_redundant_btn['state'] = 'disabled'

    def _trc_qm_method(self, *args):
        value = self.gui.var_qm_method.get()
        if value == 'DFT':
            self.gui.ui_qm_functional_type.configure(menubutton_state='normal')
            self.gui.ui_qm_functionals.configure(menubutton_state='normal')
        else:
            self.gui.ui_qm_functional_type.configure(menubutton_state='disabled')
            self.gui.ui_qm_functionals.configure(menubutton_state='disabled')

    def _trc_qm_functional_type(self, *args):
        value = self.gui.var_qm_functional_type.get()
        self.gui.ui_qm_functionals.setitems(QM_FUNCTIONALS[value], index=0)

    def _trc_qm_basis_kind(self, *args):
        basis = self.gui.var_qm_basis_kind.get()
        ext = self.gui.var_qm_basis_ext.get()
        if basis:
            self.gui.var_qm_basis_set.set('{}{}'.format(basis, ext if ext else ''))
    _trc_qm_basis_ext = _trc_qm_basis_kind
                        
    def _trg_molecule_changed(self, *args, **kwargs):
        self.model._bondorder_cache.clear()
        self.model._atoms_map.clear()
        self.gui._layers.clear()

class Model(object):

    """
    Light interface between UCSF Chimera and pygaussian.GaussianInputFile
    """

    def __init__(self, gui, *args, **kwargs):
        self.gui = gui
        self._atoms_map = {}
        self._bondorder_cache = {}

    def build_model_from_current_state(self, with_atoms=True, with_replicas=True):
        state = self.state
        self._atoms_map.clear()
        kwargs = dict(
            title=state['title'],
            processors=state['nproc'] or None,
            memory=(state['memory'] or None, state['memory_units']),
            checkpoint=state['checkpoint_path'] if state['checkpoint'] else None,
            connectivity=state['connectivity'],
            mm_external=state['mm_external'],
        )
        infile = GaussianInputFile(**kwargs)
        infile.job = state['job']
        infile.job_options = state['job_options']
        infile.qm_method = state['qm_method']
        infile.qm_basis_set = state['qm_basis_set']
        infile.charge = state['charge_qm']
        infile.multiplicity = state['multiplicity_qm']
        if state['qm_method'] == 'DFT':
            infile.qm_functional = state['qm_functional']
        if state['qm_basis_set_extra']:
            for element, bs in state['qm_basis_set_extra'].items():
                infile.add_extra_basis_set(bs, element)
        if state['calculation'] == 'ONIOM':  # enter MM details
            infile.mm_forcefield = state['mm_forcefield']
            infile.mm_water_forcefield = state['mm_water_forcefield']
            if state['mm_frcmod']:
                infile.add_mm_forcefield(state['mm_frcmod'], state['mm_types'])
            if state['charge_mm']:
                infile.mm_charge = int(state['charge_mm'])
            if state['multiplicity_mm']:
                infile.mm_multiplicity = int(state['multiplicity_mm'])
        if state['restraints']:
            for restraint in state['restraints']:
                infile.add_restraint(restraint)
        if state['qm_keywords']:
            tokens = state['qm_keywords'].split()
            for token in tokens:
                key, value = token.split('=')
                infile.add_route_option(token)

        replicas = [infile]
        
        if with_atoms:
            infile.atoms = self.process_atoms(state)
            if with_replicas and state['replicas']:
                for index, coordset in state['molecule'].coordSets.items():
                    if coordset is state['molecule'].activeCoordSet:
                        continue
                    replica = copy(infile)
                    replica.title += ' - Replica {}'.format(index)
                    for atom, coord in zip(replica.atoms, coordset.xyzArray()):
                        atom.coordinates = tuple(coord)
                    replicas.append(replica)
        
        return replicas

    @property
    def state(self):
        state = {}
        for attr in dir(self.gui):
            if attr.startswith('var_'):
                state[attr[4:]] = getattr(self.gui, attr).get()

        state['qm_basis_set_extra'] = self.gui._qm_basis_extra.copy()
        state['restraints'] = self.gui._restraints[:]
        state['molecule'] = self.gui.ui_molecules.getvalue()
        state['layers_flex'] = self.gui._layers.copy()
        state['replicas'] = self.gui.var_molecule_replicas.get()
        state['mm_types'] = self.gui._mmtypes.copy()
        state['mm_frcmod'] = self.gui._mm_frcmod[:]

        return state

    def patch_residue_names(self, state=None):
        if state is None:
            state = self.state
        for residue in state['molecule'].residues:
            #Waters are 'WAT'
            if residue.type  == 'HOH':
                residue.type = 'WAT'
            #Histidines cannot be 'HIS', have to be 'HIP'/'HIN'/'HID'/'HIE'
            elif residue.type == 'HIS':
                hd1, he2 = False, False
                for atom in residue.atoms:
                    if atom.name.upper() == "HD1":
                        hd1 = True
                    elif atom.name.upper() == "HE2":
                        he2 = True
                if hd1 and he2:
                    residue.type = 'HIP'
                elif hd1:
                    residue.type = 'HID'
                elif he2:
                    residue.type = 'HIE'
                else:
                    residue.type = 'HIN'
            #GLU/ASP with COOH have to be GLH/ASH
            elif residue.type == 'GLU' or residue.type == 'ASP':
                for atom in residue.atoms:
                    if atom.name.upper() == "HD1" or atom.name.upper() == 'HD2':
                        residue.type = residue.type[:2] + 'H'
                        break
            #Deprotonated TYR has to be TYD
            elif residue.type == 'TYR':
                tyd = True
                for atom in residue.atoms:
                    if atom.name.upper() == "HH" or atom.name.upper() == 'HO':
                        tyd = False
                        break
                if tyd:
                    residue.type = 'TYD'

            #Valorate if use len(atom.primaryNeighbors()) to discern protonated/deprotonated
            #Consult http://archive.ambermd.org/201406/0113.html

            #Terminal residues
            if len(residue.bondedResidues()) == 1:
                pass

    def gaussian_atom(self, atom, n, oniom=True, layer=None, link=None, frozen=False):
        """
        Creates a GaussianAtom instance from a chimera.Atom object.

        Parameters
        ----------
        atom : chimera.atom
        oniom : bool, optional, default=False
            Parse additional fields to make it ONIOM-compatible.

        Returns
        -------
        GaussianAtom
        """
        element = atom.element.name
        coordinates = atom.coord().data()
        gatom = GaussianAtom(element, coordinates, n)
        if oniom:
            if layer is None:
                raise ValueError('layer must be set if oniom is True')
            gatom.pdb_name = atom.name
            gatom.atom_type = getattr(atom, 'mmType',
                                          chimera2sybyl.get(atom.idatmType, atom.idatmType))
            geometry = chimera.idatm.typeInfo.get(atom.idatmType)
            if geometry is not None:
                gatom.geometry = geometry.geometry
            gatom.residue_name = atom.residue.type
            gatom.residue_number = atom.residue.id.position
            gatom.oniom_layer = layer
            gatom.oniom_link = link
            gatom.freeze_code = frozen
            gatom.charge = getattr(atom, 'charge', None)
        self._atoms_map[atom] = gatom
        return gatom

    def process_atoms(self, state=None):
        if state is None:
            state = self.state
        gaussian_atoms = []
        chimera_atoms = state['molecule'].atoms
        mapping = {}
        oniom, kw, layers_flex, frozen = False, {}, {}, {}

        if state['calculation'] == 'ONIOM':  # we have layers to deal with!
            oniom = True
            layers_flex = state['layers_flex']
            mm_types = state['mm_types']
            if not layers_flex:
                raise chimera.UserError('ONIOM layers have not been defined!')
            if not mm_types:
            	raise chimera.UserError('MM types have not been defined')

        if state['mm_residues']:
            self.patch_residue_names(state)

        for n, catom in enumerate(chimera_atoms):
            layer, frozen = layers_flex.get(catom, (None, 0))

            kw = dict(oniom=oniom, layer=layer, frozen=int(frozen))
            gatom = self.gaussian_atom(catom, n=n+1, **kw)
            gaussian_atoms.append(gatom)
            mapping[catom] = gatom

        show_warning = False
        """
        try:
            assign_bond_orders(state['molecule'], engine='openbabel')
        except:
            try:
                assign_bond_orders(state['molecule'], engine='rdkit')
            except:
                pass
        """
        for catom, gatom in zip(chimera_atoms, gaussian_atoms):
            for cneighbor, bond in catom.bondsMap.items():
                order = getattr(bond, 'order', None)
                if not order:
                    order = 1.0
                gatom.add_neighbor(mapping[cneighbor], order)
        if show_warning:
            errormsg = ('Some bonds did not specify bond order, so a default of 1.0 '
                        'was used. If you want compute them or edit them manually, '
                        'please use Tangram BondOrder extension.')
            d = NotifyDialog(errormsg, icon='warning')
            d.OK = d.Close
            d.enter()

        return gaussian_atoms
