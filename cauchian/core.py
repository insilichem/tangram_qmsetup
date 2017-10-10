#!/usr/bin/env python
# encoding: utf-8


from __future__ import print_function, division
# Python stdlib
from tkFileDialog import askopenfilename, asksaveasfilename
from copy import deepcopy
from traceback import print_exc
# Chimera stuff
import chimera
# Additional 3rd parties
# Own
from pygaussian import GaussianAtom, GaussianInputFile
from pygaussian import JOB_OPTIONS, QM_FUNCTIONALS


class Controller(object):

    def __init__(self, gui=None, model=None, *args, **kwargs):
        self.gui = gui
        self.model = model

        # Flags
        self._basis_set_dialog = None
        self._layers_dialog = None
        
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
        buttons = ('ui_molecules', 'ui_molecules_replicas', 'ui_layers',
                   'ui_solvent_cfg', 'ui_qm_basis_per_atom', 'ui_flex_btn',
                   'ui_redundant_btn', 'ui_charges_auto', 'ui_charges_manual',
                   'ui_mm_frcmod_btn', 'ui_checkpoint_btn')
        for btn in buttons:
            button = getattr(self.gui, btn)
            command = getattr(self, '_cmd' + btn[2:], None)
            if command:
                button.configure(command=command)

        # Event callbacks & variable tracing
        with_callback = ('ui_job', 'ui_job_options', 'ui_calculation',
                         'ui_solvent', 'ui_qm_methods', 'ui_qm_functional_type',
                         'ui_qm_functionals', 'ui_qm_basis_kind', 'ui_qm_basis_ext',
                         'ui_mm_forcefields', 'ui_mm_water_forcefield', 'ui_flex_policy',
                         'ui_memory_units')
        for name in with_callback:
            item = getattr(self.gui, name)
            command = getattr(self, '_cb' + name[2:], None)
            if command:
                item.configure(command=command)

        variables = ('var_molecule_replicas', 'var_job', 'var_job_options', 'var_frequencies',
                     'var_calculation', 'var_solvent', 'var_qm_method',
                     'var_qm_functional', 'var_qm_functional_type', 'var_qm_basis_set',
                     'var_qm_basis_kind', 'var_qm_basis_ext',
                     'var_mm_forcefield', 'var_mm_water_forcefield',
                     'var_mm_frcmod', 'var_charge_qm', 'var_charge_mm',
                     'var_multiplicity_qm', 'var_multiplicity_mm', 'var_flex_policy',
                     'var_flex_lbl', 'var_title', 'var_checkpoint', 'var_checkpoint_path',
                     'var_nproc', 'var_memory', 'var_memory_units')
        for name in variables:
            var = getattr(self.gui, name)
            command = getattr(self, '_trc' + name[3:], None)
            if command:
                var.trace('w', command)
                command()

    # Button actions start with _cmd
    def _cmd_Preview(self, *args):
        contents = ''
        self.gui.ui_preview.configure(text_state='normal')
        self.gui.ui_preview.clear()
        try:
            current_file = self.model.build_model_from_current_state()
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
        contents = self._cmd_Preview()
        if contents:
            path = asksaveasfilename(title='Choose destination (.com)',
                                     filetypes=[('Gaussian input', '*.com'), ('All', '*')],
                                     defaultextension='.com')
            if path:
                with open(path, 'w') as f:
                    f.write(contents)
                self.gui.status('Saved to {}!'.format(path), blankAfter=5)
    
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

    def _cmd_mm_frcmod_btn(self, *args):
        path = askopenfilename(filetypes=[('Frcmod', '*.frcmod'), ('All files', '*')])
        if path:
            self.gui.var_mm_frcmod.set(path)

    def _cmd_checkpoint_btn(self, *args):
        path = asksaveasfilename()
        if path:
            self.gui.var_mm_frcmod.set(path)

    # Event callbacks start with _cb

    # Variables are traced with _trc methods
    def _trc_molecule_replicas(self, *args):
        if self.gui.var_molecule_replicas.get():
            self.gui.ui_molecules_replicas.configure(listbox_state='normal')
        else:
            self.gui.ui_molecules_replicas.configure(listbox_state='disabled')

    def _trc_calculation(self, *args):
        value = self.gui.var_calculation.get()
        if value == 'ONIOM':
            self.gui.ui_layers.configure(state='normal')
            self.gui.ui_mm_forcefields.configure(menubutton_state='normal')
            self.gui.ui_mm_water_forcefield.configure(menubutton_state='normal')
            self.gui.ui_mm_frcmod.configure(state='normal')
            self.gui.ui_mm_frcmod_btn.configure(state='normal')
            self.gui.ui_charges_mm.configure(state='normal')
            self.gui.ui_multiplicity_mm.configure(state='normal')
        else:  # == QM
            self.gui.ui_layers.configure(state='disabled')
            self.gui.ui_mm_forcefields.configure(menubutton_state='disabled')
            self.gui.ui_mm_water_forcefield.configure(menubutton_state='disabled')
            self.gui.ui_mm_frcmod.configure(state='disabled')
            self.gui.ui_mm_frcmod_btn.configure(state='disabled')
            self.gui.ui_charges_mm.configure(state='disabled')
            self.gui.ui_multiplicity_mm.configure(state='disabled')

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
        if options:
            size = len(options) * 20
            size = size if size <= 200 else 200
            self.gui.ui_job_options.component('scrolledlist').setlist(options)
            self.gui.ui_job_options.configure(scrolledlist_hull_height=len(options)*20)
        else:
            self.gui.ui_job_options.clear()
            self.gui.ui_job_options.configure(scrolledlist_hull_height=0)

        if value in ('SP', 'Opt'):
            self.gui.ui_frequencies.configure(state='normal')
        else:
            self.gui.var_frequencies.set(0)
            self.gui.ui_frequencies.configure(state='disabled')

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


class Model(object):

    """
    Light interface between UCSF Chimera and pygaussian.GaussianInputFile
    """

    def __init__(self, gui, *args, **kwargs):
        self.gui = gui
        self._atoms_map = {}

    def build_model_from_current_state(self, with_atoms=True, with_replicas=False):
        state = self.state
        self._atoms_map.clear()
        kwargs = dict(
            title=state['title'],
            processors=state['nproc'] or None,
            memory=(state['memory'] or None, state['memory_units']),
            checkpoint=state['checkpoint_path'] if state['checkpoint'] else None,
            connectivity=state['connectivity'],
        )
        infile = GaussianInputFile(**kwargs)
        infile.job = state['job']
        infile.job_options = state['job_options']
        infile.frequencies = bool(state['frequencies'])
        infile.qm_method = state['qm_method']
        infile.qm_basis_set = state['qm_basis_set']
        if state['qm_method'] == 'DFT':
            infile.qm_functional = state['qm_functional']
        if state['qm_basis_set_extra']:
            for element, bs in state['qm_basis_set_extra'].items():
                infile.add_extra_basis_set(bs, element)
        if state['calculation'] == 'ONIOM':  # enter MM details
            infile.charge = state['charge_qm'] + state['charge_mm']
            infile.multiplicity = state['multiplicity_qm'] + state['multiplicity_mm']
            infile.mm_forcefield = state['mm_forcefield']
            infile.mm_water_forcefield = state['mm_water_forcefield']
            if state['mm_frcmod']:
                infile.add_mm_forcefield(state['mm_frcmod'])
        else:  # QM only details
            infile.charge = state['charge_qm']
            infile.multiplicity = state['multiplicity_qm']
        if state['restraints']:
            for restraint in state['restraints']:
                infile.add_restraint(restraint)
        if state['qm_keywords']:
            tokens = state['qm_keywords'].split()
            for token in tokens:
                key, value = token.split('=')
                infile.add_route_option(token)
        
        if with_atoms:
            infile.atoms = self.process_atoms(state, connectivity=state['connectivity'])
        if with_replicas:
            infile_replicas = []
            for replica in state['replicas']:
                replica_atoms = [self.gaussian_atom(a, n=i+1) 
                                 for i, a in enumerate(replica.atoms)]
                if len(replica_atoms) != len(infile.atoms):
                    raise ValueError('Replica {} has different number of atoms.'.format(replica))
                infile_replica = deepcopy(infile)
                infile_replica.atoms = replica_atoms
                infile_replicas.append(infile_replica)
            return infile_replicas
        
        return infile

    @property
    def state(self):
        state = {}
        for attr in dir(self.gui):
            if attr.startswith('var_'):
                state[attr[4:]] = getattr(self.gui, attr).get()

        state['qm_basis_set_extra'] = self.gui._qm_basis_extra.copy()
        state['restraints'] = self.gui._restraints.copy()
        state['molecule'] = self.gui.ui_molecules.getvalue()
        state['layers'] = self.gui._layers.copy()
        if self.gui.var_molecule_replicas.get():
            state['replicas'] = self.gui.ui_molecules_replicas.getvalue()
                
        return state

    def gaussian_atom(self, atom, n, oniom=True, layer=None, link=None):
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
            gatom.atom_type = atom.idatmType
            gatom.geometry = chimera.idatm[gatom.idatmType].geometry
            gatom.residue_name = atom.residue.type
            gatom.residue_number = atom.residue.id.position
            gatom.oniom_layer = layer
            gatom.oniom_link = link
            gatom.charge = getattr(atom, 'charge', None)
        self._atoms_map[atom] = gatom
        return gatom
    
    def process_atoms(self, state=None, connectivity=True):
        if state is None:
            state = self.state
        gaussian_atoms = []
        chimera_atoms = state['molecule'].atoms
        mapping = {}
        oniom, kw = False, {}

        if state['calculation'] == 'ONIOM':  # we have layers to deal with!
            oniom = True
            layers = state['layers']
            if not layers:
                raise chimera.UserError('ONIOM layers have not been defined!') 
        for n, catom in enumerate(chimera_atoms):
            if oniom:
                kw = dict(layer=layers[catom], link=None)
            gatom = self.gaussian_atom(catom, n=n+1, oniom=oniom, **kw)
            gaussian_atoms.append(gatom)
            mapping[catom] = gatom
        
        if connectivity:
            for catom, gatom in zip(chimera_atoms, gaussian_atoms):
                for neighbor in catom.neighbors:
                    gatom.neighbors.append(mapping[neighbor])
        
        return gaussian_atoms