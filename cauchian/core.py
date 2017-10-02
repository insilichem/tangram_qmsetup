#!/usr/bin/env python
# encoding: utf-8


from __future__ import print_function, division
# Python stdlib
from tkFileDialog import askopenfilename, asksaveasfilename
# Chimera stuff
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
        
        # Tie everything up
        self.set_mvc()

    def set_mvc(self):
        """
        Tie GUI to specific actions and events in the controller
        """
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
                         'ui_qm_functionals', 'ui_qm_basis', 'ui_qm_basis_ext',
                         'ui_mm_forcefields', 'ui_mm_water_forcefield', 'ui_flex_policy',
                         'ui_memory_units')
        for name in with_callback:
            item = getattr(self.gui, name)
            command = getattr(self, '_cb' + name[2:], None)
            if command:
                item.configure(command=command)

        variables = ('var_molecule_replicas', 'var_job', 'var_job_options', 'var_frequencies',
                     'var_calculation', 'var_solvent', 'var_qm_method',
                     'var_qm_functional', 'var_qm_functional_type', 'var_qm_basis',
                     'var_qm_basis_ext', 'var_qm_basis_custom', 'var_qm_basis',
                     'var_qm_basis_ext', 'var_mm_forcefield', 'var_mm_water_forcefield',
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
    def _cmd_qm_basis_per_atom(self, *args):
        if self._basis_set_dialog is None:
            from gui import BasisSetDialog
            self._basis_set_dialog = BasisSetDialog(self.gui._qm_basis_extra,
                                                    parent=self.gui)
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
            self.gui.ui_job_options.configure(menubutton_state='normal')
            self.gui.ui_job_options.setitems(options, index=0)
        else:
            self.gui.ui_job_options.setitems([])
            self.gui.ui_job_options.configure(menubutton_state='disabled')

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

    def _trc_qm_basis(self, *args):
        basis = self.gui.var_qm_basis.get()
        ext = self.gui.var_qm_basis_ext.get()
        if basis:
            self.gui.var_qm_basis_custom.set('{}{}'.format(basis, ext if ext else ''))
    _trc_qm_basis_ext = _trc_qm_basis


class Model(object):

    """
    Light interface between UCSF Chimera and pygaussian.GaussianInputFile
    """

    def __init__(self, gui, *args, **kwargs):
        self.gui = gui
        self._atoms_map = {}

    @property
    def route(self):
        state = self.state(with_atoms=False)
        i = GaussianInputFile(**state)
        return i.route

    def build_model_from_current_state(self):
        state = self.state()
        return GaussianInputFile(**state)

    def state(self, with_atoms=True):
        self._atoms_map.clear()
        state = {}
        for attr in dir(self.gui):
            if attr.startswith('var_'):
                state[attr[4:]] = getattr(self.gui, attr).get()

        state['qm_basis_extra'] = self.gui._qm_basis_extra.copy()
        print(state)
        # System
        if with_atoms:
            molecules = self.gui.ui_molecules.getvalue()
            state['atoms'] = [self.chimera_atom_to_gaussian(a)
                            for m in molecules for a in m.atoms]
        return state

    def chimera_atom_to_gaussian(self, atom, oniom=False):
        element = atom.element.name
        coordinates = tuple(*atom.coord())
        charge = getattr(atom, 'charge', None)
        gaussian_atom = GaussianAtom(element, coordinates, charge=charge)
        if oniom:
            pass # Get & set mooooer attrs
        self._atoms_map[atom] = gaussian_atom
        return gaussian_atom