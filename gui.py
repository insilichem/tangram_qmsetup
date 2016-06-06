#!/usr/bin/env python
# encoding: utf-8

# Get used to importing this in your Py27 projects!
from __future__ import print_function, division
# Python stdlib
import Tkinter as tk
import Pmw
import tkFont
# Chimera stuff
import chimera
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MoleculeScrolledListBox, SortableTable
# Additional 3rd parties

# Own
from core import Controller, Model

"""
The gui.py module contains the interface code, and only that. 
It should only 'draw' the window, and should NOT contain any
business logic like parsing files or applying modifications
to the opened molecules. That belongs to core.py.
"""

# This is a Chimera thing. Do it, and deal with it.
ui = None
def showUI(callback=None, *args, **kwargs):
    """
    Requested by Chimera way-of-doing-things
    """
    if chimera.nogui:
        tk.Tk().withdraw()
    global ui
    if not ui:  # Edit this to reflect the name of the class!
        ui = CauchianDialog(*args, **kwargs)
    model = Model()
    controller = Controller(gui=ui, model=model)
    ui.enter()
    if callback:
        ui.addCallback(callback)


STYLES = {
    tk.Entry: {
        'background': 'white',
        'borderwidth': 1,
        'highlightthickness': 0,
        'insertwidth': 1,
    },
    tk.Button: {
        'borderwidth': 1,
        'highlightthickness': 0,
    },
    Pmw.OptionMenu: {
        'menubutton_borderwidth': 1
    }
}


class CauchianDialog(ModelessDialog):

    """
    To display a new dialog on the interface, you will normally inherit from
    ModelessDialog class of chimera.baseDialog module. Being modeless means
    you can have this dialog open while using other parts of the interface.
    If you don't want this behaviour and instead you want your extension to 
    claim exclusive usage, use ModalDialog.
    """

    buttons = ('Save', 'Close')
    default = None
    help = 'https://www.insilichem.com'

    QM_METHODS = ['DFT', 'HF', 'MP2']
    QM_FUNCTIONALS = {
        'Pure': ['VSXC', 'HCTH', 'HCTH93', 'HCTH147', 'HCTH407', 'tHCTH',
                 'M06L', 'B97D', 'B97D3', 'SOGGA11', 'M11L', 'N12', 'MN12L'],
        'Hybrid': ['B3LYP', 'B3P86', 'B3PW91', 'B1B95', 'mPW1PW91', 'mPW1LYP',
                   'mPW1PBE', 'mPW3PBE', 'B98', 'B971', 'B972', 'PBE1PBE', 'B1LYP',
                   'O3LYP', 'BHandH', 'BHandHLYP', 'BMK', 'M06', 'M06HF', 'M062X',
                   'tHCTHhyb', 'APFD', 'APF', 'SOGGA11X', 'PBEh1PBE', 'TPSSh', 'X3LYP'],
        'RS hybrid': ['HSEH1PBE', 'OHSE2PBE', 'OHSE1PBE', 'wB97XD', 'wB97', 'wB97X',
                      'LC-wPBE', 'CAM-B3LYP', 'HISSbPBE', 'M11', 'N12SX', 'MN12SX']
    }
    QM_BASIS_SETS = ['STO-3G', '3-21G', '6-21G', '4-31G', '6-31G', "6-31G(d')",
                     "6-31G(d',p')", '6-311G', 'D95V', 'D95', 'SHC', 'CEP-4G',
                     'CEP-31G', 'CEP-121G', 'LanL2MB', 'LanL2DZ', 'SDD', 'SDDAll',
                     'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z', 'cc-pV6Z']
    QM_BASIS_SETS_EXT = ['', '+', '++', '*', '**']
    # QM_BASIS_SETS_EXT = ['d,f-Diffuse (+)', 'p,d,f-Diffuse (++)', 'd,f-Polarization (*)', 'p,d,f-Polarization (**)']
    MM_FORCEFIELDS = {
        'General': ['Amber', 'Dreiding', 'UFF'],
        'Water': ['TIP3P']
    }
    MEM_UNITS = ['TB', 'GB', 'MB']

    def __init__(self, *args, **kwarg):
        # GUI init
        self.title = 'Plume Cauchian'

        # Job variables
        self.var_optimization = tk.StringVar()
        self.var_frequencies = tk.IntVar()
        self.var_calculation = tk.StringVar()
        self.var_solvent = tk.StringVar()

        # QM variables
        self.var_qm_method = tk.StringVar()
        self.var_qm_functional = tk.StringVar()
        self.var_qm_functional_type = tk.StringVar()
        self.var_qm_basis = tk.StringVar()
        self.var_qm_basis_ext = tk.StringVar()
        self.var_qm_basis_custom = tk.StringVar()
        self.var_qm_basis.trace('w', self._basis_sets_custom_build)
        self.var_qm_basis_ext.trace('w', self._basis_sets_custom_build)
        self.var_qm_keywords = tk.StringVar()

        # MM variables
        self.var_mm_forcefield = tk.StringVar()
        self.var_mm_water_forcefield = tk.StringVar()
        self.var_mm_frcmod = tk.StringVar()

        # Charges & Multiplicity
        self.var_charge_qm = tk.DoubleVar()
        self.var_charge_mm = tk.DoubleVar()
        self.var_multiplicity_qm = tk.IntVar()
        self.var_multiplicity_mm = tk.IntVar()

        # Flexibility & restraints
        self.var_flex_policy = tk.StringVar()
        self.var_redundant = tk.IntVar()

        # Hardware & Output variables
        self.var_title = tk.StringVar()
        self.var_checkpoint = tk.StringVar()
        self.var_nproc = tk.IntVar()
        self.var_memory = tk.IntVar()
        self.var_memory_units = tk.StringVar()

        # Fire up
        ModelessDialog.__init__(self)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

        # Fix styles
        self._fix_styles(*self.buttonWidgets.values())

    def _basis_sets_custom_build(self, *args):
        basis = self.var_qm_basis.get()
        ext = self.var_qm_basis_ext.get()
        self.var_qm_basis_custom.set('{}{}'.format(basis, ext if ext else ''))

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def _fix_styles(self, *widgets):
        for widget in widgets:
            try:
                widget.configure(**STYLES[widget.__class__])
            except:
                pass

    def fillInUI(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

        # Select molecules
        self.ui_molecule_frame = tk.LabelFrame(self.canvas, text='Select molecules')
        self.ui_molecules = MoleculeScrolledListBox(self.ui_molecule_frame)
        self.ui_molecules.pack(expand=True, fill='both', padx=5, pady=5)

        # Modelization
        self.ui_model_frame = tk.LabelFrame(self.canvas, text='Modelization')
        self.ui_optimization = Pmw.OptionMenu(self.canvas,
                                              menubutton_textvariable=self.var_optimization,
                                              items=['No', 'Min', 'TS'])
        self.ui_frequencies = tk.Checkbutton(self.canvas, variable=self.var_frequencies,
                                             text='Get frequencies')
        self.ui_calculation = Pmw.OptionMenu(self.canvas,
                                             menubutton_textvariable=self.var_calculation,
                                             items=['QM', 'ONIOM'])
        self.ui_layers = tk.Button(self.canvas, text='Define layers')
        self.ui_solvent = Pmw.OptionMenu(self.canvas,
                                         menubutton_textvariable=self.var_solvent,
                                         items=['Implicit', 'Explicit'])
        self.ui_solvent_cfg = tk.Button(self.canvas, text='Configure')

        model_grid = [['Optimize', self.ui_optimization, self.ui_frequencies],
                      ['Calculation', self.ui_calculation, self.ui_layers],
                      ['Solvation', self.ui_solvent, self.ui_solvent_cfg]]
        self.auto_grid(self.ui_model_frame, model_grid)

        # QM configuration
        self.ui_qm_frame = tk.LabelFrame(self.canvas, text='QM Settings')
        self.ui_qm_methods = Pmw.OptionMenu(self.canvas,
                                            menubutton_textvariable=self.var_qm_method,
                                            items=self.QM_METHODS)
        self.ui_qm_functional_type = Pmw.OptionMenu(self.canvas,
                                                    menubutton_textvariable=self.var_qm_functional_type,
                                                    items=sorted(self.QM_FUNCTIONALS.keys()))
        self.ui_qm_functionals = Pmw.OptionMenu(self.canvas,
                                                menubutton_textvariable=self.var_qm_functional,
                                                items=self.QM_FUNCTIONALS['Pure'])
        self.ui_qm_basis = Pmw.OptionMenu(self.canvas,
                                          menubutton_textvariable=self.var_qm_basis,
                                          items=self.QM_BASIS_SETS)
        self.ui_qm_basis_ext = Pmw.OptionMenu(self.canvas,
                                              menubutton_textvariable=self.var_qm_basis_ext,
                                              items=self.QM_BASIS_SETS_EXT)
        self.ui_qm_basis_per_atom = tk.Button(self.canvas, text='Per-element')
        self.ui_qm_basis_custom = tk.Entry(self.canvas, textvariable=self.var_qm_basis_custom)
        self.ui_qm_keywords = tk.Entry(self.canvas, textvariable=self.var_qm_keywords)

        qm_grid = [['Method', (self.ui_qm_methods, 'Functional', self.ui_qm_functional_type, self.ui_qm_functionals)],
                   ['Basis set', (self.ui_qm_basis, self.ui_qm_basis_ext, self.ui_qm_basis_custom, self.ui_qm_basis_per_atom)],
                   ['Keywords', self.ui_qm_keywords]]
        self.auto_grid(self.ui_qm_frame, qm_grid)

        # MM Configuration
        self.ui_mm_frame = tk.LabelFrame(self.canvas, text='MM Settings')
        self.ui_mm_forcefields = Pmw.OptionMenu(self.canvas,
                                                menubutton_textvariable=self.var_mm_forcefield,
                                                items=self.MM_FORCEFIELDS['General'])
        self.ui_mm_water_forcefield = Pmw.OptionMenu(self.canvas,
                                                menubutton_textvariable=self.var_mm_water_forcefield,
                                                items=self.MM_FORCEFIELDS['Water'])
        self.ui_mm_frcmod = tk.Entry(self.canvas, textvariable=self.var_mm_frcmod)

        mm_grid = [['Forcefield', self.ui_mm_forcefields],
                   ['Waters', self.ui_mm_water_forcefield],
                   ['Frcmod', self.ui_mm_frcmod]]
        self.auto_grid(self.ui_mm_frame, mm_grid)

        # Flexibility & Restraints
        self.ui_flex_frame = tk.LabelFrame(self.canvas, text='Flexibility & Restraints')
        self.ui_flex_policy = Pmw.OptionMenu(self.canvas,
                                             menubutton_textvariable=self.var_flex_policy,
                                             items=['flexible', 'fixed'])
        self.ui_flex_btn = tk.Button(self.canvas, text='Configure')
        self.ui_redundant = tk.Checkbutton(self.canvas, text='Enable redundant coordinates',
                                           variable=self.var_redundant)
        self.ui_redundant_btn = tk.Button(self.canvas, text='Configure')

        flex_grid = [[('All atoms are', self.ui_flex_policy, 'except'), self.ui_flex_btn],
                     [self.ui_redundant, self.ui_redundant_btn]]
        self.auto_grid(self.ui_flex_frame, flex_grid, resize_columns=(0,), label_sep='...')

        # Charges & multiplicity
        self.ui_charges_frame = tk.LabelFrame(self.canvas, text='Charges & Multiplicity')
        self.ui_charges_auto = tk.Button(self.canvas, text='Automatic')
        self.ui_charges_manual = tk.Button(self.canvas, text='Manual')
        self.ui_charges_qm = tk.Entry(self.canvas, textvariable=self.var_charge_qm, width=5)
        self.ui_charges_mm = tk.Entry(self.canvas, textvariable=self.var_charge_mm, width=5)
        self.ui_multiplicity_qm = tk.Entry(self.canvas, textvariable=self.var_multiplicity_qm, width=5)
        self.ui_multiplicity_mm = tk.Entry(self.canvas, textvariable=self.var_multiplicity_mm, width=5)

        charges_grid = [['Set charges:', self.ui_charges_auto, self.ui_charges_manual],
                        ['Total charge:', (self.ui_charges_qm, '(QM)'), (self.ui_charges_mm, '(MM)')],
                        ['Multiplicity:', (self.ui_multiplicity_qm, '(QM)'), (self.ui_multiplicity_mm, '(MM)')]]
        self.auto_grid(self.ui_charges_frame, charges_grid, resize_columns=(1,2), label_sep='')

        # Hardware
        self.ui_hw_frame = tk.LabelFrame(self.canvas, text='Output')
        self.ui_title = tk.Entry(self.canvas, textvariable=self.var_title)
        self.ui_title_btn = tk.Button(self.canvas, text='Random')
        self.ui_checkpoint = tk.Entry(self.canvas, textvariable=self.var_checkpoint)
        self.ui_checkpoint_btn = tk.Button(self.canvas, text='Browse')
        self.ui_nproc = tk.Entry(self.canvas, textvariable=self.var_nproc, width=5)
        self.ui_memory = tk.Entry(self.canvas, textvariable=self.var_memory, width=5)
        self.ui_memory_units = Pmw.OptionMenu(self.canvas,
                                              menubutton_textvariable=self.var_memory_units,
                                              items=self.MEM_UNITS)
        hw_grid = [['Title', self.ui_title, self.ui_title_btn,
                    ('# CPUs', self.ui_nproc)],
                   ['Checkpoint', self.ui_checkpoint, self.ui_checkpoint_btn, 
                    ('Memory', self.ui_memory, self.ui_memory_units)]]
        self.auto_grid(self.ui_hw_frame, hw_grid, sticky='news')

        # Live output
        self.ui_preview_frame = tk.LabelFrame(self.canvas, text='Preview output')
        self.ui_preview = Pmw.ScrolledText(self.canvas, text_state='disabled',
                                           text_padx=4, text_pady=4, usehullsize=True,
                                           hull_width=300, hull_height=200,)
        self.ui_preview.pack(in_=self.ui_preview_frame, expand=True, fill='both', padx=5, pady=5)

        frames = [[self.ui_molecule_frame, self.ui_model_frame],
                  [self.ui_qm_frame, self.ui_mm_frame],
                  [self.ui_flex_frame, self.ui_charges_frame] ]
        self.auto_grid(self.canvas, frames, resize_columns=(0,1), sticky='news')
        self.ui_hw_frame.grid(row=len(frames), columnspan=2, sticky='ew', padx=5, pady=5)
        self.ui_preview_frame.grid(row=len(frames)+1, columnspan=2, sticky='ew', padx=5, pady=5)

    def Save(self):
        """
        Default! Triggered action if you click on an Apply button
        """
        pass

    def Close(self):
        """
        Default! Triggered action if you click on the Close button
        """
        global ui
        ui = None
        ModelessDialog.Close(self)
        self.destroy()

    # Below this line, implement all your custom methods for the GUI.
    def auto_grid(self, parent, grid, resize_columns=(1,), label_sep=':', **options):
        """
        Auto grid an ordered matrix of Tkinter widgets.

        Parameters
        ----------
        parent : tk.Widget
            The widget that will host the widgets on the grid
        grid : list of list of tk.Widget
            A row x columns matrix of widgets. It is built on lists.
            Each list in the toplevel list represents a row. Each row
            contains widgets, tuples or strings, in column order.  
            If it's a widget, it will be grid at the row i (index of first level
            list) and column j (index of second level list).
            If a tuple of widgets is found instead of a naked widget,
            they will be packed in a frame, and grid'ed as a single cell.
            If it's a string, a Label will be created with that text, and grid'ed. 

            For example:
            >>> grid = [['A custom label', widget_0_1, widget_0_2], # first row
            >>>         [widget_1_0, widget_1_1, widget_1_2],       # second row
            >>>         [widget_2_0, widget_2_1, (widgets @ 2_2)]]  # third row

        """
        for column in resize_columns:
            parent.columnconfigure(column, weight=int(100/len(resize_columns)))
        _kwargs = {'padx': 2, 'pady': 2, 'ipadx': 2, 'ipady': 2}
        _kwargs.update(options)
        for i, row in enumerate(grid):
            for j, item in enumerate(row):
                kwargs = _kwargs.copy()
                sticky = 'ew'
                if isinstance(item, tuple):
                    frame = tk.Frame(parent)
                    self.auto_pack(frame, item, side='left', padx=2, pady=2, expand=True, fill='x',
                                   label_sep=label_sep)
                    item = frame
                elif isinstance(item, basestring):
                    sticky = 'e'
                    item = tk.Label(parent, text=item + label_sep if item else '')
                elif isinstance(item, tk.Checkbutton):
                    sticky = 'w'
                if 'sticky' not in kwargs:
                    kwargs['sticky'] = sticky
                item.grid(in_=parent, row=i, column=j, **kwargs)
                self._fix_styles(item)

    def auto_pack(self, parent, widgets, label_sep=':', **kwargs):
        for widget in widgets:
            options = kwargs.copy()
            if isinstance(widget, basestring):
                widget = tk.Label(parent, text=widget + label_sep if widget else '')
            if isinstance(widget, (tk.Button, tk.Label)):
                options['expand'] = False
            widget.pack(in_=parent, **options)
            self._fix_styles(widget)

if __name__ == '__main__':
    master = tk.Tk()
    master.mainloop()
    showUI(parent=master)