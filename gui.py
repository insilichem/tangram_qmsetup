#!/usr/bin/env python
# encoding: utf-8

# Get used to importing this in your Py27 projects!
from __future__ import print_function, division 
# Python stdlib
import Tkinter as tk
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
def showUI(callback=None):
    """
    Requested by Chimera way-of-doing-things
    """
    if chimera.nogui:
        tk.Tk().withdraw()
    global ui
    if not ui: # Edit this to reflect the name of the class!
        ui = CauchianDialog()
    model = Model()
    controller = Controller(gui=ui, model=model)
    ui.enter()
    if callback:
        ui.addCallback(callback)

ENTRY_STYLE = {
    'background': 'white',
    'borderwidth': 1,
    'highlightthickness': 0,
    'insertwidth': 1,
}
BUTTON_STYLE = {
    'borderwidth': 1,
    'highlightthickness': 0,
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

    QM_METHODS = ['']
    QM_FUNCTIONALS = ['']
    QM_BASIS_SETS = ['']
    MM_FORCEFIELDS = ['']
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
        self.var_qm_basis = tk.StringVar()

        # MM variables
        self.var_mm_forcefield = tk.StringVar()
        self.var_mm_frcmod = tk.StringVar()
        self.var_mm_waters = tk.StringVar()

        # Charges & Multiplicity
        self.var_charges_method = tk.DoubleVar()
        self.var_total_charge = tk.DoubleVar()
        self.var_multiplicity = tk.IntVar()

        # Hardware & Output variables
        self.var_title = tk.StringVar()
        self.var_nproc = tk.IntVar()
        self.var_memory = tk.IntVar()
        self.var_memory_units = tk.StringVar()

        # Fire up
        ModelessDialog.__init__(self)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

        # Fix styles
        self._fix_styles()

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def _fix_styles(self):
        for name, btn in self.buttonWidgets.items():
            btn.configure(**BUTTON_STYLE)

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

        # Job type
        self.ui_job_frame = tk.LabelFrame(self.canvas, text='Choose job modes')
        self.ui_optimization = tk.OptionMenu(self.ui_job_frame, self.var_optimization, "Don't optimize", 'Min', 'TS')
        self.ui_frequencies = tk.Checkbutton(self.ui_job_frame, variable=self.var_frequencies)

        job_grid = [['Optimize', self.ui_optimization],
                    ['Frequencies', self.ui_frequencies]]
        self.fill_grid(self.ui_job_frame, job_grid)

        # Modelization
        self.ui_model_frame = tk.LabelFrame(self.canvas, text='Modelization')
        self.ui_calculation = tk.OptionMenu(self.ui_model_frame, self.var_calculation, 'QM', 'ONIOM')
        self.ui_layers = tk.Button(self.ui_model_frame, text='Define layers')
        self.ui_solvent = tk.OptionMenu(self.ui_model_frame, self.var_solvent, 'Implicit', 'Explicit')
        self.ui_solvent_cfg = tk.Button(self.ui_model_frame, text='Configure')

        model_grid = [['Calculation', self.ui_calculation, self.ui_layers],
                      ['Solvation', self.ui_solvent, self.ui_solvent_cfg]]
        self.fill_grid(self.ui_model_frame, model_grid)
        
        # QM configuration
        self.ui_qm_frame = tk.LabelFrame(self.canvas, text='QM Settings')
        self.ui_qm_methods = tk.OptionMenu(self.ui_qm_frame, self.var_qm_method, *self.QM_METHODS)
        self.ui_qm_functionals = tk.OptionMenu(self.ui_qm_frame, self.var_qm_functional, *self.QM_FUNCTIONALS)
        self.ui_qm_basis = tk.OptionMenu(self.ui_qm_frame, self.var_qm_basis, *self.QM_BASIS_SETS)
        self.ui_qm_basis_per_atom = tk.Button(self.ui_qm_frame, text='Per-element basis sets')

        qm_grid = [['Method', self.ui_qm_methods],
                   ['Functional', self.ui_qm_functionals],
                   ['Basis set', [self.ui_qm_basis, self.ui_qm_basis_per_atom]]]
        self.fill_grid(self.ui_qm_frame, qm_grid)

        # MM Configuration
        self.ui_mm_frame = tk.LabelFrame(self.canvas, text='MM Settings')
        self.ui_mm_forcefields = tk.OptionMenu(self.ui_mm_frame, self.var_mm_forcefield, *self.MM_FORCEFIELDS)
        self.ui_mm_frcmod = tk.Entry(self.ui_mm_frame, textvariable=self.var_mm_forcefield)
        self.ui_mm_waters = tk.Entry(self.ui_mm_frame, textvariable=self.var_mm_waters)

        mm_grid = [['Forcefield', self.ui_mm_forcefields],
                   ['Frcmod', self.ui_mm_frcmod],
                   ['Waters', self.ui_mm_waters]]
        self.fill_grid(self.ui_mm_frame, mm_grid)
        
        # Flexibility
        self.ui_flex_frame = tk.LabelFrame(self.canvas, text='Atom flexibility')
        self.ui_flex_btn = tk.Button(self.ui_flex_frame, text='Designate selected atoms as flexible')

        flex_grid = [[self.ui_flex_btn]]
        self.fill_grid(self.ui_flex_frame, flex_grid)

        # Charges & multiplicity
        self.ui_charges_frame = tk.LabelFrame(self.canvas, text='Charges & Multiplicity')
        self.ui_charges_method = tk.OptionMenu(self.ui_charges_frame, self.var_charges_method, 'Auto', 'Advanced')
        self.ui_charges_btn = tk.Button(self.ui_charges_frame, text='Proceed')
        self.ui_charges_total = tk.Entry(self.ui_charges_frame, textvariable=self.var_total_charge)
        self.ui_multiplicity = tk.Entry(self.ui_charges_frame, textvariable=self.var_multiplicity)

        charges_grid = [['Set charges', [self.ui_charges_method, self.ui_charges_btn]],
                        ['Total charge', self.ui_charges_total],
                        ['Multiplicity', self.ui_multiplicity]]
        self.fill_grid(self.ui_charges_frame, charges_grid)

        # Hardware
        self.ui_hw_frame = tk.LabelFrame(self.canvas, text='Hardware & Output')
        self.ui_title = tk.Entry(self.ui_hw_frame, textvariable=self.var_title)
        self.ui_nproc = tk.Entry(self.ui_hw_frame, textvariable=self.var_nproc)
        self.ui_memory = tk.Entry(self.ui_hw_frame, textvariable=self.var_memory)
        self.ui_memory_units = tk.OptionMenu(self.ui_hw_frame, self.var_memory_units, *self.MEM_UNITS)

        hw_grid = [['Title', self.ui_title],
                   ['# CPUs', self.ui_nproc],
                   ['Memory', [self.ui_memory, self.ui_memory_units]]]
        self.fill_grid(self.ui_hw_frame, hw_grid)

        self.ui_job_frame.pack(padx=5, pady=5, fill='x', expand=True)
        self.ui_model_frame.pack(padx=5, pady=5, fill='x', expand=True)
        self.ui_qm_frame.pack(padx=5, pady=5, fill='x', expand=True)
        self.ui_mm_frame.pack(padx=5, pady=5, fill='x', expand=True)
        self.ui_flex_frame.pack(padx=5, pady=5, fill='x', expand=True)
        self.ui_charges_frame.pack(padx=5, pady=5, fill='x', expand=True)
        self.ui_hw_frame.pack(padx=5, pady=5, fill='x', expand=True)        

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
    def fill_grid(self, parent, grid):
        for i, row in enumerate(grid):
            for j, item in enumerate(row):
                sticky = 'we'
                if isinstance(item, basestring):
                    sticky = 'e'
                    item = tk.Label(parent, text=item + ':')
                elif isinstance(item, list):
                    frame = tk.Frame(parent)
                    for widget in item:
                        widget.pack(in_=frame, side='left', padx=5)
                    item = frame
                item.grid(row=i, column=j, sticky=sticky, padx=2, pady=2)
                try:
                    item.configure(**BUTTON_STYLE)
                    item.configure(**ENTRY_STYLE)
                except:
                    pass
