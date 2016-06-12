#!/usr/bin/env python
# encoding: utf-8

# Get used to importing this in your Py27 projects!
from __future__ import print_function, division
# Python stdlib
import Tkinter as tk
import Pmw
# Chimera stuff
import chimera
from chimera.baseDialog import ModelessDialog
# 3rd party
from ebsel import EMSL_local


elements = [["H",  "",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "",    "",    "",   "",    "He" ],
            ["Li", "Be", "",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "B",  "C",   "N",   "O",  "F",   "Ne" ],
            ["Na", "Mg", "",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "Al", "Si",  "P",   "S",  "Cl",  "Ar" ],
            ["K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",  "Ge", "As",  "Se", "Br",  "Kr" ],
            ["Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",  "Sn", "Sb",  "Te", "I",   "Xe" ],
            ["Cs", "Ba", "",   "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",  "Pb", "Bi",  "Po", "At",  "Rn" ],
            ["Fr", "Ra", "",   "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"],
            ["",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "",   "",    "",    "",   "",    ""   ],
            ["",   "",   "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",  "Er", "Tm",  "Yb", "Lu",  ""   ],
            ["",   "",   "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",  "Fm", "Md",  "No", "Lr",  ""   ]]
all_elements = [element for row in elements for element in row if element]


class BasisSetDialog(ModelessDialog):

    buttons = ('OK', 'Close')
    default = None
    help = 'https://www.insilichem.com'

    def __init__(self, saved_basis, parent=None, *args, **kwarg):
        # GUI init
        self.parent = parent
        self.title = 'BasisSet Database'
        self._saved_basis = saved_basis
        self.saved_basis = saved_basis.copy()

        # Variables
        self.var_elements = {e: tk.IntVar() for e in all_elements}
        # Constants
        self.db = EMSL_local(fmt="g94")
        self.db_basissets = sorted([b for (b, d) in self.db.get_available_basis_sets()])

        # Fire up
        ModelessDialog.__init__(self)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def OK(self):
        self._saved_basis.clear()
        self._saved_basis.update(self.saved_basis)
        self.Close()

    def Close(self):
        self.parent._basis_set_dialog = None
        ModelessDialog.Close(self)
        self.destroy()

    def fillInUI(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Create main window
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')
        self.canvas.columnconfigure(1, weight=1)
        self.ui_basis_set_frame = tk.LabelFrame(self.canvas, text='Choose a basis set')
        self.ui_basis_set_frame.grid(rowspan=2, row=0, column=0, sticky='news', pady=5, padx=5)
        self.ui_basis_sets = Pmw.ScrolledListBox(self.ui_basis_set_frame,
                                                 items=self.db_basissets,
                                                 selectioncommand=self._cb_basissets_changed)
        self.ui_basis_sets.pack(fill='y', expand=True, padx=2, pady=5)
        self.ui_basis_set_restore = tk.Button(self.ui_basis_set_frame, text='Reset all',
                                              command=self._reset_all)
        self.ui_basis_set_restore.pack(fill='x', padx=2, pady=5)

        self.ui_periodic_table = tk.LabelFrame(self.canvas, text='Choose elements')
        self.ui_periodic_table.grid(row=0, column=1, columnspan=5, sticky='news', pady=5, padx=5)
        self.ui_elements = {}
        for i, row in enumerate(elements):
            for j, element in enumerate(row):
                if element:
                    w = tk.Checkbutton(self.ui_periodic_table, text=element,
                                       variable=self.var_elements[element],
                                       command=self._cb_elements_changed)
                    self.ui_elements[element] = w
                else:
                    w = tk.Label(self.ui_periodic_table)
                w.grid(row=i, column=j, sticky='w')

        self.ui_output_frame = tk.LabelFrame(self.canvas, text='Basis set')
        self.ui_output_frame.grid(row=1, column=1, columnspan=4, sticky='news', pady=5, padx=5)
        self.ui_output = Pmw.ScrolledText(self.ui_output_frame, text_state='disabled',
                                          text_padx=4, text_pady=4, usehullsize=True,
                                          hull_width=300, hull_height=200, text_font='Courier')
        self.ui_output.pack(expand=True, fill='both')

        self.ui_saved_basis_frame = tk.LabelFrame(self.canvas, text='Your saved basis sets')
        self.ui_saved_basis = Pmw.ScrolledListBox(self.ui_saved_basis_frame,
                                                  items=sorted(self.saved_basis.keys()))
        self.ui_saved_basis_add = tk.Button(self.ui_saved_basis_frame, text='Add current',
                                            command=self._cb_saved_basis_add)
        self.ui_saved_basis_del = tk.Button(self.ui_saved_basis_frame, text='Delete',
                                            command=self._cb_saved_basis_del)
        self.ui_saved_basis_frame.grid(row=1, column=5, sticky='news', pady=5, padx=5)
        self.ui_saved_basis.grid(row=0, column=0, columnspan=2, sticky='news', pady=2, padx=2)
        self.ui_saved_basis_add.grid(row=1, column=0, sticky='we', pady=2, padx=2)
        self.ui_saved_basis_del.grid(row=1, column=1, sticky='we', pady=2, padx=2)

    # Callbacks
    def _cb_basissets_changed(self):
        self._cb_selection_changed()

    def _cb_elements_changed(self):
        self._refresh_basis_sets()
        self._cb_selection_changed()

    def _cb_selection_changed(self):
        self.ui_output.settext("")
        basis = self._selected_basis_set()
        if not basis:
            return
        selected_elem = self._selected_elements()
        supported_elem = self._cb_supported_elements(basis)
        elements = [e for e in selected_elem if e in supported_elem]
        text = self.get_basis_set(basis, elements)
        self.ui_output.settext(text)

    def _cb_supported_elements(self, basis_set=None):
        if basis_set is None:
            basis_set = self._selected_basis_set()
        elements = self.db.get_available_elements(basis_set)
        self._restore_periodic_table()
        for e in elements:
            if e:
                self.ui_elements[e]['fg'] = 'blue'
        return elements

    def _cb_saved_basis_add(self):
        basis_text = self.ui_output.getvalue()
        if basis_text:
            elements = tuple(sorted(self._selected_elements()))
            if not elements:
                elements = ('*',)
            self.saved_basis[elements] = basis_text
            self.ui_saved_basis.setlist(sorted(self.saved_basis.keys()))
        
    def _cb_saved_basis_del(self):
        item = self.ui_saved_basis.getvalue()
        for i in item:
            try:
                del self.saved_basis[tuple(i)]
            except KeyError:
                pass
        self.ui_saved_basis.setlist(sorted(self.saved_basis.keys()))
        
    # Helpers
    def get_basis_set(self, basis_set, elements=()):
        try:
            basis = self.db.get_basis(basis_set, elements=elements)
        except UnboundLocalError:
            return ""
        return '\n'.join([b.replace('****\n', '****\n-') for b in basis])

    def _selected_basis_set(self):
        try:
            basis_set = self.ui_basis_sets.getvalue()[0]
        except IndexError:
            return
        if basis_set != '--None--':
            return basis_set

    def _selected_elements(self):
        return [name for name, var in self.var_elements.iteritems() if var.get()]

    def _restore_periodic_table(self):
        for wid in self.ui_elements.itervalues():
            wid['fg'] = 'black'

    def _refresh_basis_sets(self):
        current = self._selected_basis_set()
        elements = self._selected_elements()
        basis_sets = self.db.get_available_basis_sets(elements=elements)
        basis_sets_names = [b for (b, _) in basis_sets]
        self.ui_basis_sets.setlist(basis_sets_names)
        if current and current in basis_sets_names:
            self.ui_basis_sets.setvalue([current])

    def _reset_all(self):
        for var in self.var_elements.itervalues():
            var.set(0)
        for wid in self.ui_elements.itervalues():
            wid['fg'] = 'black'
        self.ui_basis_sets.setlist(self.db_basissets)
        self.ui_output.settext("")

