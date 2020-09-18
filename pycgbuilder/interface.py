from matplotlib.backend_bases import MouseButton
from matplotlib.figure import Figure
from matplotlib.cm import get_cmap
from matplotlib.backends.qt_compat import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)

from vermouth.processors import MakeBonds
from vermouth.system import System
from vermouth.pdb import read_pdb
from pysmiles import read_smiles, remove_explicit_hydrogens
import networkx as nx

from .embed_molecule import vsepr_layout
from .draw_mol import draw_molecule

class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.molecule = nx.Graph()
        self.mapping = {}

        self._main = QtWidgets.QWidget(self)
        self.setCentralWidget(self._main)

        self.page_builders = [self.molecule_selector, self.molecule_mapper]
        self.validators = [self.validate_mol, self.validate_mapping]
        self.pages = QtWidgets.QStackedWidget(self)
        for builder in self.page_builders:
            self.pages.addWidget(builder())

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.pages)

        button_bar = QtWidgets.QHBoxLayout()
        self.back_button = QtWidgets.QPushButton('Back')
        self.next_button = QtWidgets.QPushButton('Next')
        self.back_button.clicked.connect(self._prev_page)
        self.next_button.clicked.connect(self._next_page)
        # self.back_button.setEnabled(False)
        # self.next_button.setEnabled(False)
        button_bar.addWidget(self.back_button)
        button_bar.addWidget(self.next_button)

        layout.addLayout(button_bar)
        self._main.setLayout(layout)
        self.show()

    def _next_page(self):
        idx = self.pages.currentIndex()
        if self.validators[idx]():
            self.pages.setCurrentIndex(idx + 1)

    def _prev_page(self):
        idx = self.pages.currentIndex()
        self.pages.setCurrentIndex(idx - 1)

    def molecule_selector(self):
        widg = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widg)

        file_layout = QtWidgets.QHBoxLayout()
        self._pth_widget = QtWidgets.QLineEdit()
        browse_button = QtWidgets.QPushButton('Browse')
        browse_button.clicked.connect(self._select_file)
        file_layout.addWidget(QtWidgets.QLabel('PDB file: '))
        file_layout.addWidget(self._pth_widget)
        file_layout.addWidget(browse_button)

        smiles_layout = QtWidgets.QHBoxLayout()
        smiles_layout.addWidget(QtWidgets.QLabel('SMILES: '))
        self._smiles_widget = QtWidgets.QLineEdit()
        smiles_layout.addWidget(self._smiles_widget)

        layout.addLayout(file_layout)
        layout.addLayout(smiles_layout)
        self._pth_widget.setText('/home/peterkroon/python/molecules/genistein.pdb')
        self._smiles_widget.setText('C1=CC(=CC=C1C2=COC3=CC(=CC(=C3C2=O)O)O)O')

        return widg

    def _select_file(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(filter="PDB file (*.pdb)")
        filename = filename[0]
        self._pth_widget.setText(filename)

    def validate_mol(self):
        filename = self._pth_widget.text()
        if filename:
            try:
                pdb_mol = read_pdb(filename)
            except Exception as err:
                self._pth_widget.setText('')
                print(err)
                return False
            pdb_mol = pdb_mol[0]
            if not pdb_mol.edges:
                system = System()
                system.add_molecule(pdb_mol)
                MakeBonds(allow_name=False).run_system(system)
        else:
            pdb_mol = None
        smiles = self._smiles_widget.text()
        if smiles:
            try:
                smiles_mol = read_smiles(smiles)
            except Exception as err:
                # self._smiles_widget.setText('')
                print(err)
                return False
        else:
            smiles_mol = None

        if pdb_mol and smiles_mol:
            remove_explicit_hydrogens(pdb_mol)
            gm = nx.isomorphism.GraphMatcher(pdb_mol, smiles_mol,
                                             nx.isomorphism.categorical_node_match('element', None))
            match = next(gm.isomorphisms_iter(), {})
            if not match:
                print('Smiles and PDB molecule are not isomorphic!')
                return False
            for pdb_idx, smi_idx in match.items():
                smiles_mol.nodes[smi_idx].update(pdb_mol.nodes[pdb_idx])
                self.molecule = smiles_mol
        self.molecule = pdb_mol or smiles_mol
        if not self.molecule:
            self.molecule = None
            return False
        return True

    def molecule_mapper(self):
        widg = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widg)

        self.figure = Figure()
        canvas = FigureCanvas(self.figure)
        canvas.mpl_connect('button_press_event', self._click_canvas)
        ax = canvas.figure.subplots()
        # self.addToolBar(NavigationToolbar(canvas, self))
        layout.addWidget(canvas)
        print(repr(self.molecule))
        pos = vsepr_layout(self.molecule)
        draw_molecule(self.molecule, pos=pos, ax=ax)
        ax.autoscale(True)
        return widg

    def validate_mapping(self):
        return True

    def _click_canvas(self, mpl_event):
        pass

if __name__ == '__main__':
    import sys
    qapp = QtWidgets.QApplication(sys.argv)
    app = ApplicationWindow()
    app.show()
    qapp.exec_()
