from pathlib import Path

import networkx as nx

from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

from pysmiles import read_smiles, remove_explicit_hydrogens, add_explicit_hydrogens

from vermouth.processors import MakeBonds
from vermouth.system import System
from vermouth.pdb import read_pdb

class MoleculeWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        layout = QVBoxLayout(self)

        file_layout = QHBoxLayout()
        self._pth_widget = QLineEdit()
        browse_button = QPushButton('Browse')
        browse_button.clicked.connect(self._select_file)
        file_layout.addWidget(QLabel('PDB file: '))
        file_layout.addWidget(self._pth_widget)
        file_layout.addWidget(browse_button)

        smiles_layout = QHBoxLayout()
        smiles_layout.addWidget(QLabel('SMILES: '))
        self._smiles_widget = QLineEdit()
        smiles_layout.addWidget(self._smiles_widget)

        self.hydrogen_checkbox = QCheckBox('Keep hydrogen atoms')

        layout.addLayout(file_layout)
        layout.addLayout(smiles_layout)
        layout.addWidget(self.hydrogen_checkbox)

        self._pth_widget.setText('/home/peterkroon/python/molecules/genistein.pdb')
        self._smiles_widget.setText('C1=CC(=CC=C1C2=COC3=CC(=CC(=C3C2=O)O)O)O')

    def _select_file(self):
        filename = QFileDialog.getOpenFileName(filter="PDB file (*.pdb)")
        filename = filename[0]
        self._pth_widget.setText(filename)

    def get_value(self):
        filename = self._pth_widget.text()
        if filename:
            try:
                pdb_mol = read_pdb(filename)
            except Exception as err:
                self._pth_widget.setText('')
                dialog = QErrorMessage()
                dialog.showMessage(str(err))
                dialog.exec_()
                return False
            pdb_mol = pdb_mol[0]
            pdb_mol.graph['name'] = Path(filename).stem
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
                dialog = QErrorMessage()
                dialog.showMessage(str(err))
                dialog.exec_()
                self._smiles_widget.setText('')
                return False
            else:
                smiles_mol.graph['smiles'] = smiles
                smiles_mol.graph['name'] = smiles
        else:
            smiles_mol = None

        if pdb_mol and smiles_mol:
            if self.hydrogen_checkbox.checkState():
                add_explicit_hydrogens(smiles_mol)
            else:
                remove_explicit_hydrogens(pdb_mol)
            gm = nx.isomorphism.GraphMatcher(pdb_mol, smiles_mol,
                                             nx.isomorphism.categorical_node_match('element', None))
            match = next(gm.isomorphisms_iter(), {})
            if not match:
                dialog = QErrorMessage()
                dialog.showMessage('Smiles and PDB molecule are not isomorphic!')
                dialog.exec_()
                return False
            for pdb_idx, smi_idx in match.items():
                smiles_mol.nodes[smi_idx].update(pdb_mol.nodes[pdb_idx])
            smiles_mol.graph.update(pdb_mol.graph)

        molecule = smiles_mol or pdb_mol
        if not molecule:
            return False
        return molecule
