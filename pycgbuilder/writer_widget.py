from collections import defaultdict
from functools import partial
from itertools import product
from pathlib import Path

import networkx as nx

import numpy as np

from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

from vermouth.molecule import Molecule
from vermouth.system import System
from vermouth.pdb.pdb import write_pdb_string
from vermouth.gmx import write_molecule_itp
from vermouth.file_writer import open, DeferredFileWriter


def make_cg_mol(aa_mol, mapping, bead_names, bead_types):
    molname = aa_mol.graph['name']
    cg_mol = Molecule(nrexcl=1, meta=dict(moltype=molname))
    mapdict = defaultdict(list)
    for bd_idx, at_idxs in enumerate(mapping):
        name = bead_names[bd_idx]
        for member in at_idxs:
            mapdict[member].append(bd_idx)
        subgraph = aa_mol.subgraph(at_idxs)
        charge = sum(nx.get_node_attributes(subgraph, 'charge').values())
        position = np.mean([
            subgraph.nodes[idx].get('position', (np.nan, np.nan, np.nan))
            for idx in subgraph
        ], axis=0)
        cg_mol.add_node(bd_idx, atomname=name, resname=molname, resid=1,
                        atype=bead_types[bd_idx], charge_group=bd_idx+1, graph=subgraph,
                        charge=charge, position=position)
    for aa_idx, aa_jdx in aa_mol.edges:
        cg_idxs = mapdict[aa_idx]
        cg_jdxs = mapdict[aa_jdx]
        for cg_idx, cg_jdx in product(cg_idxs, cg_jdxs):
            if cg_idx != cg_jdx:
                cg_mol.add_edge(cg_idx, cg_jdx)
    for idx, jdx in cg_mol.edges:
        cg_mol.add_interaction('bonds', [idx, jdx], [])
    return cg_mol


def write_ndx(filename, cg_mol, stepsize=10):
    with open(filename, 'w') as file_out:
        for bead_idx in cg_mol:
            node = cg_mol.nodes[bead_idx]
            at_idxs = list(node.get('graph', []))
            file_out.write('[ {} ]\n'.format(node['atomname']))
            for idx in range(0, len(at_idxs), stepsize):
                idxs = (at_idx + 1 for at_idx in at_idxs[idx:idx+stepsize])
                file_out.write(' '.join(map(str, idxs)) + '\n')
            file_out.write('\n')


def write_map(filename, cg_mol):
    aa_nodes = {}
    aa_to_cg = defaultdict(list)
    for cg_idx in cg_mol:
        bead = cg_mol.nodes[cg_idx]
        aa_graph = bead['graph']
        for aa_idx in aa_graph:
            atom = aa_graph.nodes[aa_idx]
            aa_nodes[aa_idx] = atom
            aa_to_cg[aa_idx].append(cg_idx)
    with open(filename, 'w') as file_out:
        molname = cg_mol.meta['moltype']
        file_out.write('[ molecule ]\n')
        file_out.write(molname + '\n')
        file_out.write('[ martini ]\n')
        file_out.write(' '.join(cg_mol.nodes[idx]['atomname'] for idx in cg_mol) + '\n')
        file_out.write('[ mapping ]\n')
        file_out.write('<FORCE FIELD NAMES>\n')
        file_out.write('[ atoms ]\n')
        for aa_idx, atom in sorted(aa_nodes.items(), key=lambda i: i[1]['atomid']):
            cg_idxs = aa_to_cg[aa_idx]
            file_out.write('{} {} {}\n'.format(
                atom['atomid'],
                atom['atomname'],
                ' '.join(cg_mol.nodes[cg_idx]['atomname'] for cg_idx in cg_idxs)
            ))


def write_itp(path, cg_mol):
    with open(path, 'w') as out:
        write_molecule_itp(cg_mol, out)


def write_pdb(path, cg_mol):
    system = System()
    system.add_molecule(cg_mol)
    with open(path, 'w') as out:
        out.write(write_pdb_string(system))


WRITERS = {
    'ndx': write_ndx,
    'pdb': write_pdb,
    'itp': write_itp,
    'map': write_map,
}


class WriterWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        layout = QVBoxLayout(self)
        self.widgets = {}

        for ext in WRITERS:
            line = QHBoxLayout()
            check = QCheckBox()
            check.setCheckState(Qt.Checked)
            label = QLabel(ext.upper())
            pth = QLineEdit()
            pth.textChanged.connect(self._update_pths)
            button = QPushButton('Browse')
            button.clicked.connect(partial(self._browse, ext=ext, line=pth))
            widgets = [check, label, pth, button]
            for widget in widgets:
                line.addWidget(widget)
            self.widgets[ext] = widgets
            layout.addLayout(line)

        # TODO: Uncheck PDB writer if no positions

    def set_value(self, value):
        names, types, mapping, mol = value
        self.mapping = mapping
        self.bead_names = names
        self.bead_types = types
        self.aa_molecule = mol
        self.nativeParentWidget().set_final_page()
        self.nativeParentWidget().next.setEnabled(False)

    def _browse(self, ext, line):
        cur = Path(line.text())
        stem = cur.stem or self.aa_molecule.graph['name']
        filename = QFileDialog.getSaveFileName(
            directory='{}/{}.{}'.format(cur.parent, stem, ext),
            filter="{} file (*.{})".format(ext.upper(), ext)
        )
        filename = filename[0]
        line.setText(filename)

    def _update_pths(self, val):
        new_val = Path(val)
        if not new_val:
            return
        self.nativeParentWidget().next.setEnabled(True)
        for ext, widgs in self.widgets.items():
            line = widgs[2]
            if not line.text():
                line.setText('{}/{}.{}'.format(new_val.parent, new_val.stem, ext))

    def do_write(self):
        cg_mol = make_cg_mol(self.aa_molecule, self.mapping, self.bead_names, self.bead_types)
        for ext, widgs in self.widgets.items():
            check = widgs[0]
            path = widgs[2].text()
            if path and check.checkState() and WRITERS[ext]:
                WRITERS[ext](path, cg_mol)
        DeferredFileWriter().write()

    def get_value(self):
        self.do_write()
