from collections import defaultdict, Counter

from matplotlib.backend_bases import MouseButton
from matplotlib.figure import Figure
from matplotlib.cm import get_cmap
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

from .embed_molecule import vsepr_layout
from .draw_mol import draw_molecule

import networkx as nx

EMBEDDINGS = {
    'Kamada Kawai': nx.kamada_kawai_layout,
    'VSEPR': vsepr_layout,
    'Spring': nx.spring_layout,
    'Spectral': nx.spectral_layout,
    'Planar': nx.planar_layout,
}


class MappingView(FigureCanvas):
    def __init__(self, figure, atom_radius=0.05):
        super().__init__(figure)
        self.ax = figure.subplots()
        self.mpl_connect('button_press_event', self._click_canvas)
        self.atom_radius = atom_radius
        self._drawn_mapping = None
        self._embeddings = {}
        self._current_embedding = ''
        self._molecule = nx.Graph()
        self._model = MappingModel(self._molecule)
        self._selectionmodel = QItemSelectionModel(self._model)

    @property
    def mapping(self):
        members = defaultdict(list)
        for bd_idx in range(self.model.rowCount(self.model.index(0, 0))):
            atoms = self.model.index(bd_idx, 2).data(role=Qt.UserRole)
            if not atoms:
                continue
            for atom in atoms:
                members[atom].append(bd_idx)
        members = dict(members)
        return members

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, qtmodel):
        self._model = qtmodel
        self.molecule = self._model.molecule
        self._model.dataChanged.connect(self.redraw)

    def setModel(self, qtmodel):
        self.model = qtmodel

    def setSelectionModel(self, qtselectionmodel):
        self._selectionmodel = qtselectionmodel
        self._selectionmodel.selectionChanged.connect(self.redraw)

    def selectionModel(self):
        return self._selectionmodel

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, new_mol):
        self._molecule = new_mol.copy()
        self._embeddings.clear()
        self.draw_molecule()
        self.draw_mapping()

    @property
    def current_embedding(self):
        return self._current_embedding

    @current_embedding.setter
    def current_embedding(self, embedding_name):
        self._current_embedding = embedding_name
        self.draw_molecule()
        self.draw_mapping()

    def _set_embedding(self, name):
        self.current_embedding = name

    def _make_embedding(self):
        name = self.current_embedding
        if name in self._embeddings:
            return self._embeddings[name]
        mol = self._molecule
        method = EMBEDDINGS[name]
        new_pos = method(mol)
        self._embeddings[name] = new_pos
        return new_pos

    @property
    def embedding(self):
        return self._make_embedding()

    def redraw(self, *args):
        # changed = bottomright.data(role=Qt.UserRole)
        # self.draw_mapping({k: v for k, v in self.mapping.items() if k in changed})
        self.draw_molecule()
        self.draw_mapping()

    def draw_molecule(self):
        self.ax.clear()
        draw_molecule(self.molecule, pos=self.embedding, ax=self.ax,
                      labels=nx.get_node_attributes(self.molecule, 'atomname'))
        self.ax.set_aspect(1)
        self.ax.autoscale(True)
        self.ax.set_axis_off()
        self.draw_idle()

    def draw_mapping(self, mapping=None):
        pos = self.embedding
        mapping = mapping or self.mapping

        for idx in mapping:
            position = pos[idx]
            member = mapping[idx]
            counts = Counter(member)
            sizes = [counts[idx] for idx in range(max(member) + 1)]
            colors = [self.model.index(bd_idx, 0).data(role=Qt.DecorationRole).getRgbF()
                      for bd_idx in range(len(sizes))]
            pie = self.ax.pie(sizes, center=position, radius=self.atom_radius, normalize=True,
                              colors=colors)
            selected = self._selected_bead()
            for bd_idx, wedge in enumerate(pie[0]):
                if bd_idx == selected:
                    wedge.set_linewidth(1)
                    wedge.set_linestyle('-')
                    wedge.set_edgecolor('black')

        self.ax.autoscale(True)
        self.draw_idle()

    def _click_canvas(self, mpl_event):
        try:
            xpress = mpl_event.xdata
            ypress = mpl_event.ydata
        except AttributeError:
            return
        if not xpress or not ypress:
            return
        for n_idx, pos in self.embedding.items():
            if (pos[0] - xpress) ** 2 + (pos[1] - ypress) ** 2 <= self.atom_radius ** 2:
                # Clicked node n_idx
                if mpl_event.button == MouseButton.RIGHT:
                    self._select_atom(n_idx)
                elif mpl_event.button == MouseButton.LEFT:
                    self._map(n_idx)

    def _select_atom(self, n_idx):
        rev_mapping = self.model.reverse_mapping
        members = rev_mapping.get(n_idx, [])
        if not members:
            # Select new bead
            bd_idx = self.model.rowCount(0) - 1
        elif len(members) == 1:
            bd_idx = members[0]
        else:
            current = self._selected_bead()
            try:
                cur_idx = members.index(current)
            except ValueError:
                bd_idx = members[0]
            else:
                bd_idx = members[(cur_idx + 1) % len(members)]

        idx = self.model.index(bd_idx, 2)
        self._selectionmodel.select(idx, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)

    def _map(self, n_idx):
        bd_idx = self._selected_bead()
        if bd_idx == -1:
            return
        index = self.model.index(bd_idx, 2)
        self.model.setData(index, n_idx, Qt.UserRole)

    def _selected_bead(self):
        bd_idxs = set(idx.row() for idx in self._selectionmodel.selection().indexes())
        if len(bd_idxs) != 1:
            return -1
        return bd_idxs.pop()


class MappingModel(QAbstractTableModel):
    def __init__(self, molecule, mapping=None, names=None, types=None):
        super().__init__()
        self.molecule = molecule

        self._name_to_idx = defaultdict(list)
        for idx in self.molecule:
            name = self.molecule.nodes[idx]['atomname']
            self._name_to_idx[name].append(idx)
        self._name_to_idx = dict(self._name_to_idx)
        if any(len(names) != 1 for names in self._name_to_idx.values()):
            self._atom_flags = Qt.NoItemFlags
        else:
            self._atom_flags = Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable

        self.mapping = mapping or []
        self.names = names or []
        self.types = types or []

    @property
    def reverse_mapping(self):
        out = defaultdict(list)
        for bd_idx, at_idxs in enumerate(self.mapping):
            for at_idx in at_idxs:
                out[at_idx].append(bd_idx)
        return dict(out)

    def data(self, index, role):
        row = index.row()
        col = index.column()
        if role == Qt.DisplayRole or role == Qt.EditRole:
            if row >= len(self.mapping):
                return ''
            if col == 0:
                # Bead name
                return self.names[row]
            elif col == 1:
                # Bead type
                return self.types[row]
            elif col == 2:
                # atom members
                return ' '.join(self.molecule.nodes[idx]['atomname'] for idx in self.mapping[row])
        elif role == Qt.UserRole:
            if row >= len(self.mapping):
                return None
            if col == 0:
                # Bead name
                return self.names[row]
            elif col == 1:
                # Bead type
                return self.types[row]
            elif col == 2:
                # atom members
                return self.mapping[row]
        elif role == Qt.BackgroundRole:
            # row to color
            cmap = get_cmap('tab20')
            c_idx = row*2 + 1  # Row 0 gets color 1, row 1 gets color 3, etc.
            return QColor.fromRgbF(*cmap.colors[c_idx % cmap.N])
        elif role == Qt.DecorationRole and col == 0:
            cmap = get_cmap('tab10')
            return QColor.fromRgbF(*cmap.colors[row % cmap.N])

    def rowCount(self, index):
        return len(self.mapping)+1

    def columnCount(self, index):
        return 3

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            names = ['Name', 'Type', 'Atoms']
            return names[section]

    def setData(self, index, value, role):
        row = index.row()
        col = index.column()
        if row == len(self.names) and value != '':
            self.layoutAboutToBeChanged.emit([QPersistentModelIndex(index.siblingAtRow(row+1))])
            self.names.append('BD{}'.format(row))
            self.types.append('__')
            self.mapping.append([])
            self.layoutChanged.emit([QPersistentModelIndex(index.siblingAtRow(row+1))])

        if role == Qt.EditRole:
            value = value.strip()

            if col == 0:
                # Bead name
                self.names[row] = value
            elif col == 1:
                # Bead type
                self.types[row] = value
            elif col == 2:
                # Atom names
                idxs = []
                for atname in value.split():
                    if atname not in self._name_to_idx:
                        dialog = QErrorMessage()
                        dialog.showMessage('Atom with name {} not found'.format(atname))
                        dialog.exec_()
                        return False
                    idxs.extend(self._name_to_idx[atname])
                self.mapping[row] = sorted(idxs)
            self.dataChanged.emit(index, index, [role])
        elif role == Qt.UserRole and col == 2:
            if value in self.mapping[row]:
                self.mapping[row].remove(value)
            else:
                self.mapping[row] = sorted(self.mapping[row] + [value])
            self.dataChanged.emit(index, index, [role])
        return True

    def flags(self, index):
        if index.column() == 2:
            return self._atom_flags
        else:
            return Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable

class MappingWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._molecule = nx.Graph()
        self._mapping = MappingModel(self._molecule)

        layout = QHBoxLayout()
        canvas_layout = QVBoxLayout()
        self.embeddings_box = QComboBox()
        self.embeddings_box.addItems(EMBEDDINGS.keys())
        self.embeddings_box.setEditable(False)
        self.embeddings_box.currentTextChanged.connect(self._set_embedding)

        canvas_layout.addWidget(self.embeddings_box)

        self.figure = Figure()
        self.canvas = MappingView(self.figure)

        canvas_layout.addWidget(self.canvas)

        layout.addLayout(canvas_layout)

        self._table = QTableView()
        self._table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self._table)

        self.setLayout(layout)
        self._set_embedding(self.embeddings_box.currentText())

    def _set_embedding(self, name):
        self.canvas.current_embedding = name

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, new_mol):
        self._molecule = new_mol.copy()
        for idx in self._molecule:
            node = self._molecule.nodes[idx]
            node['atomname'] = node.get('atomname', '{}{}'.format(node['element'], idx))
        self._mapping.molecule = self._molecule
        self._mapping = MappingModel(self._molecule)
        self._table.setModel(self._mapping)
        self.canvas.setModel(self._mapping)
        self.canvas.setSelectionModel(self._table.selectionModel())
        self._set_embedding(self.embeddings_box.currentText())

    def set_value(self, value):
        self.molecule = value

    def get_value(self):
        return self._mapping.names, self._mapping.types, self._mapping.mapping, self.molecule
