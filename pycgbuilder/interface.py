from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

from pysmiles import read_smiles

from .mapping_widget import MappingWidget
from .molecule_widget import MoleculeWidget
from .writer_widget import WriterWidget


class PagedWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        widg = QWidget()
        layout = QVBoxLayout(widg)
        self.pages = QStackedWidget()

        layout.addWidget(self.pages)

        button_layout = QToolBar()
        self.back = QAction('Back', self, icon=self.style().standardIcon(QStyle.SP_ArrowBack))
        self.next = QAction('Next', self, icon=self.style().standardIcon(QStyle.SP_ArrowForward))
        self.back.triggered.connect(self._prev_page)
        self.next.triggered.connect(self._next_page)
        button_layout.addAction(self.back)
        sep = QWidget()
        sep.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        button_layout.addWidget(sep)
        button_layout.addAction(self.next)
        button_layout.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        self.toolbar = button_layout
        # button_layout.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Fixed)
        layout.addWidget(button_layout, alignment=Qt.AlignBottom)

        self._toggle_buttons()
        widg.setLayout(layout)
        self.setCentralWidget(widg)
        # self.setStatusBar(QStatusBar(self))

    def _next_page(self):
        self.pages.setCurrentIndex(self.pages.currentIndex() + 1)
        self._toggle_buttons()

    def _prev_page(self):
        self.pages.setCurrentIndex(self.pages.currentIndex() - 1)
        self._toggle_buttons()

    def _toggle_buttons(self):
        current_idx = self.pages.currentIndex()
        self.back.setEnabled(current_idx != 0)
        self.next.setEnabled(current_idx != self.pages.count() - 1)

    def set_final_page(self):
        self.back.setEnabled(False)
        self.next.setIcon(self.style().standardIcon(QStyle.SP_DialogApplyButton))
        self.next.triggered.connect(self.close)
        self.next.setText('Finish')
        self.next.setEnabled(True)


class CGBuilder(PagedWindow):
    def __init__(self):
        super().__init__()
        mol_picker = MoleculeWidget(self)
        self.pages.addWidget(mol_picker)

        mapper = MappingWidget(self)
        self.pages.addWidget(mapper)

        writer = WriterWidget(self)
        self.pages.addWidget(writer)

        self.pages.setCurrentIndex(0)
        self._toggle_buttons()
        mapper.molecule = read_smiles('OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O')

    def _next_page(self):
        value = self.pages.currentWidget().get_value()
        if value:
            super()._next_page()
            self.pages.currentWidget().set_value(value)


if __name__ == '__main__':
    import sys
    qapp = QApplication(sys.argv)
    app = CGBuilder()
    app.show()
    qapp.exec_()
