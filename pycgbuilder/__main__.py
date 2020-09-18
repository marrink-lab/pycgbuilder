from PyQt5.QtWidgets import QApplication
from .interface import CGBuilder
import sys

qapp = QApplication(sys.argv)
app = CGBuilder()
app.show()
qapp.exec_()
