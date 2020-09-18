from PyQt5.QtWidgets import QApplication
from .interface import CGBuilder
import sys


def main():
    qapp = QApplication(sys.argv)
    app = CGBuilder()
    app.show()
    qapp.exec_()

if __name__ == '__main__':
    main()
