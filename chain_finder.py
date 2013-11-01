# -*- coding: utf-8 -*-

__author__ = 'mactep'

import os
import sys
import logging
from PyQt4.QtGui import *
from PyQt4.Qt import pyqtSlot

from vh_finder.vh_find import process


class WidgetSanger(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.lay = QVBoxLayout()

        self.lay1 = QHBoxLayout()
        self.labin = QLabel("In: ")
        self.linein = QLineEdit()
        self.buttin = QPushButton("File...")
        self.lay1.addWidget(self.labin)
        self.lay1.addWidget(self.linein)
        self.lay1.addWidget(self.buttin)

        self.lay2 = QHBoxLayout()
        self.labout = QLabel("Out: ")
        self.lineout = QLineEdit()
        self.buttout = QPushButton("File...")
        self.lay2.addWidget(self.labout)
        self.lay2.addWidget(self.lineout)
        self.lay2.addWidget(self.buttout)

        self.lay3 = QVBoxLayout()
        self.lay3.addWidget(QLabel("Min protein length:"))
        self.spinmin = QSpinBox(self)
        self.spinmin.setMinimum(0)
        self.spinmin.setMaximum(500)
        self.spinmin.setValue(100)
        self.lay3.addWidget(self.spinmin)

        self.buttok = QPushButton("Ok")

        self.label = QLabel(" ")

        self.lay.addLayout(self.lay1)
        self.lay.addLayout(self.lay2)
        self.lay.addLayout(self.lay3)
        self.lay.addWidget(self.buttok)
        self.lay.addWidget(self.label)

        self.setLayout(self.lay)

        self.buttin.clicked.connect(self.getFilename)
        self.buttout.clicked.connect(self.getFilename)

        self.buttok.clicked.connect(self.run)

    @pyqtSlot(name="getFilename")
    def getFilename(self):
        directory = QFileDialog.getExistingDirectory(self, "Open Directory", "",
                                                     QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
        if self.sender() == self.buttin:
            self.linein.setText(directory)
        else:
            self.lineout.setText(directory)
        self.label.setText(" ")

    @pyqtSlot(name="run")
    def run(self):
        process(str(self.linein.text()), str(self.lineout.text()), int(self.spinmin.value()))
        self.label.setText("READY!")


def main():
    logging.basicConfig(format=u'%(filename)s [LINE:%(lineno)d]# %(levelname)-8s [%(asctime)s]  %(message)s',
                        level=logging.DEBUG, filename="finder.log")
    app = QApplication(sys.argv)

    ws = WidgetSanger()
    ws.setVisible(True)
    logging.debug("Window opened!")
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()