# -*- coding: utf-8 -*-

__author__ = 'mactep'

import os
import sys
import logging
from PyQt4.QtGui import *
from PyQt4.QtCore import *

from sanger_processing.sanger_processing import run_on_directory
from sanger_processing.cluster_analysis import *


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

        self.lay3 = QHBoxLayout()
        self.lay31 = QVBoxLayout()
        self.lay32 = QVBoxLayout()
        self.lay31.addWidget(QLabel("Forward:"))
        self.lay32.addWidget(QLabel("Backward:"))
        self.forward = QLineEdit("SeqR")
        self.backward = QLineEdit("H3b")
        self.lay31.addWidget(self.forward)
        self.lay32.addWidget(self.backward)
        self.lay3.addLayout(self.lay31)
        self.lay3.addLayout(self.lay32)

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
        sanger_run(str(self.linein.text()), str(self.lineout.text()),
                   str(self.forward.text()), str(self.backward.text()))
        self.label.setText("READY!")


def sanger_step1(in_dir, out_dir, fm, bm):
    run_on_directory(in_dir, out_dir, fm, bm)


def sanger_step2(out_dir):
    vl, vh, vln, vhn = load_directory(out_dir)
    save_as_alignment(vl, os.path.join(out_dir, "light-chains-alignment.txt"))
    save_as_alignment(vh, os.path.join(out_dir, "heavy-chains-alignment.txt"))


def sanger_run(in_dir, out_dir, fm, bm):
    logging.debug("Step 1: Data preprocessing")
    sanger_step1(in_dir, out_dir, fm, bm)
    logging.debug("Step 2: Data alignment")
    sanger_step2(out_dir)


def main():
    logging.basicConfig(format=u'%(filename)s [LINE:%(lineno)d]# %(levelname)-8s [%(asctime)s]  %(message)s',
                        level=logging.DEBUG, filename="sanger.log")
    app = QApplication(sys.argv)

    ws = WidgetSanger()
    ws.setVisible(True)
    logging.debug("Window opened!")
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()