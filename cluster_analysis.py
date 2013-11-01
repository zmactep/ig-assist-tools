__author__ = 'yakovlev'

import os
import sys
import logging
from PyQt4.QtGui import *
from PyQt4.Qt import pyqtSlot

from ighumanizer3.extra.share.fasta_tools import read_fasta

from cluster_processing.cluster_info import process, process_one
from cluster_processing.heavy_light_processing import hl_process

class WidgetCluster(QWidget):
    def __init__(self):
        QWidget.__init__(self)

        # Switches
        radio_group = QButtonGroup(self)
        self.radio_sanger = QRadioButton("Sanger_run results", self)
        self.radio_standalone = QRadioButton("Standalone run", self)

        radio_group.addButton(self.radio_sanger)
        radio_group.addButton(self.radio_standalone)

        self.radio_sanger.setChecked(True)

        # Group type
        type_group = QGroupBox("Type", self)
        type_lay = QVBoxLayout()
        type_lay.addWidget(self.radio_sanger)
        type_lay.addWidget(self.radio_standalone)
        type_group.setLayout(type_lay)

        # Group options
        options_group = QGroupBox("Options", self)
        opt_lay = QVBoxLayout()
        opt_llay = QHBoxLayout()
        self.labin = QLabel("SangerRun results: ")
        self.linein = QLineEdit(self)
        self.buttin = QPushButton("File...")
        opt_llay.addWidget(self.labin)
        opt_llay.addWidget(self.linein)
        opt_llay.addWidget(self.buttin)
        self.radio_confidence = QCheckBox("Merge similar", self)
        self.radio_confidence.setChecked(False)
        opt_lay.addLayout(opt_llay)
        opt_lay.addWidget(self.radio_confidence)
        options_group.setLayout(opt_lay)

        # Group additional
        add_lay = QHBoxLayout()
        self.add_group = QGroupBox("Addtional")
        labadd = QLabel("Similarity: ", self)
        self.spinadd = QSpinBox(self)
        self.spinadd.setMinimum(0)
        self.spinadd.setMaximum(100)
        self.spinadd.setSingleStep(5)
        self.spinadd.setValue(100)
        add_lay.addWidget(labadd)
        add_lay.addWidget(self.spinadd)
        self.add_group.setEnabled(False)
        self.add_group.setLayout(add_lay)

        # Submit
        self.bok = QPushButton("Ok")
        self.lok = QLabel(" ")

        # Construct
        mlay = QVBoxLayout()
        mlay.addWidget(type_group)
        mlay.addWidget(options_group)
        mlay.addWidget(self.add_group)
        mlay.addWidget(self.bok)
        mlay.addWidget(self.lok)
        self.setLayout(mlay)

        # Connections
        radio_group.buttonClicked.connect(self.radio_changed)
        self.buttin.clicked.connect(self.get_filename)
        self.bok.clicked.connect(self.run)
        self.radio_confidence.clicked.connect(self.additional_activate)

    @pyqtSlot(name="radio_changed")
    def radio_changed(self):
        if self.radio_sanger.isChecked():
            self.labin.setText("SangerRun results: ")
        else:
            self.labin.setText("File to analyse: ")

    @pyqtSlot(name="get_filename")
    def get_filename(self):
        if self.radio_sanger.isChecked():
            d = QFileDialog.getExistingDirectory(self, "Open Directory", "",
                                                 QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
        else:
            d = QFileDialog.getOpenFileName(self, "Open File", os.curdir,
                                            "FASTA files (*.fa *.fasta)")
        if self.sender() == self.buttin:
            self.linein.setText(d)
        else:
            self.lineout.setText(d)
        self.lok.setText(" ")

    @pyqtSlot(name="additional_activate")
    def additional_activate(self):
        if self.radio_confidence.isChecked():
            self.add_group.setEnabled(True)
        else:
            self.add_group.setEnabled(False)

    @pyqtSlot(name="run")
    def run(self):
        prct = False if not self.radio_confidence.isChecked() else float(self.spinadd.value()) / 100
        out = str(self.linein.text())
        if self.radio_sanger.isChecked():
            analysis_dir = "analysis"
            process(out, analysis_dir, prct)
            hl_process(os.path.join(out, analysis_dir), prct)
        else:
            analysis_dir = os.path.dirname(out)
            chain = os.path.basename(out)
            chain = chain[:chain.rfind(".")] + "-analysis"
            process_one(analysis_dir, chain, read_fasta(out, False), prct)
        self.lok.setText("READY")


def main():
    logging.basicConfig(format=u'%(filename)s [LINE:%(lineno)d]# %(levelname)-8s [%(asctime)s]  %(message)s',
                        level=logging.DEBUG, filename="cluster.log")
    app = QApplication(sys.argv)

    ws = WidgetCluster()
    ws.setVisible(True)
    logging.debug("Window opened!")
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()