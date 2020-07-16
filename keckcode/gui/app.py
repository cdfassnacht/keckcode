from select_data import dataSelectWindow, spectraFileBrowser
from PyQt5.QtWidgets import QApplication
from keckcode.deimos import deimosmask1d
import matplotlib.pyplot as plt

class spectraApp(QApplication):
    def __init__(self):
        super().__init__([])
        self.startupWindow = dataSelectWindow()
        self.connectSignals()
    
        self.startupWindow.show()

    def connectSignals(self):
        self.startupWindow.fileOpened.connect(lambda x: self.handleFileOpen(x))
    
    def handleFileOpen(self, fname):
        plt.ion()
        mask1 = deimosmask1d.DeimosMask1d(fname)
        mask1.plot(0)
