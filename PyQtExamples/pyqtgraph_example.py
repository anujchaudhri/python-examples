#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 11:03:04 2018

@author: anujchaudhri
"""
from PyQt5 import QtCore, QtGui
import numpy as np
import pyqtgraph as pg

def rand(n):
    data = np.random.random(n)
    data[int(n*0.1):int(n*0.13)] += .5
    data[int(n*0.18)] += 2
    data[int(n*0.1):int(n*0.13)] *= 5
    data[int(n*0.18)] *= 20
    data *= 1e-12
    return data, np.arange(n, n+len(data)) / float(n)

def updateData(self,p1):
        yd, xd = rand(10000)
        p1.setData(y=yd, x=xd)

class MainWindow(QtGui.QMainWindow):
    def __init__(self,parent=None):
#        Can use either def below to call the constructor of the parent class
#        since MainWindow is derived from QtGui.QMainWindow()
#        QtGui.QMainWindow.__init__(self,parent=parent)
        super(MainWindow,self).__init__(parent)
        self.setWindowTitle('Plots from Raspberry Pi')
        self.resize(1000,800)
        __widget = QtGui.QWidget()
        self.setCentralWidget(__widget)
        layout = QtGui.QVBoxLayout()
        __widget.setLayout(layout)
        
        layout.addWidget(Widget1())
        layout.addWidget(Widget2())
        layout.addWidget(Widget3())
        
class Widget1(pg.PlotWidget):
    def __init__(self,parent=None):
        pg.PlotWidget.__init__(self,parent=parent)
        p1 = self.plot()
        p1.setPen((200,200,100))

        rect = QtGui.QGraphicsRectItem(QtCore.QRectF(0,0,1,5e-11))
        rect.setPen(pg.mkPen(100,200,100))
        self.addItem(rect)
        
        self.setLabel('left','Value',unit='V')
        self.setLabel('bottom','Time',units='s')
        self.setXRange(0,2)
        self.setYRange(0,1e-10)
        
        t = QtCore.QTimer()
        t.timeout.connect(lambda: updateData(p1))
        t.start(50)
        

class Widget2(pg.PlotWidget):
    def __init__(self,parent=None):
        #pg.PlotWidget().__init__(self,parent)
        for i in range(0, 5):
            for j in range(0, 3):
                yd, xd = rand(10000)
                self.plot(y=yd*(j+1), x=xd, params={'iter': i, 'val': j})

class Widget3(pg.PlotWidget):
    def __init__(self,parent=None):
        #pg.PlotWidget().__init__(self,parent)
        curve = self.plot(np.random.normal(size=100)*1e0, clickable=True)
        curve.curve.setClickable(True)
        curve.setPen('w')  ## white pen
        curve.setShadowPen(pg.mkPen((70,70,30), width=6, cosmetic=True))
        curve.sigClicked.connect(self.clicked)

        lr = pg.LinearRegionItem([1, 30], bounds=[0,100], movable=True)
        self.addItem(lr)
        line = pg.InfiniteLine(angle=90, movable=True)
        self.addItem(line)
        line.setBounds([0,200])

    def clicked():
        print("curve clicked")


if __name__ == '__main__':
    import sys
    app = QtGui.QApplication([])
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())
