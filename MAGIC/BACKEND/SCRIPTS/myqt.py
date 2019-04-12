#/usr/bin/env python2.7

import sys

try:

    QTVERSION = 5
    from PyQt5 import QtGui, QtCore, QtWidgets, uic
    from PyQt5.Qt import QPainter
    from PyQt5.QtGui import QPixmap, QFont, QIcon, QBrush, QStandardItem, QStandardItemModel, QColor, QLinearGradient, QPen
    from PyQt5.QtCore import QSettings, QObject, QLineF, QPointF, QRect, QFileInfo, qDebug, pyqtSlot, pyqtSignal
    from PyQt5.QtWidgets import QApplication, \
                                QSplashScreen, \
                                QProgressBar, \
                                QGraphicsItem, \
                                QGraphicsLineItem, \
                                QMainWindow, \
                                QWidget, \
                                QDialog, \
                                QTreeWidget, \
                                QTreeWidgetItem, \
                                QTreeWidgetItemIterator, \
                                QStackedLayout, \
                                QDesktopWidget, \
                                QLabel, \
                                QTextEdit, \
                                QPushButton, \
                                QRadioButton, \
                                QProgressDialog, \
                                QVBoxLayout, \
                                QHBoxLayout, \
                                QSpinBox, \
                                QFrame, \
                                QCheckBox, \
                                QLineEdit, \
                                QComboBox, \
                                QDoubleSpinBox, \
                                QGridLayout, \
                                QScrollArea, \
                                QTabWidget, \
                                QSplitter, \
                                QMenu, \
                                QFileDialog, \
                                QMessageBox, \
                                QDockWidget, \
                                QTreeView, \
                                QTableView, \
                                QSlider, \
                                QToolButton, \
                                QGraphicsView, \
                                QGraphicsScene, \
                                QActionGroup, \
                                QAction, \
                                QToolBar, \
                                QFileIconProvider, \
                                QColorDialog, \
                                qApp
    
    from QVTKRenderWindowInteractor import *
    # import QVTKRenderWindowInteractor

except ImportError:

    #QTVERSION = 4
    #from PyQt4 import QtGui, QtCore, uic
    #from PyQt4.QtCore import *   # PyQt core libs
    #from PyQt4.QtGui import *    # PyQy GUI components
    
    print ('No PyQT5 support in module: ',__name__)
    sys.exit(1)
