#/usr/bin/env python2.7

try:
    import vtk
    from vtk import *
    from vtk.util import numpy_support 

except:
    print ('No VTK support in module: ' + __name__)

class ErrorObserver:

    """
    Usage on client side:

       err = ErrorObserver()
       filt = vtkMassProperties()
       filt = ...
       filt.Update()
       if err.ErrorOccurred():
          ...catch exception
       else:
          move on
    """

    def __init__(self):
        self.__ErrorOccurred = False
        self.__ErrorMessage = None
        self.CallDataType = 'string0'

    def __call__(self, obj, event, message):
        self.__ErrorOccurred = True
        self.__ErrorMessage = message

    def ErrorOccurred(self):
        occ = self.__ErrorOccurred
        self.__ErrorOccurred = False
        return occ

