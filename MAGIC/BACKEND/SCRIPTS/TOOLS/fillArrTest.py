import ctypes
kmin        = ctypes.c_int(111) 
kmax        = ctypes.c_int(222)
jmin        = ctypes.c_int(333) 
jmax        = ctypes.c_int(444)
imin        = ctypes.c_int(555) 
imax        = ctypes.c_int(666)
stencilIDim = ctypes.c_int(777) 
stencilJDim = ctypes.c_int(888)
imageIDim   = ctypes.c_int(999)  
imageJDim   = ctypes.c_int(101)  


#stenExtent = ctypes.c_int*6
#imgExtent  = ctypes.c_int*6
class Extent(ctypes.Structure):
	_fields_ = [ ("xlo", ctypes.c_int),
	             ("xhi", ctypes.c_int),
	             ("ylo", ctypes.c_int),
	             ("yhi", ctypes.c_int),
	             ("zlo", ctypes.c_int),
	             ("zhi", ctypes.c_int)]

stenExtent = Extent(1,2,3,4,5,6)
imgExtent = Extent(11,22,33,44,55,66)

print stenExtent

#extern void fillInFromStencil(int kmin,int kmax,
#		               int jmin,int jmax,
#                              int imin,int imax,
#		               void *pstencilExtent, int stencilIDim, int stencilJDim, void *pstencilArr,
#                              void *pimageExtent,   int imageIDim,   int imageJDim,   void *pimageArr)
#{

myshortArr = ctypes.c_uint*100;
#
myLib = ctypes.CDLL("./fillArr.so")
myLib.fillInFromStencil(kmin,kmax,
			jmin,jmax,
			imin,imax,
			ctypes.POINTER(Extent)(stenExtent), stencilIDim, stencilJDim, ctypes.POINTER(Extent)(stenExtent) ,
			ctypes.POINTER(Extent)(imgExtent),  imageIDim,   imageJDim,   ctypes.POINTER(Extent)(imgExtent) )  

