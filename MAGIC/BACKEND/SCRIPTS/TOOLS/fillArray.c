#include <stdio.h>
typedef unsigned short int myShort;

extern void removeFromStencil(int kmin,int kmax,
		       int jmin,int jmax,
                       int imin,int imax,
		       void *pstencilExtent, int stencilIDim, int stencilJDim, void *pstencilArr,
                       void *pimageExtent,   int imageIDim,   int imageJDim,   void *pimageArr)
{
	// printf("kmin %i, kmax %i, jmin %i, jmax %i, imin %i, imax %i, stencilIDim %i, stencilJDim %i, imageIDim %i,   imageJDim %i \n", 
	//        kmin   , kmax   , jmin   , jmax   , imin   , imax   , stencilIDim   ,  stencilJDim  , imageIDim   ,   imageJDim ); 

       int *stencilExtent = pstencilExtent;
       myShort *stencilArr    = pstencilArr;
       int *imageExtent   = pimageExtent;
       myShort *imageArr      = pimageArr;
        
        for(int k=kmin; k < kmax; k++)
        { 
        	for(int j=jmin; j < jmax; j++)
        	{
        		for(int i=imin; i < imax; i++)
        		{
                                int stindex = (k-stencilExtent[4])*stencilIDim*stencilJDim+
                                          (j-stencilExtent[2])*stencilIDim+
                                          (i-stencilExtent[0]);
                                int stvalue = stencilArr[stindex];
                                if( stvalue == 0)
        			{
                                        int wiindex = (k-imageExtent[4])*imageIDim*imageJDim+
                                                  (j-imageExtent[2])*imageIDim+
                                                  (i-imageExtent[0]);
                                        imageArr[wiindex] = 0;
        			}
        		}
        	}
        }
}
float sqrdist(float x1,float y1,float z1,float x2,float y2,float z2)
{
	return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);	
}
extern void sphereRemove( void *pimageExtent,  
                            void *pimageSpac,  
                            void *pimageOrig,  
			    int imageIDim,   int imageJDim,   void *pimageArr,
			    float x, float y, float z, float r)
{


	 int  *imageExtent   = pimageExtent; 
	 float *imgSpacing   = pimageSpac; 
	 float *imageOrig    = pimageOrig;  
         myShort *imageArr   = pimageArr;
	
	float rsqr = r*r;
	for( int k = imageExtent[4]; k <  imageExtent[5]; k++)
	{
		for( int j = imageExtent[2]; j <  imageExtent[3]; j++)
		{
			for( int i = imageExtent[0]; i < imageExtent[1]; i++)
			{
				int wiindex = (k-imageExtent[4])*imageIDim*imageJDim+
						  (j-imageExtent[2])*imageIDim+
						  (i-imageExtent[0]);
				
				float curx = imgSpacing[0]*i; 
				float cury = imgSpacing[1]*j;
				float curz = imgSpacing[2]*k;
				if( sqrdist(curx,cury,curz,x,y,z) < rsqr )
					imageArr[wiindex] = 0;
			}
		}
	}

}
