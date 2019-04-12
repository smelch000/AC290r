MODULE points2cells_mod

USE iso_c_binding
IMPLICIT NONE

INTEGER :: nn
INTEGER :: nx,ny,nz
INTEGER(kind=8) :: lnx,lny,lnxy

CONTAINS

    INTEGER FUNCTION binsearch(i4,itppp) RESULT(low)
        INTEGER(kind=8) :: i4
        INTEGER(kind=8) :: itppp(0:)
        INTEGER :: n

    INTEGER :: high, medium, lowb, highb

    lowb = LBOUND(itppp,1)
    highb = UBOUND(itppp,1) ! nn

    low = lowb
    high = highb + 1

    ! print *,'low,high:',low,high,nn

    DO WHILE(low.lt.high)

       medium = (high+low)/2

       if(itppp(medium) < i4) then
          low = medium+1
       else
          high = medium
       endif

    ENDDO

    IF (low < lowb .or. low > highb) then
       low = -99
    ELSE
       IF (itppp(low) /= i4) low = -99
    ENDIF

    ! print *,(i4==itppp(low)),low

    END FUNCTION binsearch

    !!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(kind=8) FUNCTION i4back(i,j,k) 
        INTEGER :: i,j,k
        i4back = k * lnxy + j * lnx + i
    END FUNCTION i4back

    !!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE points2cells(nn_l,nx_l,ny_l,nz_l, &
                            spacing_l, &
                            ijk1,ijk2,ijk3, &
                            ncell,icell) BIND(C,name='points2cells')

    INTEGER(c_int),VALUE,INTENT(in) :: nn_l,nx_l,ny_l,nz_l
    INTEGER(c_int),VALUE,INTENT(in) :: spacing_l
    INTEGER(c_int),INTENT(in) :: ijk1(0:nn_l-1), ijk2(0:nn_l-1), ijk3(0:nn_l-1)
    INTEGER(c_int),INTENT(out) :: ncell, icell(0:8*nn_l-1)

    INTEGER :: i,j,k,n
    INTEGER(KIND=8) :: i4
    INTEGER(KIND=8),ALLOCATABLE :: itppp(:)
    INTEGER :: ipjk,ijpk,ijkp,ipjpk,ipjkp,ijpkp,ipjpkp
    INTEGER :: mm

    nn = nn_l
    nx = nx_l
    ny = ny_l
    nz = nz_l
    lnx = nx_l
    lny = ny_l
    lnxy = nx_l * ny_l

    print *
    print *,'...running points2cells (from points2cells.so)'
    print *
    ! print *,'nn:',nn_l,'nx,ny,nz:',nx_l,ny_l,nz_l,lnx,lnxy

    ! ALLOCATE(itppp(nn))
    ALLOCATE(itppp(0:nn-1))

    DO n=0, nn-1
        i = ijk1(n); j = ijk2(n); k = ijk3(n)
        itppp(n) = i4back(i,j,k)
    ENDDO

    ! print *,'itppp bounds:',LBOUND(itppp,1),UBOUND(itppp,1)

    ncell = 0
    DO n=0, nn-1

        i = ijk1(n); j = ijk2(n); k = ijk3(n)

        ipjk   = binsearch(i4back(i+spacing_l, j,           k          ),itppp) 
        if (ipjk < 0) CYCLE

        ijpk   = binsearch(i4back(i,           j+spacing_l, k          ),itppp) 
        if (ijpk <0) CYCLE

        ijkp   = binsearch(i4back(i,           j,           k+spacing_l),itppp)
        if (ijkp <0) CYCLE

        ipjpk  = binsearch(i4back(i+spacing_l, j+spacing_l, k          ),itppp)
        if (ipjpk <0) CYCLE

        ipjkp  = binsearch(i4back(i+spacing_l,j,            k+spacing_l),itppp)
        if (ipjkp <0) CYCLE

        ijpkp  = binsearch(i4back(i  ,         j+spacing_l, k+spacing_l),itppp)
        if (ijpkp <0) CYCLE

        ipjpkp = binsearch(i4back(i+spacing_l, j+spacing_l, k+spacing_l),itppp)
        if (ipjpkp <0) CYCLE

        ! print *, itppp(ipjpkp), i4back(i+spacing_l,j+spacing_l,k+spacing_l)

        icell(8*ncell    ) = n
        icell(8*ncell + 1) = ipjk
        icell(8*ncell + 2) = ijpk
        icell(8*ncell + 3) = ipjpk
        icell(8*ncell + 4) = ijkp
        icell(8*ncell + 5) = ipjkp
        icell(8*ncell + 6) = ijpkp
        icell(8*ncell + 7) = ipjpkp

        ncell = ncell + 1

    ENDDO
    print *,'ncell:',ncell

    DEALLOCATE(itppp)

    END SUBROUTINE points2cells

END MODULE points2cells_mod

