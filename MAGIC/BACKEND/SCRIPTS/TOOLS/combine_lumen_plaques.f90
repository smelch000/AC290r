MODULE combine_lumen_plaques_mod

USE iso_c_binding
IMPLICIT NONE

INTEGER(kind=8) :: nx,ny,nz,nx2,nxy2

INTERFACE quick_sort
   MODULE PROCEDURE quick_sort_i8
   MODULE PROCEDURE quick_sort_i4
END INTERFACE

CONTAINS

   INTEGER FUNCTION binsearch(array, lastCindex, value) RESULT(low)

  ! array : sorted array (increasing order)
  ! lastCindex : last index of the input array in C-style (e.g. lastCindex = SIZE(array)-1)
  ! If the input value is larger than array(N-1), then binsearch returns N (aka insert at the rhs of array)

    INTEGER(kind=8),INTENT(in) :: array(0:)
    INTEGER(kind=8),INTENT(in) :: value
    INTEGER,INTENT(in) :: lastCindex
    INTEGER :: high, medium

    low = 0
    high = lastCindex+1

    do while(low < high)
       medium=(high+low)/2
       if(array(medium).lt.value) then
          low=medium+1
       else
          high=medium
       endif
    end do

    ! print *,LBOUND(array,1),UBOUND(array,1),lastCindex; STOP

    IF(low>UBOUND(array,1)) THEN
        low=-1
    ELSE IF(value /= array(low)) THEN
        low=-1
    ENDIF

  END FUNCTION binsearch

    !!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(kind=8) FUNCTION i4back(i,j,k) 
        INTEGER :: i,j,k
        i4back = k * nxy2 + j * nx2 + i
    END FUNCTION i4back

    !!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION ijk(idir, i4) BIND(C,name='ijk')
    INTEGER(kind=8),INTENT(in) :: i4 !< location on mesh
    INTEGER(kind=8) :: iz, iy, ix
    INTEGER :: idir, ijk

      iz = i4/INT(nxy2,KIND=8)
      IF(idir==3) THEN
         ijk = INT(iz,KIND=4)
         RETURN
      ENDIF

      iy = (i4 - iz*INT(nxy2,KIND=8))/INT(nx2,KIND=8)
      IF(idir==2) THEN
         ijk = INT(iy,KIND=4)
         RETURN
      ENDIF

      ijk = INT(i4 - iz*INT(nxy2,KIND=8) - iy*INT(nx2,KIND=8),KIND=4)

  END FUNCTION ijk


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !> Quick sort routine from:
        !! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
        !! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
        !! Modified by Alan Miller to include an associated integer array which gives
        !! the positions of the elements in the original order.
        RECURSIVE SUBROUTINE quick_sort_i8(list,order,ln)


          IMPLICIT NONE
          INTEGER :: ln
          INTEGER(kind=8), DIMENSION (1:ln), INTENT(inout)  :: list
          INTEGER, DIMENSION(1:ln), INTENT(out)  :: order

          ! Local variable
          INTEGER(kind=8) :: i

          DO i = 1, SIZE(list)
             order(i) = i
          END DO

          ! CALL quick_sort_1(1, SIZE(list))
          CALL quick_sort_1_i8(INT(1,kind=8), INT(SIZE(list),kind=8))

        CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          RECURSIVE SUBROUTINE quick_sort_1_i8(left_end, right_end)

            INTEGER(kind=8),INTENT(in) :: left_end, right_end

            !     Local variables
            INTEGER(kind=8)             :: i, j, itemp
            INTEGER(kind=8)                :: reference, temp
            INTEGER(kind=8), PARAMETER  :: max_simple_sort_size = 6

            IF (right_end < left_end + max_simple_sort_size) THEN
               ! Use interchange sort for small lists
               CALL interchange_sort_i8(left_end, right_end)
            ELSE
               ! Use partition ("quick") sort
               reference = list((left_end + right_end)/2)
               i = left_end - 1; j = right_end + 1

               DO
                  ! Scan list from left end until element >= reference is found
                  DO
                     i = i + 1
                     IF (list(i) >= reference) EXIT
                  END DO
                  ! Scan list from right end until element <= reference is found
                  DO
                     j = j - 1
                     IF (list(j) <= reference) EXIT
                  END DO


                  IF (i < j) THEN
                     ! Swap two out-of-order elements
                     temp = list(i); list(i) = list(j); list(j) = temp
                     itemp = order(i); order(i) = order(j); order(j) = itemp
                  ELSE IF (i == j) THEN
                     i = i + 1
                     EXIT
                  ELSE
                     EXIT
                  END IF
               END DO

               IF (left_end < j) CALL quick_sort_1_i8(left_end, j)
               IF (i < right_end) CALL quick_sort_1_i8(i, right_end)
            END IF

          END SUBROUTINE quick_sort_1_i8


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SUBROUTINE interchange_sort_i8(left_end, right_end)

            INTEGER(kind=8), INTENT(in) :: left_end, right_end

            !     Local variables
            INTEGER(kind=8)             :: i, j, itemp
            INTEGER(kind=8)                :: temp

            DO i = left_end, right_end - 1
               DO j = i+1, right_end
                  IF (list(i) > list(j)) THEN
                     temp = list(i); list(i) = list(j); list(j) = temp
                     itemp = order(i); order(i) = order(j); order(j) = itemp
                  END IF
               END DO
            END DO

          END SUBROUTINE interchange_sort_i8

        END SUBROUTINE quick_sort_i8


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !> Quick sort routine from:
        !! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
        !! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
        !! Modified by Alan Miller to include an associated integer array which gives
        !! the positions of the elements in the original order.
        RECURSIVE SUBROUTINE quick_sort_i4(list,order,ln)


          IMPLICIT NONE
          INTEGER :: ln
          INTEGER(kind=4), DIMENSION (1:ln), INTENT(inout)  :: list
          INTEGER, DIMENSION(1:ln), INTENT(out)  :: order

          ! Local variable
          INTEGER(kind=4) :: i

          DO i = 1, SIZE(list)
             order(i) = i
          END DO

          ! CALL quick_sort_1(1, SIZE(list))
          CALL quick_sort_1_i4(INT(1,kind=4), INT(SIZE(list),kind=4))

        CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          RECURSIVE SUBROUTINE quick_sort_1_i4(left_end, right_end)

            INTEGER(kind=4),INTENT(in) :: left_end, right_end

            !     Local variables
            INTEGER(kind=4)             :: i, j, itemp
            INTEGER(kind=4)                :: reference, temp
            INTEGER(kind=4), PARAMETER  :: max_simple_sort_size = 6

            IF (right_end < left_end + max_simple_sort_size) THEN
               ! Use interchange sort for small lists
               CALL interchange_sort_i4(left_end, right_end)
            ELSE
               ! Use partition ("quick") sort
               reference = list((left_end + right_end)/2)
               i = left_end - 1; j = right_end + 1

               DO
                  ! Scan list from left end until element >= reference is found
                  DO
                     i = i + 1
                     IF (list(i) >= reference) EXIT
                  END DO
                  ! Scan list from right end until element <= reference is found
                  DO
                     j = j - 1
                     IF (list(j) <= reference) EXIT
                  END DO


                  IF (i < j) THEN
                     ! Swap two out-of-order elements
                     temp = list(i); list(i) = list(j); list(j) = temp
                     itemp = order(i); order(i) = order(j); order(j) = itemp
                  ELSE IF (i == j) THEN
                     i = i + 1
                     EXIT
                  ELSE
                     EXIT
                  END IF
               END DO

               IF (left_end < j) CALL quick_sort_1_i4(left_end, j)
               IF (i < right_end) CALL quick_sort_1_i4(i, right_end)
            END IF

          END SUBROUTINE quick_sort_1_i4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SUBROUTINE interchange_sort_i4(left_end, right_end)

            INTEGER(kind=4), INTENT(in) :: left_end, right_end

            !     Local variables
            INTEGER(kind=4)             :: i, j, itemp
            INTEGER(kind=4)                :: temp

            DO i = left_end, right_end - 1
               DO j = i+1, right_end
                  IF (list(i) > list(j)) THEN
                     temp = list(i); list(i) = list(j); list(j) = temp
                     itemp = order(i); order(i) = order(j); order(j) = itemp
                  END IF
               END DO
            END DO

          END SUBROUTINE interchange_sort_i4

        END SUBROUTINE quick_sort_i4


    !!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE combine(lname_lumen, name_lumen, &
                       lname_plaques, name_plaques, &
                       lname_output, name_output) BIND(C,name='combine')

    INTEGER(c_int),VALUE,INTENT(in) :: lname_lumen, lname_plaques, lname_output
    CHARACTER(LEN=1,KIND=c_char),INTENT(in) :: name_lumen(*), name_plaques(*), name_output(*)
    CHARACTER(len=256) :: flumen,fplaques,foutput

    INTEGER,PARAMETER :: npop=19
    INTEGER :: icx(0:npop-1),icy(0:npop-1),icz(0:npop-1) 

    INTEGER :: i,j,k,flg,nfl,nwl,nin,nou,nv,nv_new,n,ii,jj,kk,ip,nfl_new,nwl_new,nin_new,nou_new,nfinal,nv_pl,nn,nwl_clean
    INTEGER(KIND=8) :: i4,i4p
    INTEGER,ALLOCATABLE :: iflag_new(:),iflag_wall_clean(:),iflag_final(:),iflag_fina2(:)
    INTEGER(KIND=8),ALLOCATABLE :: itppp_new(:),itppp_pl(:),itppp_wall(:),itppp_final(:)
    INTEGER(KIND=8),ALLOCATABLE :: itppp_wall_clean(:)
    INTEGER,ALLOCATABLE :: iord(:)
    CHARACTER(len=256) :: str
    LOGICAL :: lexist

    icx(0) = 0; icy(0) = 0; icz(0) = 0
    icx(1) = 0; icy(1) = 0; icz(1) = 1
    icx(2) = 1; icy(2) = 0; icz(2) = 1
    icx(3) =-1; icy(3) = 0; icz(3) = 1
    icx(4) = 0; icy(4) = 1; icz(4) = 1
    icx(5) = 0; icy(5) =-1; icz(5) = 1
    icx(6) = 0; icy(6) = 0; icz(6) =-1
    icx(7) =-1; icy(7) = 0; icz(7) =-1
    icx(8) = 1; icy(8) = 0; icz(8) =-1
    icx(9) = 0; icy(9) =-1; icz(9) =-1
    icx(10)= 0; icy(10)= 1; icz(10)=-1
    icx(11)= 1; icy(11)= 0; icz(11)= 0
    icx(12)=-1; icy(12)= 0; icz(12)= 0
    icx(13)= 0; icy(13)= 1; icz(13)= 0
    icx(14)= 0; icy(14)=-1; icz(14)= 0
    icx(15)= 1; icy(15)= 1; icz(15)= 0
    icx(16)=-1; icy(16)= 1; icz(16)= 0
    icx(17)=-1; icy(17)=-1; icz(17)= 0
    icx(18)= 1; icy(18)=-1; icz(18)= 0

        print *
        print *,'...running (from combine_lumen_plaques.so)'
        print *
        
        flumen(:) = ''
        DO i=1, lname_lumen
          flumen(i:i) = name_lumen(i)
        ENDDO
        
        fplaques(:) = ''
        DO i=1, lname_plaques
          fplaques(i:i) = name_plaques(i)
        ENDDO
        
        foutput(:) = ''
        DO i=1, lname_output
          foutput(i:i) = name_output(i)
        ENDDO
        
        OPEN(10,FILE=TRIM(fplaques)//'.hdr')
        READ(10,*) 
        READ(10,*) i,nfl,nwl,nin,nou
        CLOSE(10)

        print *,'plaques: nfl,nwl,nin,nou:', nfl,nwl,nin,nou
        IF(nin>0 .OR. nou>0) STOP 'something wrong1'
        ALLOCATE(itppp_pl(0:nfl+nwl-1))
        
        OPEN(10,FILE=TRIM(flumen)//'.hdr')
        READ(10,*) nx,ny,nz
        READ(10,*) i,nfl,nwl,nin,nou
        CLOSE(10)

        ! nv = nfl + nwl + nin + nou
        nv = nfl + nin + nou ! NO WALLS

        nx2 = (nx+2)
        nxy2 = (nx+2) * (ny+2)

        print *,'  lumen: nfl,nwl,nin,nou:', nfl,nwl,nin,nou
        ! print *,'         nv:', nv
        
        ALLOCATE(itppp_new(0:nv-1), iflag_new(0:nv-1))
        ALLOCATE(itppp_wall(0:4*nv-1))
        
        ! read from plaques (box from lumen should be read in at this point)
        OPEN(10,FILE=TRIM(fplaques)//'.dat')

        nv_pl = 0
        DO
            READ(10,*,ERR=1,END=1) i,j,k,flg

            itppp_pl(nv_pl) = i4back(i,j,k)
            nv_pl = nv_pl + 1

        ENDDO
        1 CLOSE(10)
        ! nv_pl = nv_pl - 1

        ! read from lumen and strip off wall nodes
        OPEN(10,FILE=TRIM(flumen)//'.dat')
        nv_new = 0; nfl_new = 0; nwl_new = 0; nin_new = 0; nou_new = 0
        DO
            READ(10,*,ERR=2,END=2) i,j,k,flg
        
            IF(flg==2) CYCLE ! exclude walls
        
            i4 = i4back(i,j,k)

            IF ( binsearch(itppp_pl, nv_pl-1, i4) >= 0 ) THEN
                CYCLE ! exclude points in plaques (fluid or wall)
            ENDIF
        
            IF(flg==1) THEN

               nfl_new = nfl_new + 1

            ELSE IF(flg==3) THEN

               nin_new = nin_new + 1

            ELSE IF(flg==4) THEN

               nou_new = nou_new + 1

            ENDIF

            itppp_new(nv_new) = i4
            iflag_new(nv_new) = flg
        
            nv_new = nv_new + 1
        
        ENDDO
        2 CLOSE(10)
        
        ! surround with wall nodes
        nwl_new = 0
        DO n=0, nv_new-1

            i4 = itppp_new(n); i = ijk(1,i4); j = ijk(2,i4); k = ijk(3,i4)

            DO ip=1, npop-1

                ii = i + icx(ip); jj = j + icy(ip); kk = k + icz(ip)
                i4p = i4back(ii,jj,kk)

                IF ( binsearch(itppp_new, nv_new-1, i4p) < 0 ) THEN
                    itppp_wall(nwl_new) = i4p
                    nwl_new = nwl_new + 1
                ENDIF

            ENDDO
        ENDDO
        ! print *,'  New wall nodes:',nwl_new

        ! sort wall and remove doubles
        ALLOCATE(iord(0:nwl_new-1))

        CALL quick_sort(itppp_wall(0:), iord(0:), nwl_new)

        DEALLOCATE(iord)

        nwl_clean = 0
        ALLOCATE(itppp_wall_clean(0:nwl_new-1))
        ALLOCATE(iflag_wall_clean(0:nwl_new-1))
        DO i=0, nwl_new-1

            IF(i>0) THEN
               IF(itppp_wall(i-1) == itppp_wall(i)) CYCLE

               IF(itppp_wall(i-1) > itppp_wall(i)) THEN
                   STOP 'ERROR...'
               ENDIF
            ENDIF

            itppp_wall_clean(nwl_clean) = itppp_wall(i)
            iflag_wall_clean(nwl_clean) = 2
            nwl_clean = nwl_clean + 1
        ENDDO

        ! print *,'Clean wall nodes:',nwl_clean

        ! merge current and walls
        nfinal = nv_new + nwl_clean

        ALLOCATE(itppp_final(0:nfinal-1))
        ALLOCATE(iflag_final(0:nfinal-1))
        ALLOCATE(iflag_fina2(0:nfinal-1))

        itppp_final = [itppp_new(0:nv_new-1), itppp_wall_clean(0:nwl_clean-1)]
        iflag_final = [iflag_new(0:nv_new-1), iflag_wall_clean(0:nwl_clean-1)]

        ! print *,'sample:',itppp_new(0:1),itppp_new(nv_new-2:nv_new-1)

        DEALLOCATE(itppp_new, iflag_new)
        DEALLOCATE(itppp_wall_clean, iflag_wall_clean)

        ! sort global
        ALLOCATE(iord(0:nfinal-1))

        CALL quick_sort(itppp_final(0:), iord(0:), nfinal)

        DO i=0, nfinal-1
            iflag_fina2( i ) = iflag_final( iord(i)-1 )
        ENDDO

        DEALLOCATE(iord)

        ! write final file

        OPEN(10,FILE=TRIM(foutput)//'.hdr')
        WRITE(10,'(3(i0,1x))') nx, ny, nz
        WRITE(10,'(5(i0,1x))') 5, nfl_new, nwl_clean, nin_new, nou_new
        CLOSE(10)

        OPEN(10,FILE=TRIM(foutput)//'.dat')
        DO n=0, nfinal-1

            i4 = itppp_final(n); i=ijk(1,i4); j=ijk(2,i4); k=ijk(3,i4)

            IF(n>0) THEN
                IF(i4<=itppp_final(n-1)) THEN
                    print *,n,'ERROR in unordering:',itppp_final(n-1),i4
                    STOP
                ENDIF
            ENDIF

            ! check that fluids do not touch deads
            !DO ip=1, npop-1
            !    ii = i + icx(ip); jj = j + icy(ip); kk = k + icz(ip)
            !    i4p = i4back(ii,jj,kk)
            !    IF(iflag_fina2(n)==2) CYCLE
            !    IF ( binsearch(itppp_final, nfinal-1, i4p) < 0 ) STOP '!OWEQE'
            !ENDDO

            WRITE(10,'(4(i0,1x))') ijk(1,i4), ijk(2,i4), ijk(3,i4), iflag_fina2(n)

        ENDDO
        CLOSE(10)

        INQUIRE(FILE=TRIM(flumen)//'.ios',EXIST=lexist)

        IF(.NOT. lexist) GOTO 11

        OPEN(9,STATUS='old',FILE=TRIM(flumen)//'.ios')
        OPEN(10,            FILE=TRIM(foutput)//'.ios')

        READ(9,*,ERR=11,END=11) n
        WRITE(10,*) n
        IF(n==0) GOTO 11

        DO i=1, n
            READ(9,'(a)') str
            WRITE(10,'(a)') TRIM(str)
        ENDDO

        READ(9,*) nn
        WRITE(10,*) nn

        DO n=1, nn
            READ(9,*) i,j,k,flg
            i4 = i4back(i,j,k)
            IF ( binsearch(itppp_pl, nv_pl-1, i4) > 0 ) THEN
                print *,n,nn,'ERROR.... Inlet and Plaques interfere',i,j,k
                STOP
            ELSE
                WRITE(10,'(4(i0,1x))') i,j,k,flg
            ENDIF
        ENDDO

        READ(9,*) nn
        WRITE(10,*) nn

        DO n=1, nn
            READ(9,*) i,j,k,flg
            i4 = i4back(i,j,k)
            IF ( binsearch(itppp_pl, nv_pl-1, i4) > 0 ) THEN
                print *,n,nn,'ERROR.... Outlet and Plaques interfere',i,j,k
                STOP
            ELSE
                WRITE(10,'(4(i0,1x))') i,j,k,flg
            ENDIF
        ENDDO

11      CONTINUE
        CLOSE(10)
        CLOSE(9)


        print *,'...done.'

    END SUBROUTINE combine

END MODULE combine_lumen_plaques_mod

