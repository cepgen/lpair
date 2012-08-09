       SUBROUTINE save2(ndim,ntape)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       INTEGER nm
       COMMON/maxi/mdum,mbin,ffmax,fmax(80000),nm(80000)
       max=mbin**ndim
       WRITE (ntape,100) mbin,ffmax
100    FORMAT(i10,D24.16)
       WRITE (ntape,101) (fmax(i),i=1,max)
101    FORMAT(5D24.16)
       WRITE (ntape,102) (nm(i),i=1,max)
102    FORMAT(8i10)
       RETURN
       ENTRY restr2(ndim,ntape)
       REWIND(ntape)
       READ (ntape,100) mbin,ffmax
       max = mbin**ndim
       READ (ntape,101) (fmax(i),i=1,max)
       READ (ntape,102) (nm(i),i=1,max)
       RETURN
       END
