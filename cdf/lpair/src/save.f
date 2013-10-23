      SUBROUTINE save(ndim,ntape)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(MAXDIV=250)
      COMMON/bveg2/ndo,it,si,si2,swgt,schi,xi(MAXDIV,10),scalls
     + ,d(MAXDIV,10),di(MAXDIV,10),nxi(MAXDIV,10)
*
*   stores vegas data (unit ntape) for later re-initialization
*
      WRITE(ntape,200) ndo,it,si,si2,swgt,schi
      WRITE(ntape,201)
     +      ((xi(i,j),i=1,ndo),j=1,ndim)
     +     ,((di(i,j),i=1,ndo),j=1,ndim)
      RETURN
      ENTRY restr(ndim,ntape)
*
*   enters initialization data for vegas
*
      READ(ntape,200) ndo,it,si,si2,swgt,schi
      READ(ntape,201)
     +      ((xi(i,j),i=1,ndo),j=1,ndim)
     +     ,((di(i,j),i=1,ndo),j=1,ndim)
200   FORMAT(2i8,4D24.16)
201   FORMAT(5D24.16)
*      WRITE(6,200) ndo,it,si,si2,swgt,schi
*      WRITE(6,201)
*     +      ((xi(i,j),i=1,ndo),j=1,ndim)
*     +     ,((di(i,j),i=1,ndo),j=1,ndim)

      RETURN
      END
