c     LHE file composition
c
c     Authors :
c     ...
c     N. Schul (UCL, Louvain-la-Neuve)
c     L. Forthomme (UCL, Louvain-la-Neuve), Feb 2014
c---------------------------------------------------------------------------

      SUBROUTINE LHEBEG
      IMPLICIT none
      INTEGER nout,ilhef
      COMMON/outp/nout,ilhef

      OPEN(ilhef,file='events.lhe', status='unknown')
      WRITE(ilhef,1000) '<LesHouchesEvents version="1.0">'

 1000 FORMAT((a))
      END

c---------------------------------------------------------------------------

      SUBROUTINE LHEHDR
      IMPLICIT none
      INTEGER nout,ilhef
      COMMON/outp/nout,ilhef

      DOUBLE PRECISION xsec,err,s3,s4
      COMMON/result/xsec,err,s3,s4
      INTEGER ipar(20)
      DOUBLE PRECISION lpar(20)
      COMMON/datapar/ipar,lpar
      
      WRITE(ilhef,1000) '<header>Sample generated using LPAIR</header>'
      WRITE(ilhef,1000) '<init>'
      WRITE(ilhef,1100) ipar(8),ipar(8),lpar(3),lpar(3),
     +                  ' 0 0 10042 10042 2 1'
      WRITE(ilhef,1200) xsec,err,' 0.2673112E-03 0'
      WRITE(ilhef,1000) '</init>'

 1000 FORMAT((a))
 1100 FORMAT(2I5, 2F12.1, (A))
 1200 FORMAT(1P, 2E12.4, (A))

      END

c---------------------------------------------------------------------------

      SUBROUTINE LHEEND
      IMPLICIT none
      INTEGER nout,ilhef
      COMMON/outp/nout,ilhef

      WRITE(ilhef,1000) '</LesHouchesEvents>'
      CLOSE(ilhef)

 1000 FORMAT((a))

      END

c---------------------------------------------------------------------------

      SUBROUTINE LHEFIL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c     LUND common
      REAL p(4000,5)
      INTEGER n,k(4000,5)
      COMMON/lujets/n,k,p
      save /lujets/

      INTEGER nout,ilhef
      COMMON/outp/nout,ilhef

c     Event common
      LOGICAL accepted
      INTEGER ndim,leppdg
      DOUBLE PRECISION x(10)
      COMMON/event/accepted,ndim,x,leppdg

      INTEGER I,J

      WRITE(ilhef,*) '<event>'
      WRITE(ilhef,*) N,'  661   0.2983460E-04  0.9118800E+02',
     &     '0.7821702E-02  0.1300000E+00'
      
      DO 202 I=1,N
c         WRITE(ilhef,5300) K(I,2),K(I,1),K(I,3),K(I,4),0,0,
         WRITE(ilhef,5300) K(I,2),K(I,1),K(I,3),0,0,0,
     &        (P(I,J),J=1,5),' 0. 1.'
 202  CONTINUE
      
      WRITE(ilhef,*) '</event>'

 5300 FORMAT(1P,I8,5I6,5E18.10,A6)
      END
      
