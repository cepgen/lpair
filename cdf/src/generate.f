       SUBROUTINE generate(nevent)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       DOUBLE PRECISION me,mu,f,mxcut
       LOGICAL accepted
       INTEGER leppdg
       COMMON/inpu/me,mu,ebeam,const,sq
       COMMON/tell/nn
       COMMON/ini/xxx,yyy
       COMMON/outp/nout,ilhef
       COMMON/cuts/angcut,encut,etacut,mxcut
       COMMON/event/accepted,ndim,x,leppdg
*
       INTEGER ndim,npoin,nprin,ntreat,nevent
*     
*     --- LPAIR data COMMON block
*     
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
       DOUBLE PRECISION x(10)
*     
*     
       EXTERNAL f
*     

*     General definitions
       ndim = ipar(5)
       ntreat = ipar(13)
       npoin= ipar(14)
       nbin = ipar(15)
       nprin = 1
       nstrat = 0
       n2dim = 7
       leppdg = ipar(18)

*     First we fetch the Vegas integration grid
c       OPEN (15,file='dl2.vegas.grid',status='old')
c       CALL restr(n2dim,15)
c       CLOSE (15)

*     Setgen
c       OPEN (16,file='dl2.lattice.1',status='unknown')
       CALL setgen(f,ndim,npoin,nbin,nprin,ntreat)
c       CALL save2(n2dim,16)
c       CLOSE (16)

*     
c       OPEN (15,file='dl2.vegas.grid',status='old')
c       CALL restr(n2dim,15)
c       CLOSE (15)

c       OPEN (16,file='dl2.lattice.1',status='old')
c       CALL restr2(n2dim,16)
c       CLOSE (16)

*     
c       OPEN (17,file='dl2.lattice.2',status='unknown')
       CALL genera(f,ndim,nevent,nstrat,ntreat)
c       CALL save2(n2dim,17)
c       CLOSE (17)
c       IF(IPAR(2).EQ.2) THEN
c         CALL LHEFIL
c       ENDIF
*     
*     
       END
