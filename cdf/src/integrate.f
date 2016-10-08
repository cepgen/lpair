       SUBROUTINE integrate
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       DOUBLE PRECISION me,mu,f,mxcut
       COMMON/inpu/me,mu,ebeam,const,sq
       COMMON/tell/nn
       COMMON/ini/xxx,yyy
       COMMON/outp/nout,ilhef
       COMMON/cuts/angcut,encut,etacut,mxcut
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
       REAL ulmass ! from jetset/pythia
       EXTERNAL f,ulmass
*     
       pi = dacos(-1.d+00)
       nout=6
       
       nn=0
*     
*     ---- particle masses (GeV)
       me = lpar(1)             ! incoming particle
       mu = ulmass(ipar(18))    ! outgoing particle
*     
*     ---- beam energy (GeV)
       ebeam = lpar(3)
       sq=2.*ebeam
*     
*     ---- angle cuts
       angcut = dcos(pi*lpar(4))
*     
*     ---- rapidity cuts
       etacut = lpar(5)
*     
*     ---- energy cuts
       encut = lpar(6)

       mxcut = lpar(10)
*     
*     ---- constant (convert GeV**2 to picobarns)
       const = (19.732d+03)**2
*     
*     ---- VEGAS integration
*     
       xx = ranf(1211)
c       do 10 i=1,ipar(5)
c          x(i) = 0.4
c 10    continue
c       xx = f(x)
c       print *,'x=',(x(i),i=1,ipar(5))
c       print *,xx
       
       OPEN (15,file='dl2.vegas.grid',status='unknown')
       CALL vegas(f,0.1d-03,ipar(5),ipar(6),ipar(7),0,ipar(4))
       CALL save(ipar(5),15)
       CLOSE (15)
*
       IF(IPAR(2).eq.2) THEN
         ILHEF=15
         call LHEBEG
         call LHEHDR
       ENDIF
*     
*     
       END
