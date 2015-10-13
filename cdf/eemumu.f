       PROGRAM eemumu
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       DOUBLE PRECISION me,mu,f
       COMMON/inpu/me,mu,ebeam,const,sq
       COMMON/tell/nn
       COMMON/ini/xxx,yyy
       COMMON/outp/nout
       COMMON/cuts/angcut,encut,etacut
*
       INTEGER ndim,npoin,nprin,ntreat,nevent
*     
*     --- LPAIR data COMMON block
*     
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*     
*     --- HEPEVT COMMON block
       PARAMETER (NMXHEP=100)
       COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &      JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &      VHEP(4,NMXHEP)
       REAL PHEP,VHEP
*     
       EXTERNAL f
*     
       CALL fileini(IPAR,LPAR)  ! initialize variables
       CALL pawini              ! initialize histograms
c       CALL genzini             ! initialize genz structures
*     
       pi = dacos(-1.d+00)
       nout=6
       
       IF(ipar(4).EQ.1) THEN
          OPEN (12,file='ineemm.data',status='old')
          REWIND 12
       ENDIF
       
       IF(ipar(17).EQ.1) THEN
c          OPEN(20,file='events.ascii',status='new')
          OPEN(20,file='events.ascii',status='unknown')
       ENDIF
       
       nn=0
*     
*     ---- particle masses (GeV)
       me = lpar(1)             ! incoming particle
       mu = lpar(2)             ! outgoing particle
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
*     
*     ---- constant (convert GeV**2 to picobarns)
       const = (19.732d+03)**2
*     
*     ---- VEGAS integration
*     
       xx = ranf(1211)
c       OPEN (15,file='dl2.vegas.grid',status='new')
       OPEN (15,file='dl2.vegas.grid',status='unknown')
       CALL vegas(f,0.1d-03,ipar(5),ipar(6),ipar(7),+1,ipar(4))
       CALL save(7,15)
       PRINT*,'Wrote VEGAS grid to dl2.vegas.grid'
       CLOSE (15)
*     
*     ---- cross-section calculation only: exit PROGRAM
*     
       IF(ipar(1).EQ.1) THEN
          PRINT*,'IPAR(1) = 1: cross-section calculation complete'
          STOP
       ENDIF       
*     
*     ---- SETGEN: find local min/max
*     
       ndim   = ipar(5)
       npoin  = ipar(14)
       nbin   = ipar(15)
       nprin  = 1
       ntreat = ipar(13)
*     
       OPEN (15,file='dl2.vegas.grid',status='old')
       CALL restr(7,15)
       PRINT*,'Read dl2.vegas.grid'
       CLOSE (15)
c       OPEN (16,file='dl2.lattice.1',status='new')
       OPEN (16,file='dl2.lattice.1',status='unknown')
       CALL setgen(f,ndim,npoin,nbin,nprin,ntreat)
       PRINT*,'setgen complete'
       CALL save2(7,16)
       PRINT*,'Wrote SETGEN maxima to dl2.lattice.1'
       CLOSE (16)
*     
*     ---- GENERA: generate some events
*     
       ndim   = ipar(5)
       nevent = ipar(12)
       nstrat = 0
       ntreat = ipar(13)
*     
       OPEN (15,file='dl2.vegas.grid',status='old')
       CALL restr(7,15)
       PRINT*,'Read dl2.vegas.grid'
       CLOSE (15)
       OPEN (16,file='dl2.lattice.1',status='old')
       CALL restr2(7,16)
       CLOSE (16)
       PRINT*,'Read maxima from dl2.lattice.1'
c       OPEN (17,file='dl2.lattice.2',status='new')
       OPEN (17,file='dl2.lattice.2',status='unknown')
       xxxx = ranf(1236785)
       CALL genera(f,ndim,nevent,nstrat,ntreat)
       CALL save2(7,17)
       PRINT*,'Wrote new maxima to dl2.lattice.2'
       CLOSE (17)
*     
       CLOSE(20)                ! CLOSE events.ascii
*     
       CALL pawend              ! CLOSE histograms
c       CALL genzend             ! CLOSE genz
*     
       PRINT*,'GENZ events generated: NEVHEP = ',NEVHEP
*     
       stop
       END
