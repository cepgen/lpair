       subroutine genzfil
*
       implicit double precision (a-h,o-z)
       double precision me,mu
*
* --- HEPEVT common block
       PARAMETER (NMXHEP=10000)
       COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
       REAL PHEP,VHEP
       INTEGER I,J
*
* --- LPAIR common blocks
*
       common/inpu/me,mu,ebeam,const,sq
       common/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5
     1 ,st5,cp3,sp3,cp5,sp5
       common/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
       common/lplot/xl(10),v1(2),v2(2),av(10)
       common/extra/s1,s2,t1,t2
       common/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
*
* --- GENZ quantities common block
*


*
* --- LPAIR data common block
*
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*
* --- MYGENZ common block
*
       common/mygenz/px3,py3,pz3,px4,py4,pz4,px5,py5,pz5,
     +               px6,py6,pz6,px7,py7,pz7
*
*
       IF(IPAR(3).EQ.0) RETURN
*
       IF(NEVHEP.GE.IPAR(11)) RETURN ! max # of events to GNZ output
*
* ---- begin GENZ (LPAIR)
*
* ---- fill HEPEVT common block
*
      NHEP = 7 ! elastic production
*
* ---- store entries for:
*   1  proton       (incoming, left)
*   2  proton       (incoming,right)
*   3  proton       (outgoing, left)
*   5  proton       (outgoing,right)
*   4  muon-pair (resonance product)
*   6  muon         (outgoing, left)
*   7  muon         (outgoing,right)
*
* ---- initialize all HEPEVT entries to zero
      do i=1,nhep
        isthep(i)=0.
         idhep(i)=0.
        do j=1,2
          jmohep(j,i)=0.
          jdahep(j,i)=0.
        enddo
        do j=1,5
          phep(j,i)=0.
        enddo
        do j=1,4
          vhep(j,i)=0.
        enddo
      enddo
*
* ---- proton (incoming, 1)
      isthep(1)=1.
      idhep(1)=2212.
      jmohep(1,1)=0.
      jmohep(2,1)=0.
      jdahep(1,1)=0.
      jdahep(2,1)=0.
      phep(1,1)=real(p1)*( 0.)
      phep(2,1)=real(p1)*( 0.)
      phep(3,1)=real(p1)*( 1.)
      phep(4,1)=real(e1)
      phep(5,1)=real(me)
*
* ---- proton (incoming, 2)
      isthep(2)=1.
      idhep(2)=2212.
      jmohep(1,2)=0.
      jmohep(2,2)=0.
      jdahep(1,2)=0.
      jdahep(2,2)=0.
      phep(1,2)=real(p2)*( 0.)
      phep(2,2)=real(p2)*( 0.)
      phep(3,2)=real(p2)*(-1.)
      phep(4,2)=real(e2)
      phep(5,2)=real(me)
*
* ---- proton (outgoing, 3)
      isthep(3)=1.
      idhep(3)=2212.
      jmohep(1,3)=1.
      jmohep(2,3)=1.
      jdahep(1,3)=0.
      jdahep(2,3)=0.
      phep(1,3)=real(px3)
      phep(2,3)=real(py3)
      phep(3,3)=real(pz3)
      phep(4,3)=real(e3)
      phep(5,3)=real(me)
*
* ---- proton (outgoing, 5)
      isthep(5)=1.
      idhep(5)=2212.
      jmohep(1,5)=2.
      jmohep(2,5)=2.
      jdahep(1,5)=0.
      jdahep(2,5)=0.
      phep(1,5)=real(px5)
      phep(2,5)=real(py5)
      phep(3,5)=real(pz5)
      phep(4,5)=real(e5)
      phep(5,5)=real(me)
*
* ---- lepton pair (resonance, 4)
      isthep(4)=2.
      idhep(4)=93.
      jmohep(1,4)=1.
      jmohep(2,4)=2.
      jdahep(1,4)=6.
      jdahep(2,4)=7.
      phep(1,4)=real(px4)
      phep(2,4)=real(py4)
      phep(3,4)=real(pz4)
      phep(4,4)=real(e4)
      phep(5,4)=real(sqrt(w4))
*
* ---- muon - (outgoing, 6)
      isthep(6)=1.
      idhep(6)=13.
      jmohep(1,6)=4.
      jmohep(2,6)=4.
      jdahep(1,6)=0.
      jdahep(2,6)=0.
      phep(1,6)=real(px6)
      phep(2,6)=real(py6)
      phep(3,6)=real(pz6)
      phep(4,6)=real(e6)
      phep(5,6)=real(mu)
*
* ---- muon + (outgoing, 7)
      isthep(7)=1.
      idhep(7)=-13.
      jmohep(1,7)=4.
      jmohep(2,7)=4.
      jdahep(1,7)=0.
      jdahep(2,7)=0.
      phep(1,7)=real(px7)
      phep(2,7)=real(py7)
      phep(3,7)=real(pz7)
      phep(4,7)=real(e7)
      phep(5,7)=real(mu)
*
*
* ... increment counter of lepton pair events produced
*
      NEVHEP = NEVHEP + 1
*
* ... print some diagnostics
*
      WRITE(20,*) NEVHEP, NHEP
      DO JJ=1,NHEP
        WRITE(20,*) JJ, ISTHEP(JJ),IDHEP(JJ),JMOHEP(1,JJ),JMOHEP(2,JJ),
     &             JDAHEP(1,JJ),JDAHEP(2,JJ)
        WRITE(20,*) PHEP(1,JJ),PHEP(2,JJ),PHEP(3,JJ),PHEP(4,JJ),
     &             PHEP(5,JJ)
        WRITE(20,*) VHEP(1,JJ),VHEP(2,JJ),VHEP(3,JJ),VHEP(4,JJ)
      ENDDO
*
*      translate HEPEVT to GENZ
*
c        CALL GNZFRHC(IRET)
*
c        IF(IRET.NE.0) THEN
c           WRITE(6,*) 'IRET = ',IRET,' on return from GNZFRHC. STOP'
c           CALL GNZEND
c        ENDIF
*
*      write out event
*
c        CALL GNZWRIT(IRET)
c        IF(IRET.NE.0) THEN
c          WRITE(6,*) 'GNZWRIT: ERROR ...EXIT GNZEND'
c          CALL GNZEND
c        ENDIF
*
*      GENZ printing of first ten events
c        IF(NEVHEP.LE.10) THEN
c           CALL GNZPRIN(1,4)
c        ENDIF
*
*  --- end GENZ (LPAIR)
*
       return
       end
