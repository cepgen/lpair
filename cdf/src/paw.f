      subroutine pawfil1
*
* ---- LPAIR common blocks
*
      implicit double precision (a-h,o-z)
*
      double precision me,mu,mxcut
      common/inpu/me,mu,ebeam,const,sq
      common/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5
     1                                         ,st5,cp3,sp3,cp5,sp5
      common/variac/al3,al4,be4,be5,de3,de5,pp3,pp4,pp5
      common/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
      common/lplot/xl(10),v1(2),v2(2),av(10)
      common/extra/s1,s2,t1,t2
      common/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
      common/levi/gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
      common/civita/epsi,g5,g6,a5,a6,bb
      common/dotps/q1dq,q1dq2,w6
      common/tell/nn
      common/cuts/angcut,encut,etacut,mxcut
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
     +              px6,py6,pz6,px7,py7,pz7

      REAL PT3,PT4,PT5,PT6,PT7

      DIMENSION P12(4)
      REAL PI,XTRA1(23)
*     
      PI = ACOS(-1.0)
      DPI = DACOS(-1.0D00)
*
* ---- fill MYGENZ common ntuple entries (before reflection/rotation)
*
      PX3 = PP3*CP3
      PY3 = PP3*SP3
      PZ3 = P3*CT3
      PT3 = SQRT(PX3**2+PY3**2)
      PX5 = PP5*CP5
      PY5 = PP5*SP5
      PZ5 = P5*CT5
      PT5 = SQRT(PX5**2+PY5**2)
      PX6 = PP6*CP6
      PY6 = PP6*SP6
      PZ6 = P6*CT6
      PT6 = SQRT(PX6**2+PY6**2)
      PX7 = PP7*CP7
      PY7 = PP7*SP7
      PZ7 = P7*CT7
      PT7 = SQRT(PX7**2+PY7**2)
*     
      PX4 = PP4*CP4
      PY4 = 0.
      PZ4 = P4*CT4
      PT4 = SQRT(PX4**2+PY4**2)
*
      if(ipar(16).eq.0) goto 100
*
*****************************************************************
*
* --- do we wish to random reflect/rotate?
* --- random reflection in x-z plane (VERMAS10)

      if(ranf(dummy) .ge. 0.5) then
         sp3 = -sp3
         sp5 = -sp5
         sp6 = -sp6
         sp7 = -sp7
      endif
* --- random rotation around z-axis (VERMAS10)

      ranphi=2.0*acos(-1.)*ranf(dummy)
      sinphi=sin(ranphi)
      cosphi=cos(ranphi)

* --- rotate, reflect and transform values (VERMAS10)

      PEMX=PP3*(CP3*COSPHI-SP3*SINPHI)
      PEMY=PP3*(CP3*SINPHI+SP3*COSPHI)
      PEMZ=P3*CT3
      PEPX=PP5*(CP5*COSPHI-SP5*SINPHI)
      PEPY=PP5*(CP5*SINPHI+SP5*COSPHI)
      PEPZ=P5*CT5
      PP6=P6*ST6
      PMMX=PP6*(CP6*COSPHI-SP6*SINPHI)
      PMMY=PP6*(CP6*SINPHI+SP6*COSPHI)
      PMMZ=P6*CT6
      PMPX=-PEMX-PEPX-PMMX
      PMPY=-PEMY-PEPY-PMMY
      PMPZ=-PEMZ-PEPZ-PMMZ
*
      PRESX = PP4*(CP4*COSPHI)
      PRESY = PP4*(CP4*SINPHI)
      PRESZ = P4*CT4
*
      px3 = pemx
      py3 = pemy
      pz3 = pemz
      pt3 = sqrt(px3**2+py3**2)
      px5 = pepx
      py5 = pepy
      pz5 = pepz
      pt5 = sqrt(px5**2+py5**2)
      px6 = pmmx
      py6 = pmmy
      pz6 = pmmz
      pt6 = sqrt(px6**2+py6**2)
      px7 = pmpx
      py7 = pmpy
      pz7 = pmpz
      pt7 = sqrt(px7**2+py7**2)
*
      px4 = presx
      py4 = presy
      pz4 = presz
      pt4 = sqrt(px4**2+py4**2)
*
c      px3 = p3*st3*cp3
c      py3 = p3*st3*sp3
c      pz3 = p3*ct3
c      pt3 = p3*st3
c      px5 = p5*st5*cp5
c      py5 = p5*st5*sp5
c      pz5 = p5*ct5
c      pt5 = p5*st5
c      px6 = p6*st6*cp6
c      py6 = p6*st6*sp6
c      pz6 = p6*ct6
c      pt6 = p6*st6
c      px7 = p7*st7*cp7
c      py7 = p7*st7*sp7
c      pz7 = p7*ct7
c      pt7 = p7*st7
*
*****************************************************************
*
 100  continue
*
      IF(IPAR(2).NE.1) RETURN
*
      XTRA1(1)  = real(dsign(dlog((dsqrt(pt6**2+pz6**2)+
     &                              dabs(pz6))/pt6),pz6))
      XTRA1(2)  = real(dsign(dlog((dsqrt(pt7**2+pz7**2)+
     &                              dabs(pz7))/pt7),pz7))
      XTRA1(3)  = real(dsqrt(w4))
      
      XTRA1(4)  = real(px6)
      XTRA1(5)  = real(py6)
      XTRA1(6)  = real(pz6)
      XTRA1(7)  = real(pt6)
      XTRA1(8)  = real(px7)
      XTRA1(9)  = real(py7)
      XTRA1(10) = real(pz7)
      XTRA1(11) = real(pt7)
      XTRA1(12) = real(p4*st4)
      
* ---- acoplanarity angle (deviation from PI)

      PHI6 = PYANGL(PX6,PY6)
      PHI7 = PYANGL(PX7,PY7)
      IF(PHI7.LT.0.0D0) THEN
         PHI7 = PHI7 + 2.*DPI
      ENDIF
      IF(PHI6.LT.0.0D0) THEN
         PHI6 = PHI6 + 2.*DPI
      ENDIF
      IF(PHI6.LT.PHI7) THEN
         DPHI = PHI7 - PHI6
         ACOPL = DPI - DPHI
      ELSEIF(PHI6.GT.PHI7) THEN
         DPHI = PHI6 - PHI7
         ACOPL = DPHI - DPI
      ENDIF
      XTRA1(13) = real(phi6)
      XTRA1(14) = real(phi7)
      XTRA1(15) = real(acopl)

* ---- coplanarity angle (transverse plane, 2D)

      P12(1) = PX6*PX7
      P12(2) = PY6*PY7
      P12(3) = PZ6*PZ7
      P12(4) = P12(1) + P12(2)
      ARG    = (P12(4)/(PT6*PT7))
      IF(ARG.LT.(-1.0D00)) ARG=-1.0D00
      THETOPN2 = (180.0D00/DPI)*DACOS(ARG)

* ---- collinearity angle ((x,y,z) plane, 3D)

      P12(1) = PX6*PX7
      P12(2) = PY6*PY7
      P12(3) = PZ6*PZ7
      P12(4) = P12(1) + P12(2) + P12(3)
      ARG    = (P12(4)/(P6*P7))
      IF(ARG.LT.(-1.0D0)) ARG=-1.0D00
      THETOPN3 = (180.0D00/DPI)*DACOS(ARG)
      
      XTRA1(16) = REAL(THETOPN2)
      XTRA1(17) = REAL(THETOPN3)

* ---- calculate deflection of outgoing protons

      IF(PZ3.LT.0.) THEN
         XTRA1(18) = PI - REAL(ATAN(PT3/DABS(PZ3)))
      ELSE
         XTRA1(18) = REAL(ATAN(PT3/PZ3))
      ENDIF
      IF(PZ5.LT.0.) THEN
         XTRA1(19) = PI - REAL(ATAN(PT5/DABS(PZ5)))
      ELSE
         XTRA1(19) = REAL(ATAN(PT5/PZ5))
      ENDIF
*     
      XTRA1(20) = PT3
      XTRA1(21) = PZ3
      XTRA1(22) = PT5
      XTRA1(23) = PZ5
*
      return
      end
