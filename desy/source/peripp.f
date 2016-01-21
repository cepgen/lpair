*-- Author :    O. Duenger   17/12/91
      DOUBLE PRECISION FUNCTION PERIPP(NUP,NDOWN)
*      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*      
      COMMON /PICKZZ/ W1,W2,W3,W4,W5,W31,W52,W12,TAU,SL1
      COMMON /CIVITA/ EPSI,G5,G6,A5,A6,BB
      COMMON /EXTRA/  S1,S2,T1,T2
      COMMON /DOTPS/  Q1DQ,Q1DQ2,W6
      COMMON /LEVI/   GRAM,D1,D2,D3,D4,D5,DELTA,G4,A1,A2
      COMMON /PERIC/  U1,U2,V1,V2,T11,T12,T21,T22
*
      DATA RHO /.585D+00/,Cc1/0.86926/,Cc2/2.23422/,Dd1/0.12549/
      DATA Cp/0.96/,Bp/0.63/
!     DATA RHO /1.05/,   Cc1/0.6303/, Cc2/2.2049/, Dd1/0.0468/
!     DATA Cp/1.23/,Bp/0.61/

      REAL*8 DUMMY,PSFW1,PSFW2,M1
      REAL*8 W1STRFUN,W2STRFUN

c      print *,'peripp','nup=',nup,'ndown=',ndown
      IF(NUP.EQ.1) THEN
         U1 = 1.
         U2 = 1.
      ELSEIF(NUP.EQ.2) THEN
         XT  = 1. - T1 / .71
         XT  = XT * XT
         XU  = 2.79 / XT
         U1  = XU * XU
         TAU = T1 / (4.*W1)
         U2  = (1./(XT*XT) - XU * XU * TAU) / (1.-TAU)
      ELSEIF (NUP .EQ. 4) THEN
         M1=DSQRT(W1)
         CALL PSF(T1,W3,DUMMY,PSFW1,PSFW2)
         U1=-PSFW1*2D0*M1/T1
         U2= PSFW2/2D0/M1
      ELSEIF (NUP .EQ. 5) THEN
         M1=DSQRT(W1)
         CALL W1W2F2(T1,W3,W1STRFUN,W2STRFUN)
         U1=-W1STRFUN*2.*M1/T1
         U2= W2STRFUN/2./M1
      ELSE
         X    = T1  / (T1-W3)
         EN   = W31 - T1
         TAU  = T1  / (4.*W1)
         RHOT = RHO - T1
         U1   = (-Cc1*RHO*RHO*W31/RHOT/RHOT-Cc2*W1*(1.-X)**4
     &        / (X*(X*Cp-2*Bp)+1.)) / T1
         U2   = (-TAU*U1-Dd1*W31*T1*RHO/RHOT/RHOT*W31*W31/EN/EN/W1)
     &        / (1.-EN*EN/(4.*W1*T1))
      ENDIF
      IF(NDOWN.EQ.1) THEN
         V1   = 1.
         V2   = 1.
      ELSEIF(NDOWN.EQ.2) THEN
         XT  = 1. - T2 / .71
         XT  = XT * XT
         XU  = 2.79 / XT
         V1  = XU * XU
         TAU = T2 / (4.*W2)
         V2  = (1./(XT*XT) - XU * XU * TAU) / (1.-TAU)
      ELSE
         X    = T2  / (T2-W5)
         EN   = W52 - T2
         TAU  = T2  / (4.*W2)
         RHOT = RHO - T2
         V1   = (-Cc1*RHO*RHO*W52/RHOT/RHOT-Cc2*W2*(1.-X)**4
     &        / (X*(X*Cp-2*Bp)+1.))/T2
         V2   = (-TAU*V1-Dd1*W52*T2*RHO/RHOT/RHOT*W52*W52/EN/EN/W2)
     &        / (1.-EN*EN/(4.*W2*T2))
      ENDIF
      QQQ = Q1DQ * Q1DQ
      QDQ = 4. * W6 - W4
      T22 = 512. * (BB*(DELTA*DELTA-GRAM)-(EPSI-DELTA*(QDQ+Q1DQ2))**2
     &     -A1*A6*A6-A2*A5*A5-A1*A2*QQQ)
      T12 = 128. * (-BB*(D2+G6)-2.*(T1+2.*W6)*(A2*QQQ+A6*A6))*T1
      T21 = 128. * (-BB*(D4+G5)-2.*(T2+2.*W6)*(A1*QQQ+A5*A5))*T2
      T11 = 64.*(BB*(QQQ-G4-QDQ*(T1+T2+2.*W6))-2.*(T1+2.*W6)*(T2+2.*W6)
     &     *QQQ)*T1*T2
      PERIPP = (((U1*V1*T11+U2*V1*T21+U1*V2*T12+U2*V2*T22)/(T1*T2*BB))
     &     / (T1*T2*BB))*.25
      
      RETURN
      END
