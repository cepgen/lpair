*-- Author :    O. Duenger   17/12/91
      SUBROUTINE VEGAS(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C
C  SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG*N
C    - BY G.P.  LEPAGE  SEPT 1976/ (REV) APR 1978
C
C    -FTN5 VERSION 21-8-1984
C
C  AUTHOR                                       : G. P. LEPAGE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FXN
      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
      DIMENSION XL(10),XU(10),QRAN(10),X(10)
      COMMON/VGASIO/NINP,NOUTP,NSGOUT
      COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     + ,D(50,10),DI(50,10)
      COMMON/VGRES/S1,S2,S3,S4
      REAL S1,S2,S3,S4
      real ran2
      integer idum
C
C
      DATA XL,XU/10*0.,10*1./
      DATA NDMX/50/,ALPH/1.5/,ONE/1./,MDS/1/
      data idum/-1/
C
      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'VEGAS CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
      IF(NPRN.GT.0)THEN
         IPR=0
      ELSE
         IPR=1
      ENDIF
      NDO=1
      DO 1 J=1,NDIM
         XI(1,J)=ONE
1     CONTINUE
C
      ENTRY VEGAS1(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C...INITIALIZES CUMMULATIVE VARIABLES,BUT NOT GRID
      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'VEGAS1 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
C
      IT=0
      SI=0.
      SI2=SI
      SWGT=SI
      SCHI=SI
      SCALLS=SI
C
      ENTRY VEGAS2(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C...NO INITIALIZATION
      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'VEGAS2 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
      ND=NDMX
      NG=1
      IF(MDS.NE.0)THEN
         NG=(NCALL/2.)**(1./NDIM)
         MDS=1
         IF((2*NG-NDMX).GE.0)THEN
            MDS=-1
            NPG=NG/NDMX+1
            ND=NG/NPG
            NG=NPG*ND
         ENDIF
      ENDIF
C
      K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2) NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE
c      print *,k,ng,calls,dxg,dv2g,xnd,ndm
c      stop
      DO 3 J=1,NDIM
         DX(J)=XU(J)-XL(J)
         XJAC=XJAC*DX(J)
c         print *,j,xjac,dx(j)
3     CONTINUE
C
C  REBIN PRESERVING BIN DENSITY
C
      IF(ND.NE.NDO)THEN
c         print *,'rebinning',nd,ndo
         RC=NDO/XND
         DO 7 J=1,NDIM
            K=0
            XN=0.
            DR=XN
            I=K
4           K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(K,J)
5           IF(RC.GT.DR) GO TO 4
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
c            print *,i,j,k,dr,xin(i)
            IF(I.LT.NDM) GO TO 5
            DO 6  I=1,NDM
               XI(I,J)=XIN(I)
c               print *,i,j,xi(i,j)
6           CONTINUE
            XI(ND,J)=ONE
7        CONTINUE
         NDO=ND
      ENDIF
c      print *,nd,ndo,dr,xn
c      stop
C
       IF(NPRN.NE.0.AND.NPRN.NE.10)WRITE(NOUTP,200)NDIM,CALLS,IT,ITMX
     + ,ACC,MDS,ND
       IF(NPRN.EQ.10)WRITE(NOUTP,290)NDIM,CALLS,ITMX,ACC,MDS,ND
C
      ENTRY VEGAS3(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C     - MAIN INTEGRATION LOOP
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'VEGAS3 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
9     CONTINUE
      IT=IT+1
      TI=0.
      TSI=TI
C
      DO 10 J=1,NDIM
         KG(J)=1
         DO 10 I=1,ND
            D(I,J)=TI
            DI(I,J)=TI
10    CONTINUE
C
11    FB=0.
      F2B=FB
      K=0
C
12    CONTINUE
      K=K+1
      DO 121 J=1,NDIM
         QRAN(J)=ran2(idum)
c         QRAN(j) = 0.7
121   CONTINUE
      WGT=XJAC
      DO 15 J=1,NDIM
         XN=(KG(J)-QRAN(J))*DXG+ONE
         IA(J)=XN
         IAJ=IA(J)
         IAJ1=IAJ-1
c         print *,xn,iaj,xo
         IF(IAJ.LE.1)THEN
            XO=XI(IAJ,J)
            RC=(XN-IAJ)*XO
         ELSE
            XO=XI(IAJ,J)-XI(IAJ1,J)
            RC=XI(IAJ1,J)+(XN-IAJ)*XO
         ENDIF
         X(J)=XL(J)+RC*DX(J)
c         print *,j,x(j)
c         stop
         WGT=WGT*XO*XND
c         print *,'-->',kg(j),iaj,xn
15    CONTINUE
C
      F=FXN(X)*WGT
c      print *,'passed!',k,f
      F1=F/CALLS
      W=WGT/CALLS
C
      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
         IAJ=IA(J)
         DI(IAJ,J)=DI(IAJ,J)+F/CALLS
         IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
c         print *,iaj,D(iaj,j),Di(iaj,j)
16    CONTINUE
c      print *,wgt,xo,f,f1,w,fb,f2b
      IF(K.LT.NPG) GO TO 12
C
      F2B=F2B*NPG
      F2B=DSQRT(F2B)
      F2B=DABS((F2B-FB)*(F2B+FB))
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.LT.0)THEN
         DO 17 J=1,NDIM
            IAJ=IA(J)
            D(IAJ,J)=D(IAJ,J)+F2B
c            print *,d(iaj,j),f2b
17       CONTINUE
      ENDIF

      K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
c      print *,kg(k),'--->',ng
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
c      do 314 j=1,ndim
c         do 315 i=1,ndm
c            print *,'->',i,j,d(i,j)
c 315     continue
c 314  continue
c      stop
C
C  FINAL RESULTS FOR THIS ITERATION
C
      TI=TI/CALLS
      TSI=TSI*DV2G
      TI2=TI*TI
c      print *,ti,tsi,ti2
c      STOP
      IF(TSI .EQ. 0.)THEN
         WGT = 0.
      ELSE
         WGT=TI2/TSI
      ENDIF
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      IF(SWGT .EQ. 0.)THEN
         AVGI=TI
      ELSE
         AVGI=SI/SWGT
      ENDIF
      IF(SI2 .EQ. 0.)THEN
         SD=TSI
      ELSE
         SD=SWGT*IT/SI2
      ENDIF
      SCALLS=SCALLS+CALLS
      CHI2A=0.
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      IF(SD .NE. 0.)THEN
         SD=ONE/SD
         SD=DSQRT(SD)
      ELSE
         SD=TSI
      ENDIF
      IF(NPRN.NE.0)THEN
         TSI=DSQRT(TSI)
         IF(NPRN.NE.10)WRITE(NOUTP,201)IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.EQ.10)WRITE(NOUTP,203)IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.LT.0)THEN
            DO 20 J=1,NDIM
               WRITE(NOUTP,202)J
               WRITE(NOUTP,204)(XI(I,J),DI(I,J),D(I,J),I=1,ND)
20         CONTINUE
         ENDIF
      ENDIF
C
C   REFINE GRID
C
21    IF(SD .NE. 0.)THEN
         REL = DABS(SD/AVGI)
      ELSE
         REL = 0.
      ENDIF
      IF(REL.LE.DABS(ACC).OR.IT.GE.ITMX)NOW=2
      S1=AVGI
      S2=SD
      S3=TI
      S4=TSI
C

      DO 23 J=1,NDIM
         XO=D(1,J)
         XN=D(2,J)
         D(1,J)=(XO+XN)/2.
         DT(J)=D(1,J)
         DO 22 I=2,NDM
            D(I,J)=XO+XN
            XO=XN
            XN=D(I+1,J)
            D(I,J)=(D(I,J)+XN)/3.
            DT(J)=DT(J)+D(I,J)
c            print *,'====>',i,xo,xn,dt(j),d(i,j)
22       CONTINUE
         D(ND,J)=(XN+XO)/2.
         DT(J)=DT(J)+D(ND,J)
23    CONTINUE
C
      DO 28 J=1,NDIM
         RC=0.
         DO 24 I=1,ND
            R(I)=0.
            IF(D(I,J).GT.0.)THEN
               XO=DT(J)/D(I,J)
               R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
            ENDIF
c            print *,'-->',i,r(i),xo,dlog(xo)
            RC=RC+R(I)
24       CONTINUE
         RC=RC/XND
c         print *,rc
c         stop
         K=0
         XN=0.
         DR=XN
         I=K
25       K=K+1
         DR=DR+R(K)
         XO=XN
         XN=XI(K,J)
c         print *,'-->',k,dr,xo,xn
26       IF(RC.GT.DR) GO TO 25
         I=I+1
         DR=DR-RC
         IF(DR .EQ. 0.)THEN
            XIN(I)=XN
         ELSE
            XIN(I)=XN-(XN-XO)*DR/R(K)
         ENDIF
c         print *,i,xn,xo,dr,r(k)
c         print *,i,xin(i),dr,rc,r(k),xn,xo
c         print *,i,xin(i),dr,r(k)
         IF(I.LT.NDM) GO TO 26
c         stop
         DO 27 I=1,NDM
            XI(I,J)=XIN(I)
27       CONTINUE
         XI(ND,J)=ONE
28    CONTINUE
C
c      print *,dabs(acc),rel
      IF(IT.LT.ITMX.AND.DABS(ACC).LT.REL)GO TO 9
C
      S1=AVGI
      S2=SD
      S3=CHI2A
      RETURN
C
199   FORMAT(A)
200   FORMAT('INPUT PARAMETERS FOR VEGAS   NDIM=',I3
     +,'  NCALL=',F8.0/28X,'  IT=',I5,'  ITMX =',I5/28X
     +,'  ACC=',D9.3/28X,'  MDS=',I3,'   ND=',I4//)
290    FORMAT('VEGAS  NDIM=',I3,'  NCALL=',F8.0,'  ITMX =',I5
     + ,'  ACC=',D9.3,'  MDS=',I3,'   ND=',I4)
201   FORMAT('INTEGRATION BY VEGAS'/'ITERATION NO',I3,
     +'.   INTEGRAL =',D14.8/20X,'STD DEV  =',D10.4/
     +' ACCUMULATED RESULTS.   INTEGRAL =',D14.8/
     +24X,'STD DEV  =',D10.4 / 24X,'CHI**2 PER ITN   =',D10.4)
202   FORMAT('DATA FOR AXIS',I2 / 7X,'X',7X,'  DELT I  ',
     +2X,' CONVCE    ',11X,'X',7X,'  DELT I  ',2X,' CONVCE     '
     +,11X,'X',7X,'  DELT I  ',2X,' CONVCE     '/)
204   FORMAT(1X,3D12.4,5X,3D12.4,5X,3D12.4)
203   FORMAT(1X,I3,D20.8,D12.4,D20.8,D12.4,D12.4)
C
      END
