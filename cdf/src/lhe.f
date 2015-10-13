c     LHE file composition
c
c     Authors :
c     ...
c     N. Schul (UCL, Louvain-la-Neuve)
c     L. Forthomme (UCL, Louvain-la-Neuve), Feb 2014

      SUBROUTINE LHEFIL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c     Pair produced code
      INTEGER NLHE
      INTEGER L,NLEP 

c     LUND common
      REAL*4        P(4000,5)
      INTEGER       N,K(4000,5)
      COMMON/LUJETS/N,K,P
      save /lujets/

c     Event common
      LOGICAL ACCEPTED
      INTEGER NDIM
      DOUBLE PRECISION x(10)
      COMMON/event/accepted,ndim,x

      INTEGER I,J

      OPEN(22,file='events.out', status='unknown')

      NLHE=N                    !lf
      
clf      NLHE=0
clf      DO 203 I=1,N
clf         IF(K(I,1).EQ.1.AND.K(I,2).EQ.2212.AND.I.LT.16) THEN
clf            NLHE=NLHE+1
clf         ENDIF
clf         IF(K(I,1).EQ.1.AND.I.GT.15)        NLHE=NLHE+1
clf         IF(K(I,1).EQ.1.AND.I.GT.15.AND.K(I,2).EQ.12) NLHE=NLHE-1
clf         IF(K(I,1).EQ.1.AND.I.GT.15.AND.K(I,2).EQ.-12) NLHE=NLHE-1
clf         IF(K(I,1).EQ.1.AND.I.GT.15.AND.K(I,2).EQ.14) NLHE=NLHE-1
clf         IF(K(I,1).EQ.1.AND.I.GT.15.AND.K(I,2).EQ.-14) NLHE=NLHE-1
clf         IF(K(I,1).EQ.1.AND.I.GT.15.AND.K(I,2).EQ.16) NLHE=NLHE-1
clf         IF(K(I,1).EQ.1.AND.I.GT.15.AND.K(I,2).EQ.-16) NLHE=NLHE-1
clf 203  CONTINUE
      WRITE(22,*) '<event>'
      WRITE(22,*) NLHE,'  661   0.2983460E-04  0.9118800E+02',
     &     '0.7821702E-02  0.1300000E+00'
      
      DO 202 I=1,N
         
clf         IF(K(I,1).EQ.1.AND.K(I,2).EQ.2212.AND.I.LT.14) THEN
clf            WRITE(22,*) '2212 1 1 2 0 0 0.0 0.0 ',P(I,4), ! -P(I,4) normally quoted
clf     &           P(I,4),P(I,5),' 0. 1'
clf         ENDIF
clf         IF(I.GT.13.AND.K(I,1).EQ.1) THEN
clf            IF(K(I,2).NE.12.AND.K(I,2).NE.-12.AND.K(I,2).NE.14) THEN
clf               IF(K(I,2).NE.-14.AND.K(I,2).NE.16.
clf     &              AND.K(I,2).NE.-16) THEN
clf                  WRITE(22,*) K(I,2),' 1 1 2 0 0',
clf     &                 P(I,1),P(I,2),P(I,3), ! P(I,3) normally quoted 
clf     &                 P(I,4),P(I,5),' 0. 1'
clf               ENDIF
clf            ENDIF
clf         ENDIF
         WRITE(22,5300) K(I,2),K(I,1),K(I,3),K(I,4),0,0,
     &        (P(I,J),J=1,5),' 0. 1.'
clf         WRITE(22,*) IDUP(I),ISTUP(I),MOTHUP(1,I),
clf     &        MOTHUP(2,I),ICOLUP(1,I),ICOLUP(2,I),(PUP(J,I),J=1,5),
clf     &        ' 0. 9.'
 202  CONTINUE
c     FOR MUONS
clf      WRITE(22,*) '13  1 1 2 0 0 ',PPXM1,PPYM1,PPZM1,PPEM1,
clf     &     ' 0.1057 0. 1'
clf      WRITE(22,*) '-13 1 1 2 0 0 ',PPXM2,PPYM2,PPZM2,PPEM2,
clf     &     ' 0.1057 0. 1'
      
c     FOR TAUS
c     WRITE(22,*) '15  1 1 2 0 0 ',PPXM1,PPYM1,PPZM1,PPEM1,
c     &  ' 1.77684 0. 1'
c     WRITE(22,*) '-15 1 1 2 0 0 ',PPXM2,PPYM2,PPZM2,PPEM2,
c     &  ' 1.77684 0. 1'
      
      WRITE(22,*) '</event>'

      CLOSE(22)

 5300 FORMAT(1P,I8,5I5,5E18.10,A6)
      END
      
