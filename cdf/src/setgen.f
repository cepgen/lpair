      SUBROUTINE setgen(f,ndim,npoin,nbin,nprin,ntreat)
      INTEGER i,j,k,m,mbin,mmax,ndim,npoin,nprin,ntreat,n(10),nm
      INTEGER mdum,jj,jjj,nbin
      DOUBLE PRECISION f,ffmax,sum,sum2,sum2p,fmax,x(10),fsum,fsum2,sigp
      DOUBLE PRECISION av,av2,sig,sig2,eff1,eff2,DSQRT,ranf,treat,z,eff
      EXTERNAL f
      COMMON/maxi/mdum,mbin,ffmax,fmax(80000),nm(80000)
      DOUBLE PRECISION maxar(11,100)
      INTEGER kk,kkk
      LOGICAL hassetgen
      DATA hassetgen/.FALSE./
      SAVE hassetgen
 
*     
*     One does not want to perform setgen twice
*     (awfully time consuming...)
*
      IF ( hassetgen ) RETURN
      hassetgen = .TRUE.

      DO 50 i = 1,100
         DO 50 j = 1,11
 50   maxar(j,i) = 0

      mbin  = nbin
      ffmax = 0
      sum   = 0
      sum2  = 0
      sum2p = 0
      mmax  = mbin**ndim
      IF ( nprin .GE. 2 ) PRINT 200,mbin,mmax,npoin
      DO 5 j = 1,mmax
         nm(j)   = 0
         fmax(j) = 0
 5    CONTINUE
      DO 1 j = 1,mmax
         jj = j-1
         DO 2 k=1,ndim
            jjj  = jj/mbin
            n(k) = jj-jjj*mbin
            jj   = jjj
 2       CONTINUE
         fsum  = 0
         fsum2 = 0
         DO 3 m = 1,npoin
            DO 4 k=1,ndim
               x(k)=(ranf(0)+n(k))/mbin
 4          CONTINUE
            IF ( ntreat .GT. 0 ) THEN
               z = treat(f,x,ndim)
*      PRINT *,m,(x(i),i=1,ndim),z
            ELSE
               z = f(x)
            ENDIF
            IF ( z .GT. fmax(j) ) fmax(j) = z
            fsum  = fsum+z
            fsum2 = fsum2+z*z

            IF ( z .GT. maxar(ndim+1,100) ) THEN
               DO 51 kk = 99,1,-1
                  IF ( z .GT. maxar(ndim+1,kk) ) THEN
                     DO 52 kkk = 1,ndim+1
 52                  maxar(kkk,kk+1) = maxar(kkk,kk)
                  ELSE
                     GOTO 53
                  ENDIF
 51            CONTINUE
 53            CONTINUE
               maxar(ndim+1,kk+1) = z
               DO 54 kkk = 1,ndim
 54            maxar(kkk,kk+1) = x(kkk)
            ENDIF
            
 3       CONTINUE
         av    = fsum/npoin
         av2   = fsum2/npoin
         sig2  = av2-av*av
         sig   = DSQRT(sig2)
         sum   = sum+av
         sum2  = sum2+av2
         sum2p = sum2p+sig2
         IF ( fmax(j) .GT. ffmax ) ffmax = fmax(j)
         eff = 10000
         IF ( fmax(j) .NE. 0 ) eff = fmax(j)/av
         IF ( nprin .GE. 3 )PRINT 100,j,av,sig,fmax(j),eff,
     +        (n(i),i=1,ndim)
 1    CONTINUE
      sum   = sum/mmax
      sum2  = sum2/mmax
      sum2p = sum2p/mmax
      sig   = DSQRT(sum2-sum*sum)
      sigp  = DSQRT(sum2p)
      eff1  = 0.
      DO 6 j = 1,mmax
         eff1 = eff1+fmax(j)
 6    CONTINUE
      eff1 = eff1/(mmax*sum)
      eff2 = ffmax/sum
      IF ( nprin .GE. 1 ) PRINT 101,sum,sig,sigp,ffmax,eff1,eff2
      IF ( nprin .GE. 4 ) THEN
         PRINT *,'The 100 highest function values are:'
         DO 56 kk = 1,100
            PRINT 110,kk,maxar(ndim+1,kk),(maxar(kkk,kk),kkk=1,ndim)
 110        FORMAT(i4,g14.4,10f14.10)
 56      CONTINUE
      ENDIF
      RETURN
 100  FORMAT(i6,3x,g13.6,g12.4,g13.6,f8.2,3x,10i1)
 101  FORMAT(29h the average function value =,g14.6/
     1     29h the overall std dev        =,g14.4/
     2     29h the average std dev        =,g14.4/
     3     29h the maximum function value =,g14.6/
     4     29h the average inefficiency   =,g14.3/
     5     29h the overall inefficiency   =,g14.3/)
 200  FORMAT(25h subroutine setgen uses a,i3,15h**ndim division/
     1     17h this results in ,i7,6h cubes/
     2     17h the program put ,i5,29h points in each cube to find
     3     ,30hstarting values for the maxima//)
      END
