      SUBROUTINE genera(f,ndim,nevent,nstrat,ntreat)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION am,ami,ranf,y,treat,g,h,a,b,t,z,tf,tot
      DOUBLE PRECISION ffmax,fmax
      INTEGER mbin,max,m,ndim,ncall,mcall,nev,mev,nstrat,nn,j,n(10)
      INTEGER k,jj,jjj,ncycle,ntreat,nm,mdum,nevent
      EXTERNAL f
      
      LOGICAL accepted
      INTEGER endim,leppdg
      DOUBLE PRECISION x(10)
      COMMON/event/accepted,endim,x,leppdg

      COMMON/maxi/mdum,mbin,ffmax,fmax(80000),nm(80000)

      endim = ndim
      
      accepted = .FALSE.
      am  = mbin
      ami = 1.D0/am
      max = mbin**ndim
      nev = 0
      mev = 0
      ncall = 0
      mcall = 0
      tot = 0                   !FIXME
      m = 0
      IF ( nstrat .GT. 0 ) GO TO 10
      nn = 0
 1    nn = nn+1
      j  = ranf(0)*max+1
      y  = ranf(0)*ffmax
      nm(j) = nm(j) + 1
      IF ( y .GT. fmax(j) ) GO TO 1
      jj = j-1
      DO 2 k = 1,ndim
         jjj  = jj/mbin
         n(k) = jj-jjj*mbin
         x(k) = (ranf(0)+n(k))*ami
         jj=jjj
 2    CONTINUE
      IF ( ntreat .GT. 0 ) THEN
         g = treat(f,x,ndim)
      ELSE
         g = f(x)
      ENDIF
      ncall = ncall+1
      IF ( y .GT. g ) GO TO 1
      CALL accept(1)
      nev = nev+1
      IF ( nev .GE. nevent ) GO TO 8
      IF ( g .LE. fmax(j) ) GO TO 1
 3    IF ( nm(j) .EQ. 1 ) GO TO 7
      a = nm(j)*(g-fmax(j))/ffmax
      m = m+1
      t = a
 4    IF ( t .LT. 1 ) THEN
         IF ( ranf(0) .GT. t ) GO TO 7
         t = 1
      ENDIF
      t = t-1
      DO 6 k = 1,ndim
         x(k) = (ranf(0)+n(k))*ami
 6    CONTINUE
      IF ( ntreat .GT. 0 ) THEN
         z = treat(f,x,ndim)
      ELSE
         z = f(x)
      ENDIF
      mcall = mcall+1
      ncall = ncall+1
      IF ( z .GE. fmax(j) ) THEN
         CALL accept(1)
         mev = mev+1
         nev = nev+1
         IF ( nev .GE. nevent ) GO TO 7
         IF ( z .GT. g ) THEN
            b = nm(j)*(z-fmax(j))/ffmax
            t = t+(b-a)
            a = b
            g = z
         ENDIF
      ENDIF
      IF ( t .GT. 0 ) GO TO 4
 7    fmax(j) = g
      IF ( nstrat .GT. 0 ) GO TO 18
      IF ( g .GT. ffmax ) ffmax = g
      IF ( nev .LT. nevent ) GO TO 1
 8    CONTINUE
c      PRINT 100,nev,ncall,nn
c      PRINT 101,m,mev,mcall
      accepted = .TRUE.
      RETURN
      
 10   ncycle = 0
 11   ncycle = ncycle+1
      h  = ffmax
      tf = nstrat/ffmax
      j  = 0
 12   jj = j
      j  = j+1
      DO 13 k = 1,ndim
         jjj  = jj/mbin
         n(k) = jj-jjj*mbin
         jj   = jjj
 13   CONTINUE
      tot  = fmax(j)*tf
      nm(j)= nm(j)+nstrat
      g    = fmax(j)
 14   IF ( tot .GE. 1. ) GO TO 15
      IF ( tot .LE. 0 ) GO TO 17
      IF ( ranf(0) .GT. tot ) GO TO 17
 15   tot = tot-1.
      DO 16 k = 1,ndim
         x(k) = (ranf(0)+n(k))*ami
 16   CONTINUE
      IF ( ntreat .GT. 0 ) THEN
         z = treat(f,x,ndim)
      ELSE
         z = f(x)
      ENDIF
      ncall = ncall+1
      IF ( z .LT. ranf(0)*fmax(j) ) GO TO 14
      IF ( z .GT. g ) g = z
      nev = nev+1
      CALL accept(1)
      GO TO 14
 17   IF ( g .GT. fmax(j) ) GO TO 3
 18   IF ( g .GT. h ) h = g
      IF ( j .LT. max ) GO TO 12
      CALL accept(0)
      ffmax = h
      IF ( nev .LT. nevent ) GO TO 11
      PRINT 102,nev,ncall,ncycle,nstrat
      PRINT 101,m,mev,mcall
      CALL accept(-1)
      accepted = .TRUE.
      RETURN
      
 100  FORMAT(32h1subroutine genera has produced ,i10,8h events./
     1     11h this took ,i10,19h function calls in ,i10,10h attempts./)
 101  FORMAT(1h0,i10,29h times a maximum was changed./
     1     11h this gave ,i8,11h events in ,i9,7h calls.//)
 102  FORMAT(32h1subroutine genera has produced ,i10,8h events./
     1     11h this took ,i10,16h function calls./
     2     i10,11h cycles of ,i8,27h points per bin were needed,
     3     28h during stratified sampling.///)
      
      END
