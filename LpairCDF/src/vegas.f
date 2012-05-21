      SUBROUTINE vegas(fxn,bcc,ndim,ncall,itmx,nprn,igraph)
      IMPLICIT DOUBLE PRECISION ( a-h,o-z )
      PARAMETER(MAXDIV=250)
      COMMON/bveg2/ndo,it,si,si2,swgt,schi,xi(MAXDIV,10),scalls
     + ,d(MAXDIV,10),di(MAXDIV,10),nxi(MAXDIV,10)
      DOUBLE PRECISION xin(MAXDIV),r(MAXDIV),dx(10),dt(10)
      INTEGER ia(10),kg(10)
      DOUBLE PRECISION xl(10),xu(10),qran(10),x(10)
      COMMON/result/s1,s2,s3,s4
      EXTERNAL fxn
      DATA xl,xu/10*0.D0,10*1.D0/
      DATA ndmx/MAXDIV/,alph/1.5D0/,one/1.D0/,mds/1/
      ipr = 1
      IF ( nprn .GT. 0 ) ipr = 0
      ndo = 1
      DO 1 j=1,ndim
1     xi(1,j) = one
      ENTRY vegas1(fxn,bcc,ndim,ncall,itmx,nprn,igraph)
      now=igraph
      IF ( igraph .GT. 0 ) CALL inplot(now,f1,w)
      it     = 0
      si     = 0
      si2    = si
      swgt   = si
      schi   = si
      scalls = si
      ENTRY vegas2(fxn,bcc,ndim,ncall,itmx,nprn,igraph)
      nd = ndmx
      ng = 1
      IF ( mds .EQ. 0 ) GO TO 2
      ng = ( ncall*0.5D0 ) ** ( 1.D0/ndim )
      mds = 1
      IF ( (2*ng-ndmx) .LT. 0 ) GO TO 2
      mds = -1
      npg = ng/ndmx+1
      nd = ng/npg
      ng = npg*nd
2     k = ng**ndim
      npg = ncall/k
      IF ( npg .LT. 2 ) npg = 2
      calls = npg*k
      dxg = one/ng
      dv2g = (dxg**(2*ndim))/npg/npg/(npg-one)
      xnd = nd
      ndm = nd-1
      dxg = dxg*xnd
      xjac = one
      DO 3 j = 1,ndim
      dx(j) = xu(j)-xl(j)
3     xjac  = xjac*dx(j)
      IF ( nd .EQ. ndo ) GO TO 8
      rc = ndo/xnd
      DO 7 j=1,ndim
      k  = 0
      xn = 0
      dr = xn
      i  = k
4     k  = k+1
      dr = dr+one
      xo = xn
      xn = xi(k,j)
5     IF ( rc .GT. dr ) GO TO 4
      i  = i+1
      dr = dr-rc
      xin(i) = xn-(xn-xo)*dr
      IF ( i .LT. ndm ) GO TO 5
      DO 6 i = 1,ndm
6     xi(i,j)  = xin(i)
7     xi(nd,j) = one
      ndo = nd
      acc = bcc
8     IF ( ( nprn .NE. 0 ) .AND. ( nprn .NE. 10 ) ) PRINT 200
     +  ,ndim,calls,it,itmx,acc,mds,nd
      IF ( nprn .EQ. 10 ) PRINT 290,ndim,calls,itmx,acc,mds,nd
      ENTRY vegas3(fxn,bcc,ndim,ncall,itmx,nprn,igraph)
9     it = it+1
      ti = 0
      tsi = ti
      IF ( igraph .GT. 0 ) CALL replot(now,f1,w)
      DO 10 j = 1,ndim
      kg(j) = 1
      DO 10 i = 1,nd
      nxi(i,j) = 0
      d(i,j)  = ti
10    di(i,j) = ti
11    fb = 0
      f2b = fb
      k = 0
12    k = k+1
      DO 121 j = 1,ndim
121   qran(j) = ranf(0)
      wgt = xjac
      DO 15 j = 1,ndim
      xn = (kg(j)-qran(j))*dxg+one
      ia(j) = xn
      iaj   = ia(j)
      iaj1  = iaj-1
      IF ( iaj .GT. 1 ) GO TO 13
      xo = xi(iaj,j)
      rc = (xn-iaj)*xo
      GO TO 14
13    xo   = xi(iaj,j)-xi(iaj1,j)
      rc   = xi(iaj1,j)+(xn-iaj)*xo
14    x(j) = xl(j)+rc*dx(j)
15    wgt  = wgt*xo*xnd
      f    = fxn(x)*wgt
      f1   = f/calls
      w    = wgt/calls
      IF ( igraph .GT. 0 ) CALL xplot(now,f1,w)
      f2   = f*f
      fb   = fb+f
      f2b  = f2b+f2
      DO 16 j=1,ndim
       iaj  = ia(j)
       nxi(iaj,j) = nxi(iaj,j)+1
       di(iaj,j) = di(iaj,j)+f/calls
       IF ( mds .GE. 0 ) d(iaj,j) = d(iaj,j)+f2
16    CONTINUE
      IF ( k .LT. npg ) GO TO 12
      f2b = f2b*npg
      f2b = DSQRT(f2b)
      f2b = (f2b-fb)*(f2b+fb)
      ti  = ti+fb
      tsi = tsi+f2b
      IF ( mds .GE. 0 ) GO TO 18
      DO 17 j = 1,ndim
      iaj = ia(j)
17    d(iaj,j) = d(iaj,j)+f2b
18    k = ndim
19    kg(k) = MOD(kg(k),ng)+1
      IF ( kg(k) .NE. 1 ) GO TO 11
      k = k-1
      IF ( k .GT. 0 ) GO TO 19
      ti  = ti/calls
      tsi = tsi*dv2g
      ti2 = ti*ti
      wgt = ti2/tsi
      si  = si+ti*wgt
      si2 = si2+ti2
      swgt= swgt+wgt
      schi= schi+ti2*wgt
      scalls = scalls+calls
      avgi = si/swgt
      sd   = swgt*it/si2
      chi2a= 0
      IF ( it .GT. 1 ) chi2a = sd*(schi/swgt-avgi*avgi)/(it-1)
      sd = one/sd
      sd = DSQRT(sd)
      IF ( nprn .EQ. 0 ) GO TO 21
      tsi = DSQRT(tsi)
      IF ( nprn .NE. 10 ) PRINT 201,ipr,it,ti,tsi,avgi,sd,chi2a
      IF ( nprn .EQ. 10 ) PRINT 203,it,ti,tsi,avgi,sd,chi2a
      IF ( nprn .GE.  0 ) GO TO 21
      DO 20 j = 1,ndim
      PRINT 202,j
20    PRINT 204,(xi(i,j),di(i,j),d(i,j),i=1,nd)
21    IF ( DABS(sd/avgi) .LE. DABS(acc) .OR. it .GE. itmx ) now = 2
      s1 = avgi
      s2 = sd
      s3 = ti
      s4 = tsi
      IF ( igraph .GT. 0 ) CALL plotit(now,f1,w)
*      DO 23 j=1,ndim
*      xo=d(1,j)
*      xn=d(2,j)
*      d(1,j)=(xo+xn)*0.5D0
*      dt(j)=d(1,j)
*      DO 22 i=2,ndm
*      d(i,j)=xo+xn
*      xo=xn
*      xn=d(i+1,j)
*      d(i,j)=(d(i,j)+xn)/3
*22    dt(j)=dt(j)+d(i,j)
*      d(nd,j)=(xn+xo)/2
*23    dt(j)=dt(j)+d(nd,j)
*-----this part of the vegas-algorithm is unstable
*-----it should be replaced by
      DO 23 j = 1,ndim
      dt(j) = 0
      DO 23 i = 1,nd
      IF ( nxi(i,j) .GT. 0 ) d(i,j) = d(i,j)/nxi(i,j)
      dt(j) = dt(j)+d(i,j)
23    CONTINUE
      DO 28 j=1,ndim
      rc = 0
      DO 24 i = 1,nd
      r(i) = 0
      IF ( d(i,j) .LE. 0 ) GO TO 24
      xo = dt(j)/d(i,j)
      r(i) = ((xo-one)/xo/dlog(xo))**alph
24    rc = rc+r(i)
      rc = rc/xnd
      k  = 0
      xn = 0
      dr = xn
      i  = k
25    k  = k+1
      dr = dr+r(k)
      xo = xn
      xn = xi(k,j)
26    IF ( rc .GT. dr ) GO TO 25
      i  = i+1
      dr = dr-rc
      xin(i) = xn-(xn-xo)*dr/r(k)
      IF (i .LT. ndm ) GO TO 26
      DO 27 i = 1,ndm
27    xi(i,j) = xin(i)
28    xi(nd,j)= one
      IF ( it .LT. itmx .AND. DABS(acc) .LT. DABS(sd/avgi) ) GO TO 9
200   FORMAT(35h0input parameters for vegas   ndim=,i3
     +,8h  ncall=,f8.0/28x,5h  it=,i5,8h  itmx =,i5/28x
     +,6h  acc=,g9.3/28x,6h  mds=,i3,6h   nd=,i4//)
290   FORMAT(13h0vegas  ndim=,i3,8h  ncall=,f8.0,8h  itmx =,i5
     +,6h  acc=,g9.3,6h  mds=,i3,6h   nd=,i4)
201   FORMAT(/i1,20hintegration by vegas/13h0iteration no,i3,
     +14h.   integral =,g14.8/20x,10hstd dev  =,g10.4/
     +34h accumulated results.   integral =,g14.8/
     +24x,10hstd dev  =,g10.4 / 24x,18hchi**2 per itn   =,g10.4)
202   FORMAT(14h0data for axis,i2 / 7x,1hx,7x,10h  delt i  ,
     +2x,11h convce    ,11x,1hx,7x,10h  delt i  ,2x,11h convce
     +,11x,1hx,7x,10h  delt i  ,2x,11h convce     /)
204   FORMAT(1x,3g12.4,5x,3g12.4,5x,3g12.4)
203   FORMAT(1h ,i3,g20.8,g12.4,g20.8,g12.4,g12.4)
      s1 = avgi
      s2 = sd
      s3 = chi2a
      RETURN
      END

