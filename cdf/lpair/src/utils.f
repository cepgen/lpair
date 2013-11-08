      SUBROUTINE maps2(s2,x,smin,smax,ds)
      IMPLICIT DOUBLE PRECISION (a-z)
      y=smax/smin
      s2=smin*y**x
      ds=s2*dlog(y)
      RETURN
      END

*----------------------------------------------------------------------
      
      SUBROUTINE mapla(x,y,z,u,xm,xp,d)
      IMPLICIT DOUBLE PRECISION (a-z)
      xmb=xm-y-z
      xpb=xp-y-z
      c=-4.*y*z
      alp=dsqrt(xpb*xpb+c)
      alm=dsqrt(xmb*xmb+c)
      am=xmb+alm
      ap=xpb+alp
      yy=ap/am
      zz=yy**u
      x=y+z+(am*zz-c/(am*zz))*0.5
      ax=dsqrt((x-y-z)**2+c)
      d=ax*dlog(yy)
      RETURN
      END
      
*----------------------------------------------------------------------

      SUBROUTINE inplot(now,ff,pdx)
      IMPLICIT DOUBLE PRECISION ( a-h,o-z)
      COMMON/lplot/xl(10),v1(2),v2(2),av(10)
      DIMENSION zav(10),yav(10),zsv(10),ysv(10),ztv(10)
      DIMENSION xlmax(10),xlmin(10),nlp(10),ltop(10),text(8,10),ll(10)
      DIMENSION numb(12)
      DIMENSION xls(42,10),yls(42,10),nlsn(42,10),mlsn(42,10),dls(10)
     1,xlav(10),xlsq(10),xlava(10),sxa(10),tlim(6),top(10),xltq(10)
      DIMENSION nbin(41),nlog(41),slog(41),tlog(41),hv(12)
      DIMENSION v1max(2),v1min(2),v2max(2),v2min(2),nv1(2)
     1,nv2(2),vtext(6,2)
      DIMENSION vm(12,12,2),nvm(12,12,2),bin1(2),bin2(2),vol(2)
     1,wm(12,12,2),mvm(12,12,2)
      COMMON/result/y,si,u,v
      character*1 hmin,hplus,hblank,hstar,char(40)
      save
      data tlim/1.6d+00,2.5d+00,4.0d+00,6.666666667d+00,10.d+00,
     +16.d+00/
      data hmin,hplus,hblank,hstar/'-','+',' ','*'/
      data mls,mav,ndmax/10,10,2/
      data ngraph/0/
      igraph=now
      now=0
      kk=0
      itt=0
      IF(igraph.le.0) GO TO 800
      IF(igraph.ne.ngraph)read (12,810)nls
      IF(igraph.ne.ngraph)PRINT 814,nls
      IF(nls.LT.0) nls=0
      IF(nls.EQ.0) GO TO 802
      IF(nls.GT.mls) GO TO 807
      IF(igraph.ne.ngraph)PRINT 815
      do 801 i=1,nls
      IF(igraph.ne.ngraph)
     1read (12,811)xlmin(i),xlmax(i),nlp(i),ltop(i),ll(i),
     2(text(j,i),j=1,8)
      IF(igraph.ne.ngraph) PRINT 816,i,xlmin(i),xlmax(i)
     1,nlp(i),ltop(i),ll(i),(text(j,i),j=1,8)
      IF(nlp(i).LT.1)nlp(i)=1
      IF(nlp(i).GT.40)nlp(i)=40
      dls(i)=(xlmax(i)-xlmin(i))/nlp(i)
      nlps=nlp(i)+2
      do 300 j=1,nlps
      yls(j,i)=0
300   mlsn(j,i)=0
801   CONTINUE
802   IF(igraph.ne.ngraph) read (12,810)ndd
      IF(igraph.ne.ngraph) PRINT 817,ndd
      IF(ndd.LT.0) ndd=0
      IF(ndd.EQ.0) GO TO 804
      IF(ndd.GT.ndmax) GO TO 807
      IF(igraph.ne.ngraph) PRINT 818
      do 803 i=1,ndd
      IF(igraph.ne.ngraph)
     1read (12,812)v1min(i),v1max(i),nv1(i),v2min(i),v2max(i),nv2(i)
     1,(vtext(j,i),j=1,6)
      IF(igraph.ne.ngraph) PRINT 819,i,v1min(i),v1max(i),nv1(i)
     1,v2min(i),v2max(i),nv2(i),(vtext(j,i),j=1,6)
      IF(nv1(i).LT.1)nv1(i)=1
      IF(nv2(i).LT.1)nv2(i)=1
      IF(nv1(i).GT.10)nv1(i)=10
      IF(nv2(i).GT.10)nv2(i)=10
      bin1(i)=(v1max(i)-v1min(i))/nv1(i)
      bin2(i)=(v2max(i)-v2min(i))/nv2(i)
      vol(i)=bin1(i)*bin2(i)
803   CONTINUE
      wtow=0.
      do 805 i=1,ndd
      do 805 j=1,12
      do 805 k=1,12
      wm(k,j,i)=0.
805   mvm(k,j,i)=0
804   CONTINUE
      IF(igraph.ne.ngraph)read (12,810)nave
      IF(igraph.ne.ngraph)PRINT 820,nave
      IF(nave.LT.0)nave=0
      IF(nave.GT.mav) GO TO 807
      do 11 i=1,mav
      yav(i)=0.
11    ysv(i)=0.
      kt=0
      GO TO 808
800   nave=0
      nls=0
      ndd=0
      GO TO 808
807   PRINT 813,mav,mls,ndmax
      stop
808   ngraph=igraph
      RETURN
      entry replot(now,ff,pdx)
      IF(nave.EQ.0) GO TO 49
      do 62 i=1,nave
      zav(i)=0.
      ztv(i)=0.
62    zsv(i)=0.
49    fsqa=0.
      kt=kt+1
      IF(nls.EQ.0) GO TO 303
      do 302 i=1,nls
      nlps=nlp(i)+2
      xlav(i)=0
      xltq(i)=0.
      xlsq(i)=0
      do 302 j=1,nlps
      xls(j,i)=0
302   nlsn(j,i)=0
303   CONTINUE
      IF(ndd.EQ.0) GO TO 403
      do 402 i=1,ndd
      n1=nv1(i)+2
      n2=nv2(i)+2
      do 402 i1=1,n1
      do 402 i2=1,n2
      vm(i1,i2,i)=0
402   nvm(i1,i2,i)=0
403   CONTINUE
      RETURN
      entry xplot(now,ff,pdx)
      fsqa=fsqa+ff*ff/pdx
      itt=itt+1
      IF(nls.EQ.0) GO TO 305
      do 304  i=1,nls
      nlps=(xl(i)-xlmin(i))/dls(i)+1.
      IF(nlps.LT.0)nlps=0
      IF(nlps.GT.nlp(i))nlps=nlp(i)+1
      nlps=nlps+1
      xls(nlps,i)=xls(nlps,i)+ff/dls(i)
      nlsn(nlps,i)=nlsn(nlps,i)+1
      xlav(i)=xlav(i)+ff*xl(i)
      xltq(i)=xltq(i)+ff*ff*xl(i)/pdx
304   xlsq(i)=xlsq(i)+(ff*xl(i))**2/pdx
305   CONTINUE
      IF(ndd.EQ.0) GO TO 405
      do 404 i=1,ndd
      i1=(v1(i)-v1min(i))/bin1(i)+2
      IF(i1.LT.1) i1=1
      IF(i1.GT.nv1(i)+2) i1=nv1(i)+2
      i2=(v2(i)-v2min(i))/bin2(i)+2
      IF(i2.LT.1) i2=1
      IF(i2.GT.nv2(i)+2) i2=nv2(i)+2
      vm(i1,i2,i)=vm(i1,i2,i)+ff/vol(i)
404   nvm(i1,i2,i)=nvm(i1,i2,i)+1
405   CONTINUE
      IF(nave.EQ.0) GO TO 99
      do 22 i=1,nave
      zav(i)=zav(i)+av(i)*ff
      ztv(i)=ztv(i)+ff*ff*av(i)/pdx
22    zsv(i)=zsv(i)+(av(i)*ff)**2/pdx
99    RETURN
      entry plotit(now,ff,pdx)
      IF(nls.EQ.0) GO TO 315
      IF(kk.GT.0) GO TO 307
      do 306 i=1,nls
      nlps=nlp(i)+2
      do 306 j=1,nlps
      mlsn(j,i)=nlsn(j,i)
306   yls(j,i)=xls(j,i)
      GO TO 310
307   vbef=vtot
      vu=(v/u)**2
      do 309 i=1,nls
      nlps=nlp(i)+2
      do 309 j=1,nlps
      IF(nlsn(j,i).EQ.0) GO TO 309
      IF(mlsn(j,i).EQ.0) GO TO 308
      al1=vu/nlsn(j,i)
      al2=vbef/mlsn(j,i)
      mlsn(j,i)=mlsn(j,i)+nlsn(j,i)
      yls(j,i)=(al2*xls(j,i)+al1*yls(j,i))/(al1+al2)
      GO TO 309
308   mlsn(j,i)=nlsn(j,i)
      yls(j,i)=xls(j,i)
309   CONTINUE
310   CONTINUE
      do 311 i=1,nls
      sxf=xlsq(i)-xlav(i)*xlav(i)
      sxt=xltq(i)-xlav(i)*u
      sx2=xlsq(i)/xlav(i)**2+fsqa/u**2-2.*xltq(i)/(xlav(i)*u)
      sx2=sx2*(xlav(i)/u)**2
      IF(kt.ne.1) GO TO 312
      xlava(i)=xlav(i)/u
      sxa(i)=sx2
      GO TO 311
312   xhelp=sx2+sxa(i)
      IF(xhelp.EQ.0) GO TO 311
      xlava(i)=(xlav(i)*sxa(i)/u+xlava(i)*sx2)/xhelp
      sxa(i)=sxa(i)*sx2/xhelp
311   CONTINUE
      vtot=(si/y)**2
      IF(now.ne.2) GO TO 315
      do 341 i=1,nls
      top(i)=0.
      nlps=nlp(i)+1
      do 341 j=2,nlps
      xls(j,i)=yls(j,i)/y
      IF(xls(j,i).GT.top(i))top(i)=xls(j,i)
341   CONTINUE
      do 342 i=1,nls
      IF(ltop(i).le.0)ltop(i)=i
      lto=ltop(i)
      IF(top(i).GT.top(lto))top(lto)=top(i)
342   CONTINUE
      ylog=0.5*dlog10(y*y)
      do 314 i=1,nls
      PRINT 321,i
      nlps=nlp(i)+1
      lto=ltop(i)
      top(i)=top(lto)
      IF(top(i).EQ.0)top(i)=1.
      an1=dlog10(top(i))
      n1=an1
      IF(n1.GT.an1)n1=n1-1
      z1=top(i)*10.**(-n1)
      do 343 l=1,4
      IF(z1.LT.tlim(l)) GO TO 344
343   CONTINUE
      l=5
344   IF(top(i).LT.1.6/(xlmax(i)-xlmin(i)))l=l+1
      topm=tlim(l)*10.**n1
      do 345 j=2,nlps
      nbin(j)=xls(j,i)*40./topm+1.5
      IF(ll(i).LT.0)nbin(j)=0
      IF(xls(j,i).le.0) GO TO 346
      tlog(j)=dlog10(xls(j,i))
      slog(j)=tlog(j)+ylog
      nlog(j)=(tlog(j)-n1)*8.+33.5
      IF(ll(i).GT.0)nlog(j)=0
      GO TO 345
346   slog(j)=0
      tlog(j)=0
      nlog(j)=0
345   CONTINUE
      PRINT 322,(text(j,i),j=1,8)
      n1p1=n1+1
      n1m4=n1-4
      PRINT 323,tlim(l),n1,n1p1,n1m4
      do 347 l=1,40
      char(l)=hmin
      IF(nlog(l+1).EQ.41)char(l)=hplus
      IF(nbin(l+1).EQ.41)char(l)=hstar
347   CONTINUE
      xmin=xlmin(i)
      xmax=xmin+dls(i)
      PRINT 324,xmin,xmax,yls(2,i),slog(2),xls(2,i),tlog(2),mlsn(2,i)
     1,char
***********************
*      CALL pawfil2(   i,real(xmin),real(xmax),real(xls(2,i)))
*      CALL pawfil2(10+i,real(xmin),real(xmax),real(yls(2,i)))
***********************
      do 348 j=3,nlps
      xmin=xmax
      xmax=xmin+dls(i)
      do 349 l=1,40
      char(l)=hblank
      IF(nlog(l+1).EQ.43-j)char(l)=hplus
      IF(nbin(l+1).EQ.43-j)char(l)=hstar
349   CONTINUE
      PRINT 324,xmin,xmax,yls(j,i),slog(j),xls(j,i),tlog(j),mlsn(j,i)
     1,char
***********************
*      CALL pawfil2(   i,real(xmin),real(xmax),real(xls(j,i)))
*      CALL pawfil2(10+i,real(xmin),real(xmax),real(yls(j,i)))
***********************
348   CONTINUE
      nlps1=nlps+1
      IF(nlps.EQ.41) GO TO 352
      do 351 j=nlps1,41
      do 350 l=1,40
      char(l)=hblank
      IF(nlog(l+1).EQ.43-j)char(l)=hplus
      IF(nbin(l+1).EQ.43-j)char(l)=hstar
350   CONTINUE
351   PRINT 325,char
352   do 353 l=1,40
      char(l)=hmin
      IF(nlog(l+1).EQ.1)char(l)=hplus
      IF(nbin(l+1).EQ.1)char(l)=hstar
353   CONTINUE
      PRINT 326,char
      el1=yls(1,i)*dls(i)
      el2=el1/y
      PRINT 327,el1,el2,mlsn(1,i)
      el1=yls(42,i)*dls(i)
      el2=el1/y
      PRINT 328,el1,el2,mlsn(nlps1,i)
      sxsq=dsqrt(sxa(i)/itt)
      PRINT 329,xlava(i),sxsq
314   CONTINUE
315   CONTINUE
      IF(ndd.EQ.0) GO TO 409
      wbef=wtot
      do 500 i=1,ndd
      nx=nv1(i)+2
      ny=nv2(i)+2
      IF(kk.GT.0) GO TO 502
      do 501 j=1,nx
      do 501 k=1,ny
      wm(j,k,i)=vm(j,k,i)
501   mvm(j,k,i)=nvm(j,k,i)
      GO TO 500
502   vu=(v/u)**2
      do 503 j=1,nx
      do 503 k=1,ny
      IF(nvm(j,k,i).EQ.0) GO TO 503
      IF(mvm(j,k,i).EQ.0) GO TO 504
      al1=vu/nvm(j,k,i)
      al2=vbef/mvm(j,k,i)
      mvm(j,k,i)=mvm(j,k,i)+nvm(j,k,i)
      wm(j,k,i)=(al2*vm(j,k,i)+al1*wm(j,k,i))/(al1+al2)
      GO TO 503
504   mvm(j,k,i)=nvm(j,k,i)
      wm(j,k,i)=vm(j,k,i)
503   CONTINUE
500   CONTINUE
      wtot=(si/y)**2
      IF(now.ne.2) GO TO 409
      do 408 i=1,ndd
      PRINT 481,i,(vtext(j,i),j=1,6)
      vvv=v2max(i)
      mvv=nv1(i)+2
      nvv=nv2(i)+1
      size=vol(i)/y
      do 406 i2=1,nvv
      j2=nvv+2-i2
      do 410 i1=1,mvv
410   numb(i1)=1000.*wm(i1,j2,i)*size+.5
      PRINT 486,(numb(i1),i1=1,mvv)
      PRINT 483,(wm(i1,j2,i),i1=1,mvv)
      PRINT 486,(mvm(i1,j2,i),i1=1,mvv)
      PRINT 484,vvv
      vvv=vvv-bin2(i)
      IF(dabs(vvv/bin2(i)).LT.1.e-10)vvv=0.
406   CONTINUE
      do 411 i1=1,mvv
411   numb(i1)=1000.*wm(i1,1,i)*size+.5
      PRINT 486,(numb(i1),i1=1,mvv)
      PRINT 483,(wm(i1,1,i),i1=1,mvv)
      PRINT 486,(mvm(i1,1,i),i1=1,mvv)
      PRINT 482
      mvv=mvv-1
      do 407 i1=1,mvv
      hv(i1)=v1min(i)+(i1-1)*bin1(i)
      IF(dabs(hv(i1)/bin1(i)).LT.1.d-10)hv(i1)=0.
407   CONTINUE
      PRINT 485,(hv(i1),i1=1,mvv)
408   CONTINUE
409   CONTINUE
      IF(nave.EQ.0) GO TO 23
      IF(now.EQ.2) PRINT 26
      do 24 i=1,nave
      sxf=zsv(i)-zav(i)*zav(i)
      sxt=zsv(i)/zav(i)**2+fsqa/u**2-2.*ztv(i)/(zav(i)*u)
      sx2=sxt*(zav(i)/u)**2
      IF(kt.ne.1) GO TO 21
      yav(i)=zav(i)/u
      ysv(i)=sx2
      GO TO 30
21    xhelp=sx2+ysv(i)
      IF(xhelp.EQ.0) GO TO 30
      yav(i)=(ysv(i)*zav(i)/u+yav(i)*sx2)/xhelp
      ysv(i)=ysv(i)*sx2/xhelp
30    yssq=dsqrt(ysv(i)/itt)
      IF(now.EQ.2)PRINT 27,i,yav(i),yssq
24    CONTINUE
23    now=1
      kk=kk+1
      RETURN
27    format(12x,i2,9x,d15.5,5x,d15.3)
26    format(1h1,10x,46hthe following are averages with error estimate/)
321   format(1h1,40x,40hsingle dIFferential cross-section number,i3///)
322   format(38h single dIFferential cross section of ,8a4/)
323   format(11x,6hlimits,9x,1hi,16x,19haccumulated results,15x,1hi,24x
     1,9hupper bin,6x,9hlower bin/26x,1hi,50x,20hi * linear      plot,
     2f8.2,5h*10**,i3,8x,1h0/5x,5hlower,7x,5hupper,4x,1hi,5x,5hds/dx,4x
     3,56halog10   (ds/dx)/s  alog10  points  i + logarithmic plot,6x
     4,4h10**,i3,8x,4h10**,i3/2x,24(1h-),1hi,50(1h-),1hi)
324   format(d12.4,d12.4,3h  i,2(d12.4,f8.2),i8,3h  i,4x,1hi,40a1,1hi)
325   format(26x,1hi,50x,1hi,4x,1hi,40a1,1hi)
326   format(2x,24(1h-),1hi,50(1h-),1hi,4x,1hi,40a1,1hi)
327   format(7x,15htotal underflow,4x,1hi,d12.4,d20.4,i16,2x,1hi)
328   format(7x,15htotal  overflow,4x,1hi,d12.4,d20.4,i16,2x,1hi)
329   format(//19x,21haccumulated average =,d12.5
     1/19x,21hestimated error     =,d12.5)
481   format(1h1,45x,40hdouble dIFferential cross-section number,i3/
     1/60x,7hx-axis ,3a4/60x,7hy-axis ,3a4/)
482   format(20x,11(1hi,9x))
483   format(11x,d9.3,11(1hi,d9.3))
484   format(1x,d10.3,  9h---------,11(10hi---------))
485   format(1h0,14x,11d10.3)
486   format(11x,i8,1x,11(1hi,i8,1x))
810   format(i2)
811   format(2d12.4,3i2,8a4)
812   format(2(2d10.3,i4),6a4)
813   format(12h1***error***,10x,24htoo many plots requested//
     122x,20hthe upper limits are //19x,i2,9h averages//19x,i2,22h one d
     2mensional plots//19x,i2,22h two DIMENSIONal plots////22x,25h***exe
     3cution is halted***   )
814   format(37h1number of single dIFferential cross
     1,20hsections requested =,i3/)
815   format(30h information on the data cards//
     13h  i,10x,5hxlmin,12x,5hxlmax,7x,
     224hbins  correllation  type,19x,4htext/)
816   format(i3,2e17.4,i8,i9,i10,5x,8a4)
817   format(37h0number of double dIFferential cross
     1,20hsections requested =,i3/)
818   format(30h information on the data cards//
     13h  i,8x,5hv1min,10x,5hv1max,4x,4hbins,8x,5hv2min
     2,10x,5hv2max,4x,4hbins,8x,6htext 1,8x,6htext 2/)
819   format(i3,2d15.3,i5,1x,2d15.3,i5,6x,3a4,2x,3a4)
820   format(31h0number of averages requested =,i3)
      END

*----------------------------------------------------------------------

      SUBROUTINE mapw2(w2,x,w2min,w2max,dw)
      IMPLICIT DOUBLE PRECISION (a-z)
*     wmin=1./w2min
*     wmax=1./w2max
*     dw=wmin-wmax
*     w2=1./(wmax+dw*x)
*     dw=dw*w2*w2
      y = w2max/w2min
      w2 = w2min*y**x
      dw = w2*dlog(y)
      RETURN
      END
      
      SUBROUTINE mapt1(t,x,tmin,tmax,dt)
      IMPLICIT DOUBLE PRECISION (a-z)
      y=tmax/tmin
      t=tmin*y**x
      dt=-t*dlog(y)
*     dt = tmin-tmax
*     t = tmin - x*dt
      RETURN
      END
      
      SUBROUTINE mapt2(t,x,tmin,tmax,dt)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      y=tmax/tmin
      t=tmin*y**x
      dt=-t*dlog(y)
*     dt = tmin-tmax
*     t = tmin - x*dt
      RETURN
      END
