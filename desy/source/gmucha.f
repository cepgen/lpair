*-- Author :    ZEUS Offline Group   18/08/94

      SUBROUTINE GMUCHA

***********************************************************************
*
*    SUBROUTINE GMUCHA
*
* PURPOSE: Interpret the data cards to change the default parameters of MUPAIR.
*
* INPUT: Bos text bank GMUP
*
* OUTPUT: updated parameters in MUPAIR commons.
*
* CALLED BY: GMUPA
*
* AUTHOR: OLAF DUENGER             CREATED AT: 91/12/12
*
* CHANGED BY:                              AT:
* REASON:
*
************************************************************************
*
      implicit NONE
*
*KEEP,VEGPAR.
      INTEGER          NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND,NGEN
      COMMON /VEGPAR/  NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND,NGEN

*KEEP,BEAM.
      INTEGER          INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,IPAIR,NQUARK
      REAL*8           INPE,INPP
      COMMON /BEAM/    INPE,INPP,INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,
     &                 IPAIR,NQUARK

*KEEP,CUTS.
      INTEGER      MODCUT
      REAL*4       THMAX,THMIN,MXMN,MXMX,Q2MN,Q2MX
      REAL*8       COTTH1,COTTH2,ECUT,PTCUTMIN,PTCUTMAX,MXMIN2,MXMAX2,
     &             QP2MIN,QP2MAX
      COMMON /CUTS/COTTH1,COTTH2,ECUT,PTCUTMIN,PTCUTMAX,MXMIN2,MXMAX2,
     &             THMAX,THMIN,QP2MIN,QP2MAX,MODCUT,MXMN,MXMX,Q2MN,Q2MX

*KEND.
*
C*  End of common
*
c      intrinsic iargc,getarg
c      external iargc,getarg
c      integer iargc
      integer numarguments
      external numarguments
      integer i,lun,maxln
      character(len=32) file
      character(len=6) key
      double precision value
      logical fexst
*
*--------  Read data cards and overwrite defaults:
*
      lun=15
      maxln=20

c      print *,'==>',numarguments()
c      print *,'-->',iargc()
c      if (iargc().gt.0) then
c         call getarg(1,file)
c      else
         file='lpair.card'
c      endif

*---- Make sure the file exists
      inquire(file=file,exist=fexst)
      if (fexst.eqv..false.) then
         print *,'GMUCHA: ERROR! Input card does not exist!'
         stop
      endif

*---- Read the parameter card using key/value pairs
*
      open(lun,file=file,status='old')
      do i=1,maxln
         read(lun,1000,end=10) key,value
         if (trim(key).eq."IBEG") ibeg=value
         if (trim(key).eq."IEND") iend=value
         if (trim(key).eq."NGEN") ngen=value
         if (trim(key).eq."NTRT") ntreat=value
         if (trim(key).eq."PRVG") nprin=value
         if (trim(key).eq."NCVG") ncvg=value
         if (trim(key).eq."ITVG") itmx=value
         if (trim(key).eq."NCSG") npoin=value
         if (trim(key).eq."INPP") inpp=value
         if (trim(key).eq."PMOD") pmod=value
         if (trim(key).eq."GPDF") gpdf=value
         if (trim(key).eq."SPDF") spdf=value
         if (trim(key).eq."INPE") inpe=value
         if (trim(key).eq."EMOD") emod=value
         if (trim(key).eq."PAIR") ipair=value
         if (trim(key).eq."QPDF") nquark=value
         if (trim(key).eq."MCUT") modcut=value
         if (trim(key).eq."THMX") thmax=value
         if (trim(key).eq."THMN") thmin=value
         if (trim(key).eq."ECUT") ecut=value
         if (trim(key).eq."PTCT") ptcutmin=value
         if (trim(key).eq."PTMX") ptcutmax=value
         if (trim(key).eq."Q2MN") q2mn=value
         if (trim(key).eq."Q2MX") q2mx=value
         if (trim(key).eq."MXMN") mxmn=value
         if (trim(key).eq."MXMX") mxmx=value
      enddo
 10   continue
      close(lun)
*
*--------  Return
*
 1000 format(a4,d9.0)
      return
      end
