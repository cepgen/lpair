      subroutine prtlhe(mode)

      implicit none

      integer mode
      integer i,j
* === Run common block
      integer MAXPUP
      parameter ( MAXPUP=100 )
      integer IDBMUP, PDFGUP, PDFSUP, IDWTUP, NPRUP, LPRUP
      double precision EBMUP, XSECUP, XERRUP, XMAXUP
      common/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     +     IDWTUP, NPRUP, XSECUP(MAXPUP), XERRUP(MAXPUP),
     +     XMAXUP(MAXPUP), LPRUP(MAXPUP)
      save /HEPRUP/
* === Event information
      integer MAXNUP
      parameter ( MAXNUP=500 )
      integer NUP, IDPRUP, IDUP, ISTUP, MOTHUP, ICOLUP
      double precision XWGTUP, SCALUP, AQEDUP, AQCDUP,
     +     PUP, VTIMUP, SPINUP
      common/HEPEUP/ NUP, IDPRUP, XWGTUP, SCALUP, AQEDUP, AQCDUP,
     +     IDUP(MAXNUP), ISTUP(MAXNUP), MOTHUP(2,MAXNUP),
     +     ICOLUP(2,MAXNUP), PUP(5,MAXNUP), VTIMUP(MAXNUP),
     +     SPINUP(MAXNUP)
      save /HEPEUP/
* === Kinematic information from JETSET
      integer n,k, npad
      real p,v
      common /pyjets/ N, npad, K(4000,5), P(4000,5), V(4000,5)
*
* === Input parameters
      integer      MODCUT
      real*4       THMAX,THMIN,MXMN,MXMX,Q2MN,Q2MX
      real*8       COTTH1,COTTH2,ECUT,PTCUTMIN,PTCUTMAX,MXMIN2,MXMAX2,
     +             QP2MIN,QP2MAX
      common/CUTS/    COTTH1,COTTH2,ECUT,PTCUTMIN,PTCUTMAX,
     +                MXMIN2,MXMAX2,
     +                THMAX,THMIN,QP2MIN,QP2MAX,MODCUT,MXMN,MXMX,
     +                Q2MN,Q2MX
      integer         NDIM,NCVG,ITMX,NPRN,IGRAPH,
     +                NPOIN,NPRIN,NTREAT,IBEG,IEND,NGEN
      common/VEGPAR/  NDIM,NCVG,ITMX,NPRN,IGRAPH,
     +                NPOIN,NPRIN,NTREAT,IBEG,IEND,NGEN
      integer         INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,IPAIR,NQUARK
      real*8          INPE,INPP
      common/BEAM/    INPE,INPP,INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,
     +                IPAIR,NQUARK
* === Configuration options
      integer lhe
      data lhe/10/
      logical init
      data init/.false./
*
*
      if (init.eqv..false.) then
         open(lhe,file='events.lhe',status='unknown')
         init = .true.
      endif
*
      if (mode.eq.1) then       ! header
         write(lhe,'(A)') '<LesHouchesEvents version="1.0">'
         write(lhe,'(A)') '<!--'
         write(lhe,'(A)') 'File generated with LPAIR'
         call prtpar(lhe)
         write(lhe,'(A)') '-->'
      elseif (mode.eq.2) then   ! event
         write(lhe,'(A)') '<event>'
         do 10,i=1,n
c            print *,p(j,5)
            if (VTIMUP(I).eq.0D0) then
               write(lhe,5300) IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),
     &              ICOLUP(1,I),ICOLUP(2,I),(PUP(j,i),J=1,5),
     &              ' 0. 9.'
            else
               write(lhe,5400) IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),
     &              ICOLUP(1,I),ICOLUP(2,I),(PUP(j,i),J=1,5),
     &              VTIMUP(I),' 9.'
            endif
 10      continue
         write(lhe,'(A)') '</event>'
      elseif (mode.eq.3) then   ! footer
         write(lhe,'(A)') '</LesHouchesEvents>'
      else
         print *, 'PrintLHE: Error: Unrecognized mode:', mode
      endif
*
 5300 format(1p,i8,5i5,5e18.10,a6)
 5400 format(1p,i8,5i5,5e18.10,e12.4,a3)
*
      end
