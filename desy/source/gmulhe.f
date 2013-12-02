      subroutine gmulhe

      implicit none

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
      double precision p,v
      common /pyjets/ N, npad, K(4000,5), P(4000,5), V(4000,5)
* === 
      integer MSTP, MSTI
      double precision PARP, PARI
      common/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      save /PYPARS/
*
c      integer nstable
*
*
c      nstable = 0
*
* === Event filling
      NUP = N                   ! number of particle entries in the event
      IDPRUP = 1                ! ID of the process for the event
      XWGTUP = 1.               ! event weight
      SCALUP = -1.              ! scale of the event in GeV, as used for calculation of PDFs
      AQEDUP = -1.              ! the QED coupling used for this event
      AQCDUP = -1.              ! the QCD coupling used for this event
*
      do 10 i=1,N
         IDUP(i) = K(i,2)       ! particle ID according to Particle Data Group convention
         ISTUP(i) = 0           ! status code
         MOTHUP(1,i) = 0        ! index of first and last mother
         MOTHUP(2,i) = 0
         ICOLUP(1,I) = 0        ! integer tag for the color flow line passing through the color of the
                                ! particle
         ICOLUP(2,I) = 0        ! integer tag for the color flow line passing through the anti-color of
                                ! the particle
         do 11 j=1,5
            PUP(j,i) = P(i,j)   ! lab frame momentum (Px, Py, Pz, E, M) of particle in GeV
 11      continue
         
         VTIMUP(i) = 0.         ! invariant lifetime c*tau (distance from production to decay) in mm
         SPINUP(i) = 0.         ! cosine of the angle between the spin-vector of particle I and the 3-
                                ! momentum of the decaying particle, specified in the lab frame
 10   continue
c      print *,nstable
*
c      call pylhef
*
      end
