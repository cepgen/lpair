      subroutine prtpar(file)

      implicit none
*
* === Input parameters
      integer      MODCUT
      real*4       THMAX,THMIN,MXMN,MXMX,Q2MN,Q2MX
      real*8       COTTH1,COTTH2,ECUT,PTCUT,MXMIN2,MXMAX2,QP2MIN,QP2MAX
      common/CUTS/    COTTH1,COTTH2,ECUT,PTCUT,MXMIN2,MXMAX2,
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
      integer file
*
      write(file,'(A)') 'Input parameters for this generation:'
      write(file,1002)
      write(file,1000) 'NTRT',ntreat
      write(file,1000) 'NCVG',ncvg
      write(file,1000) 'ITVG',itmx
      write(file,1000) 'NCSG',npoin
      write(file,1001) 'INPP',inpp
      write(file,1000) 'PMOD',pmod
      write(file,1001) 'INPE',inpe
      write(file,1000) 'EMOD',emod
      write(file,1000) 'GPDF',gpdf
      write(file,1000) 'SPDF',spdf
      write(file,1000) 'QPDF',nquark
      write(file,1001) 'THMN',thmin
      write(file,1001) 'THMX',thmax
      write(file,1001) 'Q2MN',q2mn
      write(file,1001) 'Q2MX',q2mx
      write(file,1001) 'MXMN',mxmn
      write(file,1001) 'MXMX',mxmx
      write(file,1000) 'PAIR',ipair
      write(file,1000) 'MCUT',modcut
      write(file,1001) 'ECUT',ecut
      write(file,1001) 'PTCT',ptcut
      write(file,1001) 'NGEN',ngen
      write(file,1002)
*
 1000 format((a),1x,i10)
 1001 format((a),1x,f10.2)
 1002 format(15('='))
*
      end
