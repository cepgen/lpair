      subroutine W1W2f2(t1,MX,W1strfun,W2strfun)
c     
c      ================================================
c      program calculates W1 and W2 structure functions
c      from the GRV95 LO parametrization
c      ================================================
c
c      standard parameters
c

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       COMMON /PICKZZ/ W1,W2,W3,W4,W5,W31,W52,W12,TAU,SL1

       am_p = 0.93827203
       Q02 = 0.8

       Q2 = -t1

       x = Q2/(W3+Q2+am_p*am_p)

       amu2 = Q2+Q02       ! scale is shifted
C       write(*,*) W3, x, amu2

       call grv95lo(x,amu2,xuv,xdv,xus,xds,xss,xg)

       F2_aux = 4./9.*(xuv+2.*xus)
     2       + 1./9.*(xdv+2.*xds)
     3       + 1./9.*2.*xss

c
c      F2 corrected for low Q^2 behaviour
c
       F2_corr = Q2/(Q2+Q02)*F2_aux

       F1 = F2_corr/(2.*x)              ! Callan-Gross relation

C       write(*,*) F1, F2_corr

       W2strfun = 2.*am_p*x/Q2*F2_corr
       W1strfun = F1/am_p

       end
