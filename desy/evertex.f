*****************************************************************
* 		EVENT VERTEX ver. 1.0 (CERN)			*
* Subroutine for a simulation of the vertex and beam divergence	*
* Created 22 June 2005, last update 29 June 2005              	*
* Author: Dariusz Bocian@cern.ch                              	*
*								*
*****************************************************************
*
      SUBROUTINE EVERTEX(phiv)
*
      REAL phiv
*
*- global parameters
*
      REAL PI,SQ2
*
      parameter (PI=3.14159265)
      parameter (SQ2=1.4142135624)
*
*- global common
*
      INTEGER iv
      REAL vx,vy,vz,dthetax,dthetay
      COMMON /PVERTEX/ iv,vx,vy,vz,dthetax,dthetay ! vx-z - production vertex coordinates
                                                   ! dthetax, dthetay - vertex smearing due to the beam divergence
*
*- normal distribution random number generation
*
      DIMENSION RVEC(10)
      REAL x_gen,y_gen,z_gen
      REAL sigma_x,sigma_y,sigma_z,sigma_alfa
      REAL x_av,y_av,z_av
      REAL dtheta_x,dtheta_y
*
      parameter (sigma_x=15.9E-6)
      parameter (sigma_y=15.9E-6)
      parameter (sigma_z=75.0E-3)
      parameter (sigma_alfa=16.0E-6)
*
*--------------------------------------------------------------
*
*- Normal distribution random number generation
*      
      LEN=1
      DO J=1,1
*
*- coordinates of IP(e+e-) (x,y,z) and its modifications
*
       CALL RNORML (RVEC,LEN)
       x_gen=RVEC(J)
       x_gen=x_gen*(sigma_x/SQ2)
       call hf1(101,x_gen,1.)

       CALL RNORML (RVEC,LEN)
       y_gen=RVEC(J)
       y_gen=y_gen*(sigma_y/SQ2)
       call hf1(102,y_gen,1.)

       CALL RNORML (RVEC,LEN)
       z_gen=RVEC(J)
       z_gen=z_gen*(sigma_z/SQ2)
       call hf1(103,z_gen,1.)
*
*- modification due to the beam divergence
* 
       CALL RNORML (RVEC,LEN)
       dtheta_x=RVEC(J)
       dtheta_x=dtheta_x*16.E-6*cos(phiv)
       call hf1(105,dtheta_x,1.)
*
       CALL RNORML (RVEC,LEN)
       dtheta_y=RVEC(J)
       dtheta_y=dtheta_y*16.E-6*sin(phiv)
       call hf1(106,dtheta_y,1.)
*
      ENDDO

      vx=x_gen
      vy=y_gen
      vz=z_gen
      dthetax=dtheta_x
      dthetay=dtheta_y

      END
