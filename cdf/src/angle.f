 
      FUNCTION MYULANGL(X,Y) 
C...Purpose: to reconstruct an angle from given x and y coordinates. 

      REAL ULANGL,R,PI
      PI = ACOS(-1.)
 
      ULANGL=0. 
      R=SQRT(X**2+Y**2) 
      IF(R.LT.1E-20) RETURN 
      IF(ABS(X)/R.LT.0.8) THEN 
        ULANGL=SIGN(ACOS(X/R),Y) 
      ELSE 
        ULANGL=ASIN(Y/R) 
        IF(X.LT.0..AND.ULANGL.GE.0.) THEN 
          ULANGL=PI-ULANGL 
        ELSEIF(X.LT.0.) THEN 
          ULANGL=-PI-ULANGL 
        ENDIF 
      ENDIF 
 
      RETURN 
      END 
