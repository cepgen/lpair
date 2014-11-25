C...PYANGL
C...Reconstructs an angle from given x and y coordinates.
 
      FUNCTION PYANGL(X,Y)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      PI = ACOS(-1.0D0) 
      PYANGL=0D0
      R=SQRT(X**2+Y**2)
      IF(R.LT.1D-20) RETURN
      IF(ABS(X)/R.LT.0.8D0) THEN
        PYANGL=SIGN(ACOS(X/R),Y)
      ELSE
        PYANGL=ASIN(Y/R)
        IF(X.LT.0D0.AND.PYANGL.GE.0D0) THEN
          PYANGL=PI-PYANGL
        ELSEIF(X.LT.0D0) THEN
          PYANGL=-PI-PYANGL
        ENDIF
      ENDIF
 
      RETURN
      END
