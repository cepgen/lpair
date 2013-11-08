       DOUBLE PRECISION FUNCTION ranf(dummy)
       DOUBLE PRECISION rvec(100)
       SAVE rvec,ncall,mcall
       INTEGER dummy
       DATA ncall,mcall/0,0/
       IF ( ncall .EQ. 0 ) THEN
           IF ( dummy .NE. 0 ) CALL rcargo(dummy)
           ncall = 1
       ENDIF
       IF ( mcall .EQ. 0 ) THEN
           CALL rcarry(rvec,100)
           mcall = 100
       ENDIF
       ranf = rvec(mcall)
       mcall = mcall - 1
       END

      SUBROUTINE rcarry(rvec,lenv)
*
*     This version is identical to that in CPC software library
*
*         Add-and-carry random number generator proposed by
*         Marsaglia and Zaman in SIAM J. Scientific and Statistical
*             Computing, to appear probably 1990.
*         modified with enhanced initialization by F. James, 1990
*
*
*  NB! Recently, F. James informed us that there is a slight mistake
*      in this implementation on line 13 and suggested the following
*      change:
*
*      Original line:
*
*      uni = seeds(i24) - seeds(j24) - carry
*
*      Suggested modification:
*
*      uni = seeds(j24) - seeds(i24) - carry
*
*      We have also tested this new versionand the test results
*      are available in the hep-lat preprint 9306008.
*
*!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*!!!  Calling sequences for rcarry:                                  ++
*!!!      CALL rcarry (rvec, LEN)   returns a vector rvec of LEN     ++
*!!!                   32-bit random floating point numbers between  ++
*!!!                   zero and one.                                 ++
*!!!      CALL rcargo(INT)     initializes the generator from one    ++
*!!!                   32-bit integer INT                            ++
*!!!      CALL RCARIN(ivec)    restarts the generator from vector    ++
*!!!                   ivec of 25 32-bit integers (see rcarut)       ++
*!!!      CALL rcarut(ivec)    outputs the current values of the 25  ++
*!!!                 32-bit integer seeds, to be used for restarting ++
*!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION rvec(lenv)
      DIMENSION seeds(24), iseeds(24), isdext(25)
      PARAMETER (TWOP12=4096.D0)
      PARAMETER (ITWO24=2**24, icons=2147483563)
      SAVE notyet, i24, j24, carry, seeds, twom24
      LOGICAL notyet
      DATA notyet/.TRUE./
      DATA i24,j24,carry/24,10,0.D0/
C
C              Default Initialization by Multiplicative Congruential
      IF (notyet) THEN
         notyet = .FALSE.
         jseed = 314159265
*         WRITE(6,'(A,I12)') ' rcarry default initialization: ',jseed
            twom24 = 1
         DO 25 i= 1, 24
            twom24 = twom24 * 0.5
         k = jseed/53668
         jseed = 40014*(jseed-k*53668) -k*12211
         IF ( jseed .LT. 0 ) jseed = jseed + icons
         iseeds(i) = MOD(jseed,ITWO24)
   25    CONTINUE
         DO 50 i= 1,24
         seeds(i) = iseeds(i)*twom24
   50    CONTINUE
         i24 = 24
         j24 = 10
         carry = 0.
         IF (seeds(24) .LT. seeds(14)) carry = twom24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 ivec= 1, lenv
      uni = seeds(j24) - seeds(i24) - carry
      IF (uni .LT. 0.)  THEN
         uni = uni + 1.0
         carry = twom24
      ELSE
         carry = 0.
      ENDIF
      seeds(i24) = uni
      i24 = i24 - 1
      IF ( i24 .EQ. 0 )  i24 = 24
      j24 = j24 - 1
      IF (j24 .EQ. 0)  j24 = 24
      rvec(ivec) = uni
  100 CONTINUE
      RETURN
C           Entry to input and float integer seeds from previous run
      ENTRY RCARIN(isdext)
         twom24 = 1
         DO 195 i= 1, 24
  195    twom24 = twom24 * 0.5D0
      WRITE(6,'(A)') ' FULL INITIALIZATION OF rcarry WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') isdext
      DO 200 i= 1, 24
      seeds(i) = isdext(i)*twom24
  200 CONTINUE
      carry = MOD(isdext(25),10)*twom24
      isd = isdext(25)/10
      i24 = MOD(isd,100)
      isd = isd/100
      j24 = isd
      RETURN
C                    Entry to ouput seeds as integers
      ENTRY rcarut(isdext)
      DO 300 i= 1, 24
         isdext(i) = seeds(i)*TWOP12*TWOP12
  300 CONTINUE
      icarry = 0
      IF (carry .GT. 0.)  icarry = 1
      isdext(25) = 1000*j24 + 10*i24 + icarry
      RETURN
C                    Entry to initialize from one integer
      ENTRY rcargo(inseed)
      jseed = inseed
      notyet = .FALSE.
*      WRITE(6,'(A,I12)') ' rcarry INITIALIZED FROM SEED ',inseed
      twom24 = 1
         DO 325 i= 1, 24
           twom24 = twom24 * 0.5D0
         k = jseed/53668
         jseed = 40014*(jseed-k*53668) -k*12211
         IF ( jseed .LT. 0 ) jseed = jseed+icons
         iseeds(i) = MOD(jseed,ITWO24)
  325    CONTINUE
         DO 350 i= 1,24
         seeds(i) = iseeds(i)*twom24
  350    CONTINUE
         i24 = 24
         j24 = 10
         carry = 0.
         IF ( seeds(24) .LT. seeds(14) ) carry = twom24
      RETURN
*###] rcarry:
      END

