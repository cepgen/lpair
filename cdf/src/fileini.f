      SUBROUTINE fileini
*
      IMPLICIT none
*
* --- define data common
*
      INTEGER IPAR(20)
      REAL*8  LPAR(20)
      COMMON/DATAPAR/IPAR,LPAR
*
      INTEGER LUN
      CHARACTER(len=32) file
      LOGICAL fexst
*
      LUN = 15
      CALL LUGIVE("MSTU(21)=1")
*
* LF workaround to pass the input card as an argument
      IF (iargc().gt.0) then
         CALL getarg(1,file)
      ELSE
         file='lpair.dat'
      ENDIF
      INQUIRE(file=file,exist=fexst) ! ensure the input card exists
      IF (fexst.eqv..false.) THEN
         PRINT *,'FILEINI: ERROR! Input card does not exist!'
         STOP
      ENDIF
*
      OPEN(LUN,file=file,status='old')
*
      READ(LUN,'(I10)') IPAR(1)  ! event or cross-section calculation
      READ(LUN,'(I10)') IPAR(2)  ! paw histograms/ntuple
      READ(LUN,'(I10)') IPAR(3)  ! gnz histograms/ntuple
      READ(LUN,'(I10)') IPAR(4)  ! read 'ineemm.data' for graph
      READ(LUN,'(I10)') IPAR(5)  ! VEGAS dimensions
      READ(LUN,'(I10)') IPAR(6)  ! VEGAS points per iteration
      READ(LUN,'(I10)') IPAR(7)  ! VEGAS iterations
      READ(LUN,'(I10)') IPAR(8)  ! incoming particle type
      READ(LUN,'(I10)') IPAR(9)  ! vermaseren or caron cuts
      READ(LUN,'(I10)') IPAR(10) ! vermaseren form factor selection
      READ(LUN,'(I10)') IPAR(11) ! max # of events to GNZ output
      READ(LUN,'(I10)') IPAR(12) ! GENERA: # of events to generate
      READ(LUN,'(I10)') IPAR(13) ! TREAT: non-smooth or smooth F
      READ(LUN,'(I10)') IPAR(14) ! SETGEN: points per bin
      READ(LUN,'(I10)') IPAR(15) ! SETGEN: nbin
      READ(LUN,'(I10)') IPAR(16) ! vector rotation flag
      READ(LUN,'(I10)') IPAR(17) ! hepevt common block to ascii output
*
      READ(LUN,'(D10.4)') LPAR(1)  ! incoming particle mass (GeV)
      READ(LUN,'(D10.4)') LPAR(2)  ! outgoing particle mass (GeV)
      READ(LUN,'(D10.4)') LPAR(3)  ! beam energy (GeV)
      READ(LUN,'(D10.4)') LPAR(4)  ! maximum angle (rad)
      READ(LUN,'(D10.4)') LPAR(5)  ! maximum rapidity
      READ(LUN,'(D10.4)') LPAR(6)  ! energy cut (GeV)
      READ(LUN,'(D10.4)') LPAR(7)  ! minimum P_t
      READ(LUN,'(D10.4)') LPAR(8)  ! minimum invariant mass
      READ(LUN,'(D10.4)') LPAR(9)  ! maximum invariant mass
      READ(LUN,'(D10.4)') LPAR(10) ! maximum MX
*
      close(LUN)
*
      PRINT*,'***** INPUT LPAIR PARAMETERS *****'
      PRINT*,'IPAR(1):  event generation  ',IPAR(1)
      PRINT*,'IPAR(2):  paw ntuple        ',IPAR(2)
      PRINT*,'IPAR(3):  gnz output        ',IPAR(3)
      PRINT*,'IPAR(4):  lpair plot        ',IPAR(4)
      PRINT*,'IPAR(5):  VEGAS dim         ',IPAR(5)
      PRINT*,'IPAR(6):  VEGAS points      ',IPAR(6)
      PRINT*,'IPAR(7):  VEGAS iterations  ',IPAR(7)
      PRINT*,'IPAR(8):  beam particle     ',IPAR(8)
      PRINT*,'IPAR(9):  cuts              ',IPAR(9)
      PRINT*,'IPAR(10): form factor       ',IPAR(10)
      PRINT*,'IPAR(11): max GNZ events    ',IPAR(11)
      PRINT*,'IPAR(12): GENERA events     ',IPAR(12)
      PRINT*,'IPAR(13): TREAT smooth      ',IPAR(13)
      PRINT*,'IPAR(14): SETGEN points/bin ',IPAR(14)
      PRINT*,'IPAR(15): SETGEN nbin       ',IPAR(15)
      PRINT*,'IPAR(16): reflection flag   ',IPAR(16)
      PRINT*,'IPAR(17): hepevt to ascii   ',IPAR(17)
*
      PRINT*,'LPAR(1):  incoming mass  ',LPAR(1)
      PRINT*,'LPAR(2):  outgoing mass  ',LPAR(2)
      PRINT*,'LPAR(3):  beam energy    ',LPAR(3)
      PRINT*,'LPAR(4):  max angle      ',LPAR(4)
      PRINT*,'LPAR(5):  max rapidity   ',LPAR(5)
      PRINT*,'LPAR(6):  energy cut     ',LPAR(6)
      PRINT*,'LPAR(7):  min P_t        ',LPAR(7)
      PRINT*,'LPAR(8):  min invmass    ',LPAR(8)
      PRINT*,'LPAR(9):  max invmass    ',LPAR(9)
      PRINT*,'LPAR(10): max dissoc.mass',LPAR(10)
*
      RETURN
      END

