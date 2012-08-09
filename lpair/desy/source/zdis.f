*-- Author :
        Program ZDIS
* ==================
*
*    Main steering programme for the Zeus DIS simulation programme
*
        Implicit NONE
*       =============
*
*KEEP,PARTAP.
C
C        system-defined parameters for the Table Package Programmer
C
      INTEGER MINC,MAXC,NEXT,INULL,IANY,INS,REP,ORD,UNO,AND,OR,DIF,HOR,
     +        VER,ALL,ID,ALLCOL,
     +        COUTAB,COUSEL,GETIND,GETSEL,GETPRO,GETDFL,GETTDF,SPATAB,
     +        CHKTAP,MAKTAB

      REAL    RNULL,RANY

      CHARACTER*4 CNULL,CANY

      LOGICAL BELSEL,BELTAB,CHKREL,CHKTAB,CHKWIN
C
C       MINC and MAXC for cursors operation
C
      PARAMETER (MINC=1, MAXC=2147483647)
C
C       null values
C
      PARAMETER (INULL=2147483647, RNULL= 699050*16.0**26, CNULL='====')
C
C       ANY values
C
      PARAMETER (IANY=-INULL, RANY=-RNULL, CANY='!@)(')
C
C       NEXT for insertion
C
      PARAMETER (NEXT=INULL)
C
C       Modes
C
      PARAMETER (INS = 1, REP = 2,
     +           ORD = 1, UNO = 2,
     +           AND = 1, OR  = 2, DIF=3,
     +           HOR = 1, VER = 2)
C
C       Indices
C
      PARAMETER (ALL = 1-INULL , ID = INULL-1, ALLCOL = 1-INULL)
*KEEP,ZDSKEY.
      INTEGER ZDSKEY,ZDSKEY_9999
      INTEGER ZDSKEY_ID,ZDSKEY_Nr1,ZDSKEY_Nr2,ZDSKEY_TStam11,
     +   ZDSKEY_TStam12,ZDSKEY_TStam21,ZDSKEY_TStam22
      CHARACTER*4 ZDSKEY_GAFTyp
      CHARACTER*32 ZDSKEY_DflNAM

      COMMON/ZDSKEY/ZDSKEY,ZDSKEY_ID,ZDSKEY_GAFTyp,ZDSKEY_Nr1,
     +   ZDSKEY_Nr2,ZDSKEY_TStam11,ZDSKEY_TStam12,ZDSKEY_TStam21,
     +   ZDSKEY_TStam22,ZDSKEY_DflNAM,ZDSKEY_9999

*KEEP,GCUNIT.
      COMMON/GCUNIT/LIN,LOUT,NUNITS,LUNITS(5)
      INTEGER LIN,LOUT,NUNITS,LUNITS
      COMMON/GCMAIL/CHMAIL
      CHARACTER*132 CHMAIL
C
*KEEP,OUTSTR.
C -------------- Output steering switsches -------------------------------
C
      COMMON/ZOUTST/LOUSWT,LOUDFL,LOURCD,LOUNET,SKIPEV
C
C  LOUSWT :== Switch ON/OFF o/p reading
C  LOUDFL :== Switch ON/OFF o/p dataflow selection
C  LOURCD :== Switch ON/OFF o/p GAFTYP selection
C  OUMAX  :== Maximum # of output files to be opened in parallel.
C  LOUNET :== Switch ON/OFF o/p reading via network ( for Number Cruncher )
C  SKIPEV :== Switch T/F skipping of event writing ( always clear dfls )
C
      Integer OUMAX
      PARAMETER ( OUMAX = 10 )
C
      LOGICAL LOUSWT,LOUDFL,LOURCD,LOUNET,SKIPEV
C -------------------------------------------------------------------------
*KEEP,ZDSEVT.
*
*..  Event parameters for the DIS generator interface
*
      INTEGER NEVT,RUNNO,GRNDMQSEED
C
      COMMON /ZDSEVT/ NEVT,RUNNO,GRNDMQSEED(2)
*
*KEEP,ZRLIMI.
C
C -- Requested Time / Event Limit / To be checked in ZINREC
C
      COMMON /ZRLIMI/ NEVT_TOT, Time_Out
C
      Integer NEVT_TOT
      Real    Time_Out
C
*KEEP,DISSTE.
*
*..  Common to hold steering parameter(s) for the DIS generator interface,
*        DISGenerator ... Name of generator user wants to run (No default).
*
        Character*8   DISGenerator, Institute
        Integer  NDISList ,IFList
        Parameter (NDISList=1)
*
      COMMON /DISSTE/ DISGenerator, Institute, IFList(NDISList)
*
*KEEP,ZREVT.
      INTEGER ZREVT,ZREVT_9999
      INTEGER ZREVT_ID,ZREVT_RunNr,ZREVT_EvtNr,ZREVT_Time,ZREVT_TrgMsk,
     +   ZREVT_SelMsk

      COMMON/ZREVT/ZREVT,ZREVT_ID,ZREVT_RunNr,ZREVT_EvtNr(3),
     +   ZREVT_Time(2),ZREVT_TrgMsk(3),ZREVT_SelMsk(3),ZREVT_9999

*KEEP,ZRRUN.
      INTEGER ZRRUN,ZRRUN_9999
      INTEGER ZRRUN_ID,ZRRUN_TExp,ZRRUN_RunNr,ZRRUN_TRun,ZRRUN_BoR,
     +   ZRRUN_EoR,ZRRUN_NEvtR,ZRRUN_NEvtA
      CHARACTER*16 ZRRUN_OpMess

      COMMON/ZRRUN/ZRRUN,ZRRUN_ID,ZRRUN_TExp,ZRRUN_RunNr,ZRRUN_TRun,
     +   ZRRUN_BoR(2),ZRRUN_EoR(2),ZRRUN_OpMess(5),ZRRUN_NEvtR,
     +   ZRRUN_NEvtA,ZRRUN_9999

*KEEP,MOSTCK.
      INTEGER mostck,mostck_9999
      INTEGER mostck_ID,mostck_cdepth,mostck_module
      CHARACTER*8 mostck_modid,mostck_caldby
      REAL mostck_sttime
      LOGICAL mostck_actflg

      COMMON/mostck/mostck,mostck_ID,mostck_modid,mostck_caldby,
     +   mostck_cdepth,mostck_actflg,mostck_sttime,mostck_module,
     +   mostck_9999

*KEEP,MOCOND.
      INTEGER mocond,mocond_9999
      INTEGER mocond_ID,mocond_sever,mocond_mostck
      CHARACTER*32 mocond_cnd
      REAL mocond_weight

      COMMON/mocond/mocond,mocond_ID,mocond_cnd,mocond_weight,
     +   mocond_sever,mocond_mostck,mocond_9999

*KEND.
*
*
        Integer IFL,Ierr, Nevt_Strt, Want, ErrMax, MoAbPr,
     + End_Of_Run
*
      INTEGER IWANT
        Real TLeft
*
        Character*4 RECTYP
*
        External MoAbPr
*
*--------  Initialisation.
*
*.. ZAINIT performs the initialisation of the programme run.
*
        Call DISIni( Ierr)
        If (Ierr.ne.0) Then
            Write(CHMAIL,9020) Ierr
            Call GMAIL(4,0)
          GOTO 7900
        EndIf
        ErrMax = MoAbPr(0)
      IF(ERRMAX.EQ.3) THEN
        WRITE(6,'(A,I3)') ' **** ERROR level from MOMO =',errmax
        CALL DTAB(MOSTCK)
        CALL DTAB(MOCOND)
        goto 7900
      endif
*
*--------  Open event LOOP  --------  --------  --------
*
        Nevt_Strt = 0
10      Continue
*
*..  Clear the union dataflow:
        Call ZRCLEA
*
*.. Count:
        Nevt_Strt = Nevt_Strt + 1
        If (NEVT_TOT.NE.-1) Then
            If (NEvt.GE.NEVT_TOT) Then
                End_Of_Run = 1
                 GOTO 7000
            EndIf
        EndIf
        Call TIMEL(TLeft)
        If (TLeft.LT.Time_Out) Then
            End_Of_Run = 2
                 GOTO 7100
        EndIf
        CALL GRNDMQ(GRNDMQSEED(1), GRNDMQSEED(2), 0,'G')
        if (mod(nevt_strt,100).eq.0) then
          WRITE(CHMAIL, 9040) NEVT,NEVT_STRT,GRNDMQSEED(1),GRNDMQSEED(2)
          Call GMail(0,0)
        endif
*
*-------- Process this event:
*
*.. Call the event generation administration routine:
*
        Call DISGEN(IErr)
C
      IF(IERR.GT.1.AND.IERR.LT.10) THEN
C
C ---   DO NOT WRITE OUT "ILL" EVENTS
          GOTO 6900
      ELSEIF(IERR.GT.10) THEN
C --- A FATAL ERROR OCCURED : TERMINTATE RUN
         WRITE(CHMAIL,'(1X,A)')
     + ' *** FATAL ERROR IN GENERATOR : RUN MUST BE STOPPED ***'
C
         GOTO 7900
C
      ENDIF
        ErrMax = MoAbPr(0)
        If (ErrMax.Eq.3) Goto 7900
        If (ErrMax.Eq.2) Goto 6900
*
*--------  Call the user event loop routine
*
C
C ---  users selection flag ( default IWANT = 1 i.e. take it )
      IWANT = 1
C
      CALL ZDUEVT(IWANT)
C
C --- enable user's selection
C
      IF(IWANT .EQ.0) GOTO 6900
C
*
* ----  finally count the accepted events
*
      NEVT = NEVT + 1
*
*--------  Final stage: => Select the event and write it out.
*
*.. Write them out:
*
5000    Continue
        IF (LOUswt) then
          ZDSKEY_GAFTyp    = 'EVT '
          ZDSKEY_NR1= RUNNO
          ZDSKEY_Nr2       = NEVT
          ZDSKEY_TStam11 = 0
          ZDSKEY_TStam12 = 0
          ZDSKEY_TStam21 = 0
          ZDSKEY_TStam22 = 0
*
C --- fill ZREVT
         CALL NULWIN(ZREVT)
         IF(RUNNO.EQ.0) THEN
           ZREVT_RUNNR = 1
         ELSE
           ZREVT_RUNNR = RUNNO
         ENDIF
         ZREVT_EVTNR(1) = NEVT
         ZREVT_EVTNR(2) = NEVT
         ZREVT_EVTNR(3) = NEVT
      ZREVT_ID = NEXT
         CALL INSTAB(ZREVT)
C
          Call ZOUREC(Ierr)
          ErrMax = MoAbPr(0)
          If (ErrMax.Eq.3) Goto 7900
          If (ErrMax.Eq.2) Goto 6900
        END if
*
*--------  Take next event
*
6900    Continue
*.. Reset MOmo tables
        Call MOmors
*
        Goto 10
*
*--------  End of event loop
*
7000    Continue
*
*.. Event limit reached
*
        Write (CHMAIL, 9050)
        Call GMAIL(0,0)
        Goto 7900
*
*.. Out of time
*
7100    Continue
        Write (CHMAIL, 9060)
        Call GMAIL(0,0)
        Goto 7900
*
*.. Eof
*
7200    Continue
        Write (CHMAIL, 9070)
        Call GMAIL(0,0)
        Goto 7900
*
*--------  End of programme run processing:
*
7900    Continue
        Call DISEnd(End_Of_Run, IErr)
        If (Ierr.Ne.0) Goto 8030
        Goto 8900
*
*--------  Error codes:
*
8030    Continue
        Write (CHMAIL, 9030) Ierr
        Call GMAIL(4,0)
        Goto 8900
*
*
8900    Continue
*
9020    Format(
     + 8x,'ZDIS. Error from DISIni code =',i4)
*
9030    Format(
     + 8x,'ZDIS. Error from DISEND code =',i4,'  programme stops.')
*
9040    Format(
     + 1X,' ZDIS. Start event no. trial=',I7,'  so far output :',I7
     + ,' Random number seeds:',2(1x,I12))
*
9050    Format(
     + 8x,'ZDIS. Event limit reached.')
*
9060    Format(
     + 8x,'ZDIS. Time limit reached.')
*
9070    Format(
     + 8x,'ZDIS. End of file reached.')
*
        End
