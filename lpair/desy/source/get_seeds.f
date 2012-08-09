      subroutine get_seeds(iprint)

      implicit none

* get geant random number generator seeds

      integer is1,is2,iseq,iprint
      logical exists
*
* init grndm
* 
      inquire(file='./grndm.seeds',exist=exists)
      if(exists) then
       open(98,file='./grndm.seeds',form='formatted')
       read(98,*) is1,is2
       close(98)
       if(iprint.gt.0) then
        write(*,'(a)') ' ************'
        write(*,'(a)') ' ************'
        write(*,'(a)') 
     +     ' Geant random number seeds read from grndm.seeds'
        write(*,'(a)') ' ************'
        write(*,'(a)') ' ************'
       endif
       CALL GRNDMQ(IS1,IS2,1,' ')
       IF(iS1.GT.0.AND.IS2.GT.0) THEN
        CALL GRNDMQ(IS1,IS2,1,'S')
       ENDIF
       IF(IS1.GT.0.AND.IS2.EQ.0) THEN
        ISEQ=IS1
        CALL GRNDMQ(IS1,IS2,ISEQ,'Q')
        if(iprint.gt.0) WRITE(6,'(A,I4,A/A,2(2X,I12))')
     +    ' ACTUAL SEED TAKEN FROM ',ISEQ,'-TH SEED IN GRNDMQ WHICH IS'
     +   ,'           SEED1,SEED2=',IS1,IS2
        CALL GRNDMQ(IS1,IS2,ISEQ,'S')
       ENDIF
      else
       if(iprint.gt.0) then
        write(*,'(a)') ' ************'
        write(*,'(a)') ' ************'
        write(*,'(a)') 
     +   ' No random number seeds file found. Default:9876 54321 start'
        is1 = 9876
        is2 = 54321
        CALL GRNDMQ(IS1,IS2,1,' ')
        IF(iS1.GT.0.AND.IS2.GT.0) THEN
         CALL GRNDMQ(IS1,IS2,1,'S')
        ENDIF
        write(*,'(a)') ' ************'
        write(*,'(a)') ' ************'
       endif
      endif

      end
