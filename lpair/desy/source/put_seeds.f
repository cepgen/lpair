      subroutine put_seeds(iprint)

      implicit none

* store geant random number generator seeds

      integer iprint
      integer is1,is2 
*
* save grndm seeds
*
      CALL GRNDMQ(is1,iS2, 0,'G')    
      open(98,file='./grndm.seeds',form='formatted')
      write(98,*) is1,is2
      close(98)
      if(iprint.gt.0) then
       write(*,'(a)') ' ************'
       write(*,'(a)') ' ************'
       write(*,'(a,2(2x,i10))') 
     +     ' Geant random number seeds stored to grndm.seeds',is1,is2
       write(*,'(a)') ' ************'
       write(*,'(a)') ' ************'
      endif

      end
