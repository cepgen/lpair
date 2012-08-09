      subroutine kwffrd(ierr)
      common /cfread/ stor(1000)
      call ffinit(1000)
      ierr=0
      end
      subroutine kwffgo(name,ierr)
      character *(*) name
      print *,'Calling FFGO'
      call ffgo
      print *,'End of FFGO'
      ierr=0
      end
      subroutine nulwin
      entry instab
      entry pritab
      entry seltab
      entry creind
      entry getind
      end
      subroutine gmail
      end
      subroutine cletab
      integer i/0/
      i=i+1
      print *,'CLETAB call',i
      end
      subroutine vertex
      end
      subroutine structf
      print *,'STRUCTF call'
      end
      subroutine pdfsta
      end



