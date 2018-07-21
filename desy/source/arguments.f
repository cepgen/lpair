      function numarguments() ! get number of command line arguments using gfortran
      implicit none
      integer numarguments ! number of command line arguments

      print *, 'haha'
      numarguments=command_argument_count()

      return
      end

      subroutine getargument(arg,argt) ! get a command line argument using gfortran
      implicit none
      integer arg ! number of argument to get
      character*(*) argt ! output

      call get_command_argument(arg,argt)

      return
      end

