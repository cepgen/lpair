      function ran2(idum)
c     random numbers uniformly distributed between 0 and 1.
c     (code by A.T. Service, Harvard-Smithsonian Center for Astrophysics;
c     can be replaced by any other suitable generator for random numbers)
      integer idum
      common/ixtbl/ix1,ix2,ix3,rm1,rm2,r(99)
      data ia1,ic1,m1/1279,351762,1664557/
      data ia2,ic2,m2/2011,221592,1048583/
      data ia3,ic3,m3/15551,6150,29101/

      if(idum.ge.0) go to 2
      ix1=mod(-idum,m1)
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ix1,m2)
      ix1=mod(ia1*ix1+ic1,m1)
      ix3=mod(ix1,m3)
      rm1=1./float(m1)
      rm2=1./float(m2)
      do 1 j=1,99
         ix1=mod(ia1*ix1+ic1,m1)
         ix2=mod(ia2*ix2+ic2,m2)
         r(j)=(float(ix1)+float(ix2)*rm2)*rm1
 1    continue
 2    ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(99*ix3)/m3
      ran2=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
c      print *,'-->',ran2,idum
      idum=ix1
      return
      end
