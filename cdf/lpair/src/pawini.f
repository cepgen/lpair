       subroutine pawini
*
* ---- define PAWC common and open histograms
*
       INTEGER ISTAT
       CHARACTER*8 CHTAG1(23)
       PARAMETER (NWPAWC=500000)
       COMMON/PAWC/H(NWPAWC)
*
* --- LPAIR data common block
*
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*
       IF(IPAR(2).EQ.0) RETURN  ! PAW histograms/ntuple
*
       DATA CHTAG1 /'eta1    ','eta2    ','sqrtw4  ',
     &              'px6     ','py6     ','pz6     ','p6st6   ',
     &              'px7     ','py7     ','pz7     ','p7st7   ',
     &              'p4st4   ','phi1    ','phi2    ','acopl   ',
     &              'thet2d  ','thet3d  ','thetp3  ','thetp5  ',
     &              'pt3     ','pz3     ','pt5     ','pz5     '/
*
       CALL HLIMIT(NWPAWC)
*
       CALL HROPEN(31,'HISTO','histo.out','N',1024,ISTAT)
       CALL HBOOKN(100,'lpair ntuple',23,'//HISTO',4000,CHTAG1)
*
       IF(IPAR(4).EQ.0) RETURN  ! cross-section plots
*
       CALL HBOOK1( 1,'log(-t1)      ',40,-0.50E+01,0.50E+01,0.)
       CALL HBOOK1( 2,'log(-t2)      ',40,-0.50E+01,0.50E+01,0.)
       CALL HBOOK1( 3,'mu-mu mass    ',40, 0.00E+00,0.30E+02,0.)
       CALL HBOOK1( 4,'p-perp mu-    ',40, 0.00E+00,2.00E+01,0.)
       CALL HBOOK1( 5,'p-perp mu+    ',40, 0.00E+00,2.00E+01,0.)
       CALL HBOOK1( 6,'p-perp mu-pair',40, 0.00E+00,2.00E+01,0.)
*
       CALL HBOOK1(11,'log(-t1)      ',40,-0.50E+01,0.50E+01,0.)
       CALL HBOOK1(12,'log(-t2)      ',40,-0.50E+01,0.50E+01,0.)
       CALL HBOOK1(13,'mu-mu mass    ',40, 0.00E+00,0.30E+02,0.)
       CALL HBOOK1(14,'p-perp mu-    ',40, 0.00E+00,2.00E+01,0.)
       CALL HBOOK1(15,'p-perp mu+    ',40, 0.00E+00,2.00E+01,0.)
       CALL HBOOK1(16,'p-perp mu-pair',40, 0.00E+00,2.00E+01,0.)
*
       return
       end
