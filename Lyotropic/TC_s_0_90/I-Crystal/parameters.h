      implicit real*8(a-h, o-z)
      integer SEED,MX,MY,MZ,NA,LIST
      integer NM,NMZ,MAPSZE,NCELL,MAPS,HEAD, neighbor, Nb
      integer NconI4,epsCI4,grClusI4,histI4
      integer NconI2,epsCI2,grClusI2,histI2
      integer NconP2,epsCP2,grClusP2,histP2
      integer NconQ6,epsCQ6,grClusQ6,histQ6
      integer OP, OPtmp
      double precision cStcount
      integer OP_right,OP_left
      double precision cSt
      double precision r,axes,q,CL,transstep,rotstep,volstep
      double precision HN,HNI,VOLN, VU, vbox, rdis, a1, b1
      double precision aspectratio,rcutsq,rcut,kspring
      double precision rcutsqnb
      double precision press,temp,beta,cellix,aspecthalf
      double precision celliy,celliz,pr
      double precision qus,axesus,CLus,HNus
      complex*8 I4
      double precision I4star
      double precision dcI4
      complex*8 I2
      double precision I2star
      double precision dcI2
      double precision P2star
      double precision dcP2
      complex*8 q6
      double precision dcQ6
      PARAMETER (pr=(3.0d0**0.5d0))
      PARAMETER (aspectratio=1.0d0, a1=2.0d0,b1=2.6d0)
      PARAMETER (rcutsq=(2.05*0.8235607261890482)**2.0)
      PARAMETER (rcutsqnb=(1.74d0)**2.0,kspring=20.0d0)
      PARAMETER (epsCI4 = 7, dcI4=0.70d0)
      PARAMETER (epsCI2 = 5, dcI2=0.85d0)
      PARAMETER (epsCP2 = 5, dcP2=0.910d0)
      PARAMETER (epsCQ6 = 7, dcQ6=0.700d0)
      PARAMETER(NA=2744,aspecthalf=1.0d0*aspectratio)
      PARAMETER( NM=6000000,NMZ=15600000,rcut=rcutsq**0.5d0)

