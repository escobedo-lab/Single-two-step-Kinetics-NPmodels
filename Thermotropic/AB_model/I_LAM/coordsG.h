	common /coords/ r(3,NA),axes(3,3,NA),q(3,NA), VU(3,24),vbox(3,8)
	common /boxsize/ CL(3),cellix,celliy,celliz,MX,MY,MZ
	common /movsizes/ transstep(3),rotstep(3),volstep(3,3)
	common /var/ press,temp,beta,SEED,LIST(NA),NA_actual
	common /var1/ neighbor(10000,NA), Nb(NA),Nb1(NA)
	common /var2/ rdis(2000,NA),neighbor1(2000,NA)
	common /map/ MAPS(NMZ),HEAD(NM),NCELL,MAPSZE 
	common /tensor/ HN(3,3),HNI(3,3),VOLN
	common /I4val/ I4(-4:4,NA)
	common /I4val2/ I4star(NA)
        common /I4val3/ NconI4(NA),grClusI4,histI4(NA)
	common /I2val/ I2(-2:2,NA),N_ZEO,ID_TZ(NA) 
	common /I2val2/ I2star(NA),qmavg2(NA)
        common /I2val3/ NconI2(NA),grClusI2,histI2(NA) 
	common /P2val2/ P2star(NA),allClus,NA_ref2
	common /P2val3/ NconP2(NA),grClusP2,histP2(NA)
        common /Q6val/ q6(-6:6,NA),qmavg1(NA),NA_ref1
        common /Q6val1/ NconQ6(NA), grClusQ6,histQ6(NA) 
	common /Clus1/ belongsto(NA), ML,atype(NA) 
	common /varUS1/ OP, OPtmp, OP_right, OP_left
	common /clusSt/ cSt(NA+1),cStcount
        common /varUS3/ qus(3,NA),axesus(3,3,NA)
	common /boxsizeus/ CLus(3),HNus(3,3)
