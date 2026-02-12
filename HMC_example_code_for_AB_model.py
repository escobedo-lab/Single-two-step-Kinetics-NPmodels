#####################################################################################
#	HMC CODE FOR KM-MODEL: MD (PARALLEL) AND MC (SERIAL) 
#
#	WITH ROTATION MC MOVES: SMALL (~30%) AND LARGE (0 to 180 DEG) STEPS
# T-12090 N=16337
#
#####################################################################################

from mpi4py import MPI
from lammps import lammps,PyLammps
from ctypes import *
import os
import sys,random,math
#import numpy as np
import subprocess
#import time
#from oplib import *
#-------------------- Run Parameters --------------------#

# Thermodynamic state -  configured for LAMMPS' "lj" units

T = 295            # Temperature
P = 0.00               # Pressure

# MC parameters
n_sweeps  = 500000     # Number of MC Sweeps (1 Sweep = 1 HMC or Volume Move)
n_steps = 100            # Number of MD Steps per HMC trial
dt = 0.615399126265             # time step [fs]
p_hmc = 0.80             # Probability of selecting HMC move (vs Volume Move); >= 1.0 for NVT ensemble

p_rotS=0.20
rots_step=0.20


lnvol_max = 0.0176422337491        # Maximum log volume displacement
rseed = 105434             # Random number seed (fixed for testing)

lnvol_max1 = 0.15
lnvol_max2 = 0.15
lnvol_max3 = 0.15

# Biased sampling parameters

nstar_Z=0 #INITIALIZING AS ZERO BUT IT IS CHANGED BASED ON THE INITIAL CONFIG
nstar_G=0
nstar_G_restart=0
nstar_Z_restart=0


nstar_min = 70            # Minimum cluster size for the sampling window [molecules
nstar_max = int(nstar_min+20)          # Maximum cluster size for the sampling window [molecules]
freq_hist_update = 10    # Frequency to accumulate histogram data [MC Sweeps]

print(nstar_min,nstar_max,'min and max')

# File output frequency
freq_thermo = 100         # Thermodynamic output
freq_flush = 500        # Flush files
freq_hist_print=5000

#-------------------- Physical Constants and Conversion Factors --------------------#

comm=MPI.COMM_WORLD
rank=comm.Get_rank()

kB = 1.380648520000e-23  # Boltzmann's constant [m^2 kg s^-2 K^-1]
Nav = 6.02214090000e23   # Avogadro's number [molecules mol^-1]
R = kB * Nav / 1000.0    # Gas constant [kJ mol^-1 K^-1]

# Constants.
#kB = 1.00  # Boltzmann's constant

# Thermal energy
# Thermal energy 
kTs = R*T                # SI units: [kJ mol^-1]
KT = kTs/4.184          # LAMMPS units [real units, kcal mol^-1] 

#KT=kB*T

# Pressure prefactor for volume change move
#Pf=P/(KT)           # Pressure
# Pressure prefactor for volume change move
Pb = P*1.01325           # Pressure [bar]
Pc = kB*1.0e30*1.0e-5    # Conversion factor [bar A^3 K^-1]
Pf = Pb/(Pc*T)           # Prefactor for Metropolis criterion [A^-3]

#-------------------- Initialization --------------------#
print('aaaaaa')
# Seed random number generator
#random.seed(10127379)


# Initialize LAMMPS and GENERATE INITIAL RANDOM VELOCITY FOR RIGID BODY
lmp = lammps(name="",cmdargs=["-log","none","-screen","none"])

#lmp = lammps(name="",cmdargs=["-screen","none"]) #UNCOMMENT THIS AFTER DEBUGGING

#lmp = lammps(name="")
lmp.command("variable dt equal %f" % dt)
lmp.command('variable seed equal %s'%(random.randint(1, 10000000)))
lmp.file("input_ini.inp") ### NOTE THIS GENERATES INITITAL RANDOM VELICITIES (THERE ARE TWO RUN 0) 
#lmp.command("compute thermo_ke all ke")




print('REACHED HERE')
#sys.exit(-1)

# HERE WE CREATE A NEW LAMMPS DATA FILE- IT WILL BE SAME AS INITIAL DATA FILE
lmp.command('write_data NPT_init.lmp')

#sys.exit(1)

Nmol=13824



#print(L.atoms[0].position,'id of 1 part')

# Get initial system properties
natoms = lmp.extract_global("natoms",0)
mass = lmp.extract_atom("mass",2)
atomid = lmp.gather_atoms("id",0,1)
atype = lmp.gather_atoms("type",0,1)
molid= lmp.gather_atoms("molecule",0,1)
#mol= L.atoms[:].mol

# Allocate coordinate and velocity arrays 
x=(3*natoms*c_double)()
x_new = (3*natoms*c_double)()
#v=(3*natoms*c_double)()

# Initialize properties
pe = 0.0
ke = 0.0
etot = 0.0
box = 0.0
vol = 0.0

# Initialize counters
n_acc_hmc = 0.0
n_try_hmc = 0.0
n_acc_vol = 0.0
n_try_vol = 0.0

n_try_rotS= 0.0
n_acc_rotS= 0.0

n_try_rotL= 0.0
n_acc_rotL= 0.0

n_acc_vol1 = 0.0
n_try_vol1 = 0.0

n_acc_vol2 = 0.0
n_try_vol2 = 0.0

n_acc_vol3 = 0.0
n_try_vol3 = 0.0



n_acc_vol_w = 0.0
n_acc_hmc_w = 0.0

dH=0
hmc_acc=0
vol_acc=0
hmc_acc_w=0
vol_acc_w=0

#n_acc_hmc_w = 0.0
#n_acc_vol_w = 0.0

# Get initial coordinates and velocities
x = lmp.gather_atoms("x",1,3)
x_restart=lmp.gather_atoms("x",1,3)
#v = lmp.gather_atoms("v",1,3)

#print(x[0],x[1],x[2],'position 1')
#sys.exit(-1)

pe = (float(lmp.extract_compute("thermo_pe",0,0)))/(KT)	#float(lmp.get_thermo("pe"))/(KT)#float(L.eval("pe"))/(KT)	#float(lmp.extract_compute("thermo_pe",None,0))/float(KT)
#print(pe,'pe')
ke = (float(lmp.extract_compute("thermo_ke",0,0)))/(KT)
print(pe*KT,'pe*KT')
print(ke, float(lmp.extract_compute("thermo_ke",0,0)),'keeeeeeeeeeeeee')
#sys.exit(-1)	#ke = float(lmp.get_thermo("ke"))/(KT)#float(L.eval("ke"))/(KT)
#density= float(lmp.get_thermo("density")) 
temp= float(lmp.extract_compute("thermo_temp",0,0))
virial= float(lmp.extract_compute("thermo_press",0,0))
#xlo = lmp.extract_global("boxxlo",dtype=None)
pe_restart=pe
# Compute box edge length and volume (assumes cubic!)
#boxlo = lmp.extract_global("boxxlo",dtype=None)
#boxhi = lmp.extract_global("boxxhi",dtype=None)
#print(boxlo,boxhi,'box val')
#box = float(boxhi)-float(boxlo)
#vol = math.pow(box,3)

boxxlo = lmp.extract_global("boxxlo",1)
boxxhi = lmp.extract_global("boxxhi",1)

boxylo = lmp.extract_global("boxylo",1)
boxyhi = lmp.extract_global("boxyhi",1)

boxzlo = lmp.extract_global("boxzlo",1)
boxzhi = lmp.extract_global("boxzhi",1)



print(boxxlo,boxxhi,boxylo,boxyhi,'box val')
box1 = float(boxxhi)-float(boxxlo)
box2 = float(boxyhi)-float(boxylo)
box3 = float(boxzhi)-float(boxzlo)
vol = box1*box2*box3
print(vol,'1')
box=box1
vol = math.pow(box,3)
vol_restart=vol
box_restart=box
print(vol,'2')

#for i in range(len(x)):
#	T_1=float(x[i])-(float(box)*round(float(x[i])/float(box)))
#	x[i]=T_1

if comm.rank == 0:

	for i in range(len(x)):
		T_1=float(x[i])-(float(box)*round(float(x[i])/float(box)))
		x[i]=T_1


#print(box,vol)
	#print(boxlo,boxhi,box,vol,'boxlo and boxhi')
#	out1=open('coordT.inp','w')
#	out1.write('%s\t0.0\t0.0\n'%(box))
#	out1.write('0.0\t%s\t0.0\n'%(box))
#	out1.write('0.0\t0.0\t%s\n'%(box))
#        out2=open('coordS.inp','w')
#        out2.write('%s\t0.0\t0.0\n'%(box))
#        out2.write('0.0\t%s\t0.0\n'%(box))
#        out2.write('0.0\t0.0\t%s\n'%(box))
        out3=open('coord.inp','w')
        out3.write('%s\t0.0\t0.0\n'%(box))
        out3.write('0.0\t%s\t0.0\n'%(box))
        out3.write('0.0\t0.0\t%s\n'%(box))


	for i in range(int(natoms)):
		id1=int(3*i)
		out3.write('%s\t%s\t%s\t1\n'%(float(x[id1]),float(x[id1+1]),float(x[id1+2])))
#		if(int(atype[i])==1):
#			id1=int(3*i)
#			out1.write('%s\t%s\t%s\t1\n'%(float(x[id1]),float(x[id1+1]),float(x[id1+2])))

#                if(int(atype[i])==2):
#                        id1=int(3*i)
#                        out2.write('%s\t%s\t%s\t1\n'%(float(x[id1]),float(x[id1+1]),float(x[id1+2])))



	out3.flush()
	out3.close()
	print('HERE')
	#sys.exit(-1)
#HERE WE CALCULATE NOT 4 NEIGH
#	subprocess.call(["./a.out"])
	subprocess.call(["./q6.out"])
	#sys.exit(1)
	subprocess.call(["./q2A.out"])
	subprocess.call(["./q2B.out"])
	subprocess.call(["./clus.out"])

#        subprocess.call(["./q8.out"])
#        subprocess.call(["./q10.out"])
#        subprocess.call(["./clusG.out"])
	subprocess.call(["./DIA.out"])

#READ LARGEST CLUSTER FROM FILE
	out1=open('Lclus.dat','r')
#        cmd1=("cp Nclus.dat Nclus_ZEO.dat")
#        os.system(cmd1)
	nstar_Z=int(out1.readlines()[0])
	out1.close()
        cmd1=("cp Nclus.dat Nclus_ZEO.dat")
        os.system(cmd1)


	print(nstar_Z,'LAM')
#        subprocess.call(["./q8.out"])
#        subprocess.call(["./q10.out"])
#        subprocess.call(["./clusG1.out"])
#	subprocess.call(["./DIA.out"])

#        out1=open('Lclus.dat','r')
#        nstar_G=int(out1.readlines()[0])
#        out1.close()
#        cmd1=("cp Nclus.dat Nclus_GYR.dat")
#        os.system(cmd1)

#	print(nstar_G,'GYROID')

	print(nstar_Z,'Lclus and Dual Lclus1 and 2')
	#sys.exit(-1)
#OP.write('%s\t%s\n'%(0,nstar))
#if MPI.COMM_WORLD.rank == 0:
	thermo=open('thermo.dat','w')
	traj= open('traj.lammpstrj','w')
	#local_gyr=open('local_clus_gyr.dat','w')
	local_zeo=open('local_clus_lam.dat','w')
#	thermo.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(0,nstar,temp,pe*(KT),box,0,0,0))
#thermo.write('%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'%(0,nstar,vol,pe*(KT),hmc_acc,vol_acc,hmc_acc_w,vol_acc_w,dt,lnvol_max))


	f10=open('Nclus_ZEO.dat','r')

	for line in f10:
		CL5=line.split()[0]
		CL6=line.split()[1]
		local_zeo.write('%s\t%s\t%s\n'%((0),(CL5),(CL6)))

	f10.close()

	cmd1=("cp Nclus_ZEO.dat Nclus_ZEO_restart.dat")
        os.system(cmd1)

	nstar_Z_restart=nstar_Z

	thermo.write('%d\t%d\t%s\t%f\t%f\t%f\n'%(0,nstar_Z,pe*(KT),vol,dt,lnvol_max))
	thermo.flush()
	local_zeo.flush()

        cmd1=("cp NPT1.lmp NPT_restart.lmp")
        os.system(cmd1)

#sys.exit(-1)
#	thermo.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(0,nstar,temp,pe,density,press))

	#WRITING INITIAL CONFIG TO TRAJ FILE
	traj.write('ITEM: TIMESTEP\n')
	traj.write('0\n')
	traj.write('ITEM: NUMBER OF ATOMS\n')
	traj.write('%s\n'%(natoms))
	traj.write('ITEM: BOX BOUNDS pp pp pp\n')
	traj.write('%s\t%s\n'%(boxxlo,boxxhi))
	traj.write('%s\t%s\n'%(boxylo,boxyhi))
	traj.write('%s\t%s\n'%(boxzlo,boxzhi))
	traj.write('ITEM: ATOMS id mol type q x y z\n')
	for i in range(int(natoms)):
		id1=int(3*i)
		traj.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(int(atomid[i]),int(molid[i]),int(atype[i]),0,float(x[id1]),float(x[id1+1]),float(x[id1+2])))
#thermo.flush()
	traj.flush()

#sys.exit(-1)
comm.Barrier()
nstar_G=comm.bcast(nstar_G, root=0)
nstar_Z=comm.bcast(nstar_Z, root=0)
nstar_G_restart=comm.bcast(nstar_G_restart, root=0)
nstar_Z_restart=comm.bcast(nstar_Z_restart, root=0)


#Nmol=512
#hist = []
#for i in range(int(Nmol)+1):
#	hist.append(0.0)	#p.zeros(int(Nmol)+1)
#Nrecord_hist=0
print('REACHED isweep')

cross_count=0
count_FE=0
#lmp1 = lammps(name="")
flag=0
#---------------------- MC -----------------------------------#

#print(ke*KT,pe,ke,KT,T,'ke and pe val')

for isweep in range(n_sweeps):
#	print(isweep+1,'isweep')

#OP STORING FOR OP CALCULATION EVERY 2 MOVES
#	if(int):
        if comm.rank == 0:
                pickmove=random.random()
        else:
                pickmove=None

	comm.Barrier()
        pickmove=comm.bcast(pickmove, root=0)




	# HMC Move
	if(pickmove<=p_hmc):
#		print('NVT MOVE')
		if(abs((ke*KT)-12155.102509739254)>0.0000001):
			print(ke*KT,pe,'ke and pe val')
			print('HMC KE FAIL')
			sys.exit(-1)
		#n_try_hmc += 1.0
		etot=pe+ke #KE HERE HAS TO BE 6444.9
		lmp= lammps(name="",cmdargs=["-log","none","-screen","none"])
		#lmp = lammps(name="")
		lmp.command('variable seed equal %s'%(random.randint(1, 10000000)))
		lmp.command("variable dt equal %f" % dt)
		lmp.command("variable nsteps equal %s" % n_steps)
#		lmp= lammps(name="",cmdargs=["-log","none","-screen","none"])
		comm.Barrier()
		lmp.file("input2.inp")
		comm.Barrier()
		#x = lmp1.gather_atoms("x",1,3)
#		sys.exit(-1)
		pe_new = (float(lmp.extract_compute("thermo_pe",0,0)))/(KT)
		ke_new = (float(lmp.extract_compute("thermo_ke",0,0)))/(KT)	#(float(L.eval("ke")))/(KT)
#		virial_new=(float(lmp.extract_compute("thermo_press",0,0)))
		#temp_new=float(lmp.extract_compute("thermo_temp",0,0))
		x_new=lmp.gather_atoms("x",1,3)
#		print(box,'BOX LENGTH')

		flag=0

		etot_new=pe_new+ke_new
		dH = (etot_new-etot)
		comm.Barrier()
		if comm.rank == 0:
			n_try_hmc += 1.0

	                for i in range(len(x_new)):
				T_1=float(x_new[i])-(float(box)*round(float(x_new[i])/float(box)))
				x_new[i]=T_1
	
			if(random.random() <= math.exp(-(dH))):
#ACCEPTED
#SINCE IT GOT ACCEPTED CHECK FOR OP WINDOW
#				if((isweep+1)%2==0):
#				print('METRO ACCEPT')
				n_acc_hmc += 1
				pe=pe_new

				for i in range(3*natoms):
                                	x[i] = float(x_new[i])
                                flag=1

				cmd1=("cp NPT_NVE.lmp NPT1.lmp")
				os.system(cmd1)


		

#		comm.Barrier()
#		pe=comm.bcast(pe, root=0)
#		nstar_G=comm.bcast(nstar_G, root=0)
#		nstar_Z=comm.bcast(nstar_Z, root=0)
#		flag=comm.bcast(flag, root=0)
#		if(flag==1):
#			comm.Barrier()
#			lmp.command('write_data NPT1.lmp')
#RESTART FOR THER NEXT STEP SINCE IT GOT ACCEPT
#			cmd1=("cp NPT1.lmp NPT_restart.lmp")
#			os.system(cmd1)
		comm.Barrier()

		if comm.rank == 0:
			if((isweep+1)%5==0):
#				print('OPcheck now')
                        	out3=open('coord.inp','w')
                        	out3.write('%s\t0.0\t0.0\n'%(box))
                        	out3.write('0.0\t%s\t0.0\n'%(box))
                        	out3.write('0.0\t0.0\t%s\n'%(box))


                        	for i in range(int(natoms)):
                        		id1=int(3*i)
                                	out3.write('%s\t%s\t%s\t1\n'%(float(x[id1]),float(x[id1+1]),float(x[id1+2])))





				out3.flush()
                        	out3.close()

#HERE WE CALCULATE NOT 4 NEIGH
#                        	subprocess.call(["./a.out"])
#                        	subprocess.call(["./q6.out"])
#                        	subprocess.call(["./clus.out"])
#                        	subprocess.call(["./DIA.out"])
#				subprocess.call(["./a.out"])
        			subprocess.call(["./q6.out"])
        			subprocess.call(["./q2A.out"])
				subprocess.call(["./q2B.out"])
        			subprocess.call(["./clus.out"])

#        			subprocess.call(["./q8.out"])
#        			subprocess.call(["./q10.out"])
#        			subprocess.call(["./clusG.out"])
        			subprocess.call(["./DIA.out"])

#READ LARGEST CLUSTER FROM FILE
                        	out1=open('Lclus.dat','r')
                        	nstarnew_Z=int(out1.readlines()[0])
                        	out1.close()
#				print('done with ZEO')

                        	cmd1=("cp Nclus.dat Nclus_ZEO.dat")
                        	os.system(cmd1)


#                        	subprocess.call(["./q8.out"])
#                        	subprocess.call(["./q10.out"])
#                        	subprocess.call(["./clusG1.out"])
#                        	subprocess.call(["./DIA.out"])
#                        	out1=open('Lclus.dat','r')
#                        	nstarnew_G=int(out1.readlines()[0])
#                        	out1.close()
#
#                        	cmd1=("cp Nclus.dat Nclus_GYR.dat")
#                        	os.system(cmd1)


#                        	print(nstarnew_G,nstarnew_Z,'over')

				if(nstarnew_Z>= nstar_min and nstarnew_Z <= nstar_max):

#                        		print('ACCEPTED')
	                        	n_acc_hmc_w += 1
                                #pe = pe_restart
                                #ke=ke_new
#                               virial=virial_new
                                #temp=temp_new
					pe_restart=pe
                                	nstar_G = 0
                                	nstar_Z = nstarnew_Z

                                        nstar_G_restart = 0
                                        nstar_Z_restart = nstar_Z


                                        box_restart=box
                                        vol_restart=vol

                                        for i in range(3*natoms):
                                                x_restart[i] = float(x[i])


					cmd1=("cp Nclus_ZEO.dat Nclus_ZEO_restart.dat")
					os.system(cmd1)
#					cmd1=("cp Nclus_GYR.dat Nclus_GYR_restart.dat")
#					os.system(cmd1)


                        #lmp.command('write_data NPT1.lmp')
#                                	for i in range(3*natoms):
#                                		x[i] = float(x_new[i])


                                #flag=1
                                	cmd1=("cp NPT1.lmp NPT_restart.lmp")
                                	os.system(cmd1)
				else:
#                        		print('rejected going back n steps')
					pe=pe_restart

                                        box=box_restart
                                        vol=vol_restart
                                        nstar_G = 0
                                        nstar_Z = nstar_Z_restart

                                        cmd1=("cp Nclus_ZEO_restart.dat Nclus_ZEO.dat")
                                        os.system(cmd1)
#                                        cmd1=("cp Nclus_GYR_restart.dat Nclus_GYR.dat")
#                                        os.system(cmd1)


					for i in range(3*natoms):
						x[i] = float(x_restart[i])

                                	cmd1=("cp NPT_restart.lmp NPT1.lmp")
                                	os.system(cmd1)
					#sys.exit(-1)


                comm.Barrier()
                pe=comm.bcast(pe, root=0)
                pe_restart=comm.bcast(pe_restart, root=0)
                box_restart=comm.bcast(box_restart, root=0)
                vol_restart=comm.bcast(vol_restart, root=0)
                nstar_G=comm.bcast(nstar_G, root=0)
                nstar_Z=comm.bcast(nstar_Z, root=0)
                nstar_G_restart=comm.bcast(nstar_G_restart, root=0)
                nstar_Z_restart=comm.bcast(nstar_Z_restart, root=0)

                vol=comm.bcast(vol, root=0)
                box=comm.bcast(box, root=0)
                comm.Barrier()

	#Log-volume MC Move
#p_rotS=0.20
#rots_step=0.20


	else:
#		print('VOL MOVE')

		if comm.rank == 0:
			n_try_vol += 1.0
	

			lnvol= math.log(vol) + (random.random() - 0.5)*lnvol_max	
		# Calculate new boxx volume, size and scale factor
			vol_new = float(math.exp(lnvol))
			box_new = float(math.pow(vol_new,1.0/3.0))
			#print(vol_new,vol,'VOL')
#		vol_new = float(math.exp(lnvol))
#		box_new = float(math.pow(vol_new,1.0/3.0))
		# Scale the coordinates and send to LAMMPS

#		for i in range(3*natoms):
#			x_new[i]=scalef*x[i]

                	scalef = box_new/box

                	for i in range((3*natoms)):
                        	x_new[i]=(x[i])*float(scalef)


					

		## NEED TO CORRECT FOR PATCHES: TRANSLATE IN DX,DY,DZ IN SAME AMOUNT AS CENTER
		## CANNOT WRAP THE CORRECTED PATCH COORDINATES AS IT GIVES ERROR IN FLAGS

		#lmp.scatter_atoms("x",1,3,x_new)
		#lmp.command('write_data NPT_new_NPT.lmp')
		#sys.exit(-1)
		#WRITING A NEW NPT1.lmp FILE: NOTE in 2023 version you can set image flags to 0
			config= open('NPT1_new.lmp','w')
			config.write('LAMMPS VOL MOVE\n')
			config.write('\n')
			config.write('%s atoms\n'%(natoms))
			config.write('2 atom types\n')
			config.write('\n')
			config.write('%s\t%s\txlo\txhi\n'%(-(box_new)/2.0,box_new/2.0))
			config.write('%s\t%s\tylo\tyhi\n'%(-(box_new)/2.0,box_new/2.0))
			config.write('%s\t%s\tzlo\tzhi\n'%(-(box_new)/2.0,box_new/2.0))
			config.write('\n')
			config.write('Masses\n')
			config.write('\n')
			config.write('1 1.0\n')
			config.write('2 1.0\n')
			config.write('\n')
			config.write('Atoms #full\n')
			config.write('\n')
			for i in range(int(natoms)):
		        	id1=int(3*i)
		        	config.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(int(atomid[i]),int(molid[i]),int(atype[i]),0,float(x_new[id1]),float(x_new[id1+1]),float(x_new[id1+2])))
			config.flush()
			config.close()
#		print(pe,'PE OLD')
		lmp= lammps(name="",cmdargs=["-log","none","-screen","none"])		
		#lmp=lammps(name="")
		#lmp.command("set atom * image 0 0 0")
		comm.Barrier()
		lmp.file("input_vol.inp")
		comm.Barrier()
#		lmp.command("set atom * image 0 0 0")
		# Compute the new PE
		pe_new= (float(lmp.extract_compute("thermo_pe",0,0)))/(KT) #(float(L.eval("pe")))/(KT)
#		virial_new=(float(lmp.extract_compute("thermo_press",0,0)))
#		print(pe_new,'print_new')



		#Calculate argument for the acceptance criterion
		comm.Barrier()
		if comm.rank == 0:
			arg=(pe_new-pe) + Pf*(vol_new-vol) - ((float(Nmol)+1.0)*math.log(vol_new/vol))
		# Perform Metropolis acceptance-rejection test
#			print(arg, 'MC ARGUMENT')
			if(random.random() <= math.exp(-(arg))):
#				print('ACCEPTED')
				n_acc_vol += 1.0
				pe = pe_new

                                vol = vol_new
#                               virial=virial_new
				box = box_new

				for i in range((3*natoms)):
					x[i]=float(x_new[i])-(float(box)*round(float(x_new[i])/float(box)))
#
				cmd=("cp NPT1_new.lmp NPT1.lmp")
				os.system(cmd)



#				if((isweep+1)%1==0):
#	                                out3=open('coord.inp','w')
#	                                out3.write('%s\t0.0\t0.0\n'%(box_new))
#	                                out3.write('0.0\t%s\t0.0\n'%(box_new))
#	                                out3.write('0.0\t0.0\t%s\n'%(box_new))


#	                                for i in range(int(natoms)):
#	                                        id1=int(3*i)
#	                                        out3.write('%s\t%s\t%s\t1\n'%(float(x_new[id1]),float(x_new[id1+1]),float(x_new[id1+2])))





#	                                out3.flush()
#	                                out3.close()



#HERE WE CALCULATE NOT 4 NEIGH
#	                                subprocess.call(["./a.out"])
#	                                subprocess.call(["./q6.out"])
#	                                subprocess.call(["./clus.out"])
#READ LARGEST CLUSTER FROM FILE
#	                                out1=open('Lclus.dat','r')
#	                                nstarnew_Z=int(out1.readlines()[0])
#	                                out1.close()

#	                                subprocess.call(["./q8.out"])
#	                                subprocess.call(["./q10.out"])
#	                                subprocess.call(["./clusG.out"])

#	                                out1=open('LclusG.dat','r')
#	                                nstarnew_G=int(out1.readlines()[0])
#	                                out1.close()


#					if(nstarnew_G >= nstar_min and nstarnew_G <= nstar_max):
#						n_acc_vol_w += 1.0
#						pe = pe_new
#						vol = vol_new
#				virial=virial_new
#				box = box_new

#	                                	box = box_new
#
#	                                        nstar_G = nstarnew_G
#	                                        nstar_Z = nstarnew_Z



#	                                	for i in range((3*natoms)):
#							x[i]=float(x_new[i])-(float(box)*round(float(x_new[i])/float(box)))
#
#						cmd=("cp NPT1_new.lmp NPT1.lmp")
#						os.system(cmd)
#				else:
#					print('CNNOT GO HERE')
#					sys.exit(-1)
#
#					pe = pe_new
#					vol = vol_new
#
#                                       box = box_new

#					for i in range((3*natoms)):
#						x[i]=float(x_new[i])-(float(box)*round(float(x_new[i])/float(box)))
#
#					cmd=("cp NPT1_new.lmp NPT1.lmp")
#					os.system(cmd)

#		else: #Reject

#		comm.Barrier()
#		pe=comm.bcast(pe, root=0)
#		nstar_G=comm.bcast(nstar_G, root=0)
#		nstar_Z=comm.bcast(nstar_Z, root=0)
#		vol=comm.bcast(vol, root=0)
#		box=comm.bcast(box, root=0)
		comm.Barrier()
#		print('REACHED OPs')
                if comm.rank == 0:
                        if((isweep+1)%5==0):
#                                print('OPcheck now')
                                out3=open('coord.inp','w')
                                out3.write('%s\t0.0\t0.0\n'%(box))
                                out3.write('0.0\t%s\t0.0\n'%(box))
                                out3.write('0.0\t0.0\t%s\n'%(box))


                                for i in range(int(natoms)):
                                        id1=int(3*i)
                                        out3.write('%s\t%s\t%s\t1\n'%(float(x[id1]),float(x[id1+1]),float(x[id1+2])))





                                out3.flush()
                                out3.close()

#HERE WE CALCULATE NOT 4 NEIGH
#                                subprocess.call(["./a.out"])
                                subprocess.call(["./q6.out"])
                                subprocess.call(["./q2A.out"])
				subprocess.call(["./q2B.out"])
                                subprocess.call(["./clus.out"])

#                                subprocess.call(["./q8.out"])
#                                subprocess.call(["./q10.out"])
#                                subprocess.call(["./clusG.out"])
                                subprocess.call(["./DIA.out"])

#READ LARGEST CLUSTER FROM FILE
                                out1=open('Lclus.dat','r')
                                nstarnew_Z=int(out1.readlines()[0])
                                out1.close()
#                                print('done with ZEO')

                                cmd1=("cp Nclus.dat Nclus_ZEO.dat")
                                os.system(cmd1)


                                #subprocess.call(["./q8.out"])
                                #subprocess.call(["./q10.out"])
#                                subprocess.call(["./clusG1.out"])
#                                subprocess.call(["./DIA.out"])
#                                out1=open('Lclus.dat','r')
#                                nstarnew_G=int(out1.readlines()[0])
#                                out1.close()
#
#                                cmd1=("cp Nclus.dat Nclus_GYR.dat")
#                                os.system(cmd1)



#                                print(nstarnew_G,nstarnew_Z,'over')

#                                print(nstarnew_G,nstarnew_Z,'over')

                                if(nstarnew_Z>= nstar_min and nstarnew_Z <= nstar_max):

#                                        print('ACCEPTED')
                               		n_acc_vol_w += 1
                                #pe = pe_restart
                                #ke=ke_new
#                               virial=virial_new
                                #temp=temp_new
                                        pe_restart=pe
					box_restart=box
					vol_restart=vol
                                        nstar_G = 0
                                        nstar_Z = nstarnew_Z

                                        nstar_G_restart = 0
                                        nstar_Z_restart = nstar_Z

                                        cmd1=("cp Nclus_ZEO.dat Nclus_ZEO_restart.dat")
                                        os.system(cmd1)
#                                        cmd1=("cp Nclus_GYR.dat Nclus_GYR_restart.dat")
#                                        os.system(cmd1)


                                        for i in range(3*natoms):
                                                x_restart[i] = float(x[i])


                        #lmp.command('write_data NPT1.lmp')
#                                       for i in range(3*natoms):
#                                               x[i] = float(x_new[i])


                                #flag=1
                                        cmd1=("cp NPT1.lmp NPT_restart.lmp")
                                        os.system(cmd1)
                                else:
#                                        print('rejected going back n steps')
                                        pe=pe_restart
					box=box_restart
					vol=vol_restart
                                        nstar_G = 0
                                        nstar_Z = nstar_Z_restart

                                        cmd1=("cp Nclus_ZEO_restart.dat Nclus_ZEO.dat")
                                        os.system(cmd1)
#                                        cmd1=("cp Nclus_GYR_restart.dat Nclus_GYR.dat")
#                                        os.system(cmd1)


                                        for i in range(3*natoms):
                                                x[i] = float(x_restart[i])

                                        cmd1=("cp NPT_restart.lmp NPT1.lmp")
                                        os.system(cmd1)
                                        #sys.exit(-1)


                comm.Barrier()

#                pe=comm.bcast(pe, root=0)
#                nstar_G=comm.bcast(nstar_G, root=0)
#                nstar_Z=comm.bcast(nstar_Z, root=0)
                vol=comm.bcast(vol, root=0)
                box=comm.bcast(box, root=0)
#                comm.Barrier()
                pe=comm.bcast(pe, root=0)
                pe_restart=comm.bcast(pe_restart, root=0)
		box_restart=comm.bcast(box_restart, root=0)
		vol_restart=comm.bcast(vol_restart, root=0)
                nstar_G=comm.bcast(nstar_G, root=0)
                nstar_Z=comm.bcast(nstar_Z, root=0)

                nstar_G_restart=comm.bcast(nstar_G_restart, root=0)
                nstar_Z_restart=comm.bcast(nstar_Z_restart, root=0)

		comm.Barrier()
#	if((isweep+1)<500):





	#Adjust vol step size

#	print(nstar,'nstar')
	if comm.rank == 0:
#		print(nstar,'nstar')
		if((isweep+1)>5000 and (isweep+1)%5==0 and nstar_Z>=nstar_min+2 and cross_count==0):
			cross_count=1
			if(cross_count==1):

				cmd=("cp NPT1.lmp NPT-Rnstar%s.lmp"%(nstar_max))
				os.system(cmd)
				#sys.exit(-1)


	#FE CALCULATION EVERY 10 STEPS


		if((isweep+1)%freq_thermo ==0): #THERMO DATA





			#f2.close()

			f10=open('Nclus_ZEO.dat','r')

			for line in f10:
				CL5=line.split()[0]
				CL6=line.split()[1]
				local_zeo.write('%s\t%s\t%s\n'%((isweep+1),(CL5),(CL6)))

			f10.close()



                        if(n_try_hmc==0.0):
                                hmc_acc=0.0
                                hmc_acc_w=0.0
                        else:
                                hmc_acc = n_acc_hmc/n_try_hmc   #acc_ratio(n_acc_hmc, n_try_hmc)
                                hmc_acc_w = n_acc_hmc_w/n_try_hmc       #acc_ratio(n_acc_hmc_w, n_try_hmc)

                        if(n_try_vol==0.0):
                                vol_acc=0.0
                                vol_acc_w=0.0
                        else:
                                vol_acc = n_acc_vol/n_try_vol   #acc_ratio(n_acc_hmc, n_try_hmc)
                                vol_acc_w = n_acc_vol_w/n_try_vol       #acc_ratio(n_acc_hmc_w, n_try_hmc)



#		thermo.write('%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n'%(0,nstar,hist[0],hist[1],nstar_dual,pe*(KT),dt,lnvol_max))
			thermo.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%((isweep+1),(nstar_Z),pe*(KT),vol,dt,lnvol_max,hmc_acc_w,hmc_acc,vol_acc,vol_acc_w))
#		thermo.flush()
#		hmc_acc = acc_ratio(n_acc_hmc, n_try_hmc)
#		vol_acc = acc_ratio(n_acc_vol, n_try_vol)
#               hmc_acc_w = acc_ratio(n_acc_hmc_w, n_try_hmc)
#                vol_acc_w = acc_ratio(n_acc_vol_w, n_try_vol)


        	if((isweep+1)%5000==0):
#                vratio = acc_ratio(n_acc_vol, n_try_vol)
                        if(n_try_vol==0.0):
                                vratio=0.0
                        else:
                                vratio = n_acc_vol/n_try_vol    #acc_ratio(n_acc_vol, n_try_vol)

	                if(vratio < 0.30):
	                	lnvol_max=lnvol_max*0.95
	                else:
	                        lnvol_max=lnvol_max*1.05
	                n_acc_vol=0
	                n_try_vol=0
	                n_acc_vol_w=0

        	if((isweep+1)%1000==0):
#                	tratio = n_acc_hmc/n_try_hmc
                        if(n_try_hmc==0.0):
                                tratio=0.0
                        else:
                                tratio = n_acc_hmc/n_try_hmc

                	if(tratio < 0.60):
                        	dt=dt*0.95
                	else:
                        	dt=dt*1.05
                	n_acc_hmc=0
                	n_try_hmc=0
                	n_acc_hmc_w=0





		if((isweep + 1) % freq_flush == 0):
			thermo.flush()
#			local_gyr.flush()
			local_zeo.flush()

#			local_DIA.flush()
#			local_BCC.flush()
#		traj.flush()

	comm.Barrier()
	dt=comm.bcast(dt, root=0)
	lnvol_max=comm.bcast(lnvol_max, root=0)
#	rots_step=comm.bcast(rots_step, root=0)	
        if(n_steps!=100):
                print('n_steps problem')
                sys.exit(-1)
	comm.Barrier()

#pe = L.eval("pe")
# Close files
thermo.close()
traj.close()
# Compute box edge length and volume (assumes cubic!)
