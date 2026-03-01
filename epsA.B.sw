# DATE: 2008-01-29 CONTRIBUTOR: Unknown CITATION: Bere and Serra, Phil Mag 86, 2159 (2006)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# GaN potential: A. Bere, and A. Serra, Phil. Mag. 86, 2159(2006)
# note that the parameters for this literature potential are pairwise
# so that there are some flexibility in the way the 
# parameters can be entered. As one way, we assume that 
# lambda_ijk is equal to lambda_ik and eps_ijk is 
# equal to sqrt(lambda_ij*eps_ij*lambda_ik*eps_ik)/lambda_ik, 
# and all other parameters in the ijk line are for ik.

# These entries are in LAMMPS "real" units:
#   epsilon = kcal/mol; sigma = Angstroms
#   other quantities are unitless
#
#         eps   sigma 	a 	lambda 	gamma  cos(theta)     A    	B    		p   q  tol
#
A A A 	1.0 	1.0 	1.8 	0.0 	0.0 	0.0	 7.049556277 	0.6022245584 	4.0 0.0 0.0
A A B  	0.0 	0.0	0.0 	0.0	0.0	0.0	 0.0		0.0		0.0 0.0 0.0
A B A   0.0   	0.0 	0.0	0.0 	0.0	0.0   	 0.0  		0.0 		0.0 0.0 0.0
B A A   1.45 	1.08 	1.8 	0.0	0.0	0.0 	 7.049556277	0.6022245584 	4.0 0.0 0.0
A B B 	1.45	1.08 	1.8 	0.0 	0.0 	0.0 	 7.049556277 	0.6022245584 	4.0 0.0 0.0
B A B  	0.0	0.0   	0.0 	0.0 	0.0 	0.0 	 0.0   		0.0  		0.0 0.0 0.0
B B A 	0.0	0.0   	0.0 	0.0 	0.0 	0.0 	 0.0   		0.0  		0.0 0.0 0.0
B B B 	1.0 	1.0   	1.8 	0.0 	0.0 	0.0	 7.049556277   	0.6022245584  	4.0 0.0 0.0



#pair_style sw
#pair_coeff * * path A B
#
