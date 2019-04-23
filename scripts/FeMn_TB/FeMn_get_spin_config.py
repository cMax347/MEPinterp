import numpy as np

deg_to_rad	=	np.pi / 180.




#	phi values are fixed
#		#1		phi = 45 deg
phi_deg		=	np.array([	45.,	-135.,	135.	, -45.	])

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
#      #4 o-------+					theta_{#3,#4} = 180 - theta_order
#		 /	  #3 / |
#		+-------o  |			--------------------------
#		|		|  o #2
#		| /		| /					theta_{#1,#2} =		theta_order
#	 #1 o------	+					
#
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# SPIN configs (radiants)# (theta,phi) are given in radians, which correspond to (beta,alpha) in Fleur.
theta_3q = [ 0.0000000000,  1.9106332362,  1.9106332362, 1.9106332362]
phi_3q   = [ 0.0000000000, -2.0943951024,  2.0943951024, 0.0000000000]

theta_2q = [ 0.6154797087,  2.5261129449,  1.5707963268,  1.5707963268]
phi_2q   = [ 1.0471975512, -2.0943951024,  2.6179938780, -0.5235987756]

theta_1q = [ 0.9553166181,  0.9553166181,  2.1862760355,  2.1862760355]
phi_1q   = [-2.0943951024, -2.0943951024,  1.0471975512,  1.0471975512]

theta_0q = [ 0.0000000000,  1.5707963268,  1.5707963268,  1.5707963268] # For test
phi_0q   = [ 0.0000000000, -2.0943951024,  2.0943951024,  0.0000000000]





def cli_print_lattice(theta_order_deg):	
	print("^^^^ ***FeMn	- SPIN CONFIG SETUP*** ^^^")
	print("order parameter theta="+str(theta_order_deg)+" deg")
	print("#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	print("#																			")
	print("#   #4 o-------+					theta_{#3,#4} = "+str(180 - theta_order_deg))
	print("#     /    #3 / |      														")
	print("#    +-------o  |      			--------------------------					")
	print("#    |  +----|- o #2	    													")
	print("#    | /     | /    					theta_{#1,#2} = "+str(theta_order_deg))
	print("# #1 o------ + ")
	print("#																			")
	print("#																			")
	print("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~											")
	print("")
	



def get_spiral(	theta_order_deg, verbose=False	):
	if verbose:
		cli_print_lattice(theta_order_deg)
	#	convert to rad
	theta_order	=	deg_to_rad * theta_order_deg
	#
	theta_dw	=	np.pi - theta_order
	theta_up	=	theta_order
	#
	theta_sp	=	np.array([	theta_up, theta_up,		theta_dw, theta_dw		])

	#
	#	constant phi
	phi_sp		=	deg_to_rad * phi_deg
	#
	#
	return theta_sp, phi_sp




def test():
	for theta in np.linspace(0.,90.,5):
		theta_sp, phi_sp =	get_spiral(	theta )
		print('theta_sp:'	, theta_sp)
		print('phi_sp:'	, phi_sp)
		print("")
		print("")
		print("")
		print("")
#test()




