from Modules.dependencies import *
from Modules.units import *
z = 0.022 #defined here to use global variable and avoid calling as parameter. Might change that since still have to call global inside of the functions


###############################################################################################
################################# UNUSED FOR NOW ##############################################
###############################################################################################



def Polarization_Temporal_Evolution(t,a, alpha):
	""" Temporal Evolution Polarization Calculation"""
	xi_max = t - 1
	x1 = (t-1)/(1+ a) 
	xi_min = max(0 , x1)
	epislon = alpha + 1
	def Pitch(xi,phi):
		""" Pitch factor calculation"""
		x1 = numpy.power((1-xi)/(1+xi),2) * numpy.power(numpy.cos(phi),2)
		x2 = numpy.power(numpy.sin(phi),2)
		return numpy.power(x1+x2,epislon/2)
	def Pos_Angle(xi,phi):
		"""Polarization angle calculation"""
		x1 = numpy.arctan(((1-xi)/(1+xi))/numpy.tan(phi))
		return phi + x1
	def InteA(xi):
		"""Integration on xi of the dividend"""
		div = 1/numpy.power((1+xi),3+alpha)
		def Inte2(phi):
			"""Integration on phi of the dividend"""
			x1 = Pitch(xi,phi)*numpy.cos(2*Pos_Angle(xi,phi))
			return x1
		I1= scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return div * I1
	def InteB(xi):
		"""Integration on xi of the divisor"""
		div = 1/numpy.power(1+xi,3+alpha)
		def Inte2(phi):
			"""Integration on phi of the divisor"""
			x1 = Pitch(xi,phi)
			return x1
		I1= scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return div * I1

	I1 = scipy.integrate.quad(InteA, xi_min, xi_max)[0]
	I2 = scipy.integrate.quad(InteB, xi_min, xi_max)[0]
	return -I1*numpy.power(I2,-1)


def Polarization_Ord_Pulse(alpha,ximax):
	""" Time integrated Polarization for a Ordered Field"""
	epislon = alpha + 1
	def Pitch(xi,phi):
		""" Pitch factor calculation"""
		x1 = numpy.power((1-xi)/(1+xi),2) * numpy.power(numpy.cos(phi),2)
		x2 = numpy.power(numpy.sin(phi),2)
		return numpy.power(x1+x2,epislon/2)
	def Pos_Angle(xi,phi):
		"""Polarization angle calculation"""
		x1 = numpy.arctan(((1-xi)/(1+xi))/numpy.tan(phi))
		return phi + x1
	def InteA(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		"""Integration on xi of the dividend"""
		def Inte2(phi):
			"""Integration on phi of the dividend"""
			x1 = Pitch(xi,phi)*numpy.cos(2*Pos_Angle(xi,phi))
			return x1
		I1= scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return div * I1
	def InteB(xi):
		"""Integration on xi of the divisor"""
		div = 1/numpy.power(1+xi,2+alpha)
		def Inte2(phi):
			"""Integration on phi of the divisor"""
			x1 = Pitch(xi,phi)
			return x1
		I1= scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return div * I1

	I1 = scipy.integrate.quad(InteA, 0, ximax)[0]
	I2 = scipy.integrate.quad(InteB, 0, ximax)[0]
	return abs(I1*numpy.power(I2,-1))



###############################################################################################
####################################### OFF-AXIS ##############################################
######################################  TOP-HAT   #############################################
##################################### JET STRUCTURE ###########################################
###############################################################################################


def Polarization_Toroidal_TopHat(xi_j, q, alpha):
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1 + alpha
	def Phi_angle(xi):
		x1 = (1-numpy.power(q,2))*xi_j - xi
		x2 = 2*q*numpy.sqrt(xi*xi_j)
		return numpy.arccos(x1/x2)
	def Pos_Angle(xi,phi):
		a = numpy.sqrt(xi/xi_j)/q
		x1 = ((1-xi)/(1+xi))
		x2 = numpy.sin(phi)/(a+numpy.cos(phi))
		arctan = numpy.arctan(x1*x2)
		return phi - arctan
	def Pitch(xi, phi):
		a = numpy.sqrt(xi/xi_j)/q
		x1 = numpy.power((1-xi)/(1+xi),2)
		x2 = 4*xi/numpy.power(1+xi,2)
		x3 = numpy.power(a+numpy.cos(phi),2)/(1+numpy.power(a,2)+2*a*numpy.cos(phi))
		return numpy.power(x1 + x2*x3, epislon/2)
	def InteA(xi):
		heav1 = numpy.heaviside(1-q,0)
		div = 1./numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch(xi,phi)*numpy.cos(2*Pos_Angle(xi,phi))
			return x1
		I1 = scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return heav1*div*I1
	def InteB(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch(xi,phi)*numpy.cos(2*Pos_Angle(xi,phi))
			return x1
		I1 = scipy.integrate.quad(Inte2, Phi_angle(xi), 2*scipy.pi-Phi_angle(xi))[0]
		return div*I1
	def InteC(xi):
		heav1 = numpy.heaviside(1-q,0)
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch(xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return heav1*div*I1
	def InteD(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch(xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte2, Phi_angle(xi), 2*scipy.pi-Phi_angle(xi))[0]
		return div*I1

	I1a = scipy.integrate.quad(InteA, 0, ximin)[0]
	I1b = scipy.integrate.quad(InteB, ximin, ximax)[0]
	I1 = I1a + I1b
	I2a = scipy.integrate.quad(InteC, 0, ximin)[0]
	I2b = scipy.integrate.quad(InteD, ximin, ximax)[0]
	I2 = I2a + I2b
	return I1*numpy.power(I2,-1)


#Parallel Polarization
def Polarization_Parallel_TopHat(xi_j,q, alpha):
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1+alpha
	def Pitch(xi):
		x1 = numpy.sqrt(4*xi)/(1+xi)
		return numpy.power(x1,epislon)
	def Phi_angle(xi):
		x1 = (1-numpy.power(q,2))*xi_j - xi
		x2 = 2*q*numpy.sqrt(xi*xi_j)
		return numpy.arccos(x1/x2)
	def InteA(xi):
		x1 = 1./(2*scipy.pi)
		x2 = 1/numpy.power((1+xi),2+alpha)
		x3 = Pitch(xi)*numpy.sin(2*Phi_angle(xi))
		return x1*x2*x3
	def InteB(xi):
		heav1 = numpy.heaviside(1-q,0)
		x2 = Pitch(xi)/numpy.power((1+xi),2+alpha)
		return heav1*x2
	def InteB2(xi):
		x2 = (scipy.pi-Phi_angle(xi))/(scipy.pi*numpy.power((1+xi),2+alpha))
		return x2*Pitch(xi)

	I1 = scipy.integrate.quad(InteA,ximin,ximax)[0]
	I2a = scipy.integrate.quad(InteB,0,ximin)[0]
	I2b = scipy.integrate.quad(InteB2, ximin, ximax)[0]

	return -I1*numpy.power(I2a+I2b,-1) 

#Perpendicular Polarization
def Polarization_Perpendicular_TopHat(xi_j,q, alpha):
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1+alpha
	def Pitch_g(xi):
		def Inte1(phi):
			x1 = numpy.power((1-xi)/(1+xi),2)*numpy.power(numpy.cos(phi),2)
			x2 = numpy.power(numpy.sin(phi),2)
			x3 = numpy.power(1-4*xi*numpy.power(numpy.cos(phi),2)*numpy.power(1+xi,-2),(epislon-2)/2)
			return (x2-x1)*x3
		I1 = scipy.integrate.quad(Inte1, 0, scipy.pi)[0]
		return I1
	def Pitch_f(xi):
		def Inte1(phi):
			x3 = 1-4*xi*numpy.power(numpy.cos(phi),2)*numpy.power(1+xi,-2)
			return numpy.power(x3,epislon/2)
		I1 = scipy.integrate.quad(Inte1, 0, scipy.pi)[0]
		return I1
	def Phi_angle(xi):
		x1 = (1-numpy.power(q,2))*xi_j - xi
		x2 = 2*q*numpy.sqrt(xi*xi_j)
		return numpy.arccos(x1/x2)
	def InteA(xi):
		x1 = 1./(2*scipy.pi)
		x2 = 1/numpy.power((1+xi),2+alpha)
		x3 = Pitch_g(xi)*numpy.sin(2*Phi_angle(xi))
		return x1*x2*x3
	def InteB(xi):
		heav1 = numpy.heaviside(1-q,0)
		x2 = Pitch_f(xi)/numpy.power((1+xi),2+alpha)
		return heav1*x2
	def InteB2(xi):
		x2 = (scipy.pi-Phi_angle(xi))/(scipy.pi*numpy.power((1+xi),2+alpha))
		return x2*Pitch_f(xi)

	I1 = scipy.integrate.quad(InteA,ximin,ximax)[0]
	I2a = scipy.integrate.quad(InteB,0,ximin)[0]
	I2b = scipy.integrate.quad(InteB2, ximin, ximax)[0]

	return I1*numpy.power(I2a+I2b,-1)

###############################################################################################
#######################################  GAUSSIAN  ############################################
##################################### JET STRUCTURE ###########################################
###############################################################################################

def Polarization_Structured_Par(q,alpha,xi_c, b, a):
	epislon = 1 + alpha
	def Theta(xi,phi):
		til_xi = xi + numpy.power(q,2)*xi_c + 2*q*numpy.sqrt(xi*xi_c)*numpy.cos(phi) 
		x1 = numpy.sqrt(1+(til_xi/xi_c))
		return x1	
	def Pitch(xi,phi):
		til_xi = xi + numpy.power(q,2)*xi_c + 2*q*numpy.sqrt(xi*xi_c)*numpy.cos(phi) 
		hat_xi = numpy.power(Theta(xi,phi),-2*b)*til_xi
		x1 = numpy.sqrt(4*hat_xi)/(1+hat_xi)
		return numpy.power(x1,epislon)
	def Doppler(xi,phi):
		return 2*numpy.power(Theta(xi,phi),-b)/(1+numpy.power(Theta(xi,phi),-2*b)*xi)
	def InteA(xi):
		def Inte1(phi):
			x1 = numpy.power(Doppler(xi,phi), 2+alpha)*numpy.power(Theta(xi,phi),b)*numpy.power(Theta(xi,phi),-a)*Pitch(xi,phi)*numpy.cos(2*phi)
			return x1
		I1 = scipy.integrate.quad(Inte1, 0, 2*scipy.pi)[0]
		return I1
	def InteB(xi):
		def Inte1(phi):
			x1 = numpy.power(Doppler(xi,phi), 2+alpha)*numpy.power(Theta(xi,phi),b)*numpy.power(Theta(xi,phi),-a)*Pitch(xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte1, 0, 2*scipy.pi)[0]
		return I1
	xi_max = numpy.power(1+q,2)*xi_c
	I1 = scipy.integrate.quad(InteA, 0, xi_max)[0]
	I2 = scipy.integrate.quad(InteB, 0, xi_max)[0]

	return I1*numpy.power(I2,-1)


def Polarization_Structured_Per(q,alpha,xi_c, b, a):
	epislon = 1 + alpha
	def Theta(xi,phi):
		til_xi = xi + numpy.power(q,2)*xi_c + 2*q*numpy.sqrt(xi*xi_c)*numpy.cos(phi) 
		x1 = numpy.sqrt(1+(til_xi/xi_c))
		return x1	
	def Pitch_g(xi):
		def Inte1(phi):
			til_xi = xi + numpy.power(q,2)*xi_c + 2*q*numpy.sqrt(xi*xi_c)*numpy.cos(phi) 
			hat_xi = numpy.power(Theta(xi,phi),-2*b)*til_xi
			x1 = numpy.power((1-hat_xi)/(1+hat_xi),2)*numpy.power(numpy.cos(phi),2)
			x2 = numpy.power(numpy.sin(phi),2)
			x3 = numpy.power(1-4*hat_xi*numpy.power(numpy.cos(phi),2)*numpy.power(1+hat_xi,-2),(epislon-2)/2)
			return (x2-x1)*x3
		I1 = scipy.integrate.quad(Inte1, 0, scipy.pi)[0]
		return I1
	def Pitch_f(xi):
		def Inte1(phi):
			til_xi = xi + numpy.power(q,2)*xi_c + 2*q*numpy.sqrt(xi*xi_c)*numpy.cos(phi) 
			hat_xi = numpy.power(Theta(xi,phi),-2*b)*til_xi
			x3 = 1-4*hat_xi*numpy.power(numpy.cos(phi),2)*numpy.power(1+hat_xi,-2)
			return numpy.power(x3,epislon/2)
		I1 = scipy.integrate.quad(Inte1, 0, scipy.pi)[0]
		return I1
	def Doppler(xi,phi):
		return 2*numpy.power(Theta(xi,phi),-b)/(1+numpy.power(Theta(xi,phi),-2*b)*xi)
	def InteA(xi):
		def Inte1(phi):
			x1 = numpy.power(Doppler(xi,phi), 2+alpha)*numpy.power(Theta(xi,phi),b)*numpy.power(Theta(xi,phi),-a)*Pitch_g(xi)*numpy.cos(2*phi)
			return x1
		I1 = scipy.integrate.quad(Inte1, 0, 2*scipy.pi)[0]
		return I1
	def InteB(xi):
		def Inte1(phi):
			x1 = numpy.power(Doppler(xi,phi), 2+alpha)*numpy.power(Theta(xi,phi),b)*numpy.power(Theta(xi,phi),-a)*Pitch_f(xi)
			return x1
		I1 = scipy.integrate.quad(Inte1, 0, 2*scipy.pi)[0]
		return I1
	xi_max = numpy.power(1+q,2)*xi_c
	I1 = scipy.integrate.quad(InteA, 0, xi_max)[0]
	I2 = scipy.integrate.quad(InteB, 0, xi_max)[0]

	return I1*numpy.power(I2,-1)


def Polarization_Structured_Tor(q,alpha,xi_c, b, a):
	epislon = 1 + alpha
	def Theta(xi,phi):
		til_xi = xi + numpy.power(q,2)*xi_c + 2*q*numpy.sqrt(xi*xi_c)*numpy.cos(phi) 
		x1 = numpy.sqrt(1+(til_xi/xi_c))
		return x1	
	def Pitch(xi, phi):
		til_xi = xi + numpy.power(q,2)*xi_c + 2*q*numpy.sqrt(xi*xi_c)*numpy.cos(phi) 
		hat_xi = numpy.power(Theta(xi,phi),-2*b)*til_xi
		a = numpy.sqrt(hat_xi/xi_j)/q
		x1 = numpy.power((1-hat_xi)/(1+hat_xi),2)
		x2 = 4*xi/numpy.power(1+hat_xi,2)
		x3 = numpy.power(a+numpy.cos(phi),2)/(1+numpy.power(a,2)+2*a*numpy.cos(phi))
		return numpy.power(x1 + x2*x3, epislon/2)
	def Doppler(xi,phi):
		return 2*numpy.power(Theta(xi,phi),-b)/(1+numpy.power(Theta(xi,phi),-2*b)*xi)
	def InteA(xi):
		def Inte1(phi):
			x1 = numpy.power(Doppler(xi,phi), 2+alpha)*numpy.power(Theta(xi,phi),b)*numpy.power(Theta(xi,phi),-a)*Pitch_g(xi,phi)*numpy.cos(2*phi)
			return x1
		I1 = scipy.integrate.quad(Inte1, 0, 2*scipy.pi)[0]
		return I1
	def InteB(xi):
		def Inte1(phi):
			x1 = numpy.power(Doppler(xi,phi), 2+alpha)*numpy.power(Theta(xi,phi),b)*numpy.power(Theta(xi,phi),-a)*Pitch_f(xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte1, 0, 2*scipy.pi)[0]
		return I1
	xi_max = numpy.power(1+q,2)*xi_c
	I1 = scipy.integrate.quad(InteA, 0, xi_max)[0]
	I2 = scipy.integrate.quad(InteB, 0, xi_max)[0]

	return I1*numpy.power(I2,-1)


###############################################################################################
#######################################  Fluence  ############################################
###############################################################################################

def F_iso(q,alpha,xi_j, Delta):
	epislon = 1 + alpha
	def L_theta(xi):
		if(xi<=xi_j):
			x1 = 1
		else:
			x1 = numpy.exp((numpy.sqrt(xi_j)-numpy.sqrt(xi))/Delta)
		return x1
	def Pitch(xi, phi):
		a = numpy.sqrt(xi/xi_j)/q
		x1 = numpy.power((1-xi)/(1+xi),2)
		x2 = 4*xi/numpy.power(1+xi,2)
		x3 = numpy.power(a+numpy.cos(phi),2)/(1+numpy.power(a,2)+2*a*numpy.cos(phi))
		return numpy.power(x1 + x2*x3, epislon/2)
	def InteA(xi):
		x1 = numpy.power(1+xi, -(2+alpha))*L_theta(xi)
		def Inte1(phi):
			x1 = Pitch(xi,phi)*numpy.cos(2*phi)
			return x1
		I1 = scipy.integrate.quad(Inte1, 0, 2*scipy.pi)[0]
		return x1*I1
	def InteB(xi):
		x1 = numpy.power(1+xi, -(2+alpha))*L_theta(xi)
		def Inte1(phi):
			x1 = Pitch(xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte1, 0, 2*scipy.pi)[0]
		return x1*I1
	def xi_max(q):
		return numpy.power(1+q,2)*xi_j
	I1 = scipy.integrate.quad(InteA, 0, xi_max(q))[0]
	I2 = scipy.integrate.quad(InteB, 0, xi_j)[0]

	return I1*numpy.power(I2,-1)


def Polarization_F_iso(q,alpha,xi_j,Delta):
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1+alpha
	def L_theta(xi):
		if(xi<=xi_j):
			x1 = 1
		else:
			x1 = numpy.exp((numpy.sqrt(xi_j)-numpy.sqrt(xi))/Delta)
		return x1
	def Pitch(xi):
		x1 = numpy.sqrt(4*xi)/(1+xi)
		return numpy.power(x1,epislon)
	def Phi_angle(xi):
		x1 = (1-numpy.power(q,2))*xi_j - xi
		x2 = 2*q*numpy.sqrt(xi*xi_j)
		return numpy.arccos(x1/x2)
	def InteA(xi):
		x1 = 1./(2*scipy.pi)
		x2 = 1/numpy.power((1+xi),2+alpha)*L_theta(xi)
		x3 = Pitch(xi)*numpy.sin(2*Phi_angle(xi))
		return x1*x2*x3
	def InteB(xi):
		heav1 = numpy.heaviside(1-q,0)
		x2 = Pitch(xi)*L_theta(xi)/numpy.power((1+xi),2+alpha)
		return heav1*x2
	def InteB2(xi):
		x2 = L_theta(xi)*(scipy.pi-Phi_angle(xi))/(scipy.pi*numpy.power((1+xi),2+alpha))
		return x2*Pitch(xi)

	I1 = scipy.integrate.quad(InteA,ximin,ximax)[0]
	I2a = scipy.integrate.quad(InteB,0,ximin)[0]
	I2b = scipy.integrate.quad(InteB2, ximin, ximax)[0]

	return -I1*numpy.power(I2a+I2b,-1) 



###############################################################################################
#######################################  With Change  #########################################
###############################################################################################
def run(Toggle_Geom, debug):
	"""Runtime parameters for the code. Specifies geometry, whether or not to show calculations and allows easy expansion of conditions"""
	if(Toggle_Geom==0):
		Geom = "Off_Axis"
	elif(Toggle_Geom==1):		
		Geom = "Spherical"
	else:
		print("Please input  0 for Off Axis Jet or 1 for Spherical Outflow.")
	return Geom, debug

def timebreak_sphe(alpha_s, n, E, theta_j, theta):
	"""Function to calculate the time where regime changes (differs between off axis and cocoon with the value of the dtheta exponential)"""
	global z
	A1=numpy.power(3./(2.*math.pi*mp), 1./3.)
	kk=0.5
	dtheta = theta - theta_j
	return kk *A1 *numpy.power(1+z, 1.) *numpy.power(n, -1./3.) *numpy.power(E, 1./3.) *numpy.power(dtheta, (alpha_s+6.)/3.)*1/day

def timebreak_off(alpha_s, n, E, theta_j, theta):
	"""Function to calculate the time where regime changes (differs between off axis and cocoon with the value of the dtheta exponential)"""
	global z
	A1=numpy.power(3./(2.*math.pi*mp), 1./3.)
	kk=0.5
	dtheta = theta - theta_j
	return kk *A1 *numpy.power(1+z, 1.) *numpy.power(n, -1./3.) *numpy.power(E, 1./3.) *numpy.power(dtheta, 2.)*1/day

def Polarization_Ordered_Off_Changed(alpha_s,n, E, theta_j, theta, alpha, t, debug, Geom):
	""" Time evolution Polarization for a Ordered Field (same as perp)"""
	global z
	beta = 1
	def Gamma_Off():
		if(t/day<=timebreak_off(alpha_s, n, E, theta_j, theta)):
		#Gamma for off axis relativistic
	 		Gamma = 321.1 *numpy.power((1+z/1.022),3/2) *numpy.power(n,-1/2) *pow(E,1/2) *numpy.power(theta_j,-1) *numpy.power(theta-theta_j,3) *numpy.power(t,-3/2) #numpy.power(erg,-1/2) #*numpy.power(cm,-3/2) *numpy.power(math.pi/180.,-2)* numpy.power(day,3/2) 
		else:
		#Gamma for Latexal Expansion
			Gamma = 321.1 *numpy.power((1+z/1.022),1/2) *numpy.power(n,-1/6) *pow(E,1/6) *numpy.power(t,-1/2) #*numpy.power(erg,-1/6) *numpy.power(cm,3/6) * numpy.power(day,1/2)
		return Gamma
	def Gamma_Sphe():
		if(t/day<=timebreak_sphe(alpha_s, n, E, theta_j, theta)):
		#Gamma for Spherical Outflow relativistic
			Gamma = 3.8  *numpy.power((1+z/1.022),3/(alpha_s+8)) *numpy.power(n,-1/(alpha_s+8)) * pow(E,1/(alpha_s+8)) * numpy.power(t,-3/(alpha_s+8)) #* numpy.power(cm, -3/(alpha_s+8)) * numpy.power(day, 3/(alpha_s+8)) * numpy.power(erg, -1/(alpha_s+8))  
		else:
		#Gamma for Spherical Outflow expansion
			Gamma = 2.1 *numpy.power((1+z/1.022),3/(alpha_s+6)) *numpy.power(n,-1/(alpha_s+6)) *numpy.power(beta,-alpha_s/(alpha_s+6))  *pow(E,1/(alpha_s+6)) * numpy.power(t,-3/(alpha_s+6)) 
		return Gamma

	if(Geom == "Off_Axis"):
		Gamma = Gamma_Off()
		t_break = timebreak_off(alpha_s, n, E, theta_j, theta)
	elif(Geom == "Spherical"):
		Gamma = Gamma_Sphe()
		t_break = timebreak_sphe(alpha_s, n, E, theta_j, theta)
	else:
		return "Please, choose Spherical or Off as a string value to the variable Geom."

	xi_j = numpy.power(Gamma*theta_j,2)
	q = theta/theta_j #(could also use our delta theta?)
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1 + alpha

	if(debug==True):
		print("q = {:5.5f}, \tximin = {:5.5e}, \tximax = {:5.5e}".format(q, ximin, ximax))
		print("E = {:5e}, \te = {:5.5f}, \tn = {:5.5f}".format(E, pow(E,1/(alpha_s+8)), numpy.power(n,-1/(alpha_s+8))))
		print("Gamma = {:5.5f}, \tz = {:5.5f}".format(Gamma, numpy.power((1+z/1.022),3/(alpha_s+8))))
		print("t = {:3e}, \tt_break = {:3e}".format(t/day, t_break))

	def Pitch(xi,phi):
		""" Pitch factor calculation"""
		x1 = numpy.power((1-xi)/(1+xi),2) * numpy.power(numpy.cos(phi),2)
		x2 = numpy.power(numpy.sin(phi),2)
		return numpy.power(x1+x2,epislon/2)
	def Pos_Angle(xi,phi):
		"""Polarization angle calculation"""
		x1 = numpy.arctan(((1-xi)/(1+xi))/numpy.tan(phi))
		return phi + x1
	def InteA(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		"""Integration on xi of the dividend"""
		def Inte2(phi):
			"""Integration on phi of the dividend"""
			x1 = Pitch(xi,phi)*numpy.cos(2*Pos_Angle(xi,phi))
			return x1
		I1= scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return div * I1
	def InteB(xi):
		"""Integration on xi of the divisor"""
		div = 1/numpy.power(1+xi,2+alpha)
		def Inte2(phi):
			"""Integration on phi of the divisor"""
			x1 = Pitch(xi,phi)
			return x1
		I1= scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return div * I1

	I1 = scipy.integrate.quad(InteA, 0, ximax)[0]
	I2 = scipy.integrate.quad(InteB, 0, ximax)[0]
	return abs(I1*numpy.power(I2,-1))

def Polarization_Perpendicular_Off_Changed(alpha_s,n, E, theta_j, theta, alpha, t, debug, Geom):
	""" Time evolution Polarization for a Perpendicular Field (same as ord?)"""
	global z
	beta = 1
	def Gamma_Off():
		if(t/day<=timebreak_off(alpha_s, n, E, theta_j, theta)):
		#Gamma for off axis relativistic
	 		Gamma = 321.1 *numpy.power((1+z/1.022),3/2) *numpy.power(n,-1/2) *pow(E,1/2) *numpy.power(theta_j,-1) *numpy.power(theta-theta_j,3) *numpy.power(t,-3/2) #numpy.power(erg,-1/2) #*numpy.power(cm,-3/2) *numpy.power(math.pi/180.,-2)* numpy.power(day,3/2) 
		else:
		#Gamma for Latexal Expansion
			Gamma = 321.1 *numpy.power((1+z/1.022),1/2) *numpy.power(n,-1/6) *pow(E,1/6) *numpy.power(t,-1/2) #*numpy.power(erg,-1/6) *numpy.power(cm,3/6) * numpy.power(day,1/2)
		return Gamma
	def Gamma_Sphe():
		if(t/day<=timebreak_sphe(alpha_s, n, E, theta_j, theta)):
		#Gamma for Spherical Outflow relativistic
			Gamma = 3.8  *numpy.power((1+z/1.022),3/(alpha_s+8)) *numpy.power(n,-1/(alpha_s+8)) * pow(E,1/(alpha_s+8)) * numpy.power(t,-3/(alpha_s+8)) #* numpy.power(cm, -3/(alpha_s+8)) * numpy.power(day, 3/(alpha_s+8)) * numpy.power(erg, -1/(alpha_s+8))  
		else:
		#Gamma for Spherical Outflow expansion
			Gamma = 2.1 *numpy.power((1+z/1.022),3/(alpha_s+6)) *numpy.power(n,-1/(alpha_s+6)) *numpy.power(beta,-alpha_s/(alpha_s+6))  *pow(E,1/(alpha_s+6)) * numpy.power(t,-3/(alpha_s+6)) 
		return Gamma

	if(Geom == "Off_Axis"):
		Gamma = Gamma_Off()
		t_break = timebreak_off(alpha_s, n, E, theta_j, theta)
	elif(Geom == "Spherical"):
		Gamma = Gamma_Sphe()
		t_break = timebreak_sphe(alpha_s, n, E, theta_j, theta)
	else:
		return "Please, choose Spherical or Off as a string value to the variable Geom."

	xi_j = numpy.power(Gamma*theta_j,2)
	q = theta/theta_j #(could also use our delta theta?)
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1 + alpha

	if(debug==True):
		print("q = {:5.5f}, \tximin = {:5.5e}, \tximax = {:5.5e}".format(q, ximin, ximax))
		print("E = {:5e}, \te = {:5.5f}, \tn = {:5.5f}".format(E, pow(E,1/(alpha_s+8)), numpy.power(n,-1/(alpha_s+8))))
		print("Gamma = {:5.5f}, \tz = {:5.5f}".format(Gamma, numpy.power((1+z/1.022),3/(alpha_s+8))))
		print("t = {:3e}, \tt_break = {:3e}".format(t/day, t_break))

	def Pitch_g(xi):
		def Inte1(phi):
			x1 = numpy.power((1-xi)/(1+xi),2)*numpy.power(numpy.cos(phi),2)
			x2 = numpy.power(numpy.sin(phi),2)
			x3 = numpy.power(1-4*xi*numpy.power(numpy.cos(phi),2)*numpy.power(1+xi,-2),(epislon-2)/2)
			return (x2-x1)*x3
		I1 = scipy.integrate.quad(Inte1, 0, scipy.pi)[0]
		return I1
	def Pitch_f(xi):
		def Inte1(phi):
			x3 = 1-4*xi*numpy.power(numpy.cos(phi),2)*numpy.power(1+xi,-2)
			return numpy.power(x3,epislon/2)
		I1 = scipy.integrate.quad(Inte1, 0, scipy.pi)[0]
		return I1
	def Phi_angle(xi):
		x1 = (1-numpy.power(q,2))*xi_j - xi
		x2 = 2*q*numpy.sqrt(xi*xi_j)
		return numpy.arccos(x1/x2)
	def InteA(xi):
		x1 = 1./(2*scipy.pi)
		x2 = 1/numpy.power((1+xi),2+alpha)
		x3 = Pitch_g(xi)*numpy.sin(2*Phi_angle(xi))
		return x1*x2*x3
	def InteB(xi):
		heav1 = numpy.heaviside(1-q,0)
		x2 = Pitch_f(xi)/numpy.power((1+xi),2+alpha)
		return heav1*x2
	def InteB2(xi):
		x2 = (scipy.pi-Phi_angle(xi))/(scipy.pi*numpy.power((1+xi),2+alpha))
		return x2*Pitch_f(xi)

	I1 = scipy.integrate.quad(InteA,ximin,ximax)[0]
	I2a = scipy.integrate.quad(InteB,0,ximin)[0]
	I2b = scipy.integrate.quad(InteB2, ximin, ximax)[0]
	return I1*numpy.power(I2a+I2b,-1)

def Polarization_Parallel_Off_Changed(alpha_s,n, E,theta_j, theta, alpha,t, debug, Geom):
	""" Time evolution Polarization for a Parallel Field"""
	beta = 1
	def Gamma_Off():
		if(t/day<=timebreak_off(alpha_s, n, E, theta_j, theta)):
		#Gamma for off axis relativistic
	 		Gamma = 321.1 *numpy.power((1+z/1.022),3/2) *numpy.power(n,-1/2) *pow(E,1/2) *numpy.power(theta_j,-1) *numpy.power(theta-theta_j,3) *numpy.power(t,-3/2) #numpy.power(erg,-1/2) #*numpy.power(cm,-3/2) *numpy.power(math.pi/180.,-2)* numpy.power(day,3/2) 
		else:
		#Gamma for Latexal Expansion
			Gamma = 321.1 *numpy.power((1+z/1.022),1/2) *numpy.power(n,-1/6) *pow(E,1/6) *numpy.power(t,-1/2) #*numpy.power(erg,-1/6) *numpy.power(cm,3/6) * numpy.power(day,1/2)
		return Gamma
	def Gamma_Sphe():
		if(t/day<=timebreak_sphe(alpha_s, n, E, theta_j, theta)):
		#Gamma for Spherical Outflow relativistic
			Gamma = 3.8  *numpy.power((1+z/1.022),3/(alpha_s+8)) *numpy.power(n,-1/(alpha_s+8)) * pow(E,1/(alpha_s+8)) * numpy.power(t,-3/(alpha_s+8)) #* numpy.power(cm, -3/(alpha_s+8)) * numpy.power(day, 3/(alpha_s+8)) * numpy.power(erg, -1/(alpha_s+8))  
		else:
		#Gamma for Spherical Outflow expansion
			Gamma = 2.1 *numpy.power((1+z/1.022),3/(alpha_s+6)) *numpy.power(n,-1/(alpha_s+6)) *numpy.power(beta,-alpha_s/(alpha_s+6))  *pow(E,1/(alpha_s+6)) * numpy.power(t,-3/(alpha_s+6)) 
		return Gamma

	if(Geom == "Off_Axis"):
		Gamma = Gamma_Off()
		t_break = timebreak_off(alpha_s, n, E, theta_j, theta)
	elif(Geom == "Spherical"):
		Gamma = Gamma_Sphe()
		t_break = timebreak_sphe(alpha_s, n, E, theta_j, theta)
	else:
		return "Please, choose Spherical or Off as a string value to the variable Geom."

	xi_j = numpy.power(Gamma*theta_j,2)
	q = theta/theta_j #(could also use our delta theta?)
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1 + alpha

	if(debug==True):
		print("q = {:5.5f}, \tximin = {:5.5e}, \tximax = {:5.5e}".format(q, ximin, ximax))
		print("E = {:5e}, \te = {:5.5f}, \tn = {:5.5f}".format(E, pow(E,1/(alpha_s+8)), numpy.power(n,-1/(alpha_s+8))))
		print("Gamma = {:5.5f}, \tz = {:5.5f}".format(Gamma, numpy.power((1+z/1.022),3/(alpha_s+8))))
		print("t = {:3e}, \tt_break = {:3e}".format(t/day, t_break))

	def Pitch(xi):
		x1 = numpy.sqrt(4*xi)/(1+xi)
		return numpy.power(x1,epislon)
	def Phi_angle(xi):
		x1 = (1-numpy.power(q,2))*xi_j - xi
		x2 = 2*q*numpy.sqrt(xi*xi_j)
		return numpy.arccos(x1/x2)
	def InteA(xi):
		x1 = 1./(2*scipy.pi)
		x2 = 1/numpy.power((1+xi),2+alpha)
		x3 = Pitch(xi)*numpy.sin(2*Phi_angle(xi))
		return x1*x2*x3
	def InteB(xi):
		x2 = Pitch(xi)/numpy.power((1+xi),2+alpha)
		return x2
	def InteB2(xi):
		x2 = (scipy.pi-Phi_angle(xi))/(scipy.pi*numpy.power((1+xi),2+alpha))
		return x2*Pitch(xi)

	I1 = scipy.integrate.quad(InteA,ximin,ximax)[0]
	I2a = numpy.heaviside(1-q,0)*scipy.integrate.quad(InteB,0,ximin)[0]
	I2b = scipy.integrate.quad(InteB2, ximin, ximax)[0]
	return -I1*numpy.power(I2a+I2b,-1) 

def Polarization_Toroidal_Off_Changed(alpha_s,n, E, theta_j, theta, alpha, t, debug, Geom):
	""" Time evolution Polarization for a Toroidal Field"""
	global z
	beta = 1
	def Gamma_Off():
		if(t/day<=timebreak_off(alpha_s, n, E, theta_j, theta)):
		#Gamma for off axis relativistic
	 		Gamma = 321.1 *numpy.power((1+z/1.022),3/2) *numpy.power(n,-1/2) *pow(E,1/2) *numpy.power(theta_j,-1) *numpy.power(theta-theta_j,3) *numpy.power(t,-3/2) #numpy.power(erg,-1/2) #*numpy.power(cm,-3/2) *numpy.power(math.pi/180.,-2)* numpy.power(day,3/2) 
		else:
		#Gamma for Latexal Expansion
			Gamma = 321.1 *numpy.power((1+z/1.022),1/2) *numpy.power(n,-1/6) *pow(E,1/6) *numpy.power(t,-1/2) #*numpy.power(erg,-1/6) *numpy.power(cm,3/6) * numpy.power(day,1/2)
		return Gamma
	def Gamma_Sphe():
		if(t/day<=timebreak_sphe(alpha_s, n, E, theta_j, theta)):
		#Gamma for Spherical Outflow relativistic
			Gamma = 3.8  *numpy.power((1+z/1.022),3/(alpha_s+8)) *numpy.power(n,-1/(alpha_s+8)) * pow(E,1/(alpha_s+8)) * numpy.power(t,-3/(alpha_s+8)) #* numpy.power(cm, -3/(alpha_s+8)) * numpy.power(day, 3/(alpha_s+8)) * numpy.power(erg, -1/(alpha_s+8))  
		else:
		#Gamma for Spherical Outflow expansion
			Gamma = 2.1 *numpy.power((1+z/1.022),3/(alpha_s+6)) *numpy.power(n,-1/(alpha_s+6)) *numpy.power(beta,-alpha_s/(alpha_s+6))  *pow(E,1/(alpha_s+6)) * numpy.power(t,-3/(alpha_s+6)) 
		return Gamma

	if(Geom == "Off_Axis"):
		Gamma = Gamma_Off()
		t_break = timebreak_off(alpha_s, n, E, theta_j, theta)
	elif(Geom == "Spherical"):
		Gamma = Gamma_Sphe()
		t_break = timebreak_sphe(alpha_s, n, E, theta_j, theta)
	else:
		return "Please, choose Spherical or Off as a string value to the variable Geom."

	xi_j = numpy.power(Gamma*theta_j,2)
	q = theta/theta_j #(could also use our delta theta?)
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1 + alpha

	if(debug==True):
		print("q = {:5.5f}, \tximin = {:5.5e}, \tximax = {:5.5e}".format(q, ximin, ximax))
		print("E = {:5e}, \te = {:5.5f}, \tn = {:5.5f}".format(E, pow(E,1/(alpha_s+8)), numpy.power(n,-1/(alpha_s+8))))
		print("Gamma = {:5.5f}, \tz = {:5.5f}".format(Gamma, numpy.power((1+z/1.022),3/(alpha_s+8))))
		print("t = {:3e}, \tt_break = {:3e}".format(t/day, t_break))

	def Phi_angle(xi):
		x1 = (1-numpy.power(q,2))*xi_j - xi
		x2 = 2*q*numpy.sqrt(xi*xi_j)
		return numpy.arccos(x1/x2)
	def Pos_Angle(xi,phi):
		a = numpy.sqrt(xi/xi_j)/q
		x1 = ((1-xi)/(1+xi))
		x2 = numpy.sin(phi)/(a+numpy.cos(phi))
		arctan = numpy.arctan(x1*x2)
		return phi - arctan
	def Pitch(xi, phi):
		a = numpy.sqrt(xi/xi_j)/q
		x1 = numpy.power((1-xi)/(1+xi),2)
		x2 = 4*xi/numpy.power(1+xi,2)
		x3 = numpy.power(a+numpy.cos(phi),2)/(1+numpy.power(a,2)+2*a*numpy.cos(phi))
		return numpy.power(x1 + x2*x3, epislon/2)
	def InteA(xi):
		heav1 = numpy.heaviside(1-q,0)
		div = 1./numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch(xi,phi)*numpy.cos(2*Pos_Angle(xi,phi))
			return x1
		I1 = scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return heav1*div*I1
	def InteB(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch(xi,phi)*numpy.cos(2*Pos_Angle(xi,phi))
			return x1
		I1 = scipy.integrate.quad(Inte2, Phi_angle(xi), 2*scipy.pi-Phi_angle(xi))[0]
		return div*I1
	def InteC(xi):
		heav1 = numpy.heaviside(1-q,0)
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch(xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return heav1*div*I1
	def InteD(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch(xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte2, Phi_angle(xi), 2*scipy.pi-Phi_angle(xi))[0]
		return div*I1

	I1a = scipy.integrate.quad(InteA, 0, ximin)[0]
	I1b = scipy.integrate.quad(InteB, ximin, ximax, limit=100)[0]
	I2a = scipy.integrate.quad(InteC, 0, ximin)[0]
	I2b = scipy.integrate.quad(InteD, ximin, ximax, limit=100)[0]
	return (I1a + I1b)*numpy.power(I2a + I2b,-1)



