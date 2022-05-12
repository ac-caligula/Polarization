from Modules.dependencies import *

#Ordered field

def Pitch_Ord_Pulse(xi,phi,alpha):
	r"""Calculates the factor associated to the pitch angle in the formula for light aberration
	for an Ordered Magnetic Field in a Pulse.


	Parameters
	----------
	xi : float
		:math:`\xi = (\Gamma\Theta)^2`, the parametrization utilized to reduce Polarization to a
		geometric dependance.
	phi : float 
		:math:`\phi` the azimuthal angle between Jet Axis and the Magnetic Field.
	alpha : float
		Power law polarization parameter, :math:`\alpha`.

	Returns
	-------
	float
		Returns a float object representing the mathematical factor associated to the pitch angle.
	
	Notes
	-----		
	This function is part of the integrand to be calculated for the polarization and 
	it's dependent on two "ghost" variables, :math:`\xi` and :math:`\phi`, who are 
	integrated out.

	.. math:: 
		\Lambda &= \left[\left(\frac{1-\xi}{1+\xi}\right)^2\cos^2(\phi) - \sin^2(\phi)\right]^{\epsilon/2}\\
		\epsilon &= 1 + \alpha

	"""
	epislon = alpha + 1
	x1 = numpy.power((1-xi)/(1+xi),2) * numpy.power(numpy.cos(phi),2)
	x2 = numpy.power(numpy.sin(phi),2)
	return numpy.power(x1+x2,epislon/2)

def Pos_Angle_Ord_Pulse(xi,phi):
	r"""Calculates the Polarization Position Angle for an Ordered Magnetic Field in a Pulse.

	Parameters
	----------
	xi : float
		:math:`\xi = (\Gamma\theta)^2`, the parametrization utilized to reduce Polarization to a
		geometric dependance. 
	phi : float 
		:math:`\phi` the azimuthal angle between Jet Axis and the Magnetic Field.
	
	Returns
	-------
	float
		Return a float object representing the Position Angle in rads.

	Notes
	-----		
	The angle is calculated by utilizing the geometrical considerations of a thin shell emitting
	area. 

	.. math:: 
		\theta_p = \phi + \arctan\left(\frac{1-\xi}{1+\xi}\cot(\phi)\right)

	"""
	x1 = numpy.arctan(((1-xi)/(1+xi))/numpy.tan(phi))
	return phi + x1

def Polarization_Ord_Pulse(alpha,ximax):
	r"""Calculates the Polarization for a Ordered Magnetic field integrated over a single pulse.
	

	Parameters
	----------
	alpha : float 
		This parameter is the value :math:`\alpha` used as index to parametrize the power laws.
	ximax : float
		This parameter is the maximum :math:`\xi` of the system.

	Returns
	-------
	float 
		Return a float object containing the result of the integration

	Notes
	-----
	The Polarization is calculated with a double integration, one over the azimuthal angle
	:math:`\phi` and another on over the :math:`\xi` parameter, which describes a relationship between
	the Lorentz factor and :math:`\theta` (the angle between the jet axis and bulk velocity).

	.. math::

		\frac{\Pi_{ord}}{\Pi_{max}} = \frac{\int\frac{d\xi}{(1+\xi)^{2+\alpha}} \int d\phi \Lambda(\xi,\phi)\cos(2\theta_p)}{\int \frac{d\xi}{(1+\xi)^{2+\alpha}}\int d\phi \Lambda(\xi,\phi)}

	"""	

	def InteA(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch_Ord_Pulse(xi,phi,alpha)*numpy.cos(2*Pos_Angle_Ord_Pulse(xi,phi))
			return x1
		I1= scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return div * I1
	def InteB(xi):
		div = 1/numpy.power(1+xi,2+alpha)
		def Inte2(phi):
			x1 = Pitch_Ord_Pulse(xi,phi,alpha)
			return x1
		I1= scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return div * I1

	I1 = scipy.integrate.quad(InteA, 0, ximax)[0]
	I2 = scipy.integrate.quad(InteB, 0, ximax)[0]
	return abs(I1*numpy.power(I2,-1))


#Toroidal Field Polarization

def Phi_angle(q,xi_j,xi):
	r"""Calculates the angle :math:`\psi`.
	
	Parameters
	----------
	q : float
		Fraction between Observation Angle, :math:`\theta_{obs}`, and Half Opening Angle of 
		the Jet, :math:`\theta_j`.
	xi_j : float
		The parametrization utilized to estabilish a relationship with the Half Opening Angle and
		Lorentz factor, :math:`\xi_j = (\Gamma\theta_j)^2`.
	xi : float
		:math:`\xi = (\Gamma\theta)^2`, the parametrization utilized to reduce Polarization to a
		geometric dependance.

	Returns
	-------
	float
		Return a float object representing the angle in rads.

	Notes
	-----
	The angle :math:`\psi` represents the relationship between the angle parameters of the system.
	And it is calculated by

	.. math::
		\psi = \arccos\left[\frac{(1-q^2)\xi_j - \xi}{2q\sqrt{\xi \xi_j}}\right]

	"""
	x1 = (1-numpy.power(q,2))*xi_j - xi
	x2 = 2*q*numpy.sqrt(xi*xi_j)
	return numpy.arccos(x1/x2)

def Pos_Angle_Toroidal(q, xi_j, xi,phi):
	r"""Calculates the Polarization Position Angle for a Toroidal Magnetic Field in a Pulse.

	Parameters
	----------
	q : float
		Fraction between Observation Angle, :math:`\theta_{obs}`, and Half Opening Angle of 
		the Jet, :math:`\theta_j`.
	xi_j : float
		The parametrization utilized to estabilish a relationship with the Half Opening Angle and
		Lorentz factor, :math:`\xi_j = (\Gamma\theta_j)^2`.
	xi : float
		:math:`\xi = (\Gamma\theta)^2`, the parametrization utilized to reduce Polarization to a
		geometric dependance.
	phi : float 
		:math:`\phi` the azimuthal angle between Jet Axis and the Magnetic Field.
	
	Returns
	-------
	float
		Return a float object representing the Polarization Angle.

	Notes
	-----		
	The angle is calculated by utilizing the geometrical considerations of a thin shell emitting
	area. 

	.. math:: 
		\theta_p &= \phi - \arctan\left(\frac{1-\xi}{1+\xi}\frac{\sin(\phi)}{a+ \cos(\phi)}\right) \\
		a &= \frac{\sqrt{\frac{\xi}{\xi_j}}}{q}

	"""
	a = numpy.sqrt(xi/xi_j)/q
	x1 = ((1-xi)/(1+xi))
	x2 = numpy.sin(phi)/(a+numpy.cos(phi))
	arctan = numpy.arctan(x1*x2)
	return phi - arctan

def Pitch_Toroidal(q, xi_j, alpha, xi, phi):
	r"""Calculates the factor associated to the pitch angle in the formula for light aberration
	for a Toroidal Magnetic Field in a Pulse.

	Parameters
	----------
	q : float
		Fraction between Observation Angle, :math:`\theta_{obs}`, and Half Opening Angle of 
		the Jet, :math:`\theta_j`.
	xi_j : float
		The parametrization utilized to estabilish a relationship with the Half Opening Angle and
		Lorentz factor, :math:`\xi_j = (\Gamma\theta_j)^2`.
	alpha : float 
		This parameter is the value :math:`\alpha` used as index to parametrize the power laws.
	xi : float
		:math:`\xi = (\Gamma\theta)^2`, the parametrization utilized to reduce Polarization to a
		geometric dependance.
	phi : float 
		:math:`\phi` the azimuthal angle between Jet Axis and the Magnetic Field.

	Returns
	-------
	float
		Return a float object during integration
	
	Notes
	-----		
	This function is part of the integrand to be calculated for the polarization and 
	it's dependent on two "ghost" variables, :math:`\xi` and :math:`\phi`, who are 
	integrated out.

	.. math::
		\Lambda &= \bigg[\left(\frac{1-\xi}{1+\xi}^2\right) + \frac{4\xi}{(1+\xi)^2} \frac{(a+\cos(\phi))^2}{(1+a+2a\cos(\phi))} \bigg]\\
		a &= \frac{\sqrt(\frac{\xi}{\xi_j})}{q}
	"""
	epislon = 1 + alpha
	a = numpy.sqrt(xi/xi_j)/q
	x1 = numpy.power((1-xi)/(1+xi),2)
	x2 = 4*xi/numpy.power(1+xi,2)
	x3 = numpy.power(a+numpy.cos(phi),2)/(1+numpy.power(a,2)+2*a*numpy.cos(phi))
	return numpy.power(x1 + x2*x3, epislon/2)

def Polarization_Toroidal(xi_j, q, alpha):
	r"""Calculates the Polarization for a Toroidal Magnetic field integrated over a single pulse.
	
	Parameters
	----------
	xi_j : float
		The parametrization utilized to estabilish a relationship with the Half Opening Angle and
		Lorentz factor, :math:`\xi_j = (\Gamma\theta_j)^2`.
	q : float
		Fraction between Observation Angle, :math:`\theta_{obs}`, and Half Opening Angle of 
		the Jet, :math:`\theta_j`.
	alpha : float 
		This parameter is the value :math:`\alpha` used as index to parametrize the power laws.

	Returns
	-------
	float 
		Return a float object containing the result of the integration

	Notes
	-----
	The Polarization is calculated with a double integration, one over the azimuthal angle
	:math:`\phi` and another on over the :math:`\xi` parameter, which describes a relationship between
	the Lorentz factor and :math:`\theta` (the angle between the jet axis and bulk velocity).

	.. math::

		\frac{\Pi_{ord}}{\Pi_{max}} = H(1-q)\int_0^{\xi_-} \frac{d\xi}{(1+\xi)^{2+\alpha}} \int_0{2\pi} d\phi \Lambda_{tor}(\xi, \phi, a)\cos(2\theta_p) 

	"""
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1 + alpha
	def InteA(xi):
		heav1 = numpy.heaviside(1-q,0)
		div = 1./numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch_Toroidal(q,xi_j,alpha,xi,phi)*numpy.cos(2*Pos_Angle_Toroidal(q, xi_j, xi,phi))
			return x1
		I1 = scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return heav1*div*I1
	def InteB(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch_Toroidal(q,xi_j,alpha,xi,phi)*numpy.cos(2*Pos_Angle_Toroidal(q, xi_j, xi,phi))
			return x1
		I1 = scipy.integrate.quad(Inte2, Phi_angle(q, xi_j, xi), 2*scipy.pi- Phi_angle(q, xi_j, xi))[0]
		return div*I1
	def InteC(xi):
		heav1 = numpy.heaviside(1-q,0)
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch_Toroidal(q,xi_j,alpha,xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte2, 0, 2*scipy.pi)[0]
		return heav1*div*I1
	def InteD(xi):
		div = 1/numpy.power((1+xi),2+alpha)
		def Inte2(phi):
			x1 = Pitch_Toroidal(q,xi_j,alpha,xi,phi)
			return x1
		I1 = scipy.integrate.quad(Inte2, Phi_angle(q, xi_j, xi), 2*scipy.pi-Phi_angle(q, xi_j, xi))[0]
		return div*I1

	I1a = scipy.integrate.quad(InteA, 0, ximin)[0]
	I1b = scipy.integrate.quad(InteB, ximin, ximax)[0]
	I1 = I1a + I1b
	I2a = scipy.integrate.quad(InteC, 0, ximin)[0]
	I2b = scipy.integrate.quad(InteD, ximin, ximax)[0]
	I2 = I2a + I2b
	return I1*numpy.power(I2,-1)


def Pitch_Parallel(xi, alpha):
	r"""Calculates the factor associated to the pitch angle in the formula for light aberration
	for a Parallel Magnetic Field in a Pulse.

	Parameters
	----------
	alpha : float 
		This parameter is the value :math:`\alpha` used as index to parametrize the power laws.
	xi : float
		:math:`\xi = (\Gamma\theta)^2`, the parametrization utilized to reduce Polarization to a
		geometric dependance.

	Returns
	-------
	float
		Return a float object during integration
	
	Notes
	-----		
	This function is part of the integrand to be calculated for the polarization and 
	it's dependent on two "ghost" variables, :math:`\xi` and :math:`\phi`, who are 
	integrated out.

	.. math::
		\Lambda &= \bigg[\frac{\sqrt{4\xi}}{1+\xi}\bigg]^\epsilon
	"""
	epislon = 1+alpha
	x1 = numpy.sqrt(4*xi)/(1+xi)
	return numpy.power(x1,epislon)


def Polarization_Parallel(xi_j,q, alpha):
	r"""Calculates the Polarization for a Parallel Magnetic field integrated over a single pulse.
	
	Parameters
	----------
	xi_j : float
		The parametrization utilized to estabilish a relationship with the Half Opening Angle and
		Lorentz factor, :math:`\xi_j = (\Gamma\theta_j)^2`.
	q : float
		Fraction between Observation Angle, :math:`\theta_{obs}`, and Half Opening Angle of 
		the Jet, :math:`\theta_j`.
	alpha : float 
		This parameter is the value :math:`\alpha` used as index to parametrize the power laws.

	Returns
	-------
	float 
		Return a float object containing the result of the integration

	Notes
	-----
	The Polarization is calculated with a double integration, one over the azimuthal angle
	:math:`\phi` and another on over the :math:`\xi` parameter, which describes a relationship between
	the Lorentz factor and :math:`\theta` (the angle between the jet axis and bulk velocity).

	.. math::

		\frac{\Pi_{ord}}{\Pi_{max}} = \frac{\frac{1}{2\pi}\int_{\epsilon_-}^{\epsilon_+}\frac{d\xi\Lambda(\xi)\sin(\psi(\xi))}{(1+\xi)^{2+alpha}}}{H(1-q)\int_{0}^{\epsilon_-}\frac{d\xi\Lambda(\xi)}{(1+\xi)^{2+alpha}}+\int_{\epsilon_-}{\epsilon_+}d\xi \Lambda(\xi)\frac{\pi-\psi(\xi)}{\pi(1+\xi)^{2+\alpha}}}

	"""
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1+alpha
	def InteA(xi):
		x1 = 1./(2*scipy.pi)
		x2 = 1/numpy.power((1+xi),2+alpha)
		x3 = Pitch_Parallel(xi)*numpy.sin(2*Phi_angle(xi))
		return x1*x2*x3
	def InteB(xi):
		heav1 = numpy.heaviside(1-q,0)
		x2 = Pitch_Parallel(xi)/numpy.power((1+xi),2+alpha)
		return heav1*x2
	def InteB2(xi):
		x2 = (scipy.pi-Phi_angle(xi))/(scipy.pi*numpy.power((1+xi),2+alpha))
		return x2*Pitch(xi)

	I1 = scipy.integrate.quad(InteA,ximin,ximax)[0]
	I2a = scipy.integrate.quad(InteB,0,ximin)[0]
	I2b = scipy.integrate.quad(InteB2, ximin, ximax)[0]

	return -I1*numpy.power(I2a+I2b,-1) 


def Polarization_Perpendicular(xi_j,q, alpha):
	ximin = numpy.power(1-q,2)*xi_j
	ximax = numpy.power(1+q,2)*xi_j
	epislon = 1+alpha
	def Pitch_g(xi):
		def Inte1(phi):
			x1 = numpy.power((1-xi)/(1+xi),2)*numpy.power(numpy.cos(phi),2)
			x2 = numpy.power(numpy.sin(phi),2)
			x3 = numpy.power(1-4*xi*numpy.power(numpy.cos(phi),2)*numpy.power(1+xi,-2),(2-epislon)/2)
			return (x1-x2)/x3
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


