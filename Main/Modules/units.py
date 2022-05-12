# Energy scales
eV=1.
KeV=pow(10,3) * eV
MeV=pow(10,6) * eV
GeV=pow(10,9) * eV
TeV=pow(10,12) * eV
PeV=pow(10,15) * eV
EeV=pow(10,18) * eV
ZeV=pow(10,21) * eV
Hz=eV/(2.418 * pow(10,14))
GHz=10**9*Hz
Kelvin=1/(300 * 38.681686) * eV
Jule=1/1.602176462 * pow(10,19) * eV
erg=pow(10,-7) * Jule


# Mass
me=0.510998902 * MeV
mp=938.271998 * MeV
md=1875.612762 * MeV
mu=931.494013 * MeV
mpion=139.57 * MeV
mdelta=1232 * MeV
mkaon=493.67 * MeV
Mw=80.419 * GeV
Mz=91.1882 * GeV
gram=1/(1.782661731) * pow(10,33) * eV
Kg=pow(10,3) * gram
Msun=1.98844 * pow(10,30) * Kg
Mearth=5.9723 * pow(10,24) * Kg

#Unity of electric charge
alph=7.297352533 * pow(10,-3) # e=(4 Pi alph)^(0.5)
e=8.542 * pow(10,-2) # Gaussian
#e=0.3028  # Lorentz - Heaviside

#Unity of magnetic field
Gauss=5.788381749 *  pow(10,-15) *  MeV * 2 * me/e
Tesla=pow(10,4) * Gauss


#Distance Scales  
fm=1/197.3269602 * 1/MeV 
cm=pow(10,13) * fm
meter=pow(10,2) * cm
inc=0.0254 * meter
Km=pow(10,3) * meter
pc=3.0856775807 * pow(10,16) * meter
Ly=0.9462 * pow(10,16) * meter
kpc=pow(10,3) * pc
Mpc=pow(10,6) * pc
Gpc=pow(10,9) * pc
barn=pow(10,-28) * pow(cm,2)
mbarn=pow(10,-3) * barn

Rsun=6.961 * pow(10,8) * meter
Rearth=6.378140 *  pow(10,6) * meter


#Time Scales


sec=299792458*meter
msec=pow(10,-3)*sec
minute= 60 * sec
hour=60 * minute
day=24 * hour
year=31556925.2 * sec

#flux unity

Jy=pow(10,-23)*erg/pow(cm,2)*1/sec*1/Hz
mJy=pow(10,-3)*Jy
micJy=pow(10,-6)*Jy

#Cross Section
sigmaT=6.65 * pow(10,-25) * pow(cm,2)
sigmaTp=pow((me/mp),2) * sigmaT
sigmapp=4.5 * pow(10,-26) * pow(cm,2)



GN=6.673 * pow(10,-11) * pow(meter,3)/Kg * 1/pow(sec,2)
GF=1.16639 * pow(10,-5)/pow(GeV,2)
s2w=0.23122 #Sin^2 theta_W
alphaz=0.1176 # at Mz scale
