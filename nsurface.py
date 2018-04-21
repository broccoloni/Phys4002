from math import pi, exp, sqrt, log, sinh, cosh, cos, acos
import matplotlib.pyplot as plt
from numpy import *

#constants
q = 1.60217646 * 10 ** -19		# Coulomb		(Charge of an electron)
g = 1.                          	# no units      	(Neutrino degeneracy)
hbar = 6.582119514*10**-22      	# MeV s         	(Planck)
c = 2.99792458*10**23           	# fm / s        	(Speed of light)
Gf = 1.16637*10**-11*(hbar*c)**3        # MeV * fm**3           (Fermi constant)
me = 0.511                              # MeV                   (Mass of an electron)
mn = 939.5                              # MeV                   (Mass of a neutron)
sigma0 = 4*Gf**2*me**2/(pi*hbar**4*c**4)# fm**2                 Note: c**4 comes from mass
gA = 1.26                               # no units              (degeneracy)
tri = 1.293                             # MeV                   (neutron proton mass difference)
sinTw = sqrt(0.23)			# no units		(Weinberg angle)
Cv = 1/2.+2*sinTw**2   			# Note: for mu and tao neutrinos Cv = -1/2.+2*sinTw**2 and Ca = -1/2.
Ca = 1/2.               		# Note: for antineutrinos, Ca -> -Ca
C1 = Ca + Cv				# for convenience
C2 = Cv - Ca				# for convenience
Na = 6.0221409*10**23			# particles/mol		(Avogadros number)
k = 1.					# unitless (MeV/MeV)	(Boltzmann constant)
G = 6.67408 * 10 **-11			# gravitational constant (m^3 kg^-1 s^-2)
conv = (G*100**3/1000)/(c*10**-13)**2   # conversion factor from grams to cm (ie. units are cm/g)
M_s = 1.98855*10**30			# Solar mass in kg
Ms = (M_s * 1000) * conv /10.**5	# Solar mass in km
Mbh = 3 * Ms				# mass of the black hole for this simulation in km (3 solar masses)
asbh = 0.8				# spin of the black hole (unitless)
abh = asbh * Mbh			# spin of the black hole (km)

#didn't really use these in np and nn functions as we approximated them to be one
#Mp = 1.00797				# g/mol			(Molar mass of protons)
#Mn = 1.00866				# g/mol			(Molar mass of neutrons)

#function that makes plotting a single function easier
def oneplot(xmin,xmax,xstep,f,xlabel,ylabel,fname,xleg,yleg,Tinc,datacollect,dataname):	
	xs = []
	ys = []
	data = open(dataname, 'w')							
	for x in range(xmin,xmax,xstep):
		xs.append(x/1.)
		if Tinc == 'n':	
			ys.append(f(x))	
		else:
			ys.append(f(x,Tinc))
		if datacollect == 'y':
			data.write(str(x)+' '+str(f(x))+'\n')
	data.close()
	plt.plot(xs, ys, label=fname)
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.legend(bbox_to_anchor=(xleg,yleg),loc=2, borderaxespad=0)
	plt.show()

#function that makes plotting two functions together easier
def twoplot(xmin,xmax,xstep,f1,f2,xlabel,ylabel,f1name,f2name,xleg,yleg,Tinc,datacollect,dataname):
	xs = []
	y1s = []
	y2s = []
	data = open(dataname, 'w')
	for x in range(xmin,xmax,xstep):
		xs.append(x/1.)
		if Tinc == 'n':
			y1s.append(f1(x))
			y2s.append(f2(x))
		else: 
			y1s.append(f1(x,Tinc))
			y2s.append(f2(x,Tinc))	
		if datacollect == 'y':
			data.write(str(x)+' '+str(f1(x))+' '+str(f2(x))+'\n')
	data.close()
	plt.plot(xs, y1s, label=f1name)
	plt.plot(xs, y2s, label=f2name)
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.legend(bbox_to_anchor=(xleg,yleg),loc=2, borderaxespad=0)
	plt.show()

######################################## 
#Fermi-Dirac flux       equation 5 from Liliana's paper
def fdflux(E,T):
	E = E/1.
	a = g*c/(2*pi**2*(hbar*c)**3)   # MeV**-3 * fm**-2 * s**-1
	phi = a*E**2/(exp(E/T)+1)       # MeV**-1 * fm**-2 * s**-1      (Flux)
	return phi

#oneplot(0,101,2,fdflux,'Energy [MeV]','Flux [MeV^-1 * fm^-2 * s^-1]','fd flux',0.75,0.95,10,'n','fd_flux.txt') 
         
         
########################################
#Cross Sections 

# Note: for antineutinos, gA -> -gA (but this doesn't affect these equations                        
# Note: these equations assume mass is multiplied by c**2

#Charged reaction for electron neutrinos (eqn 6), cross section of eqn 12 and 13:
def vn_pe(E):
	E = E/1.					# Note: E not possible on [-1.804,-0.782] eg [-me-tri,me-tri], due to the denominator, x
	x = E + tri 					
	Wm1 = 1. + 1.1 * E / mn                               #Eq. 13
	b = sigma0 * (1+3 * gA**2)/(4 * me**2)
	sigma = x**2 * sqrt(1 - (me/x)**2) * Wm1 * b		#Eq. 12
	return sigma 

#Charged reaction for electron antineutrinos (eqn 7), cross section of eqn 14 and 15
def vp_en(E):						# Note: E not possible on [0.782,1.804] eg [tri-me, tri+me]
	E = E/1.
        x = E - tri
        Wm2 - 1. - 7.1 * E/mn					#Eq. 14
	b = sigma0 * (1+3 * gA**2)/(4 * me**2)
	sigma = x**2 * sqrt(1 - (me/x)**2) * Wm2 * b		#Eq. 13
	return sigma

#Neutral current reactions for electron neutrinos, eqns 16 - 23

#Scattering of e- neutrinos/antineutrinos off of protons, eqn 16
def vp_vp(E):
	b = sigma0 / (4 * me**2)
        return b * ((Cv - 1) ** 2 + 3 * gA ** 2 * (Ca - 1) **2) * E**2 			#eq. 16 for neutrinos and antineutrinos

#Scattering of e- neutrinos/antineutrinos off of neutrons, eqn 17
def vn_vn(E):
	return sigma0 * (1 + 3 * gA ** 2) * E ** 2 / (16 * me ** 2)			#eq. 17 for neutrinos and antineutrinos

#Scattering of e- neutrinos off of electrons, eqn 18 after being analyticaly integrated
def ve_ve(E):
	x = me + 2*E
	sigma = sigma0*2*E**2/(8*me*x)*(C1**2+2*C2/x**2*(2*E**2*C2/3.-Ca*me*x))    	# eq. 18 aafter integration
	return sigma							

#Scattering e- antineutrinos off of electrons, eqn 18 where C1 -> C2 and C2 - C1 to account for antineutrinos
def ve_ve_A(E):
	x = me + 2*E
	sigma = sigma0*2*E**2/(8*me*x)*(C2**2+2*C1/x**2*(2*E**2*C1/3.+Ca*me*x))    	# eq. 18 for antineutrinos
	return sigma	

#Average energy required for eqns 20 - 23. Obtained by averaging over fermi-dirac flux
#The analytic numbers come from using Wolfram Alpha
def avgE(T):                               
	avg = 5.6822*T/1.80309							# like eq. 4 but with E replacing xsctn. 5.6822 comes from top, all else cancels but T
	return avg

#e- neutrinos annihilating with e- antineutrinos, eqn 11 for neutrinos
def vv_ee(E,T):
	K = (1+4*sinTw**2+8*sinTw**4)/(6*pi)				       	# eq. 22	
	sigma = 4/3.*K*sigma0*E*avgE(T)/me**2 				       	# eq. 20
	return sigma								# Note: me**2 in denom. is not in the paper
									        # otherwise xsctn has units of E^2 not area

#e- antineutrinos annihilating with e- neutrinos, eqn 11 for antineutrinos
def vv_ee_A(E,T):
	K = (1+4*sinTw**2+8*sinTw**4)/(6*pi)                                    # eq. 23   
	sigma = 4/3.*K*sigma0*E*avgE(T)/me**2                                   # eq. 21
	return sigma								# Note: me**2 in denom. is not in the paper
										# otherwise xsctn has units of E^2 not area


#some plots with the plotting functions
#twoplot(0,101,1,vp_vp,vn_vn,'Energy [Mev]','Cross Section [fm^2]','Eq. 16','Eq. 17',0.05,0.95,'n','n','ntrl_xsctn.txt')
#oneplot(0,101,1,ve_ve,'Energy [MeV]','Cross Section [fm^2]','Eq. 18',0.05,0.95,'n','n','ve_ve')
#twoplot(0,101,1,vv_ee,vv_ee_A,'Energy [MeV]','Cross Section [fm^2]','Eq. 20','Eq. 21',0.05,0.95,10,'n','vv')

###########################################
#Average cross secions found by averaging the cross sections over the fermi dirac distribution
# Top of equatiion 4 as bottom can be evaluated analytically with Wolfram Alpha

#The data was output for plotting and comparison with Lilianas results into the file:
#data= open('avgxsctns.txt','w')

#For the xsections that require numerical integration this function calculates it (top of eqn 4)
#This function integrates a function numerically with the fermi dirac flux between points A and B with N steps at temperature T
#It uses the trapezoid method of numerical integration
def int_w_fdflux(f,A,B,N,T):
	h = (B-A)/(float(N-1))
	mysum = (f(A)*fdflux(A,T) + f(B)*fdflux(B,T))/2.
	for i in range(2,N):
		mysum += f(A+(i-1)*h)*fdflux(A+(i-1)*h,T)
	return h*mysum

#This function calculates the average cross section numerically putting together the int_w_fdflux function which is the
#top of eqn 4, and the analytic results from Wolfram Alpha
def avgxsection(f,A,B,N,T): 					# equation 4
	a = g*c/(2*pi**2*(hbar*c)**3)   			# MeV**-3 * fm**-2 * s**-1	(from fd flux)
	rm_3 = 1.202						# riemann zeta function (3)
	return int_w_fdflux(f,A,B,N,T)/(1.5*rm_3*a*T**3)	# numerical top of eqn 4 divided by the analytic results of the bottom


#Average cross section for the e- neutrino charged reaction
#This had to be evaluated numerically
def avg_vn_pe(T):
	T = T/1.
	A = 0
	B = mn/10.						#Stopped integrating here as the weak magnetism effects become unrealistic				
	N = 200
	avg = avgxsection(vn_pe,A,B,N,T)
	return avg

#Average cross section for the e- antineutrino charged reaction
#This had to be evaluated numerically
def avg_vp_en(T):
	T = T/1.
	A = tri+me+0.0001					#lowest possible E for this cross section
	B = mn/7.1                                              #Stopped integrating here as the weak magnetism effects become unrealistic 
	N = 200
	avg = avgxsection(vp_en,A,B,N,T)
	return avg

#Average cross section for electron neutrinos scattering off of electrons
#This had to be evaluated numberically
def avg_ve_ve(T):
	T = T/1.
	A = 0
	B = mn/10.
	N = 200
	avg = avgxsection(ve_ve,A,B,N,T)
	return avg

#Average cross section for electron antineutrinos scattering off of electrons
#This had to be evaluated numberically
def avg_ve_ve_A(T):
	T = T/1.
	A = 0
	B = mn/10.
	N = 200
	avg = avgxsection(ve_ve_A,A,B,N,T)
	return avg

#Average cross section for electron neutrinos/antineutrinos scattering off of protons
#This was evaluated analytically using Wolfram Alpha
def avg_vp_vp(T):
	T = T/1.
	rm_3 = 1.202                                            # riemann zeta function (3)
	b = sigma0*((Cv-1)**2+3*gA**2*(Ca-1)**2)/(4*me**2)	# constants pulled out of integration
	integ = 23.3309                                         # result from analytic integral to inf.
	return T**2*b*integ/(1.5*rm_3)                          # a's from fd flux cancel and T**3 cancels

#Average cross section for electron neutrinos/antineutrinos scattering off of neutrons
#This was evaluated analytically using Wolfram Alpha
def avg_vn_vn(T):
	T = T/1.
	rm_3 = 1.202                                            # riemann zeta function (3)
	b = sigma0*(1+3*gA**2)/(16*me**2)			# constants pulled out of integration
	integ = 23.3309                                         # result from analytic integral to inf. (fd flux * E**2)
	return T**2*b*integ/(1.5*rm_3)                          # a's from fd flux cancel and T**3 cancels

#Average cross section for electron neutrinos annihilating with electron antineutrinos
#This was evaluated analytically myself
def avg_vv_ee(T):
	T = T/1.
	K = (1+4*sinTw**2+8*sinTw**4)/(6*pi)                    # eq. 22
	b = 4/3.*K*sigma0*avgE(T)/me**2				# constants from eq. 20
	return b*avgE(T)					# a's from fd flux cancel and all but one T cancels

#Average cross section for electron neutrinos annihilating with electron antineutrinos
#This was evaluated analytically myself
def avg_vv_ee_A(T):
	T = T/1.
	K = (1+4*sinTw**2+8*sinTw**4)/(6*pi)                    # eq. 23
	b = 4/3.*K*sigma0*avgE(T)/me**2                         # constants from eq. 20
	return b*avgE(T)	  		        	# a's from fd flux cancel and all but one T cancels


#Plot of all avg cross sections if un - commented
'''
avgxsctns = [avg_vn_pe, avg_vp_en, avg_vp_vp, avg_vn_vn, avg_ve_ve, avg_ve_ve_A, avg_vv_ee, avg_vv_ee_A]
#xsctnnames = ['vn CC','vp CC','ve','ve A','vp','vn','vv','vv A']
#xsctnnames = ['ve + n -> p + e-','ave + p -> e+ + n','ve + e- -> ve + e-','ave + e- -> ave + e-','ve + p -> ve + p','ave + n -> ave + n','ve + ave -> e- + e+','ave + ve -> e+ + e-']
xsctnnames = ['Equation (1)', 'Equation (2)','Equation (3)','Equation (4)','Equation (5)','Equation (5) Antineutrinos','Equation (6)','Equation (6) Antineutrinos']
Ts = [[],[],[],[],[],[],[],[]]
ys = [[],[],[],[],[],[],[],[]]


for i in range(0,8):
	for T in range(1,21,1):
		Ts[i].append(T/1.)
		ys[i].append(avgxsctns[i](T))
	plt.plot(Ts[i],ys[i],label = xsctnnames[i])		

#outputs the data to a file in variable data if un - commented
#for i in range(Trange):
#	data.write(str(Ts[0][i])+' '+str(ys[0][i])+' '+str(ys[1][i])+' '+str(ys[4][i])+' '+str(ys[5][i])+' '+str(ys[2][i])+' '+str(ys[3][i])+' '+str(ys[6][i])+' '+str(ys[7][i])+'\n')
#data.close()


plt.ylabel('Cross Section [fm^2]', fontsize = 12)
plt.xlabel('Temperature [MeV]', fontsize = 12)
plt.grid(True)
plt.tick_params(labelsize = 12)
plt.xticks(arange(0, 21, 10.0))
plt.yticks(arange(0*10**-14, 6*10**-14, 1*10**-14))
plt.legend(bbox_to_anchor=(0.05,0.95),loc=2, borderaxespad=0)
plt.show()
'''

################################################
# Number densities 

#proton and neutron number densities found in paragraph before equation 24

#number density for protons. It's only really a function of Ye and density
def np(T,Ye,density):							# T is for convenience in mean free path function
	np = density*Na*Ye/((10**13)**3)				# division by (10^13)^3 is converting /cm^3 to /fm^3
	return np

#number density for neutrons. It's only really a function of Ye and density
def nn(T,Ye,density):							# T is for convenience in mean free path
	nn =density*Na*(1-Ye)/((10**13)**3)				# division by (10^13)^3 is converting /cm^3 to /fm^3
	return nn

# Electron and Positron Number Density calculations:
#Process:
#The way we find the e- and e+ number densities is by finding the mu value that results in 0 = np + ne+ - ne-
#The ne+ and ne- values are found from integrating phase space numerically and so with the contraint mu_e- + mu_e+ = 0
#we can use bisection with the y = np + ne+ - ne- equation to find the mu that corresponds to a zero. 
#The e- and e+ number densities can then be found from plugging the found mu value into the phase space integral and 
#then plugging that value with the already known value of np into 0 = np + ne+ - ne-

#This function does numerical integration for a function it's given between bounds A and B with N steps
#Its numerical integration specific to the parameters of the number density integral and inputs values of mu and T and integrates over E
def num_dens_int(f,A,B,N,mu,T):
	A = A/1.
	B = B/1.
	h = (B-A)/float(N-1)
	mysum = (f(mu,T,A) + f(mu,T,B))/2.
	for i in range(2,N):
		mysum += f(mu,T,A+(i-1)*h)
	return h*mysum

#This function is the equation y =  ne- - ne+
#The number density as found from phase space is plugged into the ne+ and ne- values
#and the function has been rearranged analytically to get rid of overflow
def num_dens_f(mu,T,E):
	x = exp(-mu/T)
	w = exp(-E/T)
	bottom = 2 + 2*x + x**2 + x + x*w + x*w**2
	y = sqrt(E**2 - me**2)*E*2*(1-x**2)/bottom
	return y

#This function takes the previous function and numerically integrates it
#It then takes that value and subtracts the corresponding np value
#The value y that it gets is then y = ne- - ne+ - np of which we want to obtain the 0
def num_density(Emax,Nmax,mu,T,Ye,density):
	b = 1./(pi**2*(hbar*c)**3)						# b is constants in front of num dens fns that need to be integrated
	y = b * num_dens_int(num_dens_f,me,Emax,Nmax,mu,T) - np(T,Ye,density)	# Note: E cant be less than me because of sqrt in num dens fns
	return y						

#This is the zero finding algorithm which is applied to the equation y = np + ne+ - ne- to obtain the correct
#chemical potentials for the electron and positron
#The zero finding algorithm is a bisection method. ie., two points surrounding the zero are evaluated as well
#as thebisector of those endpoints. The bisector then replaces one of the enpoints in a way so that the two
#endpoints always surround the zero.
def bisection(f,T,Emax,xminus,xplus,Nmax,eps,Ye,density):
	T = T/1.
	found = False                                           
	for i in range(Nmax):					
		x = 0.5*(xplus+xminus)	
		fx = f(Emax,Nmax,x,T,Ye,density)				#reusing this value reduces the number of times the computer evaluates it
		if abs(fx) < eps:						#thereby saving computation time
			found = True
			break
		if f(Emax,Nmax,xplus,T,Ye,density)*fx > 0: 			#Essentially, if sign(f(Emax,Nmax,xplus,T,Ye,density)) == sign(fx):
			xplus = x						#If the bound to the left of the zero * the bisector is positive, the bisector becomes 
		else: 								#the left bound
			xminus = x						#Otherwise the bisector becomes the right bound
	if found == True:							#if a zero is found to a given precision, we return that value	
		return  x
	else: 									#Otherwise we return 0 which will propogate through the algorithm and we will be
		return 0							#able to tell


#This is the number density as found from integrating phase space
#It is used to evaluate the ne- value after the chemical potential is found
def num_dens_single(mu,T,E):
        x = exp(-mu/T)
        w = exp(-E/T)                                   # number density function 1 without constants
        y1 = sqrt(E**2 - me**2)*E*w/(x + w)             # Sqrt(E^2 - me^2)*E * fd-distribution (not fd flux)
        return y1

#This function uses the previous functions to obtain the electron number density
def n_el(T,Ye,density):	
	Emax = 5000							#These 4 values are used so that we can change the precision of the numerical 
	Nmax = 5000							#Integration when doing the bisection method to save computation time
	E = 5000							#We can then change the precision after the chemical potential is found
	N = 5000
	eps = 10**-9
	mu = bisection(num_density,T,E,0.,500*me,N,eps,Ye,density)
	b = 1./(pi**2*(hbar*c)**3)		                        # when the mu is found its put back into num_dens_int to
	n_el = b*num_dens_int(num_dens_single,me,Emax,Nmax,mu,T)  	# find the number density of electrons
	return n_el														

#Number density for positrons
#We do not use this value as neutrinos don't interact with the positrons in any of our reactions
#It still could be useful to have a function for it in case we do want to know the value
def n_pos(T,Ye,density):
	n_pos =  n_el(T,Ye,density) - np(T,Ye,density)				# uses np = ne- - ne+ to get n_pos ( or ne+ )
	return n_pos


#This function evaluates the neutrino/antineutrino number density
#It is found analytically through integrating the phase space integral.
#This can be done as we are approximating the neutrinos to be massless and therefore they have 0 chemical potential as well
#We obtained the analytic results using Wolfram Alpha
def n_neut(T,Ye,density):					# b has 1/2. factor because g_v = 1 because neutrinos have only 1 spin state
	b = 1./(2*pi**2*(hbar*c)**3)				# integral of y = E**2/(exp((E-mu)/T)+1) from 0 to inf with constants
	y = b * 3 * 1.202 * T**3 / 2.				# different than f1 or f2 bc mass = 0 (approx) for neutrinos
	return y						# so relativistic effects are different
								# 1.202 = riemann_zeta(3)


#Prints out number density values for T = 10MeV, Ye = 0.3, dens = 10^11 g/cm^3
'''
numdens = [nn,np,n_el, n_neut]
numdensnames = ['Neutron','Proton','Electron','Neutrino']
for i in range(4):
	print (numdensnames[i], numdens[i](10, 0.3, 10**11))
'''

#################################
# Neutrino mean free path

#We get the mfp by computing 1/sum(number density * average cross section) where we sum over all the reactions for each neutrino

#Ordered list of number density functions to iterate through for the electron neutrino mfp
ordered_ns1 = [nn,np,nn,n_el,n_neut]
#Ordered list of average cross sections to iterate through for the electron neutrino mfp
ordered_avg_xsctns1 = [avg_vn_pe, avg_vp_vp, avg_vn_vn, avg_ve_ve, avg_vv_ee]

#Ordered list of number density functions to iterate through for the electron antineutrino mfp
ordered_ns2 = [np,np,nn,n_el,n_neut]                                         
#Ordered list of average cross sections to iterate through for the electron antineutrino mfp
ordered_avg_xsctns2 = [avg_vp_en, avg_vp_vp, avg_vn_vn, avg_ve_ve_A, avg_vv_ee_A]   

#This function calculates the electron neutino mfp
def MFP_ele_neut(T,Ye,density):
	l_el_neut = 0.
	for n in range(0,5):									
		l_el_neut += (ordered_ns1[n](T,Ye,density)) * (ordered_avg_xsctns1[n](T))
	return 1./(l_el_neut*10**18)								# eq. 3 for electron neut.
												# /10**18 converts fm to km

#This function calculates the electron antineutrino mfp
def MFP_ele_antineut(T,Ye,density):
	l_el_antineut = 0.
	for n in range(0,5):								
		l_el_antineut += (ordered_ns2[n](T,Ye,density)) * (ordered_avg_xsctns2[n](T))
	return 1./(l_el_antineut*10**18)						# /10**18 converts fm to km   

#Test cases using data points from the simulation data
#print MFP_ele_antineut(11.776653,0.318951,4.406495E+11)
#print MFP_ele_antineut(0.786770,0.486975,1.062447E+05)

####################################
# Finding mean free path from simulation data 

#Various data files to write to
#simul_data = loadtxt('M3a0.8t20.dat')				# opens simulation data which is in the format of 
#test_data = open('simul_all2.dat','w')				# r [km], z[km] , T[MeV], density [g/cm^3], Ye


# sepereating M3a0.8t20.dat file at each r value as required by gnuplot as the data doesn't have 
# that seperation inherently
'''
test_data = open("M3a0.8t20s.dat",'w')
r = simul_data[:,0]
z = simul_data[:,1]
T = simul_data[:,2]
dens = simul_data[:,3]
Ye = simul_data[:,4]

datasize = size(r)
for i in range(datasize):
	if i%500 == 0:
		print (i)
	if i < 482960:
		if r[i] == r[i+1]:
			test_data.write(str(r[i])+' '+str(z[i])+' '+str(T[i])+' '+str(dens[i])+' '+str(Ye[i])+'\n')
		else:
			test_data.write(str(r[i])+' '+str(z[i])+' '+str(T[i])+' '+str(dens[i])+' '+str(Ye[i])+'\n\n')
	else:
		test_data.write(str(r[i])+' '+str(z[i])+' '+str(T[i])+' '+str(dens[i])+' '+str(Ye[i])+'\n')
'''



# Trying to reduce the computation time
# writing to file not using numpy - didnt turn out to be any faster
#testing density at constant radius and varying height
'''
counter = 0
xs = []
ys = []
ys2 = []
ys3 = []
ys4 = []
ys5 = []
ys6 = []
ys7 = []

for lines in simul_data:
	data = lines.split()
	r = data[0]
	z = data[1]
	T = data[2]
	dens = data[3]
	Ye = data[4]
	if float(r) == 300.00:
		xs.append(z)
		ys.append(dens)
	if float(r) == 320.00:
		ys2.append(dens)
	if float(r) == 340.00:
		ys3.append(dens)
	if float(r) == 360.00:
		ys4.append(dens)
	if float(r) == 380.00:
		ys5.append(dens)
	if float(r) == 400.00:
		ys6.append(dens)
	if float(r) == 420.00:
		ys7.append(dens)	
	print (counter)
	counter += 1

for i in range(len(xs)):
	outfile.write(str(xs[i])+' '+str(ys[i])+' '+str(ys2[i])+' '+str(ys3[i])+' '+str(ys4[i])+' '+str(ys5[i])+' '+str(ys6[i])+' '+str(ys7[i])+'\n')

plt.plot(xs,ys, label = 'r = 300')
plt.plot(xs,ys2, label = 'r = 320')
plt.plot(xs,ys3, label = 'r = 340')
plt.plot(xs,ys4, label = 'r = 360')
plt.plot(xs,ys5, label = 'r = 380')
plt.plot(xs,ys6, label = 'r = 400')
plt.plot(xs,ys7, label = 'r = 420')
plt.ylabel('Density [g/cm^3]')
plt.xlabel('Height [km]')
plt.legend(bbox_to_anchor=(0.8,0.95),loc=2, borderaxespad=0)
plt.show()
outfile.close()
'''


# using numpy to generate mfp datafile 

'''
r = simul_data[:,0]						#
z = simul_data[:,1]						#
T = simul_data[:,2]						# picks out each column
dens = simul_data[:,3]						#
Ye = simul_data[:,4]						#
				

#datasize = size(r[r <= 19])					# only using r values of 18.5 km if uncommented
								# have to comment the other datasize value then
								
#T = T[:datasize]						#
#Ye = Ye[:datasize]						# makes other values the length of number of r = 18.5 values
#dens = dens[:datasize]						# if above is uncommented uncomment these
#z = z[:datasize]						#

datasize = size(r)
print (datasize)

#Iterates through all data values and write them back into a new file but includes the mfp for e- neutrinos and antineutrinos
#seperates the values after each new r value so we can plot using gnuplot
for i in range(datasize):
	if i%500 == 0:
		print (i)
	ri = r[i]
	zi = z[i]
	Ti = T[i]
	Yei = Ye[i] 
	densi = dens[i]
	if ri == r[i+1]:
		test_data.write(str(ri)+' '+str(zi)+' '+str(Ti)+' '+str(densi)+' '+str(Yei)+' '+str(MFP_ele_neut(Ti,Yei,densi))+' '+str(MFP_ele_antineut(Ti,Yei,densi))+'\n')
	else:
		test_data.write(str(ri)+' '+str(zi)+' '+str(Ti)+' '+str(densi)+' '+str(Yei)+' '+str(MFP_ele_neut(Ti,Yei,densi))+' '+str(MFP_ele_antineut(Ti,Yei,densi))+'\n\n')

	# prints in order of r[km], z[km], T[MeV], dnesity[g/cm^3], Ye, mfp electron neutrinos[km], mfp electron antineutrinos[km]

#test_data.close()
'''

# plotting number density data file for single r value only
# for when values werent matching Liliana's
# To plot with multiple r values we need to plot wth gnuplot to make the heatmaps

'''
numdata = loadtxt('num_dens.dat')
z = numdata[:,0]
n_e = numdata[:,2]
n_p = numdata[:,3]
n_n = numdata[:,4]
n_v = numdata[:,5]

plt.semilogy(z,n_e, label = 'ne')
plt.semilogy(z,n_p, label = 'np')
plt.semilogy(z,n_n, label = 'nn')
plt.semilogy(z,n_v, label = 'nv')
plt.legend(bbox_to_anchor=(0.8,0.95),loc=2, borderaxespad=0)
plt.show()
'''

# plotting avg cross sections files when we found that some of the cross sections are slightly off
'''
xdata1 = loadtxt('crossM3a0.8t20r18.5')
xdata2 = loadtxt('simul_xsctns.dat')

#old cross sections
z1 = xdata1[:,0]    
ve1 = xdata1[:,1]
vp1 = xdata1[:,2]
vn1 = xdata1[:,3]
vp_I1 = xdata1[:,4]
vv1 = xdata1[:,5]

#new cross sections
z2 = xdata2[:,0]
ve2 = xdata2[:,1]
vp2 = xdata2[:,2]
vn2 = xdata2[:,3]
vp_I2 = xdata2[:,4]
vv2 = xdata2[:,5]

#comparisons
plt.plot(z1,ve1,label = 've1')
plt.plot(z2,ve2,label = 've2')
plt.legend(bbox_to_anchor=(0.8,0.95),loc=2, borderaxespad=0)
plt.show()

plt.plot(z1,vp1,label = 'vp1')
plt.plot(z2,vp2,label = 'vp2')
plt.legend(bbox_to_anchor=(0.8,0.95),loc=2, borderaxespad=0)
plt.show()

plt.plot(z1,vn1,label = 'vn1')
plt.plot(z2,vn2,label = 'vn2')
plt.legend(bbox_to_anchor=(0.8,0.95),loc=2, borderaxespad=0)
plt.show()

plt.plot(z1,vp_I1,label = 'vp_I1')
plt.plot(z2,vp_I2,label = 'vp_I2')
plt.legend(bbox_to_anchor=(0.8,0.95),loc=2, borderaxespad=0)
plt.show()

plt.plot(z1,vv1,label = 'vv1')
plt.plot(z2,vv2,label = 'vv2')
plt.legend(bbox_to_anchor=(0.8,0.95),loc=2, borderaxespad=0)
plt.show()
'''
####################################
# Finding opacity
# We have to integrate 1/mfp going from far awway to close for each r value

#data = loadtxt("simul_all2.dat") 			#datafile with all data including mfp values
#data2 = open('opac_all.dat','w')			#datafile to write opac values to
#v_surface = open('v_surface.dat','w')			#where we'll put the neutrino surface values
#vA_surface = open('vA_surface.dat','w')

'''
r = data[:,0]
z = data[:,1]                                  		#
T = data[:,2]						#	
dens = data[:,3]					#
Ye = data[:,4]                                        	# picks out each column
mfp = data[:,5]                                        	#
mfpA = data[:,6]                                      	#

datasize = size(z) - 1
'''

#This function calculates the opacity for both e- neutrinos and antineutrinos and writes the data to a file
def opac(datasize):
	opac = 0
	opacA = 0                     
	for i in range(datasize):
		j = datasize - i			#To iterate through each r value going from highest to lowest
		print (i)				#even though data goes loest to highest
		if i != 0:
			if z[j+1] == 0:			#reset if we start a new radius value
				opac = 0
				opacA = 0
				data2.write("\n")
			dz = 0.5			#The height change between data points
			oldopac = opac			
			oldopacA = opacA
			opac += 1./mfp[j]*dz
			opacA += 1./mfpA[j]*dz
			if opac >= 2./3:		#If the opacity is greater than 2/3 and the previous value wasnt, it is the neutrino surface
				if oldopac < 2./3:							#so it gets written to a file
					v_surface.write(str(r[j])+' '+str(z[j])+'\n')
					print('\t', r[j],z[j])	
			if opacA >= 2./3:								#same thing but for antineutrinos
				if oldopacA < 2./3:
					vA_surface.write(str(r[j])+' '+str(z[j])+'\n')
					print('\t', r[j],z[j])
			data2.write(str(r[j])+' '+str(z[j])+' '+str(opac)+' '+str(opacA)+'\n')		#writes it to another file where rest of opacty data is 
		else:
			dz = 0.5									
			opac += 1./mfp[j]*dz
			opacA += 1./mfpA[j]*dz								
			data2.write(str(r[j])+' '+str(z[j])+' '+str(opac)+' '+str(opacA)+'\n')		# writes the rest of opacity data to a file
	return 0


#######################################
#important radii calculations

#various datafiles to write to
#outfile = open("Risco.dat",'w')

def Rs(M):
	Rs = 2 * M 						# schwarzchild radius in km
	return Rs

def Rph(a, M):
	Rph = 2 * M * (1 + cos(2/3. * acos(-a/M)))		# photon orbit
	return Rph

def Rmb(a, M):
	Rmb = 2 * M - a + 2 * sqrt(M) * sqrt(M - a)	 	# marginally bound circular orbit radius
	return Rmb

def Rms(a, M):                                            
	Z1 = 1+(1 - a**2/M**2)**(1/3.)*((1 + a/M)**(1/3.)+(1 - a/M)**(1/3.))
	Z2 = (3 * a**2/M**2 + Z1**2)**(1/2.)
	Rms = M*(3 + Z2 - ((3 - Z1)*(3 + Z1 + 2*Z2))**(1/2.))	# marginally stable circular orbit/ same as R_ISCO
	return Rms


#Prints all the radii for the simulation
#print ("Rs:", Rs(Mbh),"\nRph:",Rph(abh,Mbh),"\nRmb:",Rmb(abh, Mbh),"\nRms:",Rms(abh, Mbh), '\n\n')

#Various tests with the R_ISCO function

#plots how R_ISCO changes with changing spin parameter
'''
xs = []
ys = []
xs2 = []
ys2 = []
for i in range(0,1000):
	i = (10 - 2.39) * i/1000. + 2.40
	xs.append(i)
	ys.append(Rms(abh,i * Ms))
	i = ((i - 2.40)/(10 - 2.39)) 
	xs2.append(i)
	ys2.append(Rms(i*Mbh,Mbh))


plt.plot(xs,ys, label = 'a = 0.8')
plt.xlabel("Mass [Solar Masses]")
#plt.plot(xs2,ys2, label = 'a = 0.8')
#plt.xlabel("Spin [unitless]")
plt.ylabel("Marginally stable orbit [km]")
plt.grid(True)
plt.legend(bbox_to_anchor=(0.05,0.95),loc=2, borderaxespad=0)
plt.show()
'''

#plots how R_ISCO changes with changing mass value
'''
Mxs = []
Mys1 = []
Mys2 = []
Mys3 = []
Mys4 = []
Mys5 = []
Mys6 = []
Mys = [Mys1,Mys2,Mys3,Mys4,Mys5,Mys6]
for i in range(1, 1000):
	i /= 100.
	Mxs.append(i)
	for n in range(0,6):
		Mys[n].append(Rms(n * 2 / 10. * Mbh, i * Ms).real) 

for  i in range(6):
      plt.plot(Mxs,Mys[i], label = 'a ='+str(i * 2 / 10.))
plt.xlabel("Mass [Solar Masses]")
plt.ylabel("Marginally stable orbit [km]")
plt.grid(True)
plt.legend(bbox_to_anchor=(0.75,0.40),loc=2, borderaxespad = 0)
plt.show()

for i in range(999):
	outfile3.write(str(Mxs[i])+' '+str(Mys1[i])+' '+str(Mys2[i])+' '+str(Mys3[i])+' '+str(Mys4[i])+' '+str(Mys5[i])+' '+str(Mys6[i])+'\n')
outfile3.close()
'''		












