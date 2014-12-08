# Example illustrating the application of MBAR to compute a 1D PMF from an umbrella sampling simulation.
#
# The data represents an umbrella sampling simulation for the chi torsion of a valine sidechain in lysozyme L99A with benzene bound in the cavity.
# 
# REFERENCE
# 
# D. L. Mobley, A. P. Graves, J. D. Chodera, A. C. McReynolds, B. K. Shoichet and K. A. Dill, "Predicting absolute ligand binding free energies to a simple model site," Journal of Molecular Biology 371(4):1118-1134 (2007).
# http://dx.doi.org/10.1016/j.jmb.2007.06.002

import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
# Constants.
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K

temperature = 300 # assume a single temperature -- can be overridden with data from center.dat 
# Parameters
K = 26 # number of umbrellas
N_max = 501 # maximum number of snapshots/simulation
T_k = numpy.ones(K,float)*temperature # inital temperatures are all equal 
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
chi_min = -180.0 # min for PMF
chi_max = +180.0 # max for PMF
nbins = 36 # number of bins for 1D PMF

# Allocate storage for simulation data
N_k = numpy.zeros([K], numpy.int32) # N_k[k] is the number of snapshots from umbrella simulation k
K_k = numpy.zeros([K], numpy.float64) # K_k[k] is the spring constant (in kJ/mol/deg**2) for umbrella simulation k
chi0_k = numpy.zeros([K], numpy.float64) # chi0_k[k] is the spring center location (in deg) for umbrella simulation k
chi_kn = numpy.zeros([K,N_max], numpy.float64) # chi_kn[k,n] is the torsion angle (in deg) for snapshot n from umbrella simulation k
u_kn = numpy.zeros([K,N_max], numpy.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
g_k = numpy.zeros([K],numpy.float32);

# Read in umbrella spring constants and centers.
infile = open('data/centers.dat', 'r')
lines = infile.readlines()
infile.close()
for k in range(K):
    # Parse line k.
    line = lines[k]
    tokens = line.split()
    chi0_k[k] = float(tokens[0]) # spring center locatiomn (in deg)
    K_k[k] = float(tokens[1]) * (numpy.pi/180)**2 # spring constant (read in kJ/mol/rad**2, converted to kJ/mol/deg**2)    
    if len(tokens) > 2:
        T_k[k] = float(tokens[2])  # temperature the kth simulation was run at.

beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
DifferentTemperatures = True
if (min(T_k) == max(T_k)):
    DifferentTemperatures = False            # if all the temperatures are the same, then we don't have to read in energies.
# Read the simulation data
for k in range(K):
    # Read torsion angle data.
    filename = 'data/prod%d_dihed.xvg' % k
    print "Reading %s..." % filename
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    # Parse data.
    n = 0
    for line in lines:
        if line[0] != '#' and line[0] != '@':
            tokens = line.split()
            chi = float(tokens[1]) # torsion angle
            # wrap chi_kn to be within [-180,+180)
            while(chi < -180.0):
                chi += 360.0
            while(chi >= +180.0):
                chi -= 360.0
            chi_kn[k,n] = chi
            
            n += 1
    N_k[k] = n

    if (DifferentTemperatures):  # if different temperatures are specified the metadata file, 
                                 # then we need the energies to compute the PMF
        # Read energies
        filename = 'data/prod%d_energies.xvg' % k
        print "Reading %s..." % filename
        infile = open(filename, 'r')
        lines = infile.readlines()
        infile.close()
        # Parse data.
        n = 0
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()            
                u_kn[k,n] = beta_k[k] * (float(tokens[2]) - float(tokens[1])) # reduced potential energy without umbrella restraint
                n += 1

    # Compute correlation times for potential energy and chi
    # timeseries.  If the temperatures differ, use energies to determine samples; otherwise, use the cosine of chi
            
    if (DifferentTemperatures):        
        g_k[k] = timeseries.statisticalInefficiency(u_kn[k,:], u_kn[k,0:N_k[k]])
        print "Correlation time for set %5d is %10.3f" % (k,g_k[k])
        indices = timeseries.subsampleCorrelatedData(u_kn[k,0:N_k[k]])
    else:
        chi_radians = chi_kn[k,0:N_k[k]]/(180.0/numpy.pi)
        g_cos = timeseries.statisticalInefficiency(numpy.cos(chi_radians))
        g_sin = timeseries.statisticalInefficiency(numpy.sin(chi_radians))
        print "g_cos = %.1f | g_sin = %.1f" % (g_cos, g_sin)
        g_k[k] = max(g_cos, g_sin)
        print "Correlation time for set %5d is %10.3f" % (k,g_k[k])
        indices = timeseries.subsampleCorrelatedData(chi_radians, g=g_k[k]) 
    # Subsample data.
    N_k[k] = len(indices)
    u_kn[k,0:N_k[k]] = u_kn[k,indices]
    chi_kn[k,0:N_k[k]] = chi_kn[k,indices]

N_max = numpy.max(N_k) # shorten the array size
u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l

# Set zero of u_kn -- this is arbitrary.
u_kn -= u_kn.min()

# Construct torsion bins
print "Binning data..."
delta = (chi_max - chi_min) / float(nbins)
# compute bin centers
bin_center_i = numpy.zeros([nbins], numpy.float64)
for i in range(nbins):
    bin_center_i[i] = chi_min + delta/2 + delta * i
# Bin data
bin_kn = numpy.zeros([K,N_max], numpy.int32)
for k in range(K):
    for n in range(N_k[k]):
        # Compute bin assignment.
        bin_kn[k,n] = int((chi_kn[k,n] - chi_min) / delta)

# Evaluate reduced energies in all umbrellas
print "Evaluating reduced potential energies..."
for k in range(K):
    for n in range(N_k[k]):
        # Compute minimum-image torsion deviation from umbrella center l
        dchi = chi_kn[k,n] - chi0_k
        for l in range(K):
            if (abs(dchi[l]) > 180.0):
                dchi[l] = 360.0 - abs(dchi[l])

        # Compute energy of snapshot n from simulation k in umbrella potential l
        u_kln[k,:,n] = u_kn[k,n] + beta_k[k] * (K_k/2.0) * dchi**2

# Initialize MBAR.
print "Running MBAR..."
mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive')

# Compute PMF in unbiased potential (in units of kT).
(f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins)

# Write out PMF
print "PMF (in units of kT)"
print "%8s %8s %8s" % ('bin', 'f', 'df')
for i in range(nbins):
    print "%8.1f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i])

################
# Now compute PMF assuming a cubit spline
import pdb
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import minimize
verbose = True

method = 'spline'
nspline = 10
#method = 'periodic' - not really working, try spline
#nperiod = 3

chi_n = chi_kn[mbar.indices]
# compute KL divergence to the empirical distribution for the trial distribution F

# define the bias functions
def fbias(k,x):
    dchi = x - chi0_k[k]
        # vectorize the conditional
    i = numpy.fabs(dchi) > 180.0
    dchi = i*(360.0 - numpy.fabs(dchi)) + (1-i)*dchi
    return beta_k[k] * (K_k[k]/2.0) * dchi**2

# define functions that can change each iteration

if method == 'spline':
    xstart = numpy.linspace(chi_min,chi_max,nspline)
    def trialf(t):
        f = interp1d(xstart, t, kind='cubic')
        return f
    tstart = 0*xstart
    #for i in range(nspline):
    #    tstart[i] = f_i[numpy.argmin(numpy.abs(bin_center_i-xstart[i]))]  # start with nearest PMF value

if method == 'periodic':
    # vary the magnitude, phase, and period: not clear this is really working at all
    def trialf(t):
        def interperiod(x):
            y = numpy.zeros(numpy.size(x))
            for i in range(nperiod):
                t[i+2*nperiod] = t[i+2*nperiod]%(360.0) # recenter the offsets, all in range
                y += t[i]*numpy.cos(t[i+nperiod]*x+t[i+2*nperiod])
            return y
        return interperiod

    tstart = numpy.zeros(3*nperiod)
    # initial values of ampliudes, period, and phase
    tstart[0:nperiod] = 0
    d = (chi_max - chi_min)/(nperiod+1)
    tstart[nperiod:2*nperiod] = (2*numpy.pi/(chi_max-chi_min))
    tstart[2*nperiod:3*nperiod] = numpy.linspace(chi_min+d/2,chi_max-d/2,nperiod)

# = - <ln P(t,x)> = - \sum_{x_n} w_i(x_n) ln (P(t,x_n)/\int P(t,x_n)) dx
#                  = - \sum_{x_n} w_i(x_n) ln P(t,x_n) + ln \int P(t,x_n) dx
#                          set P(t,x) = exp(-F(t,x))
#                  = [\sum_{x_n} w_i(x_n) F(t,x_n)] + ln \int exp(-F(t,x_n) dx
#
# gradient of the parameters could be difficult in general
#                  = d/dt [\sum_{x_n} w_i(x_n) F(t,x_n)] + ln \int exp(-F(t,x_n) dx
#                  = \sum_{x_n} w_i(x_n) dF(t,x_n)/dt + d/dt[\int exp(-F(t,x_n))dx ]/\int exp(-F(t,x_n) dx
#                  = \sum_{x_n} w_i(x_n) dF(t,x_n)/dt - \int dF/dt exp(-F(t,x_n)) dx ]/\int exp(-F(t,x_n) dx
#                    so we just need to calculate dF/dt.  however, for functional forms like splines, this
#                    can be rather difficult, since changes in the input parameters (value at boundaries) change
#                    throughout the function.
#

def sumkldiverge(t,ft,x_n,K,w_kn,fbias,xrange):

    # we are interested in finding the potential that minimizes the
    # sum of the KL divergence over all of the distributions

    # define the function f, the sum of the KL divergence over all of
    # the samples, based on the current parameters t

    t -= t[0] # set a reference state, may make the minization faster by removing degenerate solutions
    #t[len(t)-1] = t[0]
    feval = ft(t)
    fx = feval(x_n)  # only need to evaluate this over all points outside
    kl = 0 
    # figure out the bias 
    for k in range(K):
        # what is the biasing function for this state
        bias = lambda x: fbias(k,x)
        # define the exponential of f based on the current parameters t.
        expf = lambda x: numpy.exp(-feval(x)-bias(x))
        pE = numpy.dot(w_kn[:,k],fx+bias(x_n))
        pF = numpy.log(quad(expf,xrange[0],xrange[1])[0])  #0 is the value of quad
        kl += (pE + pF)
    if verbose:    
        print kl,t    
    return kl

w_kn = numpy.exp(mbar.Log_W_nk) # normalized weights

# inputs to kldivergence in minimize are:
# the function that we are computing the kldivergence of
# the x values we have samples at
# the number of umbrellas
# the weights at the samples
# the umbrella restraints strengths
# the umbrella restraints centers
# the domain of the function

# set minimizer options to display. Apprently does not exist for
# BFGS. Probably don't need to set eps.
options = {'disp':True, 'eps':10**(-3)}
result = minimize(sumkldiverge,tstart,args=(trialf,chi_n,K,w_kn,fbias,[chi_min,chi_max]),options=options)
pmf_final = trialf(result.x)
nplot = 1000
import matplotlib.pyplot as plt
x = numpy.linspace(chi_min,chi_max,nplot)
plt.plot(bin_center_i,f_i,'rx')
yout = pmf_final(x)
ymin = numpy.min(yout)
yout -= ymin
# Write out PMF
print "PMF (in units of kT)"
print "%8s %8s %8s" % ('bin', 'f', 'df')
for i in range(nbins):
    print "%8.1f %8.3f" % (bin_center_i[i], pmf_final(bin_center_i[i])-ymin)
plt.plot(x,yout,'k-')
plt.xlim([chi_min,chi_max])
plt.show()

