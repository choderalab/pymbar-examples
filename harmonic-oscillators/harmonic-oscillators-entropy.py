#!/usr/bin/python

#=============================================================================================
# Test MBAR by performing statistical tests on a set of of 1D harmonic oscillators, for which
# the true free energy differences can be computed analytically.
#
# A number of replications of an experiment in which i.i.d. samples are drawn from a set of
# K harmonic oscillators are produced.  For each replicate, we estimate the dimensionless free
# energy differences and mean-square displacements (an observable), as well as their uncertainties.
#
# For a 1D harmonic oscillator, the potential is given by
#   V(x;K) = (K/2) * (x-x_0)**2
# where K denotes the spring constant.
#
# The equilibrium distribution is given analytically by
#   p(x;beta,K) = sqrt[(beta K) / (2 pi)] exp[-beta K (x-x_0)**2 / 2]
# The dimensionless free energy is therefore
#   f(beta,K) = - (1/2) * ln[ (2 pi) / (beta K) ]
#
#=============================================================================================

#=============================================================================================
# TODO
#=============================================================================================
# * Generate a plota after completion, similar to the plot from WHAM paper.
#=============================================================================================

#=============================================================================================
# VERSION CONTROL INFORMATION
#=============================================================================================
__version__ = "$Revision: 186 $ $Date: 2010-05-06 08:05:26 -0700 (Thu, 06 May 2010) $"
# $Date: 2010-05-06 08:05:26 -0700 (Thu, 06 May 2010) $
# $Revision: 186 $
# $LastChangedBy: mrshirts $
# $HeadURL: https://simtk.org/svn/pymbar/trunk/examples/harmonic-oscillators/harmonic-oscillators.py $
# $Id: harmonic-oscillators.py 186 2010-05-06 15:05:26Z mrshirts $

#=============================================================================================
# IMPORTS
#=============================================================================================
import pymbar
import numpy
import confidenceintervals
import matplotlib.pyplot as plt
import pdb
#=============================================================================================
# PARAMETERS
#=============================================================================================

#K_k = numpy.array([1, 2, 4, 8, 16, 32, 64]) # spring constants for each state
#K_k = numpy.array([1, 2, 4, 8, 16, 32, 64]) # spring constants for each state
#O_k = numpy.array([0, 1, 2, 3, 4, 5, 6]) # offsets for spring constants
#N_k = numpy.array([4, 400, 40, 40, 40, 40, 40]) # number of samples from each state (can be zero for some states)

K_k = 24*numpy.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]) # 
O_k = numpy.array([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8]) # offsets for spring constants
N_k = numpy.array([40, 40,  40,  40,  40,  40,  40]) # number of samples from each state (can be zero for some states)

beta = 1.0 # inverse temperature for all simulations
nreplicates = 100 # number of replicates of experiment for testing uncertainty estimate

observe = 'position' # the observable, one of 'mean square displacement','position', or 'potential energy'
#observe = 'position^2' # the observable, one of 'mean square displacement','position', or 'potential energy'
#observe = 'RMS displacement' # the observable, one of 'RMS displacement', 'position', or 'potential energy'
#observe = 'potential energy' # the observable, one of 'RMS displacement','position', or 'potential energy'

# Uncomment the following line to seed the random number generated to produce reproducible output.
numpy.random.seed(0)

#=============================================================================================
# MAIN
#=============================================================================================

# Determine number of simulations.
K = numpy.size(N_k)
if numpy.shape(K_k) != numpy.shape(N_k): raise "K_k and N_k must have same dimensions."

# Determine maximum number of samples to be drawn for any state.
N_max = numpy.max(N_k)

# Compute widths of sampled distributions.
# For a harmonic oscillator with spring constant K,
# x ~ Normal(x_0, sigma^2), where sigma = 1/sqrt(beta K)
sigma_k = (beta * K_k)**-0.5
print "Gaussian widths:"
print sigma_k

# Compute the absolute dimensionless free energies of each oscillator analytically.
# f = - ln(sqrt((2 pi)/(beta K)) )
print 'Computing dimensionless free energies analytically...'
f_k_analytical = - numpy.log(numpy.sqrt(2 * numpy.pi) * sigma_k )

f_k_analytical -= f_k_analytical[0]
u_k_analytical = 0.5 * numpy.zeros([K], numpy.float64)
s_k_analytical = - f_k_analytical + u_k_analytical
print "f_k = "
print f_k_analytical
print "u_k = "
print u_k_analytical
print "s_k = "
print s_k_analytical

# Compute true free energy differences.
Deltaf_ij_analytical = numpy.zeros([K,K], dtype = numpy.float64)
Deltau_ij_analytical = numpy.zeros([K,K], dtype = numpy.float64)
Deltas_ij_analytical = numpy.zeros([K,K], dtype = numpy.float64)
for i in range(0,K):
  for j in range(0,K):
    Deltaf_ij_analytical[i,j] = f_k_analytical[j] - f_k_analytical[i]
    Deltas_ij_analytical[i,j] = u_k_analytical[j] - u_k_analytical[i]
    Deltas_ij_analytical[i,j] = u_k_analytical[j] - s_k_analytical[i]    

# Compute ensemble averages analytically 
if observe == 'RMS displacement':
  A_k_analytical = sigma_k              # mean square displacement
elif observe == 'potential energy':
  A_k_analytical = 1/(2*beta)*numpy.ones([K],float)  # By eqipartition
elif observe == 'position': 
  A_k_analytical  = O_k                 # observable is the position
elif observe == 'position^2': 
  A_k_analytical  = (1+ beta*K_k*O_k**2)/(beta*K_k)                 # observable is the position
else:
  raise "Observable %s not known." % observe

# DEBUG info
print "This script will perform %d replicates of an experiment where samples are drawn from %d harmonic oscillators." % (nreplicates, K)
print "The harmonic oscillators have equilibrium positions"
print O_k
print "and spring constants"
print K_k
print "and the following number of samples will be drawn from each (can be zero if no samples drawn):"
print N_k
print ""

# Conduct a number of replicates of the same experiment
replicates_observable = [] # storage for one hash for each replicate
replicates_df = [] # storage for one hash for each replicate
replicates_fdf = [] # storage for one hash for final observable
replicates_bar = [] # storage for one hash of bar free energy
replicates_ds = [] # storage for one hash for each replicate
replicates_dsw = [] # storage for one hash for each replicate

for replicate_index in range(0,nreplicates):
  print "Performing replicate %d / %d" % (replicate_index+1, nreplicates)

  # Initialize a hash to store data for this replicate.
  replicate_observable = { }
  replicate_df = { }
  replicate_fdf = { }
  replicate_bar = { }
  replicate_ds = { }
  replicate_dsw = { }
  #=============================================================================================
  # Generate independent data samples from K one-dimensional harmonic oscillators centered at q = 0.
  #=============================================================================================
  
  print 'generating samples...'

  # Initialize storage for samples.
  x_kn = numpy.zeros([K,N_max], dtype = numpy.float64) # x_kn[k,n] is the coordinate x of independent snapshot n of simulation k
  U_kln = numpy.zeros([K,K,N_max], dtype = numpy.float64) # U_kmt(k,m,t) is the value of the energy for snapshot t from simulation k at potential m

  # Generate samples.
  for k in range(0,K):
    # generate N_k[k] independent samples with spring constant K_k[k]
    x_kn[k,0:N_k[k]] = numpy.random.normal(O_k[k], sigma_k[k], [N_k[k]])
    
    # compute potential energy of all samples in all potentials
    for l in range(0,K):      
      U_kln[k,l,:] = (K_k[l]/2.0) * (x_kn[k,:]-O_k[l])**2
 
  #=============================================================================================
  # Estimate free energies and expectations.
  #=============================================================================================

  # Estimate free energies from simulation using MBAR.
  print "Estimating relative free energies from simulation (this may take a while)..."

  # Compute reduced potential energies.
  u_kln = beta * U_kln

  # Initialize the MBAR class, determining the free energies.
  mbar = pymbar.MBAR(u_kln, N_k) 
  
  O = mbar.computeOverlap(output='matrix')
  # Get matrix of dimensionless free energy differences and uncertainty estimate.
  #print "Computing covariance matrix..."
  #(Deltaf_ij_estimated, dDeltaf_ij_estimated) = mbar.getFreeEnergyDifferences()

  # Compute decomposition into entropy and enthalpy.
  (Deltaf_ij_estimated, dDeltaf_ij_estimated, Deltau_ij_estimated, dDeltau_ij_estimated, Deltas_ij_estimated, dDeltas_ij_estimated) = mbar.computeEntropyAndEnthalpy()

  # Compute entropy from Wyczalkowki estimator.
  f = mbar.f_k

  u = numpy.zeros([K], numpy.float64)
  W_nk = mbar.getWeights()
  for i in range(K):
    u_kn = u_kln[:,i,:]
    u[i] = sum(W_nk[:,i] * u_kn[mbar.indices])

  mu_ijn = numpy.zeros([K,K,N_max], numpy.float64)
  index = 0
  for i in range(K):
    for j in range(K):
      mu_ijn[i,j,0:N_k[i]] = W_nk[index:(index+N_k[i]),j]
    index += N_k[i]

  
  s = numpy.zeros(K, float)
  junk = numpy.zeros(K, float)  
  maxits = 150
  for iteration in range(maxits):
    for i in range(K):      
      s[i] = - f[i] + u[i]

      junk[i] = 0.0
      for j in range(K):
        junk[i] += N_k[j] * (mu_ijn[j,i,0:N_k[j]] * u_kln[j,j,0:N_k[j]]).mean()
        junk[i] += - N_k[j] * (mu_ijn[j,i,0:N_k[j]]).mean() * (u_kln[j,j,0:N_k[j]]).mean()
        for k in range(K):
          junk[i] += N_k[j] * N_k[k] * (f[k] + s[k]) * (mu_ijn[j,i,0:N_k[j]] * mu_ijn[j,k,0:N_k[j]]).mean()
          junk[i] += - N_k[j] * N_k[k] * (mu_ijn[j,i,0:N_k[j]] * mu_ijn[j,k,0:N_k[j]] * u_kln[j,k,0:N_k[j]]).mean()

      s[i] += junk[i]
  print s[:] - s[0]
  print Deltas_ij_estimated[0,:]
  print junk

  Deltasw_ij_estimated = numpy.zeros([K,K], numpy.float64)
  dDeltasw_ij_estimated = numpy.zeros([K,K], numpy.float64)  
  for i in range(K):
    for j in range(K):
      Deltasw_ij_estimated[i,j] = s[j] - s[i]
  Deltasw_ij_error = Deltasw_ij_estimated - Deltas_ij_analytical

  Deltas_ij_error = Deltas_ij_estimated - Deltas_ij_analytical
      
  # Compute error from analytical free energy differences.
  Deltaf_ij_error = Deltaf_ij_estimated - Deltaf_ij_analytical

  # compute the energies and variance with BAR
  df_k = numpy.zeros(K,float)
  ddf_k = numpy.zeros(K,float)

  Bar_ij = numpy.zeros([K,K],float)
  dBar_ij = numpy.zeros([K,K],float)

  for k in range(0,K-1):
    w_F = u_kln[k,k+1,0:N_k[k]] - u_kln[k,k,0:N_k[k]] 
    w_R = u_kln[k+1,k,0:N_k[k+1]] - u_kln[k+1,k+1,0:N_k[k+1]] 

    (df_k[k],ddf_k[k]) = pymbar.computeBAR(w_F,w_R)

  for i in range(0,K):
    for j in range(i+1,K):
      Bar_ij[i,j] = numpy.sum(df_k[range(i,j)])
      dBar_ij[i,j] = numpy.sqrt(numpy.sum(ddf_k[range(i,j)]**2))
      Bar_ij[j,i] = -Bar_ij[i,j]
      dBar_ij[j,i] = dBar_ij[i,j]

  Bar_ij_error = Bar_ij - Deltaf_ij_analytical    

  # Estimate the expectation of the mean-squared displacement at each condition.
  if observe == 'RMS displacement':
    A_kn = numpy.zeros([K,K,N_max], dtype = numpy.float64);
    for k in range(0,K):
      for l in range(0,K):
        A_kn[k,l,0:N_k[k]] = (x_kn[k,0:N_k[k]] - O_k[l])**2 # observable is the squared displacement

  # observable is the potential energy
  elif observe == 'potential energy':
    A_kn = U_kln

  # observable for estimation is the position
  elif observe == 'position': 
    A_kn = numpy.zeros([K,N_max], dtype = numpy.float64);
    for k in range(0,K):
      A_kn[k,0:N_k[k]] = x_kn[k,0:N_k[k]]  

  elif observe == 'position^2': 
    A_kn = numpy.zeros([K,N_max], dtype = numpy.float64);
    for k in range(0,K):
      A_kn[k,0:N_k[k]] = x_kn[k,0:N_k[k]]**2   

  (A_k_estimated, dA_k_estimated) = mbar.computeExpectations(A_kn)

  # need to additionally transform to get the square root
  if observe == 'RMS displacement':
    A_k_estimated = numpy.sqrt(A_k_estimated)

    # Compute error from analytical observable estimate.
    dA_k_estimated = dA_k_estimated/(2*A_k_estimated)

  A_k_error = A_k_estimated - A_k_analytical

  #=============================================================================================
  # Store data for this replicate.
  #=============================================================================================  
  replicate_df['estimated'] = Deltaf_ij_estimated.copy()
  replicate_df['destimated'] = dDeltaf_ij_estimated.copy()
  replicate_df['error'] = Deltaf_ij_error.copy()
  replicates_df.append(replicate_df)

  replicate_bar['estimated'] = Bar_ij.copy()
  replicate_bar['destimated'] = dBar_ij.copy()
  replicate_bar['error'] = Bar_ij_error.copy()
  replicates_bar.append(replicate_bar)

  replicate_observable['estimated'] = A_k_estimated.copy()
  replicate_observable['destimated'] = dA_k_estimated.copy()
  replicate_observable['error'] = A_k_error.copy()  
  replicates_observable.append(replicate_observable)

  replicate_ds['estimated'] = Deltas_ij_estimated.copy()
  replicate_ds['destimated'] = dDeltas_ij_estimated.copy()
  replicate_ds['error'] = Deltas_ij_error.copy()
  replicates_ds.append(replicate_ds)

  replicate_dsw['estimated'] = Deltasw_ij_estimated.copy()
  replicate_dsw['destimated'] = dDeltasw_ij_estimated.copy()
  replicate_dsw['error'] = Deltasw_ij_error.copy()
  replicate_dsw['junk'] = junk.copy()
  replicates_dsw.append(replicate_dsw)

ds_error = numpy.zeros([nreplicates], numpy.float64)
ds_stderr = numpy.zeros([nreplicates], numpy.float64)
dsw_error = numpy.zeros([nreplicates], numpy.float64)
dsw_stderr = numpy.zeros([nreplicates], numpy.float64)
junk = numpy.zeros([nreplicates], numpy.float64)
for replicate_index in range(nreplicates):
  ds_error[replicate_index] = replicates_ds[replicate_index]['error'][0,K-1]
  ds_stderr[replicate_index] = replicates_ds[replicate_index]['destimated'][0,K-1]
  dsw_error[replicate_index] = replicates_dsw[replicate_index]['error'][0,K-1]
  dsw_stderr[replicate_index] = replicates_dsw[replicate_index]['destimated'][0,K-1]
  junk[replicate_index] = replicates_dsw[replicate_index]['junk'][K-1]

print "RMS error"
print "s = -f + u : %12.6f" % ds_error.std()
print "Eq. 33     : %12.6f" % dsw_error.std()
print "junk       : %12.6f" % junk.std()
print ""
print "bias"
print "s = -f + u : %12.6f +- %12.6f" % (ds_error.mean(), ds_error.std() / numpy.sqrt(nreplicates))
print "Eq. 33     : %12.6f +- %12.6f" % (dsw_error.mean(), dsw_error.std() / numpy.sqrt(nreplicates))
print "junk       : %12.6f +- %12.6f" % (junk.mean(), junk.std() / numpy.sqrt(nreplicates))
print ""
print "mean error estimates"
print "s = -f + u : %12.6f" % ds_stderr.mean()
print "Eq. 33     : %12.6f" % dsw_stderr.mean()

# compute the probability distribution of all states
(alpha_ds,Pobs_ds,Plow_ds,Phigh_ds,dPobs_ds,Pnorm_ds) = confidenceintervals.generateConfidenceIntervals(replicates_ds,K);

override = {
  'family'              : 'sans-serif',
  'verticalalignment'   : 'bottom',
  'horizontalalignment' : 'center',
  'weight'              : 'bold',
  'size'                : 30
}

formatstrings = ['b-','g-','c-','y-', 'r-','m-']

plt.figure(1);
plt.axis([0.0, 4.0, 0.0, 1.0])
plt.plot(alpha_ds,Pnorm_ds,'k-',label="Normal")
plt.title('Cumulative Probabilty vs. Normal Distribution',size=24)
plt.xlabel('Standard Deviations',size = 18)
plt.ylabel('Cumulative Probability',size= 18)
plt.legend(loc=4)
plt.savefig('alpha_ds_harm.pdf')

plt.figure(2);
plt.axis([0.0, 4.0, 0.0, 1.0])
plt.plot(alpha_ds,Pnorm_ds,'r-',label = 'Normal Distribution')
plt.plot(alpha_ds,Plow_ds,'k_')
plt.plot(alpha_ds,Pobs_ds,'k-',label = 'Observed')
plt.plot(alpha_ds,Phigh_ds,'k_')
plt.title('Cumulative Probabilty vs. Normal Distribution',size=24)
plt.xlabel('Standard Deviations',size = 18)
plt.ylabel('Cumulative Probability',size= 18)
plt.legend(loc=4)
plt.savefig('alpha_ds_harm2.pdf')

