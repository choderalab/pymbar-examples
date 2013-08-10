#!/usr/bin/python

# Example illustrating the use of MBAR for computing the hydration free energy of OPLS 3-methylindole
# in TIP3P water through alchemical free energy simulations.

#===================================================================================================
# IMPORTS
#===================================================================================================

import numpy
import pymbar # multistate Bennett acceptance ratio estimator.
import timeseries # for timeseries analysis 
import commands
import re
import sys
import pdb

#===================================================================================================
# CONSTANTS
#===================================================================================================

convert_atmnm3_to_kJmol = 1.01325e5*(1e-09)**3 * 6.02214 * (1e23) / 1000 # Convert pV from atmospheres*nm^3 into kJ/mol
kB = 1.381*6.02214/1000.0  # Boltzmann's constant (kJ/mol/K)
relative_tolerance = 1e-10
verbose = True

#===================================================================================================
# PARAMETERS
#===================================================================================================

if (len(sys.argv) > 1):
   # directory in which datafiles are stored   
   datafile_directory = sys.argv[1]
else:
   print "assuming 38 step solvation"
   datafile_directory = 'data/3-methylindole-38steps/' 
datafile_prefix  = 'dhdl' # prefixes for datafile sets 

temperature = 298.0 # temperature (K)
pressure = 1.0 # pressure (atm)

nequil = 100; # number of samples assumed to be equilibration, and thus omitted.
#datafile_directory = 'data/3-methylindole-38steps/' # directory in which datafiles are stored
#datafile_directory = 'data/3-methylindole-11steps/' # directory in which datafiles are stored
datafile_directory = '../../../large-datasets/trp38/' # directory in which datafiles are stored
#datafile_directory = '../../../large-datasets/trp50ns/' # directory in which datafiles are stored
#datafile_prefix  = 'dhdl' # prefixes for datafile sets 

#===================================================================================================
# HELPER FUNCTIONS
#===================================================================================================
def sortbynum(item):
   vals = item.split('.')
   for v in reversed(vals):
      if v.isdigit():
         return int(v)
   print "Error: No digits found in filename, can't sort ", item


#===================================================================================================
# MAIN
#===================================================================================================

# Compute inverse temperature and thermal energy.
kT = (kB*temperature) # thermal energy, in units of kJ/mol
beta = 1./kT # inverse temperature, in units of 1/(kJ/mol)

#===================================================================================================
# Read all snapshot data
#===================================================================================================

# Generate a list of all datafiles with this prefix.
# NOTE: This uses the file ordering provided by 'ls' to determine which order the states are in.
filenames = commands.getoutput('ls %(datafile_directory)s/%(datafile_prefix)s.*' % vars()).split()

# Determine number of alchemical intermediates.
K = len(filenames)
# Determine length of files.
filenames = commands.getoutput('ls %(datafile_directory)s/%(datafile_prefix)s*' % vars()).split()   
   
nsnapshots = numpy.zeros(K, int) # nsnapshots[k] is the number of snapshots for state k; 
# no larger than number of lines starting out.

# sort the files numerically.
filenames.sort(key=sortbynum)

for k in range(K):
   # Temporarily read the file into memory.
   infile = open(filenames[k], 'r')
   lines = infile.readlines()
   infile.close()
   
   # Determine maxnumber of snapshots from quickly parsing file and ignoring header lines.
   nsnapshots[k] = 0
   for line in lines:
      if ((line[0] == '#') or (line[0] == '@')):
         continue
      nsnapshots[k] += 1   

maxn = max(nsnapshots)   

# Load all of the data
pe = numpy.zeros(K,numpy.float64);
u_klt = numpy.zeros([K,K,maxn], numpy.float64) # u_klt[k,m,t] is the reduced potential energy of snapshot t of state k evaluated at state m
for k in range(K):
   # File to be read
   filename = filenames[k]   
   
      # Read contents of file into memory.
   print "Reading %s..." % filename
   infile = open(filename, 'r')
   lines = infile.readlines()
   infile.close()

   t = 0
   n_components = 0
   n_states = 0
   bPV = False

   # Parse the file.

   for line in lines:

      # split into elements
      elements = line.split()

      # This section automatically parses the header to count the number
      # of dhdl components, and the number of states at which energies
      # are calculated, and should be modified for different file input formats.
      #                                                                          
      if ((line[0] == '#') or (line[0] == '@')):
         if (line[0] == '@'):
            # it's an xvg legend entry -- load in the information
            if (line[2] == 's') and (line[3] != 'u'):  
               # it's a xvg entry, and a lambda component, so note it down
               if (re.search('-lambda',line)):     
                  #it's a listing of the lambdas
                  n_components +=1
                  lv = numpy.zeros([K,n_components],float)
               elif (re.search("\\\\xD\\\\f\\{\\}H \\\\x\\l\\\\f\\{\\}",line)): 
               lambdaone = re.sub('[(),"]', '', elements[6])
	       lambdatwo = re.sub('[(),"]', '', elements[7])
               lambdas = [lambdaone, lambdatwo]
                  for i in range(n_components):
                     lv[n_states,i] = float(lambdas[i])
                  n_states+=1;   
               elif (re.search("pv",line)):     
                   bPV = True;   
      else:                           
         if ((t==0) and (k==0)):     # we don't know the number of components until here.
            dhdlt = numpy.zeros([K,n_components,maxn],float) 
         time = float(elements.pop(0))
            
         #
         #  If we print the energy in the dhdl file; if not, delete this line.
         #

         energy = float(elements.pop(0))

         # 
         # In this section, store the derivative with respect to lambda
         # 
      
         for nl in range(n_components):
            dhdlt[k,nl,t] = float(elements.pop(0))
         # now record the potential energy differences.   
         for l in range(K):
            pe[l] = float(elements.pop(0)) + energy
                  
         # pressure-volume contribution
         if (bPV):
            pv = float(elements.pop(0))
         else:
            pv = 0

         # compute and store reduced potential energy at each state
         for l in range(K):
            u_klt[k,l,t] = beta * (pe[l] + pv)
               
         t += 1   
   nsnapshots[k] = t

#===================================================================================================
# Preliminaries: Subsample data to obtain uncorrelated samples
#===================================================================================================   
Nequilibrated = max(nsnapshots) - nequil
u_kln = numpy.zeros([K,K,Nequilibrated], numpy.float64) # u_kln[k,m,n] is the reduced potential energy of uncorrelated sample index n from state k evaluated at state m
dhdl = numpy.zeros([K,n_components,Nequilibrated], float) #dhdl is value for dhdl for each component in the file at each time.
N_k = numpy.zeros(K, int) # N_k[k] is the number of uncorrelated samples from state k
g = numpy.zeros(K,float) # autocorrelation times for the data
for k in range(K):
   # Determine indices of uncorrelated samples from potential autocorrelation analysis at state k.
   # alternatively, could use the energy differences -- here, we will use total dhdl
   dhdl_sum = numpy.sum(dhdlt[k,:,:],axis=0)
   g[k] = timeseries.statisticalInefficiency(dhdl_sum[nequil:nsnapshots[k]],suppress_warning=True)
   indices = numpy.array(timeseries.subsampleCorrelatedData(dhdl_sum[nequil:nsnapshots[k]],g=g[k])) # indices of uncorrelated samples
   N = len(indices) # number of uncorrelated samples
   N_k[k] = N      
   indices += nequil
   for n in range(n_components):
      dhdl[k,n,0:N] = dhdlt[k,n,indices]
   for l in range(K):         
      u_kln[k,l,0:N] = u_klt[k,l,indices]
print "Correlation times:"
print g
print ""
print "number of uncorrelated samples:"
print N_k
print ""

#===================================================================================================
# Estimate free energy difference with MBAR -- all states at once
#===================================================================================================   

# Initialize MBAR (computing free energy estimates, which may take a while)
print "Solving for MBAR estimate..."
MBAR = pymbar.MBAR(u_kln, N_k, verbose = verbose, method = 'adaptive', relative_tolerance = relative_tolerance, initialize = 'BAR') # use faster Newton-Raphson solver

# Estimate unitless entropy and enthalpy contributions to free energy differences.
print "Estimating entropic and enthalpic contributions..."
(Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij) = MBAR.computeEntropyAndEnthalpy()

# estimate enthalpy by endpoint method: <U>= \sum_{i=1}^N U / N
u_i = numpy.average(u_kln[0,0,0:N_k[0]])
du_i = numpy.sqrt(numpy.var(u_kln[0,0,0:N_k[0]])/(N_k[0]-1))

u_f = numpy.average(u_kln[K-1,K-1,0:N_k[K-1]])
du_f = numpy.sqrt(numpy.var(u_kln[K-1,K-1,0:N_k[K-1]])/(N_k[K-1]-1))

dH = u_f-u_i
ddH = numpy.sqrt(du_i**2 + du_f**2)
print "Enthalpic contribution (dH) by endpoint method"
#print "%6.3f +/- %6.3f kT" % (dH,ddH)
print " dH  : %6.3f +/- %6.3f kJ*mol^-1" % (dH*kT,ddH*kT)

print "Entropic contribution (T*dS) and entropy (dS) by endpoint method"
# estimate entropy by endpoint method: TdS = <U>-<G>
TdS = dH - Delta_f_ij[0,K-1]
dTdS = numpy.sqrt(dDelta_f_ij[0,K-1]**2+ddH**2)
#print "%6.3f +/- %6.3f kT" % (TdS,dTdS)
print "T*dS : %6.3f +/- %6.3f J*mol^-1*K^-1" % (kB*temperature*TdS,kB*temperature*dTdS)
print " dS  : %6.3f +/- %6.3f J*mol^-1*K^-1" % (1000*kB*TdS,1000*kB*dTdS)

# Display free energy difference and entropic/enthalpic contributions in units of kT.
# Note that all internal calculations are in dimesionless units, and must be converted back to physical units.
print ""
print "Free energy difference (dG) between states %d and %d" % (0, K-1)
# in kT print "%6.3f +/- %6.3f kT" % (Delta_f_ij[0,K-1], dDelta_f_ij[0,K-1])
print "%6.3f +/- %6.3f kJ*mol^-1" % (Delta_f_ij[0,K-1] * kT, dDelta_f_ij[0,K-1] * kT)
print ""

print "Enthalpic contribution (dH) by MBAR"
#print "%6.3f +/- %6.3f kT" % (Delta_u_ij[0,K-1], dDelta_u_ij[0,K-1])
print " dH  : %6.3f +/- %6.3f kJ*mol^-1" % (Delta_u_ij[0,K-1] * kT, dDelta_u_ij[0,K-1] * kT)

print "Entropic contribution (T*dS) and entropy (dS) by MBAR"
#print "%6.3f +/- %6.3f kT" % (Delta_s_ij[0,K-1], dDelta_s_ij[0,K-1])
print "T*dS : %6.3f +/- %6.3f kJ*mol^-1" % (Delta_s_ij[0,K-1] * kB*temperature, dDelta_s_ij[0,K-1] * kB*temperature)
print " dS  :  %6.3f +/- %6.3f J*mol^-1*K^-1" % (Delta_s_ij[0,K-1] * 1000*kB, dDelta_s_ij[0,K-1] * 1000*kB)
print ""

print "Note that the enthalpic and entropic contributions have correlated errors, since dG = dH - T*dS."

