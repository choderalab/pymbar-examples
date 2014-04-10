#!/bin/env python

# Originally written by Michael Shirts as:
# Example illustrating the use of MBAR for computing the hydration free energy of OPLS 3-methylindole
# in TIP3P water through alchemical free energy simulations.

# Adapted by P. Klimovich and D. Mobley, March 2011, to be slightly more general.
# Additionally adapted by Michael Shirts and P. Klimovich, May 2013

#===================================================================================================
# IMPORTS
#===================================================================================================

import pickle
import numpy
from collections import Counter # for counting elements in an array
from commands import getoutput # a note for the future: in Python 3.X, 'getoutput' is moved to the 'subprocess' module
import pymbar # multistate Bennett acceptance ratio estimator
import pymbar
from pymbar import timeseries # for timeseries analysis
from pymbar import bar,exponential_averaging
from pymbar.bar import BAR
from pymbar.exponential_averaging import EXP,EXPGauss
from optparse import OptionParser # for parsing command-line options
import os # operating system dependent modules of Python
import time as ttt_time # for timing
import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib.font_manager import FontProperties as FP
import pdb

#===================================================================================================
# CONSTANTS AND PARAMETERS
#===================================================================================================

kB = 1.381*6.02214/1000.0 # Boltzmann's constant (kJ/mol/K).
relative_tolerance = 1e-10 # Convergence criterion of the energy estimates for BAR and MBAR.
initialize_with = 'BAR' # The initial MBAR free energy guess is set to BAR.
methods = ['TI','TI-CUBIC','DEXP','IEXP','BAR','MBAR'] # Free energy estimation methods.
#methods = ['IEXP','DEXP'] # Free energy estimation methods.
#methods = ['TI','TI-CUBIC','DEXP','IEXP','GDEL','GINS','BAR','UBAR','RBAR','MBAR'] # All supported free energy estimation methods.
datafile_suffix = '.xvg' # Suffix of the energy file(s).

#===================================================================================================
# INPUT OPTIONS
#===================================================================================================

parser = OptionParser()
parser.add_option('-d', '--dir', dest = 'dir', help = 'Directory in which data files are stored. Default: Current dir.', default = '.')
parser.add_option('-f', '--forwrev', dest = 'bForwrev', help = 'Plotting the free energy change as a function of time in both directions. Default: False.', default = False, action = 'store_true')
parser.add_option('-g', '--breakdown', dest = 'breakdown', help = 'Plotting the free energy differences evaluated for each pair of adjacent states for all methods. Default: True.', default = True, action = 'store_true')
parser.add_option('-p', '--prefix', dest = 'prefix', help = 'Prefix for datafile sets, i.e. prod (default).', default = 'prod')
parser.add_option('-s', '--skiptime', dest = 'equiltime', help = 'Discard data prior to this specified time as \'equilibration\' data. Units picoseconds. Default: 100 ps.', default = 100.0)
parser.add_option('-t', '--temperature', dest = 'temperature', help = "Temperature in K. Default: 298 K.", default = 298, type=float)
parser.add_option('-u', '--units', dest = 'units', help = 'Units to report energies: \'kJ\', \'kcal\', and \'kBT\'. Default: \'kJ\'', default = 'kJ')
parser.add_option('-v', '--verbose', dest = 'verbose', help = 'Verbose option for BAR and MBAR. Default: False.', default = False, action = 'store_true')
parser.add_option('-x', '--ignoreWL', dest = 'bignoreWL', help = 'Do not check whether the WL weights are equilibrated. No log file needed as an accompanying input.', default = False, action = 'store_true')

(options, args) = parser.parse_args()
datafile_directory = options.dir
datafile_prefix = options.prefix
equiltime = float(options.equiltime)
units = options.units
verbose = options.verbose

#===================================================================================================
# MAIN
#===================================================================================================

stime=ttt_time.time()
print "Started on %s" % ttt_time.asctime()

# Define units the energies will be reported in.
beta = 1./(kB*options.temperature)
if units == 'kJ':
   beta_report = beta
   units = '(kJ/mol)'
elif units == 'kcal':
   beta_report = 4.184*beta
   units = '(kcal/mol)'
elif units == 'kBT':
   beta_report = 1
   units = '(k_BT)'
else:
   parser.error('\nI don\'t understand the unit type \'%s\': the only options \'kJ\', \'kcal\', and \'kBT\'' % units)

#===================================================================================================
# Preliminaries I: List the .xvg files, count them up, sort them; read in the @-header information.
#===================================================================================================

# Generate a list of all datafiles, sort them numerically if needed, and count them up.
filenames = getoutput('ls %(datafile_directory)s/%(datafile_prefix)s*%(datafile_suffix)s' % vars()).split()
if filenames[0] == 'ls:':
   parser.error("\nNo files found within directory '%s' with prefix '%s': check your inputs." % (datafile_directory, datafile_prefix))
if len(filenames)>1:
   # determine which values are digits:
   isdigits = []
   for f in filenames:
      vals = f.split('.')
      for i,v in enumerate(vals):
         if v.isdigit():
            isdigits.append(len(vals)-i)
   wherenumber = -1*numpy.bincount(numpy.array(isdigits)).argmax()      # what part of filename is a number?
   print "I'm assuming the file numbering is in the %dth to the last part of the file name" % (numpy.abs(wherenumber))
   filenames = [i for i in filenames if i.split('.')[wherenumber].isdigit()]
   filenames = sorted(filenames, key=lambda i: int(i.split('.')[wherenumber]))
n_files = len(filenames)

# Initialize some variables.
bEnergy, bPV, bExpanded = numpy.zeros(n_files, bool), numpy.zeros(n_files, bool), numpy.zeros(n_files, bool)
header_lines, snap_size, lv, lv_names, n_states = [], [], [], [], []

# Read in the metadata.
for filename in filenames:
   print "Reading metadata from %s..." % filename
   #MRS: may eventually want to make this processing a bit smarter rather than just reading off the first 1200 lines.
   excerpt = getoutput("""head -n 1200 %s | grep '@ s[0-9]' -A 2 -n | tr '\-:/,()"\\' ' '""" % filename).split('\n')
   header_lines.append(int(excerpt[-3].split()[0]))              # Count up the header lines in the file.
   n_components = 0                                              # Initiallize the number of components tracer.
   for line in excerpt:
      elements = line.split()
      if 'legend' in elements:
         if 'dH' in elements:
            lv_names.append(elements[8])                         # Store the lambda type names.
            n_components += 1                                    # Count up the energy components.
         if 'xD' in elements:
            lv.append(elements[-n_components:])                  # Store the lambda type values.

         if 'Energy' in elements:                                # Keep track of the booleans for                               
            bEnergy[len(n_states)] = True                        #   the total energy,
         if 'pV' in elements:
            bPV[len(n_states)] = True                            #   the PV energy,
         if 'state' in elements:                                 
            bExpanded[len(n_states)] = True                      #   and file's type (regular vs. expanded ensemble).
      else:
         snap_size.append(float(elements[1]))                    # Store the time to determine the snapshot size.

   n_states.append(len(lv) - sum(n_states))                      # Keep track of the dE columns' number.

#===================================================================================================
# Preliminaries II: Analyze the information for validity; raise errors if the data is flawed.
#===================================================================================================

# Check whether the 'n_components' are the same.
try:
   lv = numpy.array(lv, float)
except ValueError:
   parser.error("\nFiles do not contain the same number of lambda gradient components; I cannot combine the files.")

# Check whether the 'lv_names' are the same.
lv_names = numpy.array(lv_names).reshape(n_files, n_components)
if not (lv_names[0] == lv_names).all():
   parser.error("\nThe lambda gradient components, though the same in number, have different names; I cannot combine the data.")
lv_names = lv_names[0].tolist()

# Figure out the unique elements of 'n_states' (the allowed scenarios are [K], [2,3], and [2,3,K], [2,K], [3,K]) and build up the proper 'lv'.
ncol = numpy.unique(n_states)

# Scenario #1: Each file has all the dE columns -- can use MBAR.
if len(ncol) == 1:
   if (lv.reshape(n_files, n_states[0], n_components) == lv[:n_states[0]]).all():
      lv = lv[:n_states[0]]
   else:
      parser.error("\nThe lambda gradient components, though the same in number, have different values; I cannot combine the data.")

elif len(ncol) <= 3:
   # Scenario #2: Have the adjacent states only; 2 dE columns for the terminal states, 3 for inner ones.
   if ncol.tolist() == [2, 3]:
      lv  = numpy.array([lv[i] for i in range(len(lv)) if i%3==0])
   # Scenario #3: Have a mixture of formats (adjacent and all): either [2,3,K], or [2,K], or [3,K].
   else:
      ndx = n_states.index(ncol[-1])
      lv = lv[sum(n_states[:ndx]) : sum(n_states[:ndx+1])]
   if 'MBAR' in methods:
      print "\nNumber of states is NOT the same for all simulations; I'm assuming that we only evaluate"
      print "nearest neighbor states, and so cannot use MBAR, removing the method."
      methods.remove('MBAR')
   print "\nStitching together the dhdl files. I am assuming that the files are numbered in order of"
   print "increasing lambda; otherwise, results will not be correct."

else:
   print "The files contain the number of the dE columns I cannot deal with; will terminate.\n\n%-10s %s " % ("# of dE's", "File")
   for nf in range(n_files):
      print "%6d     %s" % (n_states[nf], filenames[nf])
   parser.error("\nThere are more than 3 groups of files (%s, to be exact) each having different number of the dE columns; I cannot combine the data." % len(ncol))

# Assign the proper position of each file's dE columns in the u_klt array (will be of service in Prelim. IV).
Kmap = [numpy.arange(n_states[i]) + i-1 if (i>0 and n_states[i]<=3) else numpy.arange(n_states[i]) for i in range(len(n_states))]

# If some contain total energy, and some do not, ignore the total energies.
if (bEnergy == bEnergy[0]).all():
   bUseTotalEnergy = True
else:
   bUseTotalEnergy = False
   print "\nWARNING: Some files contain total energies, some do not; I will ignore the total energies to be able to combine all data.\n"

# We won't combine if bPV isn't all the same (wrong ensemble).
if (bPV == bPV[0]).all():
   bPV = bPV[0]
else:
   parser.error("\nSome files contain the PV energies, some do not; I cannot combine the files.")

bUseTotalEnergy = False
#===================================================================================================
# Preliminaries III: Count up the equilibrated snapshots; initialize arrays to be loaded later on.
#===================================================================================================

# Initialize an array to store the number of lines from the top that are not to be read in.
skip_lines = numpy.zeros(n_files, int)

# In the .xvg files, the snapshots are stored every 'snap_size' ps.
snap_size = numpy.diff(numpy.array(snap_size).reshape(n_files, 2).transpose(), axis=0)[0]

# Create a list that stores information on whether the last line of the file is corrupted (will be of service in Prelim. IV).
last_line = [None for i in range(n_files)]

# Initialize an array to store the number of equilibrated snapshots per each state.
K = len(lv)
nsnapshots = numpy.zeros(K, int)

for nf in range(n_files):
   
   # Check if the last line in the .xvg file is cutoff and do not count it if so.
   len1, len2 = getoutput("tail -n 2 %s | awk '{print NF}'" % filenames[nf]).split()
   if len1 != len2:
      if bExpanded[nf]:
         last_line[nf] = -1 # Ignore the last line when the file is read.
      else: 
         nsnapshots[nf] -= 1 # Update the number of equilibrated snapshots.
         last_line[nf] = -1 # Ignore the last line when the file is read.
      
   # If it is ee, find the information on when the WL weights equilibration was reached and discard to the greater of WLequiltime and equiltime.
   if bExpanded[nf]:
      logfilename = filenames[nf].replace('.xvg', '.log')
      if not(options.ignoreWL):
         if not os.path.isfile(logfilename):
            parser.error('\nExpanded ensemble also requires log files as input, and the file \'%s\' was not found.\nYou may rerun with the -x flag and the data will be discarded to \'equiltime\', not bothering\nwith the extraction of the information on when the WL weights equilibration was reached.\nOtherwise, put the proper log file into the directory which is subject to the analysis.' % logfilename)
         try:
            WLequilstep, WLequiltime = getoutput("grep equilibrated -A 8 -B 8 %s | grep Time -A 1 | tail -n 1 | awk '{print $1, $2}'" % logfilename).split()
         except ValueError:
            parser.error("\nThe Wang-Landau weights haven't equilibrated yet.\nIf it is OK with you, rerun with the -x flag and the data will be discarded to \'equiltime\'.")
         equiltime = max(float(WLequiltime), equiltime)

      # Convert equilibration time into equilibration snapshots and estimate how many lines in each file are to be skipped.
      equilsnapshots = int(equiltime/snap_size[nf])
      skip_lines[nf] = header_lines[nf] + equilsnapshots

      # Calculate how many times each state is encountered after the equilibration is reached.
      extract_states = numpy.genfromtxt(filenames[nf], dtype=int, skiprows=skip_lines[nf], skip_footer=(1 if last_line[nf] == -1 else 0), usecols=1)
      nsnapshots_dict = Counter(extract_states)
      for k in range(K):
         nsnapshots[k] += nsnapshots_dict[k]
   else:
      equilsnapshots = int(equiltime/snap_size[nf])
      skip_lines[nf] = header_lines[nf] + equilsnapshots

      nsnapshots[nf] += int(getoutput("awk 'END {print NR}' %s" % filenames[nf])) # 'wc -l' defines a line as an entry ending with "\n", while we need to count any entry
      nsnapshots[nf] -= skip_lines[nf]

   print "First %s ps (%s snapshots) will be discarded due to equilibration from file %s..." % (equiltime, equilsnapshots, filenames[nf])

# Determine maximum number of the equilibrated snapshots from any state and initialize arrays with the third dimension of 'maxn'.
maxn = max(nsnapshots)
dhdlt = numpy.zeros([K,n_components,int(maxn)], float)    # dhdlt[k,n,t] is the derivative of energy component n with respect to state k of snapshot t
u_klt = numpy.zeros([K,K,int(maxn)], numpy.float64)       # u_klt[k,m,t] is the reduced potential energy of snapshot t of state k evaluated at state m

#===================================================================================================
# Preliminaries IV: Load all the intact data, skipping snapshots prior to equilibration.
#===================================================================================================   

def equilibratedLines(filename):
   """Returns the equilibrated portion of 'filename' by discarding the header, the non-equilibrated
   (at the top) and corrupted (the very last line) snapshots if any."""
   print "Loading in data from %s..." % filename
   infile = open(filename, 'r')
   lines = infile.readlines()[skip_lines[nf] : last_line[nf]]
   infile.close()
   return lines

nss = numpy.zeros(K, int) # snapshot tracer

# Load all of the data.
for nf in range(n_files):

   # The default settings for the state, the total and PV energies.
   state, energy, pv = nf, 0, 0

   for line in equilibratedLines(filenames[nf]):
      
      # Split line into elements and get rid of the first one.
      elements = line.split()
      time = elements.pop(0)

      # Update the default settings if there is a need for it.
      if bExpanded[nf]:
         state = int(elements.pop(0))
      if bEnergy[nf]:
         if bUseTotalEnergy:
            energy = float(elements.pop(0))
         else:
            rubbish = elements.pop(0)
      if bPV:
         pv = float(elements.pop(-1))

      # Store the derivative with respect to lambda.
      for nl in range(n_components):
         dhdlt[state,nl,nss[state]] = float(elements.pop(0))

      # Compute and store reduced potential energy at each state.
      for l in range(n_states[nf]):
         u_klt[state,Kmap[nf][l],nss[state]] = beta * (float(elements.pop(0)) + energy + pv)

      # Update the snapshot index tracer.
      nss[state] += 1

#===================================================================================================
# The functions we'll be using presently.
#===================================================================================================   

def getNkandUkln(do_dhdl=False):
   """Identifies uncorrelated samples and updates the arrays of the reduced potential energy and dhdlt retaining data entries of these samples only.
      Assumes that 'dhdlt' and 'u_klt' are in memory, as well as proper values for 'sta' and 'fin', i.e. the starting and
      final snapshot positions to be read, both are arrays of dimension K."""
   u_kln = numpy.zeros([K,K,max(fin-sta)], numpy.float64) # u_kln[k,m,n] is the reduced potential energy of uncorrelated sample index n from state k evaluated at state m
   N_k = numpy.zeros(K, int) # N_k[k] is the number of uncorrelated samples from state k
   g = numpy.zeros(K,float) # autocorrelation times for the data
   if do_dhdl:
      dhdl = numpy.zeros([K,n_components,max(fin-sta)], float) #dhdl is value for dhdl for each component in the file at each time.
      print "\n\nNumber of correlated and uncorrelated samples:\n\n%6s %12s %12s %12s\n" % ('State', 'N', 'N_k', 'N/N_k') 
   for k in range(K):
      # Sum up over the energy components; notice, that only the relevant data is being used in the third dimension.
      dhdl_sum = numpy.sum(dhdlt[k,:,sta[k]:fin[k]], axis=0)
      # Determine indices of uncorrelated samples from potential autocorrelation analysis at state k
      # (alternatively, could use the energy differences -- here, we will use total dhdl).
      g[k] = timeseries.statisticalInefficiency(dhdl_sum)
      indices = numpy.array(timeseries.subsampleCorrelatedData(dhdl_sum, g=g[k])) # indices of uncorrelated samples
      N = len(indices) # number of uncorrelated samples
      # Handle case where we end up with too few.
      if N < 50:
         if do_dhdl:
            print "WARNING: Only %s uncorrelated samples found at lambda number %s; proceeding with analysis using correlated samples..." % (N, k)
         indices = numpy.arange(len(dhdl_sum))
         N = len(indices)
      N_k[k] = N # Store the number of uncorrelated samples from state k.
      for l in range(K):
         u_kln[k,l,0:N] = u_klt[k,l,indices]
      if do_dhdl:
         print "%6s %12s %12s %12.2f" % (k, fin[k], N_k[k], g[k])
         for n in range(n_components): 
            dhdl[k,n,0:N] = dhdlt[k,n,indices]
   if do_dhdl:
      return (dhdl, N_k, u_kln)
   return (N_k, u_kln)

def plotOverlapMatrix(O):
   """Plots the probability of observing a sample from state i (row) in state j (column).
   For convenience, the neigboring state cells are fringed in bold."""
   max_prob = O.max()
   fig = pl.figure(figsize=(K/2.,K/2.))
   fig.add_subplot(111, frameon=False, xticks=[], yticks=[])

   for i in range(K):
      if i!=0:
         pl.axvline(x=i, ls='-', lw=0.5, color='k', alpha=0.25)
         pl.axhline(y=i, ls='-', lw=0.5, color='k', alpha=0.25)
      for j in range(K):
         if O[j,i] < 0.005:
            ii = ''
         else:
            ii = ("%.2f" % O[j,i])[1:]
         alf = O[j,i]/max_prob
         pl.fill_between([i,i+1], [K-j,K-j], [K-(j+1),K-(j+1)], color='k', alpha=alf)
         pl.annotate(ii, xy=(i,j), xytext=(i+0.5,K-(j+0.5)), size=8, textcoords='data', va='center', ha='center', color=('k' if alf < 0.5 else 'w'))

   cx = sorted(2*range(K+1))
   cy = sorted(2*range(K+1), reverse=True)
   pl.plot(cx[2:-1], cy[1:-2], 'k-', lw=2.0)
   pl.plot(numpy.array(cx[2:-3])+1, cy[1:-4], 'k-', lw=2.0)
   pl.plot(cx[1:-2], numpy.array(cy[:-3])-1, 'k-', lw=2.0)
   pl.plot(cx[1:-4], numpy.array(cy[:-5])-2, 'k-', lw=2.0)

   pl.xlim(0, K)
   pl.ylim(0, K)
   pl.savefig('O_MBAR.pdf', bbox_inches='tight', pad_inches=0.0)
   pl.close(fig)
   return

def estimatewithMBAR(reltol=relative_tolerance, regular_estimate=False):
   """Computes the MBAR free energy given the reduced potential and the number of relevant entries in it.
      Assumes that proper arrays for u_kln and N_k are in memory."""
   MBAR = pymbar.MBAR(u_kln, N_k, verbose = verbose, method = 'adaptive', relative_tolerance = reltol, initialize = initialize_with)
   # Get matrix of dimensionless free energy differences and uncertainty estimate.
   (Deltaf_ij, dDeltaf_ij) = MBAR.getFreeEnergyDifferences(uncertainty_method='svd-ew')
   if (verbose): 
      print "Matrix of free energy differences\nDeltaf_ij:\n%s\ndDeltaf_ij:\n%s" % (Deltaf_ij, dDeltaf_ij)
      print "The overlap matrix is..."
      O = MBAR.computeOverlap('matrix')
      plotOverlapMatrix(O)
      for k in range(K):
         line = ''
         for l in range(K):
            line += ' %5.2f ' % O[k, l]
         print line
   if regular_estimate:
      return (Deltaf_ij, dDeltaf_ij)
   return (Deltaf_ij[0,K-1]/beta_report, dDeltaf_ij[0,K-1]/beta_report)

#===================================================================================================
# Calculate free energies with different methods.
#===================================================================================================    

# The starting and final snapshot positions to be read from dhdlt and u_klt are:
sta = numpy.zeros(K, int)
fin = nsnapshots
# Obtain uncorrelated samples, as well as dhdl (for the TI methods) and u_kln (for all others).
dhdl, N_k, u_kln = getNkandUkln(do_dhdl=True)
# Estimate free energy difference with MBAR -- all states at once.
if 'MBAR' in methods:
   print "\nEstimating the free energy change with MBAR..."
   Deltaf_ij, dDeltaf_ij = estimatewithMBAR(regular_estimate=True)

### The TI preliminaries ###
if any(m in methods for m in ['TI', 'TI-CUBIC']):
   # Compute <dhdl> and std(dhdl) for each component, for each lambda; multiply them by beta to make unitless.
   ave_dhdl = numpy.zeros([K,n_components],float)
   std_dhdl = numpy.zeros([K,n_components],float)
   for k in range(K):
      ave_dhdl[k,:] = beta*numpy.average(dhdl[k,:,0:N_k[k]],axis=1)
      std_dhdl[k,:] = beta*numpy.std(dhdl[k,:,0:N_k[k]],axis=1)/numpy.sqrt(N_k[k]-1)

   # Lambda vectors spacing.
   dlam = numpy.diff(lv, axis=0)

   lchange = numpy.zeros([K,n_components],bool)   # booleans for which lambdas are changing 
   for j in range(n_components):
      # need to identify range over which lambda doesn't change, and not interpolate over that range.
      for k in range(K-1):
         if (lv[k+1,j]-lv[k,j] > 0):
            lchange[k,j] = True
            lchange[k+1,j] = True

### The TI-CUBIC preliminaries ###
class naturalcubicspline:

   def __init__(self, x):

      # define some space
      L = len(x)
      H = numpy.zeros([L,L],float)
      M = numpy.zeros([L,L],float)
      BW = numpy.zeros([L,L],float)
      AW = numpy.zeros([L,L],float)
      DW = numpy.zeros([L,L],float)

      h = x[1:L]-x[0:L-1]
      ih = 1.0/h

      # define the H and M matrix, from p. 371 "applied numerical methods with matlab, Chapra"
      H[0,0] = 1
      H[L-1,L-1] = 1
      for i in range(1,L-1):
         H[i,i] = 2*(h[i-1]+h[i])
         H[i,i-1] = h[i-1]
         H[i,i+1] = h[i]

         M[i,i] = -3*(ih[i-1]+ih[i])
         M[i,i-1] = 3*(ih[i-1])
         M[i,i+1] = 3*(ih[i])

      CW = numpy.dot(numpy.linalg.inv(H),M)  # this is the matrix translating c to weights in f.
                                                   # each row corresponds to the weights for each c.

      # from CW, define the other coefficient matrices
      for i in range(0,L-1):
         BW[i,:]    = -(h[i]/3)*(2*CW[i,:]+CW[i+1,:])
         BW[i,i]   += -ih[i]
         BW[i,i+1] += ih[i]
         DW[i,:]    = (ih[i]/3)*(CW[i+1,:]-CW[i,:])
         AW[i,i]    = 1

      # Make copies of the arrays we'll be using in the future.
      self.x  = x.copy()
      self.AW = AW.copy()
      self.BW = BW.copy()
      self.CW = CW.copy()
      self.DW = DW.copy()

      # find the integrating weights
      self.wsum = numpy.zeros([L],float)
      self.wk = numpy.zeros([L-1,L],float)
      for k in range(0,L-1):
         w = DW[k,:]*(h[k]**4)/4.0 + CW[k,:]*(h[k]**3)/3.0 + BW[k,:]*(h[k]**2)/2.0 + AW[k,:]*(h[k])
         self.wk[k,:] = w
         self.wsum += w

   def interpolate(self,y,xnew):
      if len(self.x) != len(y):
         parser.error("\nThe length of 'y' should be consistent with that of 'self.x'. I cannot perform linear algebra operations.")
      # get the array of actual coefficients by multiplying the coefficient matrix by the values  
      a = numpy.dot(self.AW,y)
      b = numpy.dot(self.BW,y)
      c = numpy.dot(self.CW,y)
      d = numpy.dot(self.DW,y)

      N = len(xnew)
      ynew = numpy.zeros([N],float)
      for i in range(N):
         # Find the index of 'xnew[i]' it would have in 'self.x'.
         j = numpy.searchsorted(self.x, xnew[i]) - 1
         lamw = xnew[i] - self.x[j]
         ynew[i] = d[j]*lamw**3 + c[j]*lamw**2 + b[j]*lamw + a[j]
      # Preserve the terminal points.
      ynew[0] = y[0]
      ynew[-1] = y[-1]
      return ynew

if 'TI-CUBIC' in methods:
   # construct a map back to the original components
   mapl = numpy.zeros([K,n_components],int)   # map back to the original k from the components
   for j in range(n_components):
      incr = 0
      for k in range(K):
         if (lchange[k,j]):
            mapl[k,j] += incr
            incr +=1

   # put together the spline weights for the different components
   cubspl = list()
   for j in range(n_components):
      lv_lchange = lv[lchange[:,j],j]
      if len(lv_lchange) == 0: # handle the all-zero lv column
         cubspl.append(0)
      else:
         spl = naturalcubicspline(lv_lchange)
         cubspl.append(spl)
### end of the TI-CUBIC preliminaries ###

print ("Estimating the free energy change with %s..." % ', '.join(methods)).replace(', MBAR', '')
df_allk = list(); ddf_allk = list()

for k in range(K-1):
   df = dict(); ddf = dict()

   for name in methods:

      if name == 'TI':
         #===================================================================================================
         # Estimate free energy difference with TI; interpolating with the trapezoidal rule.
         #===================================================================================================   
         df['TI'] = 0.5*numpy.dot(dlam[k],(ave_dhdl[k]+ave_dhdl[k+1]))        
         ddf['TI'] = 0.5*numpy.sqrt(numpy.dot(dlam[k]**2,std_dhdl[k]**2+std_dhdl[k+1]**2))               

      if name == 'TI-CUBIC':
         #===================================================================================================
         # Estimate free energy difference with TI; interpolating with the natural cubic splines.
         #===================================================================================================   
         df['TI-CUBIC'], ddf['TI-CUBIC'] = 0, 0
         for j in range(n_components):
            if dlam[k,j] > 0:
               lj = lchange[:,j]
               df['TI-CUBIC'] += numpy.dot(cubspl[j].wk[mapl[k,j]],ave_dhdl[lj,j])
               ddf['TI-CUBIC'] += numpy.dot(cubspl[j].wk[mapl[k,j]]**2,std_dhdl[lj,j]**2)
         ddf['TI-CUBIC'] = numpy.sqrt(ddf['TI-CUBIC'])

      if any(name == m for m in ['DEXP', 'GDEL', 'BAR', 'UBAR', 'RBAR']):
         w_F = u_kln[k,k+1,0:N_k[k]] - u_kln[k,k,0:N_k[k]] 

      if name == 'DEXP':
         #===================================================================================================
         # Estimate free energy difference with Forward-direction EXP (in this case, Deletion from solvent).
         #===================================================================================================   
         (df['DEXP'], ddf['DEXP']) = EXP(w_F)

      if name == 'GDEL':
         #===================================================================================================
         # Estimate free energy difference with a Gaussian estimate of EXP (in this case, deletion from solvent)
         #===================================================================================================   
         (df['GDEL'], ddf['GDEL']) = EXPGauss(w_F)

      if any(name == m for m in ['IEXP', 'GINS', 'BAR', 'UBAR', 'RBAR']):
         w_R = u_kln[k+1,k,0:N_k[k+1]] - u_kln[k+1,k+1,0:N_k[k+1]] 

      if name == 'IEXP':
         #===================================================================================================
         # Estimate free energy difference with Reverse-direction EXP (in this case, insertion into solvent).
         #===================================================================================================   
         (rdf,rddf) = EXP(w_R)
         (df['IEXP'], ddf['IEXP']) = (-rdf,rddf)

      if name == 'GINS':
         #===================================================================================================
         # Estimate free energy difference with a Gaussian estimate of EXP (in this case, insertion into solvent)
         #===================================================================================================   
         (rdf,rddf) = EXPGauss(w_R)
         (df['GINS'], ddf['GINS']) = (-rdf,rddf)

      if name == 'BAR':
         #===================================================================================================
         # Estimate free energy difference with BAR; use w_F and w_R computed above.
         #===================================================================================================   
         (df['BAR'], ddf['BAR']) = pymbar.bar.BAR(w_F, w_R, relative_tolerance=relative_tolerance, verbose = verbose)      

      if name == 'UBAR':
         #===================================================================================================
         # Estimate free energy difference with unoptimized BAR -- assume dF is zero, and just do one evaluation
         #===================================================================================================   
         (df['UBAR'], ddf['UBAR']) = pymbar.bar.BAR(w_F, w_R, verbose = verbose,iterated_solution=False)

      if name == 'RBAR':
         #===================================================================================================
         # Estimate free energy difference with Unoptimized BAR over range of free energy values, and choose the one 
         # that is self consistently best.
         #===================================================================================================   
         min_diff = 1E6
         best_udf = 0
         for trial_udf in range(-10,10,1):
            (udf, uddf) = pymbar.bar.BAR(w_F, w_R, DeltaF=trial_udf, iterated_solution=False, verbose=verbose)
            diff = numpy.abs(udf - trial_udf)
            if (diff < min_diff):
               best_udf = udf
               best_uddf = uddf
               min_diff = diff
         (df['RBAR'], ddf['RBAR']) = (best_udf,best_uddf)

      if name == 'MBAR':
         #===================================================================================================
         # Store the MBAR free energy difference (already estimated above) properly, i.e. by state.
         #===================================================================================================   
         (df['MBAR'], ddf['MBAR']) =  Deltaf_ij[k,k+1], dDeltaf_ij[k,k+1]

   df_allk = numpy.append(df_allk,df)
   ddf_allk = numpy.append(ddf_allk,ddf)

#===================================================================================================
# If the proper flag was used, plot the forward and reverse dF(t).
#===================================================================================================   

def plotdFvsTime(f_ts, r_ts, F_df, R_df, F_ddf, R_ddf):
   """Plots the free energy change computed using the equilibrated snapshots between the proper target time frames (f_ts and r_ts)
   in both forward (data points are stored in F_df and F_ddf) and reverse (data points are stored in R_df and R_ddf) directions."""
   fig = pl.figure(figsize=(8,6))
   ax = fig.add_subplot(111)
   pl.setp(ax.spines['bottom'], color='#D2B9D3', lw=3, zorder=-2)
   for dir in ['left', 'top', 'right']:
      ax.spines[dir].set_color('none')
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('left')
   ax.spines['bottom'].set_position(('outward',16))
   ax.spines['left'].set_position(('outward',-16))

   line0  = pl.fill_between([r_ts[0], f_ts[-1]], R_df[0]-R_ddf[0], R_df[0]+R_ddf[0], color='#D2B9D3', zorder=-1)
   for i in range(len(f_ts)):
      line1 = pl.plot([f_ts[i]]*2, [F_df[i]-F_ddf[i], F_df[i]+F_ddf[i]], color='#736AFF', ls='-', lw=3, solid_capstyle='round', zorder=1)
   line11 = pl.plot(f_ts, F_df, color='#736AFF', ls='-', lw=3, marker='o', mfc='w', mew=2.5, mec='#736AFF', ms=12, zorder=2)
   for i in range(len(r_ts)):
      line2 = pl.plot([r_ts[i]]*2, [R_df[i]-R_ddf[i], R_df[i]+R_ddf[i]], color='#C11B17', ls='-', lw=3, solid_capstyle='round', zorder=3)
   line22 = pl.plot(r_ts, R_df, color='#C11B17', ls='-', lw=3, marker='o', mfc='w', mew=2.5, mec='#C11B17', ms=12, zorder=4)

   miny, maxy = ax.get_ylim()
   linev1 = pl.plot([r_ts[0], r_ts[0]], [miny, maxy], color='#736AFF', ls='-', lw=3, zorder=-4)
   linev2 = pl.plot([f_ts[-1], f_ts[-1]], [miny, maxy], color='#C11B17', ls='-', lw=3, zorder=-3)

   xshift = (f_ts[-1] - r_ts[0])/25.
   pl.xlim(r_ts[0]-xshift, f_ts[-1]+xshift)

   pl.xticks(r_ts[::2] + f_ts[-1:], fontsize=8)
   pl.yticks(fontsize=8)

   pl.title(r'$\mathrm{The\ free\ energy\ change\ %s\ as\ a\ function\ of\ time\ (ps)}$' % units, fontsize=12)
   leg = pl.legend((line1[0], line2[0]), (r'$Forward$', r'$Reverse$'), loc=9, prop=FP(size=10), frameon=False)
   pl.tick_params(axis='x', color='#D2B9D3')
   pl.tick_params(axis='y', color='#736AFF')
   pl.savefig(os.path.join(datafile_directory, 'dF_t.pdf'))
   pl.close(fig)
   return

if options.bForwrev:
   if not 'MBAR' in methods:
      parser.error("\nCurrent version of the dF(t) analysis works with MBAR only and the method is not found in the list.")
   if not (snap_size[0] == snap_size).all():
      parser.error("\nThe snapshot size isn't the same for all the files; cannot perform the dF(t) analysis.")
   # Define a list of 10 equidistant time frames at which the free energy is to be estimated; count up the snapshots embounded between the time frames.
   n_tf = 11
   nss_tf = numpy.zeros([n_tf, K], int)
   if bExpanded[nf]:
      tf = numpy.arange(0,1.1,0.1)*(numpy.sum(nsnapshots)-1)+1
      tf[0] = 0
      for i in range(n_tf-1):
         nss = Counter(extract_states[tf[i]:tf[i+1]])       
         nss_tf[i+1] = numpy.array([nss[j] for j in range(K)])
   else:
      tf = numpy.arange(0,1.1,0.1)*(max(nsnapshots)-1)+1
      tf[0] = 0
      for i in range(n_tf-1):
         nss_tf[i+1] = numpy.array([min(j, tf[i+1]) for j in nsnapshots]) - numpy.sum(nss_tf[:i+1],axis=0)

   # Define the real time scale (in ps) rather than a snapshot sequence.
   ts = ["%.1f" % ((i-1)*snap_size[0] + equiltime) if i!=tf[0] else "%.1f" % (i*snap_size[0] + equiltime) for i in tf]
   # Initialize arrays to store data points to be plotted.
   F_df  = numpy.zeros(n_tf-1, float)
   F_ddf = numpy.zeros(n_tf-1, float)
   R_df  = numpy.zeros(n_tf-1, float)
   R_ddf = numpy.zeros(n_tf-1, float)
   # Store the MBAR energy that accounts for all the equilibrated snapshots (has already been computed in the previous section).
   F_df[-1], F_ddf[-1] = (Deltaf_ij[0,K-1]/beta_report, dDeltaf_ij[0,K-1]/beta_report)
   R_df[0], R_ddf[0]   = (Deltaf_ij[0,K-1]/beta_report, dDeltaf_ij[0,K-1]/beta_report)
   # Do the forward analysis.
   print "Forward dF(t) analysis...\nEstimating the free energy change using the data up to"
   sta = nss_tf[0]
   for i in range(n_tf-2):
      print "%60s ps..." % ts[i+1]
      fin = numpy.sum(nss_tf[:i+2],axis=0)
      N_k, u_kln = getNkandUkln()
      F_df[i], F_ddf[i] = estimatewithMBAR()
   # Do the reverse analysis.
   print "Reverse dF(t) analysis...\nUsing the data starting from"
   fin = numpy.sum(nss_tf[:],axis=0)
   for i in range(n_tf-2):
      print "%34s ps..." % ts[i+1]
      sta = numpy.sum(nss_tf[:i+2],axis=0)
      N_k, u_kln = getNkandUkln()
      R_df[i+1], R_ddf[i+1] = estimatewithMBAR()

   print """\n   The free energies %s evaluated by using the trajectory
   snaphots corresponding to various time intervals for both the
   reverse and forward (in parentheses) direction.\n""" % units
   print "%s\n %20s %19s %20s\n%s" % (70*'-', 'Time interval, ps','Reverse', 'Forward', 70*'-')
   print "%10s -- %s\n%10s -- %-10s %11.3f +- %5.3f %16s\n" % (ts[0], ts[-1], '('+ts[0], ts[0]+')', R_df[0], R_ddf[0], 'XXXXXX')
   for i in range(1, len(ts)-1):
      print "%10s -- %s\n%10s -- %-10s %11.3f +- %5.3f %11.3f +- %5.3f\n" % (ts[i], ts[-1], '('+ts[0], ts[i]+')', R_df[i], R_ddf[i], F_df[i-1], F_ddf[i-1])
   print "%10s -- %s\n%10s -- %-10s %16s %15.3f +- %5.3f\n%s" % (ts[-1], ts[-1], '('+ts[0], ts[-1]+')', 'XXXXXX', F_df[-1], F_ddf[-1], 70*'-')

   # Plot the forward and reverse dF(t); store the data points in the text file.
   print "Plotting data to the file dF_t.pdf...\n\n"
   plotdFvsTime([float(i) for i in ts[1:]], [float(i) for i in ts[:-1]], F_df, R_df, F_ddf, R_ddf)
   outtext = ["%12s %10s %-10s %17s %10s %s\n" % ('Time (ps)', 'Forward', units, 'Time (ps)', 'Reverse', units)]
   outtext+= ["%10s %11.3f +- %5.3f %18s %11.3f +- %5.3f\n" % (ts[1:][i], F_df[i], F_ddf[i], ts[:-1][i], R_df[i], R_ddf[i]) for i in range(len(F_df))]
   file = open(os.path.join(datafile_directory, 'dF_t.txt'), 'w'); file.writelines(outtext); file.close()

#===================================================================================================
# All done with calculations, now summarize and print stats.
#===================================================================================================   

# Count up the charging states.
numcharging = 0
for lv_n in ['coul', 'fep']:
   if lv_n in lv_names:
      ndx_char = lv_names.index(lv_n)
      lv_char = lv[:, ndx_char]
      if not (lv_char == lv_char[0]).all():
         numcharging = (lv_char != 1).sum()
         break

if numcharging == K:
   numcharging = K-1

# Split the total energies into segments; initialize lists to store them.
segments      = ['Coulomb'  , 'vdWaals'  , 'TOTAL']
segmentstarts = [0          , numcharging, 0      ]
segmentends   = [numcharging, K-1        , K-1    ]
dFs = []; ddFs = []

# Perform the energy segmentation; be pedantic about the TI cumulative ddF's (see the explanation in Section 3.1 of the paper).
for i in range(len(segments)):
   segment = segments[i]; segstart = segmentstarts[i]; segend = segmentends[i]
   dF  = dict.fromkeys(methods, 0)
   ddF = dict.fromkeys(methods, 0)

   for name in methods:
      if name == 'MBAR':
         dF['MBAR']  =  Deltaf_ij[segstart, segend]
         ddF['MBAR'] = dDeltaf_ij[segstart, segend]

      elif name[0:2] == 'TI':
         for k in range(segstart, segend):
            dF[name] += df_allk[k][name]

         if segment == 'Coulomb':
            jlist = [ndx_char] if numcharging>0 else []
         elif segment == 'vdWaals':
            jlist = []
         elif segment == 'TOTAL':
            jlist = range(n_components)

         for j in jlist:
            lj = lchange[:,j]
            if not (lj == False).all(): # handle the all-zero lv column
               if name == 'TI-CUBIC':
                  ddF[name] += numpy.dot((cubspl[j].wsum)**2,std_dhdl[lj,j]**2)
               elif name == 'TI':
                  h = numpy.trim_zeros(dlam[:,j])
                  wsum = 0.5*(numpy.append(h,0) + numpy.append(0,h))
                  ddF[name] += numpy.dot(wsum**2,std_dhdl[lj,j]**2)
         ddF[name] = numpy.sqrt(ddF[name])

      else:
         for k in range(segstart,segend):
            dF[name] += df_allk[k][name]
            ddF[name] += (ddf_allk[k][name])**2
         ddF[name] = numpy.sqrt(ddF[name])

   dFs.append(dF)
   ddFs.append(ddF)

for name in methods: # 'vdWaals' = 'TOTAL' - 'Coulomb'
   ddFs[1][name] = (ddFs[2][name]**2 - ddFs[0][name]**2)**0.5

# Display results.
def printLine(str1, str2, d1=None, d2=None):
   """Fills out the results table linewise."""
   print str1,
   text = str1
   for name in methods:
      if d1 == 'plain':
         print str2,
         text += ' ' + str2
      if d1 == 'name':
         print str2 % (name, units),
         text += ' ' + str2 % (name, units)
      if d1 and d2:
         print str2 % (d1[name]/beta_report, d2[name]/beta_report),
         text += ' ' + str2 % (d1[name]/beta_report, d2[name]/beta_report)
   print ''
   outtext.append(text + '\n')
   return

def plotTI():
   """Plots the ave_dhdl as a function of lambda.
   If (TI and TI-CUBIC in methods) -- plots the TI integration area and the TI-CUBIC interpolation curve,
   elif (only one of them in methods) -- plots the integration area of the method."""
   min_dl = dlam[dlam != 0].min()
   S = int(0.4/min_dl)
   if S>8:
      fig = pl.figure(figsize = (S,6))
   else:
      fig = pl.figure(figsize = (8,6))
   ax = fig.add_subplot(1,1,1)
   ax.spines['bottom'].set_position('zero')
   ax.spines['top'].set_color('none')
   ax.spines['right'].set_color('none')
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('left')
  
   xs, ndx, dx = [0], 0, 0.001
   colors = ['r', 'g', '#7F38EC', '#9F000F', 'b', 'y']
   min_y, max_y = 0, 0
  
   for j in range(n_components):
      y = ave_dhdl[:,j]
      if not (y == 0).all():
 
         # Get the coordinates.
         lj = lchange[:,j]
         x = lv[:,j][lj]
         y = y[lj]/beta_report
 
         if 'TI' in methods:
            # Plot the TI integration area.
            ss = 'TI'
            for i in range(len(x)-1):
               min_y = min(y.min(), min_y)
               max_y = max(y.max(), max_y)
               if i%2==0:
                  pl.fill_between(x[i:i+2]+ndx, 0, y[i:i+2], color=colors[ndx], alpha=1.0)
               else:
                  pl.fill_between(x[i:i+2]+ndx, 0, y[i:i+2], color=colors[ndx], alpha=0.5)

            if 'TI-CUBIC' in methods:
               # Plot the TI-CUBIC interpolation curve.
               ss += ' and TI-CUBIC'
               xnew = numpy.arange(0, 1+dx, dx)
               ynew = cubspl[j].interpolate(y, xnew)
               min_y = min(ynew.min(), min_y)
               max_y = max(ynew.max(), max_y)
               pl.plot(xnew+ndx, ynew, color='#B6B6B4', ls ='-', solid_capstyle='round', lw=2.0)
        
         else:
            # Plot the TI-CUBIC integration area.
            ss = 'TI-CUBIC'
            for i in range(len(x)-1):
               xnew = numpy.arange(x[i], x[i+1]+dx, dx)
               ynew = cubspl[j].interpolate(y, xnew)
               ynew[0], ynew[-1] = y[i], y[i+1]
               min_y = min(ynew.min(), min_y)
               max_y = max(ynew.max(), max_y)
               if i%2==0:
                  pl.fill_between(xnew+ndx, 0, ynew, color=colors[ndx], alpha=1.0)
               else:
                  pl.fill_between(xnew+ndx, 0, ynew, color=colors[ndx], alpha=0.5)

         # Store the abscissa values and update the subplot index.
         xs += (x+ndx).tolist()[1:]
         ndx += 1

   # Make sure the tick labels are not overcrowded.
   xt = range(K)
   if S>5:
      i = 0
      while i < len(xs)-1:
         if xs[i+1]-xs[i] < min_dl*1.0001:
            xt[i+1] = ''
            i += 1
         i += 1
   pl.xticks(xs[1:], xt[1:], fontsize=8)
   pl.yticks(fontsize=8)

   # Remove the abscissa ticks and set up the axes limits.
   for tick in ax.get_xticklines():
      tick.set_visible(False)
   pl.xlim(0, ndx)
   min_y = min_y + 0.01*min_y
   max_y = max_y + 0.01*max_y
   pl.ylim(min_y, max_y)

   for i in xs[1:]:
      pl.annotate('%.2f' % (i-1.0 if i>1.0 else i), xy=(i, 0), xytext=(i, 0.01), size=8, rotation=90, textcoords=('data', 'axes fraction'), va='bottom', ha='center', color='#CDCD90')
   if ndx>1:
      lenticks = len(ax.get_ymajorticklabels()) - 1
      if min_y<0: lenticks -= 1
      if lenticks < 5:
         from matplotlib.ticker import AutoMinorLocator as AML
         ax.yaxis.set_minor_locator(AML())
   pl.grid(which='both', color='w', lw=0.25, axis='y', zorder=12)
   s = 'as a function of state, the data interpolated by means of ' + ss
   pl.title(r"$\mathrm{The\ \ } \langle{\frac{ \partial H } { \partial \lambda }}\rangle \mathrm{\ \ %s \ %s}$" % ('\ '.join(s.split()), units), fontsize = 10)
   pl.savefig(os.path.join(datafile_directory, 'dhdl_TI.pdf'), bbox_inches='tight')
   pl.close(fig)
   return

def plotdFvsLambda():
   """Plots the free energy differences evaluated for each pair of adjacent states for all methods."""
   x = numpy.arange(len(df_allk))
   if x[-1]<8:
      fig = pl.figure(figsize = (8,6))
   else:
      fig = pl.figure(figsize = (len(x),6))
   width = 1./(len(methods)+1)
   elw = 30*width
   colors = {'TI':'#C45AEC', 'TI-CUBIC':'#33CC33', 'DEXP':'#F87431', 'IEXP':'#FF3030', 'GINS':'#EAC117', 'GDEL':'#347235', 'BAR':'#6698FF', 'UBAR':'#817339', 'RBAR':'#C11B17', 'MBAR':'#F9B7FF'}
   lines = tuple()
   for name in methods:
      y = [df_allk[i][name]/beta_report for i in x]
      ye = [ddf_allk[i][name]/beta_report for i in x]
      line = pl.bar(x+len(lines)*width, y, width, color=colors[name], yerr=ye, lw=0.1*elw, error_kw=dict(elinewidth=elw, ecolor='black', capsize=0.5*elw))
      lines += (line[0],)
   pl.xlabel('States', fontsize=12, color='#151B54')
   pl.ylabel('$\Delta F$ '+units, fontsize=12, color='#151B54')
   pl.xticks(x+0.5*width*len(methods), tuple(['%d--%d' % (i, i+1) for i in x]), fontsize=8)
   pl.yticks(fontsize=8)
   pl.xlim(x[0], x[-1]+len(lines)*width)
   ax = pl.gca()
   for dir in ['right', 'top', 'bottom']:
      ax.spines[dir].set_color('none')
   ax.yaxis.set_ticks_position('left')
   for tick in ax.get_xticklines():
      tick.set_visible(False)

   leg = pl.legend(lines, tuple(methods), loc=3, ncol=2, prop=FP(size=10), fancybox=True)
   leg.get_frame().set_alpha(0.5)
   pl.title('The free energy change breakdown', fontsize = 12)
   pl.savefig(os.path.join(datafile_directory, 'dF_state_long.pdf'), bbox_inches='tight')
   pl.close(fig)
   return

def plotdFvsLambda2(nb=10):
   """Plots the free energy differences evaluated for each pair of adjacent states for all methods.
   The layout is approximately 'nb' bars per subplot."""
   x = numpy.arange(len(df_allk))
   if len(x) < nb:
      return
   xs = numpy.array_split(x, len(x)/nb+1)
   mnb = max([len(i) for i in xs])
   fig = pl.figure(figsize = (8,6))
   width = 1./(len(methods)+1)
   elw = 30*width
   colors = {'TI':'#C45AEC', 'TI-CUBIC':'#33CC33', 'DEXP':'#F87431', 'IEXP':'#FF3030', 'GINS':'#EAC117', 'GDEL':'#347235', 'BAR':'#6698FF', 'UBAR':'#817339', 'RBAR':'#C11B17', 'MBAR':'#F9B7FF'}
   ndx = 1
   for x in xs:
      lines = tuple()
      ax = pl.subplot(len(xs), 1, ndx)
      for name in methods:
         y = [df_allk[i][name]/beta_report for i in x]
         ye = [ddf_allk[i][name]/beta_report for i in x]
         line = pl.bar(x+len(lines)*width, y, width, color=colors[name], yerr=ye, lw=0.05*elw, error_kw=dict(elinewidth=elw, ecolor='black', capsize=0.5*elw))
         lines += (line[0],)
      for dir in ['left', 'right', 'top', 'bottom']:
         if dir == 'left':
            ax.yaxis.set_ticks_position(dir)
         else:
            ax.spines[dir].set_color('none')
      pl.yticks(fontsize=8)
      ax.xaxis.set_ticks([])
      for i in x+0.5*width*len(methods):
         ax.annotate('%d--%d' % (i, i+1), xy=(i, 0), xycoords=('data', 'axes fraction'), xytext=(0, -2), size=8, textcoords='offset points', va='top', ha='center')
      pl.xlim(x[0], x[-1]+len(lines)*width + (mnb - len(x)))
      ndx += 1
   leg = ax.legend(lines, tuple(methods), loc=0, ncol=2, prop=FP(size=8), title='$\Delta F$ '+units, fancybox=True)
   leg.get_frame().set_alpha(0.5)
   pl.savefig(os.path.join(datafile_directory, 'dF_state.pdf'), bbox_inches='tight')
   pl.close(fig)
   return

if options.breakdown:
   plotdFvsLambda()
   plotdFvsLambda2()
   if ('TI' in methods or 'TI-CUBIC' in methods):
      plotTI()

outtext = []
printLine(12*'-', 21*'-', 'plain')
printLine('%-12s' % '   States', '%9s %-11s', 'name')
printLine(12*'-', 21*'-', 'plain')
for k in range(K-1):
   printLine('%4d -- %-4d' % (k, k+1), '%10.3f  +- %6.3f', df_allk[k], ddf_allk[k])
printLine(12*'-', 21*'-', 'plain')
for i in range(len(segments)):
   printLine('%9s:  ' % segments[i], '%10.3f  +- %6.3f', dFs[i], ddFs[i])
print """\n\n
                              In the table shown are the free energy differences computed by means of various methods
                              for each pair of adjacent states with 'TOTAL' being the sum of the rows, thus yielding
                              the free energy difference between the terminal states 0 and %d.\n\n
                              A remark on the energy components interpretation:
                              'vdWaals' is computed as ('TOTAL' - 'Coulomb'), where 'Coulomb' is thought of as the free
                              energy change between the states defined by the lambda vectors (0,0,...,0) and (1,0,...,0),
                              the only varying vector component being either 'coul-lambda' or 'fep-lambda'.""" % (K-1)
etime = ttt_time.time()
tminutes = int((etime-stime)/60.)
thours = int(tminutes/60.)
tseconds = '%.2f' % (etime-stime-60*(tminutes+60*thours)) 
print "\n\nTime spent: %s hours, %s minutes, and %s seconds.\nFinished on %s" % (thours, tminutes, tseconds, ttt_time.asctime())   

# Store results.
file = open(os.path.join(datafile_directory, 'results.txt'), 'w')
file.writelines(outtext)
file.close()

file = open(os.path.join(datafile_directory, 'results.pickle'), 'w')
pickle.dump((beta, dFs, ddFs), file)
file.close()

#===================================================================================================
#                                   End of the script 
#===================================================================================================
