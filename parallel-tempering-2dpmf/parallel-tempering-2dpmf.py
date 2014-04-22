#!/usr/bin/python

# Estimate 2D potential of mean force for alanine dipeptide parallel tempering data using MBAR.
#
# PROTOCOL
#
# * Potential energies and (phi, psi) torsions from parallel tempering simulation are read in by temperature
# * WHAM [2] is used to rapidly generate an initial guess for the dimensionless free energies f_k.
# * Replica trajectories of potential energies and torsions are reconstructed to reflect their true temporal
# correlation, and then subsampled to produce statistically independent samples, collecting them again by temperature
# * The MBAR class is initialized with this initial guess at dimensionless free energies f_k, reducing time for
# solution of self-consistent equations
# * The torsions are binned into sequentially labeled bins in two dimensions
# * The relative free energies and uncertainties of these torsion bins at the temperature of interest is estimated
# * The 2D PMF is written out
# 
#
# REFERENCES
#
# [1] Shirts MR and Chodera JD. Statistically optimal analysis of samples from multiple equilibrium states.
# J. Chem. Phys. 129:124105, 2008
# http://dx.doi.org/10.1063/1.2978177
#
# [2] Kumar S, Bouzida D, Swensen RH, Kollman PA, and Rosenberg JM. The weighted histogram analysis method
# for free-energy calculations on biomolecules. I. The Method. J. Comput Chem. 13:1011, 1992.

#===================================================================================================
# IMPORTS
#===================================================================================================

import numpy
from math import *
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import commands
import os
import os.path

#===================================================================================================
# CONSTANTS
#===================================================================================================

kB = 1.3806503 * 6.0221415 / 4184.0 # Boltzmann constant in kcal/mol/K

#===================================================================================================
# PARAMETERS
#===================================================================================================

data_directory = 'data/' # directory containing the parallel tempering data
temperature_list_filename = os.path.join(data_directory, 'temperatures') # file containing temperatures in K
potential_energies_filename = os.path.join(data_directory, 'energies', 'potential-energies') # file containing total energies (in kcal/mol) for each temperature and snapshot
trajectory_segment_length = 20 # number of snapshots in each contiguous trajectory segment
niterations = 500 # number of iterations to use
target_temperature = 302 # target temperature for 2D PMF (in K)
nbins_per_torsion = 10 # number of bins per torsion dimension

#===================================================================================================
# SUBROUTINES
#===================================================================================================

def read_file(filename):
   """Read contents of the specified file.
      
   ARGUMENTS
     filename (string) - the name of the file to be read
     
   RETURNS
     lines (list of strings) - the contents of the file, split by line

   """

   infile = open(filename, 'r')
   lines = infile.readlines()
   infile.close()
   
   return lines

def logSum(log_terms):
   """Compute the log of a sum of terms whose logarithms are provided.

   REQUIRED ARGUMENTS  
      log_terms is the array (possibly multidimensional) containing the logs of the terms to be summed.

   RETURN VALUES
      log_sum is the log of the sum of the terms.

   """

   # compute the maximum argument
   max_log_term = log_terms.max()

   # compute the reduced terms
   terms = numpy.exp(log_terms - max_log_term)

   # compute the log sum
   log_sum = log( terms.sum() ) + max_log_term

   # return the log sum
   return log_sum

def histogram_wham(beta_k, U_kn, N_k, nbins = 100, bin_width = None, maximum_iterations = 5000, relative_tolerance = 1.0e-8, initial_f_k = None):
   """Construct an initial guess of the f_k by histogram reweighting (specifically, WHAM [2]).

   ARGUMENTS
     beta_k (numpy K array) - inverse temperatures (in units of 1/energy)
     U_kn (numpy K x N_max array) - potential energies (in units of energy)
     N_k (numpy K array of numpy.int32) - number of samples per states, N_max = N_k.max()

   OPTIONAL ARGUMENTS
     nbins (int) - if specified, sets the number of bins to use (default: 100)
     bin_width (float) - if specified, sets the bin width (overrides nbins) (defulat: None)
     maximum_iterations (int) - maximum number of iterations to use
     relative_tolerance (floeat) - relative convergence tolerance (default: 1.0e-8)

   RETURNS
     f_k (numpy K array) - guess at initial state dimensionless free energies

   REFERENCE
     [2] Kumar S, Bouzida D, Swensen RH, Kollman PA, and Rosenberg JM. The weighted histogram analysis method
     for free-energy calculations on biomolecules. I. The Method. J. Comput Chem. 13:1011, 1992.
   
   """

   # Get sizes
   K = N_k.size
   N_max = N_k.max()

   # Create a list of indices of all configurations in kn-indexing.
   mask_kn = numpy.zeros([K,N_max], dtype=numpy.bool)
   for k in range(0,K):
      mask_kn[k,0:N_k[k]] = True
   # Create a list from this mask.
   sample_indices = numpy.where(mask_kn)

   # Construct histogram bins
   M = nbins # number of energy bins
   SMALL = 1.0e-6
   U_min = U_kn[sample_indices].min()
   U_max = U_kn[sample_indices].max()
   U_max += (U_max - U_min) * SMALL # increment by a bit
   delta_U = (U_max - U_min) / float(M)
   if (bin_width != None):
      delta_U = bin_width
      M = int(numpy.ceil((U_max - U_min) / bin_width))
      print "Using %d bins to achieve bin width of %f" % (M, delta_U)
   else:
      print "Bin width is %f energy units to achieve nbins = %d" % (delta_U, M)
   U_m = U_min + delta_U * (0.5 + numpy.arange(0,M,dtype = numpy.float64))
   H_m = numpy.zeros([M], numpy.float64)
   # assign snapshots to energy bins
   bin_kn = numpy.zeros([K,N_max], numpy.int32) # bin_kn[k,n] is bin index of U_kn[k,n]
   bin_kn[sample_indices] = numpy.array(((U_kn[sample_indices] - U_min) / delta_U), numpy.int32)
   H_m = numpy.bincount(bin_kn[sample_indices])
   
   # compute logs of various quantities
   LOG_ZERO = -1000.0 # replacement for log(0)
   log_H_m = numpy.ones([M], numpy.float64) * LOG_ZERO
   for m in range(M):
      if (H_m[m] > 0):
         log_H_m[m] = log(H_m[m])
   log_N_k = numpy.ones([K], numpy.float64) * LOG_ZERO
   for k in range(K):
      if (N_k[k] > 0):
         log_N_k[k] = log(N_k[k])

   # initialize free energies
   f_k = numpy.zeros([K], numpy.float64)
   if (initial_f_k != None):
      f_k = initial_f_k.copy()

   # iterate
   f_k_new = numpy.zeros([K], numpy.float64)
   max_delta = 0.0
   for iteration in range(maximum_iterations):
      # print "iteration %d" % iteration
      
      # Form auxiliary matrices, used in both self-consistent iteration and Newtom-Raphson.
      # log space
      log_w_mi = numpy.zeros([M,K], numpy.float64)
      for m in range(M):
         # denominator = \sum_k N_k exp[f_k - \beta_k U_m]
         exp_arg_k = log_N_k[:] + f_k[:] - beta_k[:]*U_m[m]
         log_w_mi[m,:] = exp_arg_k[:] - logSum(exp_arg_k[:])
      # real space
      w_mi = numpy.zeros([M,K], numpy.float64)
      for m in range(M):
         exp_arg_k = f_k[:] - beta_k[:]*U_m[m]
         max_arg = exp_arg_k.max()
         numerators_i = N_k[:] * numpy.exp(exp_arg_k[:] - max_arg)
         w_mi[m,:] = numerators_i[:] / numerators_i.sum()
         
      # Compute new estimates of log weights using self-consistent iteration.
      for i in range(K):
         # Compute new estimate of log weights.
         f_k_new[i] = f_k[i] + log(N_k[i]) - logSum(log_H_m[:] + log_w_mi[:,i])
      
      # shift to ensure f_k_new[0] = 0.0
      f_k_new -= f_k_new[0]
      # store difference
      Delta_f_k = f_k_new - f_k
      # update f_k
      f_k = f_k_new.copy()

      # If all f_k are zero, terminate
      if numpy.all(f_k == 0.0):
        break
      
      # Terminate when max((f - fold) / f) < relative_tolerance for all nonzero f.
      max_delta = (numpy.abs(Delta_f_k) / (numpy.abs(f_k_new)).max()).max()
      print "iteration %8d relative max_delta = %8e" % (iteration, max_delta)
      if numpy.isnan(max_delta) or (max_delta < relative_tolerance):
         break

   # if unconverged
   if (iteration == maximum_iterations):
      print "Did not converge in %d iterations (final relative tolerance %e)" % (maximum_iterations, max_delta)
      print "f_k = "
      print f_k
      return f_k
   
   # summary
   print "Converged to relative tolerance of %e (convergence tolerance %e) in %d iterations" % (max_delta, relative_tolerance, iteration)
   print "f_k = "
   print f_k

   # return the estimate of the dimensionless free energies
   return f_k

def write_free_energies(filename, f_k):
   """Write free energies to the given file.

   ARGUMENTS
     filename (string) - the name of the file to write to
     f_k (numpy array) - weights to write
   """

   # get size
   K = f_k.size

   # Write to file.
   outfile = open(filename, 'w')
   for k in range(K):
      outfile.write("%f " % f_k[k])
   outfile.write("\n")
   outfile.close()

   return

def read_free_energies(filename):
   """Read free energies from the given file.

   ARGUMENTS
     filename (string) - the name of the file to write to
   """

   print "Reading free energies from %s..." % filename

   # read file
   infile = open(filename, 'r')
   lines = infile.readlines()
   infile.close()

   # split into elements
   elements = list()
   for line in lines:
      elements += line.split()
   print elements

   # parse
   K = len(elements)
   f_k = numpy.zeros([K], numpy.float64)
   for k in range(K):
      f_k[k] = float(elements[k])

   return f_k

#===================================================================================================
# MAIN
#===================================================================================================

#===================================================================================================
# Read temperatures
#===================================================================================================

# Read list of temperatures.
lines = read_file(temperature_list_filename)
# Construct list of temperatures
temperatures = lines[0].split()
# Create numpy array of temperatures
K = len(temperatures)
temperature_k = numpy.zeros([K], numpy.float32) # temperature_k[k] is temperature of temperature index k in K
for k in range(K):
   temperature_k[k] = float(temperatures[k])
# Compute inverse temperatures
beta_k = (kB * temperature_k)**(-1) 

# Define other constants
T = trajectory_segment_length * niterations # total number of snapshots per temperature

#===================================================================================================
# Read potential eneriges
#===================================================================================================

print "Reading potential energies..."
U_kt = numpy.zeros([K,T], numpy.float32) # U_kn[k,t] is the potential energy (in kcal/mol) for snapshot t of temperature index k
lines = read_file(potential_energies_filename)
print "%d lines read, processing %d snapshots" % (len(lines), T)
for t in range(T):
   # Get line containing the energies for snapshot t of trajectory segment n
   line = lines[t]
   # Extract energy values from text
   elements = line.split()
   for k in range(K):
      U_kt[k,t] = float(elements[k])

#===================================================================================================
# Read phi, psi trajectories
#===================================================================================================

print "Reading phi, psi trajectories..."
phi_kt = numpy.zeros([K,T], numpy.float32) # phi_kt[k,n,t] is phi angle (in degrees) for snapshot t of temperature k
psi_kt = numpy.zeros([K,T], numpy.float32) # psi_kt[k,n,t] is psi angle (in degrees) for snapshot t of temperature k
for k in range(K):
   phi_filename = os.path.join(data_directory, 'backbone-torsions', '%d.phi' % (k))
   psi_filename = os.path.join(data_directory, 'backbone-torsions', '%d.psi' % (k))   
   phi_lines = read_file(phi_filename)
   psi_lines = read_file(psi_filename)
   print "k = %d, %d phi lines read, %d psi lines read" % (k, len(phi_lines), len(psi_lines))   
   for t in range(T):
      # Extract phi and psi
      phi_kt[k,t] = float(phi_lines[t])
      psi_kt[k,t] = float(psi_lines[t])

#===================================================================================================
# Read replica indices
#===================================================================================================

print "Reading replica indices..."
filename = os.path.join(data_directory, 'replica-indices')
lines = read_file(filename)
replica_ik = numpy.zeros([niterations,K], numpy.int32) # replica_ki[i,k] is the replica index of temperature k for iteration i
for i in range(niterations):
   elements = lines[i].split()
   for k in range(K):
      replica_ik[i,k] = int(elements[k])
print "Replica indices for %d iterations processed." % niterations

#===================================================================================================
# Permute data by replica and subsample to generate an uncorrelated subset of data by temperature
#===================================================================================================

assume_uncorrelated = False
if (assume_uncorrelated):
   # DEBUG - use all data, assuming it is uncorrelated
   print "Using all data, assuming it is uncorrelated..."
   U_kn = U_kt.copy()
   phi_kn = phi_kt.copy()
   psi_kn = psi_kt.copy()
   N_k = numpy.zeros([K], numpy.int32)
   N_k[:] = T
   N_max = T
else:
   # Permute data by replica
   print "Permuting data by replica..."
   U_kt_replica = U_kt.copy()
   phi_kt_replica = psi_kt.copy()
   psi_kt_replica = psi_kt.copy()
   for iteration in range(niterations):
      # Determine which snapshot indices are associated with this iteration
      snapshot_indices = iteration*trajectory_segment_length + numpy.arange(0,trajectory_segment_length)
      for k in range(K):
         # Determine which replica generated the data from temperature k at this iteration
         replica_index = replica_ik[iteration,k]
         # Reconstruct portion of replica trajectory.
         U_kt_replica[replica_index,snapshot_indices] = U_kt[k,snapshot_indices]
         phi_kt_replica[replica_index,snapshot_indices] = phi_kt[k,snapshot_indices]
         psi_kt_replica[replica_index,snapshot_indices] = psi_kt[k,snapshot_indices]
   # Estimate the statistical inefficiency of the simulation by analyzing the timeseries of interest.
   # We use the max of cos and sin of the phi and psi timeseries because they are periodic angles.
   # The 
   print "Computing statistical inefficiencies..."
   g_cosphi = timeseries.statisticalInefficiencyMultiple(numpy.cos(phi_kt_replica * numpy.pi / 180.0))
   print "g_cos(phi) = %.1f" % g_cosphi
   g_sinphi = timeseries.statisticalInefficiencyMultiple(numpy.sin(phi_kt_replica * numpy.pi / 180.0))
   print "g_sin(phi) = %.1f" % g_sinphi   
   g_cospsi = timeseries.statisticalInefficiencyMultiple(numpy.cos(psi_kt_replica * numpy.pi / 180.0))
   print "g_cos(psi) = %.1f" % g_cospsi
   g_sinpsi = timeseries.statisticalInefficiencyMultiple(numpy.sin(psi_kt_replica * numpy.pi / 180.0))
   print "g_sin(psi) = %.1f" % g_sinpsi
   # Subsample data with maximum of all correlation times.
   print "Subsampling data..."
   g = numpy.max(numpy.array([g_cosphi, g_sinphi, g_cospsi, g_sinpsi]))
   indices = timeseries.subsampleCorrelatedData(U_kt[k,:], g = g)   
   print "Using g = %.1f to obtain %d uncorrelated samples per temperature" % (g, len(indices))
   N_max = int(numpy.ceil(T / g)) # max number of samples per temperature   
   U_kn = numpy.zeros([K, N_max], numpy.float64)
   phi_kn = numpy.zeros([K, N_max], numpy.float64)
   psi_kn = numpy.zeros([K, N_max], numpy.float64)
   N_k = N_max * numpy.ones([K], numpy.int32)
   for k in range(K):
      U_kn[k,:] = U_kt[k,indices]
      phi_kn[k,:] = phi_kt[k,indices]
      psi_kn[k,:] = psi_kt[k,indices]
   print "%d uncorrelated samples per temperature" % N_max
         
#===================================================================================================
# Generate a list of indices of all configurations in kn-indexing
#===================================================================================================

# Create a list of indices of all configurations in kn-indexing.
mask_kn = numpy.zeros([K,N_max], dtype=numpy.bool)
for k in range(0,K):
   mask_kn[k,0:N_k[k]] = True
# Create a list from this mask.
indices = numpy.where(mask_kn)

#===================================================================================================
# Compute reduced potential energy of all snapshots at all temperatures
#===================================================================================================

print "Computing reduced potential energies..."
u_kln = numpy.zeros([K,K,N_max], numpy.float32) # u_kln[k,l,n] is reduced potential energy of trajectory segment n of temperature k evaluated at temperature l
for k in range(K):
   for l in range(K):
      u_kln[k,l,0:N_k[k]] = beta_k[l] * U_kn[k,0:N_k[k]]

#===================================================================================================
# Bin torsions into histogram bins for PMF calculation
#===================================================================================================

# Here, we bin the (phi,psi) samples into bins in a 2D histogram.
# We assign indices 0...(nbins-1) to the bins, even though the histograms are in two dimensions.
# All bins must have at least one sample in them.
# This strategy scales to an arbitrary number of dimensions.

print "Binning torsions..."
# Determine torsion bin size (in degrees)
torsion_min = -180.0
torsion_max = +180.0
dx = (torsion_max - torsion_min) / float(nbins_per_torsion)
# Assign torsion bins
bin_kn = numpy.zeros([K,N_max], numpy.int16) # bin_kn[k,n] is the index of which histogram bin sample n from temperature index k belongs to
nbins = 0
bin_counts = list()
bin_centers = list() # bin_centers[i] is a (phi,psi) tuple that gives the center of bin i
for i in range(nbins_per_torsion):
   for j in range(nbins_per_torsion):
      # Determine (phi,psi) of bin center.      
      phi = torsion_min + dx * (i + 0.5)
      psi = torsion_min + dx * (j + 0.5)

      # Determine which configurations lie in this bin.
      in_bin = (phi-dx/2 <= phi_kn[indices]) & (phi_kn[indices] < phi+dx/2) & (psi-dx/2 <= psi_kn[indices]) & (psi_kn[indices] < psi+dx/2)

      # Count number of configurations in this bin.
      bin_count = in_bin.sum()

      # Generate list of indices in bin.
      indices_in_bin = (indices[0][in_bin], indices[1][in_bin])

      if (bin_count > 0):
         # store bin (phi,psi)
         bin_centers.append( (phi, psi) )
         bin_counts.append( bin_count )

         # assign these conformations to the bin index         
         bin_kn[indices_in_bin] = nbins
         
         # increment number of bins
         nbins += 1

print "%d bins were populated:" % nbins
for i in range(nbins):
   print "bin %5d (%6.1f, %6.1f) %12d conformations" % (i, bin_centers[i][0], bin_centers[i][1], bin_counts[i])

#===================================================================================================
# Read free energies, if they exist.
# Otherwise, use WHAM to quickly compute an initial guess of dimensionless free energies f_k.
# Then initialize MBAR.
#===================================================================================================

free_energies_filename = 'f_k.out'
if os.path.exists(free_energies_filename):
   # Read free energies
   print "Reading free energies from %s ..." % free_energies_filename
   f_k = read_free_energies(free_energies_filename)

   # Initialize MBAR without any iterations.
   mbar = pymbar.MBAR(u_kln, N_k, maximum_iterations = 0, verbose = True, initial_f_k = f_k)
else:
   # Compute initial guess using histogram WHAM.
   print "Using WHAM to generate histogram-based initial guess of dimensionless free energies f_k..."
   f_k = histogram_wham(beta_k, U_kn, N_k)

   # Initialize MBAR with Newton-Raphson
   print "Initializing MBAR (will estimate free energy differences first time)..."
   mbar = pymbar.MBAR(u_kln, N_k, verbose = True, initial_f_k = f_k)

   # Write free energies
   write_free_energies(free_energies_filename, mbar.f_k)

#===================================================================================================
# Compute PMF at the desired temperature.
#===================================================================================================

print "Computing potential of mean force..."

# Compute reduced potential energies at the temperaure of interest
target_beta = 1.0 / (kB * target_temperature)
u_kn = target_beta * U_kn
# Compute PMF at this temperature, returning dimensionless free energies and uncertainties.
# f_i[i] is the dimensionless free energy of bin i (in kT) at the temperature of interest
# d2f_ij[i,j] is an estimate of the covariance in the estimate of (f_i[i] - f_j[j])
(f_i, d2f_ij) = mbar.computePMF_states(u_kn, bin_kn, nbins)

# Find index of bin with lowest free energy.
imin = f_i.argmin()

# Show free energy and uncertainty of each occupied bin relative to lowest free energy
print "2D PMF"
print ""
print "%8s %6s %6s %8s %10s %10s" % ('bin', 'phi', 'psi', 'N', 'f', 'df')
for i in range(nbins):
   print '%8d %6.1f %6.1f %8d %10.3f %10.3f' % (i, bin_centers[i][0], bin_centers[i][1], bin_counts[i], f_i[i], sqrt(d2f_ij[i,imin]))

