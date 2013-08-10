Parallel tempering data for alanine dipeptide in explicit solvent.
Only (phi,psi) trajectories and energies have been collected here.

Data for the alanine dipeptide data can be found here:

http://www.dillgroup.ucsf.edu/~jchodera/exports/jan-hendrik.prinz/alanine-dipeptide

there is a torsion-trajectories/ directory with one file for each
temperature.  All torsion trajectories are concatenated, with 20 ps
per trajectory and samples written every 0.1 ps (for 200 samples /
trajectory).  The torsion angles are written in scientific notation,
such as

L0 3.587457e+01 1.769885e+02 -1.573196e+02 1.750121e+02 -1.714699e+02

Where you should ignore the first 'L0' column and use the remaining
numbers, which are

psi1
omega1
phi2
psi2
omega2

You want 'phi2' and 'psi2'.

The temperatures are listed in the file 'temps', with the first
temperature starting from 0.  That means that 5.tors is for 302 K and
19.tors is for 400 K.

The first line of the files are blank.

The coordinate trajectories (which I've only copied for temperature
indices 5 and 19) are in coordinate-trajectories/.  I've provided you
with gzipped AMBER trajectories.

I've also copied a PDB file of the solute and the whole system so you
know which atoms are which, and the number of atoms in the coordinate
files.

etot.out - total energies for each trajectory

The etot.out file has 100200 lines and 40 columns.  Each column is the total energy trajectory for one of the temperatures.  You will probably want to average the energy for all 200 snapshots in each trajectory segment of 20 ps to produce the energy you will use in reweighting.

It appears that I *did* run 501 iterations of 20 ps each, because I started counting at 0 and stopped counting after 500.  I don't know why I didn't remember this, but I think my analysis may have excluded the last iteration for the paper.