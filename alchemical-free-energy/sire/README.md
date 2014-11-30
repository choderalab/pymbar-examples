To run the analysis, execute the script with python

`python alchemical_analysis.py`

within the directory with the data files.

All the flags were customized to Sire; therefore, there is no need for you -- at least at the stage of getting it to know -- to provide any. Here is a brief overview what these flags are. (The focus is on those that are relevant to Sire, or, more accurately, to the situation when the files to be analyzed contain only the dV/dLambda data.  
`python alchemical_analysis.py -h` would yield a more detailed description).

-a is the name of the software the files come from; set to 'Sire'

-p is the data file prefix; set to 'actual_grad_'

-q is the data file suffix; set to 'dat'

-d is the path to where the data files are; set to '.'

-u is the units the free energy are to be reported in; kcal/mol

-r is the number of decimal places the free energies are to be reported with

-m is the methods the free energies are to be estimated with: TI and TI-CUBIC.
If you want just TI-CUBIC, use -m ti_cubic

-g will produce graphs: the TI as a filled area under the interpolation curve (dhdl_TI.pdf) and the bar plot of ∆G's evaluated for each pair of adjacent states (dF_state.pdf). This requires matplotlib.

-s is to be used whenever you want to discard some initial snapshots

The file parser (`parser_sire.py`) is separated from the analysis proper (`alchemical_analysis.py`); make sure the former is either handy or a pythonpath is established for it. (There is `parser_gromacs.py`, as well, in case you want to run the analysis on Gromacs' files).

Also, to make it self-contained, all imports of not built-in modules needed for some non-trivial tasks have been hidden under the conditional statements. In other words, if you do not want to bother with the autocorrelation analysis (the -i flag) there is no need to checkout `timeseries.py` from the `pymbar` repository. numpy and matplotlib (for the graphs, optional) are the only prerequisites for the script.

(If you do not have python installed) Install one of its scientific distributions, like Enthought Canopy or Anaconda, and you will get it with a bunch of libraries (among which are numpy and matplotlib) rather than a stand-alone python.

Below is a list of the command the output files were obtained with.

screen_printout_1.txt:   
`python alchemical_analysis.py -d data/ -o .`
(Analysis with default settings)

screen_printout_2.txt:   
`python alchemical_analysis.py -d data/ -o . -m ti_cubic -u kJ -r 8`
(The free energies are to be estimated with TI-CUBIC and reported in kJ/mol, with 8 decimal places)

screen_printout_3.txt:   
`python alchemical_analysis.py -d data/ -o . -s 50`
(Skip first 50 ps)

The dataset contained in the data/ directory is obtained from a 16-window simulation of the gas-phase methane-to-ethane transformation run in the Michel lab at the University of Edinburgh.

Help for `alchemical_analysis.py` (obtained with `python alchemical_analysis.py -h`) is:

```Options:
  -h, --help            show this help message and exit
  -a SOFTWARE, --software=SOFTWARE
                        Package's name the data files come from. Default:
                        Sire.
                        
  -c, --cfm             The Curve-Fitting-Method-based consistency inspector.
                        Default: False.
  -d DATAFILE_DIRECTORY, --dir=DATAFILE_DIRECTORY
                        Directory in which data files are stored. Default:
                        Current directory.
  -f BFORWREV, --forwrev=BFORWREV
                        Plotting the free energy change as a function of time
                        in both directions. The number of time points (an
                        integer) is to be followed the flag. Default: 0
  -g, --breakdown       Plotting the free energy differences evaluated for
                        each pair of adjacent states for all methods. Default:
                        False.
  -i UNCORR_THRESHOLD, --threshold=UNCORR_THRESHOLD
                        Perform the analysis with rather all the data if the
                        number of uncorrelated samples is found to be less
                        than this number. If 0 is given, the time series
                        analysis will not be performed at all.
  -k BSKIPLAMBDAINDEX, --koff=BSKIPLAMBDAINDEX
                        Give a string of lambda indices separated by '-' and
                        they will be removed from the analysis. (Another
                        approach is to have only the files of interest present
                        in the directory). Default: None.
  -m METHODS, --methods=METHODS
                        A list of the methods to esitimate the free energy
                        with. Default: [TI, TI-CUBIC, DEXP, IEXP, BAR, MBAR].
                        To add/remove methods to the above list provide a
                        string formed of the method strings preceded with +/-.
                        For example, '-ti_cubic+gdel' will turn methods into
                        [TI, DEXP, IEXP, BAR, MBAR, GDEL]. 'ti_cubic+gdel', on
                        the other hand, will call [TI-CUBIC, GDEL]. 'all'
                        calls the full list of supported methods [TI, TI-
                        CUBIC, DEXP, IEXP, GINS, GDEL, BAR, UBAR, RBAR, MBAR].
  -o OUTPUT_DIRECTORY, --out=OUTPUT_DIRECTORY
                        Directory in which the output files produced by this
                        script will be stored. Default: Same as
                        datafile_directory.
  -p PREFIX, --prefix=PREFIX
                        Prefix for datafile sets, i.e.'actual_grad_'
                        (default).
  -q SUFFIX, --suffix=SUFFIX
                        Suffix for datafile sets, i.e. 'dat' (default).
  -r DECIMAL, --decimal=DECIMAL
                        The number of decimal places the free energies are to
                        be reported with. No worries, this is for the text
                        output only; the full-precision data will be stored in
                        'results.pickle'. Default: 3.
  -s EQUILTIME, --skiptime=EQUILTIME
                        Discard data prior to this specified time as
                        'equilibration' data. Units picoseconds. Default: 0
                        ps.
  -t TEMPERATURE, --temperature=TEMPERATURE
                        Temperature in K. Default: 298 K.
  -u UNITS, --units=UNITS
                        Units to report energies: 'kJ', 'kcal', and 'kBT'.
                        Default: 'kcal'
  -v, --verbose         Verbose option for BAR and MBAR. Default: False.
  -w, --overlap         Print out and plot the overlap matrix. Default: False.
  -x, --ignoreWL        Do not check whether the WL weights are equilibrated.
                        No log file needed as an accompanying input.
  -y RELATIVE_TOLERANCE, --tolerance=RELATIVE_TOLERANCE
                        Convergence criterion for the energy estimates with
                        BAR and MBAR. Default: 1e-10.
  -z INIT_WITH, --initialize=INIT_WITH
                        The initial MBAR free energy guess; either 'BAR' or
                        'zeroes'. Default: 'BAR'.
```
