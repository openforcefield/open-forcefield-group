# PYTHON script to set up free energy simulation by creating mdp files from a template (genericmdp.mdp)

import optparse
import sys
import string

# User inputted values and options
parser = optparse.OptionParser(description="Create .#.mdp files for running GROMACS simulations to determine free energy of solvation.")
parser.add_option("-f", "--fileprefix", help="the full prefix of the neccessary files (.mdp file)")
parser.add_option("-s", "--states", help="the number of alchemical states (including endpoints). Default = 16", default = 16, type = int)

(options, args) = parser.parse_args()

#Check for -f user input value (necessary)
datafile_prefix = options.fileprefix
if datafile_prefix == None:
    sys.exit("No data file prefix given, exiting.")
    
num_states = options.states

# Read in datafile_prefix.mdp
filename = datafile_prefix + ".mdp"
infile = open(filename, 'r')
lines = infile.readlines()

#Create a dictionary with parameter values for modifying generic job file
parameterdict = {'NAME' : datafile_prefix}

#Loop over the lambda states
for x in range(num_states):
	
	#Create a new datafile_prefix.#.mdp
	outfilename = "%(datafile_prefix)s.%(x)s.mdp" % locals()
	outfile = open(outfilename, 'w')
	
	#Edit init-lambda-states = #
	for line in lines:
		if 'init-lambda-state' in line:
			outfile.write(line[:27] + repr(x))
		else:
			outfile.write(line)
	
	#Close outfile
	outfile.close
