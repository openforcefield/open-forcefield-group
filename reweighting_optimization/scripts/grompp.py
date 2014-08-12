# PYTHON script to write and submit pbsjobs to queue

import optparse
import sys
import string
import os
import stat

# User inputted values and options
parser = optparse.OptionParser(description="submit GROMACS grompp and mdrun commands to cluster queue.")
parser.add_option("-f", "--fileprefix", help="the full prefix of the files of concern")
parser.add_option("-s", "--states", help="the number of alchemical states (including endpoints). Default = 16", default = 16, type = int)
parser.add_option("-G", "--gromacs", help="location of GROMACS binaries, default = /h3/n1/shirtsgroup/gromacs_46/Install_v462/bin", default = '/h3/n1/shirtsgroup/gromacs_46/Install_v462/bin')    

(options, args) = parser.parse_args()

#Check for -f user input value (necessary)
datafile_prefix = options.fileprefix
if datafile_prefix == None:
    sys.exit("No data file prefix given, exiting.")
    
num_states = options.states
gromacs_loc = options.gromacs

#Create dictionary with parameter values for modifying generic job file
parameterdict = {'LOCATION' : gromacs_loc, 'NAME' : datafile_prefix}

#Loop over lambda states
for x in range(num_states):
    os.system('%(gromacs_loc)s/grompp_d -f %(datafile_prefix)s.%(x)s.mdp -c %(datafile_prefix)s.gro -p %(datafile_prefix)s.top -o %(datafile_prefix)s.%(x)s.tpr -maxwarn 2' % locals())
print 'finished'
