# PYTHON script to write and submit pbsjobs to queue

import optparse
import sys
import string
import os
import stat

# Set as qsub when finished debugging, set to cat when debugging
qsub_command = 'qsub'

# User inputted values and options
parser = optparse.OptionParser(description="submit GROMACS mdrun commands to cluster queue.")
parser.add_option("-f", "--fileprefix", help="the full prefix of the files of concern")
parser.add_option("-s", "--states", help="the number of alchemical states (including endpoints). Default = 16", default = 16, type = int)
parser.add_option("-j", "--jobfile", help="the full name of the generic PBS job shell script. Default = /bigtmp/bsz9ur/tools/pbsfile.sh", default = '/bigtmp/bsz9ur/tools/pbsfile.sh')
parser.add_option("-G", "--gromacs", help="location of GROMACS binaries. Default = /h3/n1/shirtsgroup/gromacs_46/Install_v462/bin", default = '/h3/n1/shirtsgroup/gromacs_46/Install_v462/bin')

(options, args) = parser.parse_args()

#Check for -f user input value (necessary)
datafile_prefix = options.fileprefix
if datafile_prefix == None:
    sys.exit("No data file prefix given, exiting.")

num_states = options.states 
gromacs_loc = options.gromacs
job_script = options.jobfile

#Create dictionary with parameter values for modifying generic job file
parameterdict = {'LOCATION' : gromacs_loc, 'NAME' : datafile_prefix}

#Loop over lambda states
for x in range(num_states):
	#Create shell script to submit PBS job (single simulation at lambda = x)
	#Add parameter dictionary member for current lambda state
	parameterdict['NUM'] = x
	
	#Open generic job file for reading
	genericjobfile = open(job_script, 'r')
	
	#Create new job file and write to it
	jobfilename = "pbsjob%s.sh" % (x)
	jobfile = open(jobfilename, 'w')
	
	#Iterate through the lines of the generic file
	for line in genericjobfile:
		#Substitute keys for values in the template line and output result
		jobfile.write(string.Template(line).substitute(parameterdict))
		
	genericjobfile.close()
	jobfile.close()
	
	#Call shell script to submit PBS job (single simulation at lambda = x)
	
	#comment out os.system(qsub if running jobs on local machine
	os.system('%s %s' % (qsub_command, jobfilename))
	
	#comment out os.system(jobfilename) if running on cluster
	#st = os.state(jobfilename)
	#os.chmod(jobfilename, st.st_mode \ stat.S_IXUSR)
	#jobfilename = './' + jobfilename
	#os.system(jobfilename)
	
	#Sleep a little while to make sure that PBS has time to copy the script, and also so as not to overload qsub.
	
	os.system('sleep 1')
