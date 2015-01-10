#!/bin/bash

#SBATCH --time=01:00:00 # walltime
#SBATCH --ntasks=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1 # number of nodes
#SBATCH --mem-per-cpu=12288M # memory per CPU cores
 
# NECESSARY Input Files:
#	freqinHP - Standard output from a gaussian 09 frequency calculation *CAN BE EDITED*
#	progdyn.conf - File of configuration options for the progdyn script suite (for explanation on this, see the progdyn.conf file itself) *CAN BE EDITED*

# OPTIONAL Input Files:
#	isomernumber - Provides a start number(1+) for number runs. If no file is provided, progdyn defaults to 1
#	detour - A signal file that, by existing, signals the program to do side calculations
#	nogo - A signal file that, by existing, signals the program to stop between points
#	bypassproggen - A signal file that, by existing, signals the program to use a supplied geoPlusVel file instead of generating one for itself
#	methodfile - A file that contains lines to be added at the end of each g09.com input file, such as lines that call for an NMR calculation
#	ZMAT - An input file for the CFOUR suite of programs. When ZMAT is supplies, progdynstarterHP will run call CFOUR by making use of the progcfour script
#	cannontraj - A file containing a vector for each atom, used to fire an initial geometry in a particular direction

# NECESSARY Program Files:
#	progdynstarterHP (Bash script) - Controls the overall execution of the progdyn suite of scripts
#	proggenHP (Awk script) - Starts a trajectory, giving each mode its zero point energy (if a quasiclassical calculation)
#				 plus random additional excitations depending on the temperature
#	prog1stpoint (Awk script) - Creates the first Gaussian input file for each run
#	prog2ndpoint (Awk script) - Creates the second Gaussian input file for each run. prog2ndpoint also checks the energy of the first point to
#		       see if it fits with the desired energy and aborts the run if it does not by creating appropriate output in file Echeck
#	progdynb (Awk script) - Creates subsequent Gaussian input files until run is completed
#	proganal (Awk script) - Analyzes the latest point and see if a run is done. This program must be redone for each new system. Elaborate changes are
#				often programmed into proganal, such as the automatic changing of configuration variables, proganal creates the output to
#				dynfollowfile and NMRlist or NMRlistdis
#	randgen (C Program) - Generates random numbers between 0 and 1. These are generated all at once and stored in a file for use by proggenHP

# OPTIONAL Program Files:
#	progcfour (????) - *WARNING NOT CURRENTLY INCLUDED IN THIS SUITE* Runs CFOUR calculations (not needed for most kinds of runs)

# Output Files (Not all are necessarily output, depends on the specific run configuration)
#	isomernumber - A running tab of the trajectory number
#	runpointnumber - A running tab of the point in the trajectory
#	Echeck - Output from prog2ndpoint checks the energy of the trajectory to see if it fits with the desired energy
#	geoRecord - A record of all the geoPlusVel files
#	geoPlusVel - Created by proggen, this gives the starting positions, velocities, isotopic masses, excitations of the normal modes, and initial
#		     displacements of the normal modes for current run
#	g09.com - Created by prog1stpoint, prog2ndpoint, and progdynb, this is the latest input file for Gaussian09 for current run and latest point
#	olddynrun and olderdynrun - Contains the last two outputs from Gaussian, for creation of the next point
#	traj, traj1, traj2, traj3, etc. - Contains the geometries and energies for each trajectory, numbered by the isomernumber, in a format suitable for reading by Molden
#	dyn - A record of the Gaussian outputs
#	dynfollowfile - A short record of the runs and their results. The data desired for dynfollowfile must be programmed into the scrip proganal as needed for each system studied
#	NMRlist or NMRlistdis - output of NMR predictions at each point in a trajectory
#	skipstart - A signal file that, by existing, tells progdynstarterHP that we are in the middle of a run. For trajectories that are propagated
#		    forward and backward in time, skipstart keeps track of whether one is in the forward or reverse part
#	diagnostics - optional output that follows which subprograms are running and configuration variables, decided by variable in progdyn.conf
#	vellist - Optional output that lists the velocities of each atom, decided by variable in progdyn.conf, or lists the total kinetic energy of the system and the classical temperature

#   the NOTE in the "cleanup_scratch" signal handler

export FILES_NEEDED_DIR=REPLACEME
export TEMPORARY_DIR="/tmp/$SLURM_JOB_ID"
export GAUSS_SCRDIR="$TEMPORARY_DIR/temporary_files"
export PROG_SCRDIR="$TEMPORARY_DIR/prog_files"
export JOB_NAME=freqinHP
export PROG_HOME=$TEMPORARY_DIR
export G09_ROOT="/blueapps/chem"

#export PATH=$PATH:/blueapps/chem/gaussian/g09 #Makes it so g09.profile executes properly by giving it the location of executable file gau-machine
 
#set up function.  this isn't called/run here.  It's just used 
#   if the job is canceled via a signal
cleanup_scratch()
{
        echo "cleanup_scratch function called"

	if (test -d "$TEMPORARY_DIR") then
		cp -a /$TEMPORARY_DIR/. /$HOME/ProgdynSuite/RunOutputs/out$SLURM_JOB_ID/
		rm -rf "$TEMPORARY_DIR"
	fi

        exit 0
}
 
#Associate the function "cleanup_scratch" with the TERM signal, which is usually how jobs get killed
#trap "kill $pid; cleanup_scratch; exit 1" TERM SIGTERM KILL SIGKILL EXIT
 
#basic diagnostic output
echo "---"
echo "Beginning-of-job Diagnostic information:"
echo "---"
echo "Temporary Directory:"
echo "$TEMPORARY_DIR"
echo "---"
echo "Scratch Directory:"
echo "$GAUSS_SCRDIR"
echo "---"
echo "Job Source Directory:"
echo "$SLURM_SUBMIT_DIR"
echo "---"
echo "Current Time:"
date
echo "---"
 
#create scratch directory
echo "Creating scratch directory at $GAUSS_SCRDIR"
echo "Creating scratch directory at $PROG_SCRDIR"
mkdir -pv "$GAUSS_SCRDIR"
mkdir -pv "$PROG_SCRDIR"
echo "---"

mv $FILES_NEEDED_DIR/* $TEMPORARY_DIR

if (test -f $PROG_HOME/runpointnumber) then
	
	if [ "`cat $PROG_HOME/runpointnumber`" == "1" ]; then
		:
	elif [ "`cat $PROG_HOME/runpointnumber`" == "2" ]; then
		:
	elif [ "`cat $PROG_HOME/runpointnumber`" == "3" ]; then
		:
	else
		rm $PROG_HOME/runpointnumber
		echo 1 > $PROG_HOME/runpointnumber
	fi
else
	echo 1 > $PROG_HOME/runpointnumber
fi
#cp "$PROG_HOME/GENBAS" $PROG_SCRDIR
#cp "$PROG_HOME/progcfour" $PROG_SCRDIR
#echo "---"

echo "Changing directory to temporary dir at $TEMPORARY_DIR"
cd "$TEMPORARY_DIR"
echo "---"
 
mkdir "$HOME/ProgdynSuite/RunOutputs/out$SLURM_JOB_ID" # Where all the output will be stored after the job is finished
mkdir "$HOME/ProgdynSuite/RunOutputs/out$SLURM_JOB_ID/comlogfiles" # Where each .com and .log file is stored
export COM_LOG_FILES="$HOME/ProgdynSuite/RunOutputs/out$SLURM_JOB_ID/comlogfiles" # progdynstarterHP uses this to know where to store the .com and .log files
echo "Starting PROG Run at:"
date

. $PROG_HOME/progdynstarterHP & # Starts Progdyn

#Associate the function "cleanup_scratch" with the TERM signal, which is usually how jobs get killed
pid=$!
trap "kill $pid; cleanup_scratch; exit 1" TERM SIGTERM KILL SIGKILL EXIT
wait $pid

echo "---"
echo "PROG Run ended at:"
date
echo "---"

cleanup_scratch

exit 0
