#!/bin/bash

# walltime
#SBATCH --time=00:30:00
# number of processor cores (i.e. tasks)
#SBATCH --ntasks=1
# number of nodes
#SBATCH --nodes=1
# memory per CPU core
#SBATCH --mem-per-cpu=12288M

export JOB_NAME=REPLACEME
export FILES_NEEDED_DIR=REPLACEME
export TEMPORARY_DIR="/tmp/$SLURM_JOB_ID"
export GAUSS_SCRDIR="$TEMPORARY_DIR/gauss_scrdir"

#set up function.  this isn't called/run here.  It's just used if the job is canceled via a signal
cleanup_scratch()
{
	if (test -d "$TEMPORARY_DIR") then
		exit 0
	else
		echo "Jobbed was terminated mid-run."
		cd "$SLURM_SUBMIT_DIR"
		rm -rf "$TEMPORARY_DIR"
	fi
	exit 0
}

mkdir $TEMPORARY_DIR
cp -ar  $FILES_NEEDED_DIR/* $TEMPORARY_DIR

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
mkdir -pv "$GAUSS_SCRDIR"
echo "---"

#changing directory to $TEMPORARY_DIR
echo "Changing directory to temporary dir at $TEMPORARY_DIR"
cd "$TEMPORARY_DIR"
ls
echo "---"


echo "Starting Gaussian Run at:"
date

#the actual gaussian run starts here
/fslapps/chem/bin/rung09 < "$JOB_NAME.com" > "$JOB_NAME.log" &
pid=$!
trap "kill $pid; cleanup_scratch; exit 1" TERM SIGTERM KILL SIGKILL EXIT
wait $pid
/fslapps/chem/bin/rmipc

export PREVIOUS_JOBID=$SLURM_JOB_ID

grep "Normal termination" $TEMPORARY_DIR/"$JOB_NAME".log > $TEMPORARY_DIR/successfulrun 
if (test -s $TEMPORARY_DIR/successfulrun) then
	echo "SUCCESSFUL RUN"
	#rm $TEMPORARY_DIR/successfulrun
	#mkdir $HOME/ProgdynSuite/RunOutputs/out"$SLURM_JOB_ID"
	#cp -av $TEMPORARY_DIR/* $HOME/ProgdynSuite/RunOutputs/out"$SLURM_JOB_ID"
	cp $TEMPORARY_DIR/"$JOB_NAME".log $HOME/ProgdynSuite/ProgdynScripts/freqinHP
	echo Normal Termination. Submitting to ProgdynSuite
	$HOME/ProgdynSuite/SubmitScripts/ProgSubmitWrapper.sh
else
	echo 
fi

echo "---"
echo "Job ending time:"
date
echo "---"

rm -rf "$GAUSS_SCRDIR"
mkdir $HOME/ProgdynSuite/RunOutputs/out"$PREVIOUS_JOBID"
mv "$TEMPORARY_DIR"/* $HOME/ProgdynSuite/RunOutputs/out"$PREVIOUS_JOBID"
mv $FILES_NEEDED_DIR/slurm* $HOME/ProgdynSuite/RunOutputs/out"$PREVIOUS_JOBID"
rm -rf "$FILES_NEEDED_DIR" 

exit 0
