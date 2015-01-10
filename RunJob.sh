#!/bin/bash

com_files=(`find *.com -maxdepth 1 2>/dev/null`) # creates an array of all the com files in the current directory
if [ ${#com_files[@]} -gt 0 ]; then # test if there are any com files
	com_directory=$PWD # there are so store the location of the com files
else
	echo "No .com files were found in the current directory" # there aren't so inform the user
	exit 0 # and then quit
fi
echo Current sbatch values:
walltime=`grep "#SBATCH --time=" $HOME/ProgdynSuite/SubmitScripts/g09.sh`
echo "Walltime:                  ${walltime/"#SBATCH --time="}"
ntasks=`grep "#SBATCH --ntasks=" $HOME/ProgdynSuite/SubmitScripts/g09.sh`
echo "Number of processor cores: ${ntasks/"#SBATCH --ntasks="}"
nodes=`grep "#SBATCH --nodes=" $HOME/ProgdynSuite/SubmitScripts/g09.sh`
echo "Number of nodes:           ${nodes/"#SBATCH --nodes="}"
mem_per_cpu=`grep "#SBATCH --mem-per-cpu=" $HOME/ProgdynSuite/SubmitScripts/g09.sh`
echo "Memory per CPU:            ${mem_per_cpu/"#SBATCH --mem-per-cpu="}"
echo
echo All of the .com files in the current directory "$PWD":
ls *.com # shows the user all of the com files in the current directory
echo
for file in "${com_files[@]}" # for each com file
do
	echo Would you like to submit "$file"? # prompt if the user wants to submit the com file
	echo Enter y/n:
	read y_n_response # get the response
	if [ "$y_n_response" = "y" -o "$y_n_response" = "Y" ] ; then # case insensitive
		cd $HOME/ProgdynSuite/RunOutputs
		remove=".com"
		JOB_NAME=${file/$remove} # JOB_NAME is the name of the com file minus the .com file extension
		directories=( `find -maxdepth 1 -type d`) # creates an array of every directory name in the RunOutputs folder
		readarray -t sorted < <(for a in "${directories[@]}"; do echo "$a"; done | sort) # sorts the array of directory names
		newdirnum=${sorted[@]:(-1)} # gets the last directory name which will always be the highest numbered temp folder
		if [ "${newdirnum/.\/temp}" = "$newdirnum" ] ; then # there are no temp directories
			newdirnum=1 # start at temp1
		else
			newdirnum=${newdirnum/.\/temp} # get the number of the highest numbered temp directory
			((newdirnum++)) # increment it by one
		fi
		TEMPORARY_DIR=$HOME/ProgdynSuite/RunOutputs/"temp$newdirnum" # file path of our new temp directory
		mkdir $HOME/ProgdynSuite/RunOutputs/temp$newdirnum # make the new temp directory
		cp $HOME/ProgdynSuite/SubmitScripts/g09.sh $HOME/ProgdynSuite/RunOutputs/temp$newdirnum # copy the script to run gaussian to the temp directory
		sed -i "s~export JOB_NAME=REPLACEME~export JOB_NAME=$JOB_NAME~g" $HOME/ProgdynSuite/RunOutputs/temp"$newdirnum"/g09.sh # set the JOB_NAME varable
		sed -i "s~export FILES_NEEDED_DIR=REPLACEME~export FILES_NEEDED_DIR="$TEMPORARY_DIR"~g" $HOME/ProgdynSuite/RunOutputs/temp"$newdirnum"/g09.sh # set the TEMPORARY_DIR variable
		cp $com_directory/$file $HOME/ProgdynSuite/RunOutputs/temp$newdirnum # copy the com file to the new temp directory
		cd $HOME/ProgdynSuite/RunOutputs/temp$newdirnum # change the current directory to the new temp directory
		echo Submitting $file
		sbatch g09.sh # submit g09.sh to the scheduling queue
		rm $HOME/ProgdynSuite/RunOutputs/temp"$newdirnum"/g09.sh
	elif [ "$y_n_response" = "n" -o "$y_n_response" = "N" ] ; then
		: # do nothing
	else
		echo "$y_n_response" is an unrecognized response. # do nothing and inform the user the input was not recognized
	fi
done

exit 0
