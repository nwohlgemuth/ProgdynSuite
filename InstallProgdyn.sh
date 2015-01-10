#!/bin/bash

export INSTALL_SCRIPT_LOCATION=$PWD

cd $HOME

if grep -Fxq "PATH=\$PATH:$HOME/ProgdynSuite/SubmitScripts" .bash_profile
then
	:
else
	echo "PATH=\$PATH:$HOME/ProgdynSuite/SubmitScripts" >> $HOME/.bash_profile
fi

if [ -d "ProgdynSuite" ]
then
	line=$(head -1 $HOME/ProgdynSuite/readMe)
	echo An Installation of ProgdynSuite was detected: $line
	echo Would you like to keep the run outputs from this version?
	echo Enter y/n:
	read y_n_response
	if [ "$y_n_response" = "y" -o "$y_n_response" = "Y" ] ; then
		mkdir $HOME/ProgdynSuiteTemp
		cp -a /$HOME/ProgdynSuite/RunOutputs/. /$HOME/ProgdynSuiteTemp/	
		rm -rf "ProgdynSuite"
		mkdir ProgdynSuite
		mkdir ProgdynSuite/RunOutputs
		cp -a /$HOME/ProgdynSuiteTemp/. /$HOME/ProgdynSuite/RunOutputs/
		rm -rf $HOME/ProgdynSuiteTemp
	elif [ "$y_n_response" = "n" -o "$y_n_response" = "N" ] ; then
		rm -rf "ProgdynSuite"
		mkdir ProgdynSuite
		mkdir ProgdynSuite/RunOutputs
	else
		echo Input not recognized. Performing a clean install...
		rm -rf "ProgdynSuite"
		mkdir ProgdynSuite
		mkdir ProgdynSuite/RunOutputs
	fi
else
	mkdir ProgdynSuite
	mkdir ProgdynSuite/RunOutputs
fi

cd $HOME/ProgdynSuite

mkdir ProgdynScripts
mkdir SubmitScripts

mv $INSTALL_SCRIPT_LOCATION/progdynstarterHP $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/ProgSubmit.sh $HOME/ProgdynSuite/SubmitScripts
mv $INSTALL_SCRIPT_LOCATION/freqinHP $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/progdyn.conf $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/proggenHP $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/prog1stpoint $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/prog2ndpoint $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/progdynb $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/proganal $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/randgen.c $HOME/ProgdynSuite/ProgdynScripts
mv $INSTALL_SCRIPT_LOCATION/ProgSubmitWrapper.sh $HOME/ProgdynSuite/SubmitScripts
mv $INSTALL_SCRIPT_LOCATION/RemoveEmptyTempFiles.sh $HOME/ProgdynSuite/RunOutputs
mv $INSTALL_SCRIPT_LOCATION/readMe $HOME/ProgdynSuite
mv $INSTALL_SCRIPT_LOCATION/progdynInfo $HOME/ProgdynSuite
mv $INSTALL_SCRIPT_LOCATION/g09.sh $HOME/ProgdynSuite/SubmitScripts
mv $INSTALL_SCRIPT_LOCATION/RunJob.sh $HOME/ProgdynSuite/SubmitScripts
mv $INSTALL_SCRIPT_LOCATION/MakeXYZ.sh $HOME/ProgdynSuite/SubmitScripts
mv $INSTALL_SCRIPT_LOCATION/num_atoms_coords.awk $HOME/ProgdynSuite/SubmitScripts

chmod -w $HOME/ProgdynSuite/ProgdynScripts/progdynstarterHP
chmod -w $HOME/ProgdynSuite/ProgdynScripts/proggenHP
chmod -w $HOME/ProgdynSuite/ProgdynScripts/prog1stpoint
chmod -w $HOME/ProgdynSuite/ProgdynScripts/prog2ndpoint
chmod -w $HOME/ProgdynSuite/ProgdynScripts/progdynb
chmod -w $HOME/ProgdynSuite/ProgdynScripts/proganal
chmod -w $HOME/ProgdynSuite/readMe
chmod -w $HOME/ProgdynSuite/progdynInfo

cd $HOME/ProgdynSuite/ProgdynScripts

g++ randgen.c -o randgen
rm randgen.c
chmod -w randgen

cd $HOME/ProgdynSuite

#sed -i 's~export PROG_HOME=REPLACEME~export PROG_HOME="$HOME/ProgdynSuite/ProgdynScripts"~g' ./SubmitScripts/ProgSubmit.sh

cd $INSTALL_SCRIPT_LOCATION

rm ProgdynSuite_1_3.zip
rm InstallProgdyn.sh

PATH=$PATH:$HOME/ProgdynSuite/SubmitScripts

exit 0
