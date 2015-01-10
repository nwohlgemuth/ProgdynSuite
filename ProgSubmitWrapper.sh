#!/bin/bash

cd $HOME/ProgdynSuite/RunOutputs
directories=( `find -maxdepth 1 -type d`)
readarray -t sorted < <(for a in "${directories[@]}"; do echo "$a"; done | sort)
filename=${sorted[@]:(-1)}
if [ "${filename/.\/temp}" = "$filename" ] ; then
	export filename=1
else
	filename=${filename/.\/temp}
	((filename++))
	export filename
fi
mkdir temp$filename
cp -av $HOME/ProgdynSuite/ProgdynScripts/. $HOME/ProgdynSuite/RunOutputs/temp$filename/
cp $HOME/ProgdynSuite/SubmitScripts/ProgSubmit.sh $HOME/ProgdynSuite/RunOutputs/temp"$filename"/ProgSubmit.sh
export TEMPORARY_DIR=$HOME/ProgdynSuite/RunOutputs/temp"$filename"
sed -i "s~export FILES_NEEDED_DIR=REPLACEME~export FILES_NEEDED_DIR="$TEMPORARY_DIR"~g" $HOME/ProgdynSuite/RunOutputs/temp"$filename"/ProgSubmit.sh

cd $HOME/ProgdynSuite/RunOutputs/temp"$filename"
echo "SUBMITTING PROGSUBMIT.SH"
sbatch ProgSubmit.sh
echo "DONE SUBMITTING PROGSUBMIT.SH"

#rm $HOME/ProgdynSuite/RunOutputs/ProgSubmit.sh

exit 0
