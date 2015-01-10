#!/bin/bash
one_directory_up="$(dirname "$PWD")"
two_directories_up="$(dirname "$one_directory_up")"
xyz_file=${one_directory_up/$two_directories_up\/out}.xyz
if test -s $xyz_file
then
	rm $xyz_file
fi
all_log_files=( `find *.log` )
#echo ${all_log_files[${#all_log_files[@]}-1]}
longest_file=${all_log_files[0]}
for log_file in "${all_log_files[@]}"
do
	if (( ${#log_file} > ${#longest_file} ))
	then
		longest_file=$log_file
	fi
done
longest_file=${longest_file/g09}
for log_file in "${all_log_files[@]}"
do
	changed_value=0
	old_log_file_name=$log_file
	log_file_num=${log_file/g09}
	while (( ${#log_file_num} < ${#longest_file} ))
	do
		log_file_num="0"$log_file_num
		changed_value=1
	done
	if (($changed_value == "1"))
	then
		mv $old_log_file_name "g09"$log_file_num
	fi
done
all_log_files=( `find *.log`) 
for log_file in "${all_log_files[@]}"
do
	awk '/Standard orientation:/,/tional const/ {if ($3==0) print}' $log_file > tempstangeos
	awk -f "$HOME/ProgdynSuite/SubmitScripts/num_atoms_coords.awk" $log_file >> $xyz_file
done
rm tempstangeos
exit 0
