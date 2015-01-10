#!/bin/bash

directories=( `find -maxdepth 1 -type d` )
top_dir=$PWD
for dir in "${directories[@]}"
do
	if ( test -s $top_dir/$dir/deleteme ) then
		echo Deleting $dir
		rm -rf $top_dir/$dir
	fi
done
exit 0
