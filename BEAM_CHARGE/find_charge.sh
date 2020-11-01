#!/bin/sh
rm -f outcharge.csv
while read run
do
	echo $run
	root -l -b -q "run_charge.C($run)"
done<run_list.list
