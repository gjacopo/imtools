#!/bin/bash


echo  REPLACE FOLDER PATHS
for d in *; 
do
    if ! [ -d $d ]; then 
	continue;
    fi
    for file in */*html;
    do
	replace $file ../../$d/html/ ../$d/
    done
    break
done
exit
echo
echo  REPLACE FILENAMES
for f in */*html;
do
    fname=`basename $f .html`
    Ufname=`echo $fname | tr /a-z/ /A-Z/`
    for file in */*html;
    do
	replace $file $Ufname.html $fname.html
    done
    break
done
    
