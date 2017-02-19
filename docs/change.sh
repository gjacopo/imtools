#!/bin/bash


# echo  REPLACE FOLDER PATHS
# for d in *; 
# do
#     if ! [ -d $d ]; then 
# 	continue;
#     fi
#     echo replacements for folder $d...
#     for file in */*html;
#     do
# 	replace $file ../../$d/html/ ../$d/
#     done
# done

echo
echo  REPLACE FILENAMES
for f in */*.html;
do
    echo replacements for file $f...
    fname=`basename $f .html`
    Ufname=`echo $fname | tr /a-z/ /A-Z/`
    for file in */*html;
    do
	replace $file $Ufname.html $fname.html
    done
done
    
