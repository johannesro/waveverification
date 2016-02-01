#!/bin/bash


#export PYTHONPATH="$PYTHONPATH:/home/johannesro/.fou_python/lib/python2.7/site-packages"
#export PYTHONPATH="$PYTHONPATH:/home/johannesro/python"
#export R_LIBS_SITE="/metno/vvfelles/pakker/R-Pakker"
#export R_LIBS="/disk1/local/R-packages"

printf -v year %2.2i `date +%Y`
printf -v month %2.2i `date +%m`

cd /vol/hindcast3/waveverification/waveverification

echo 'start wave verification script' >> wv.log

# gather data
./collectdata.py 

# make plots
#./validate.py

# update webpage
#cd /vol/hindcast3/waveverification/website 

# september is not working with the above commands because the shell misinterprets the lines !?
#./makepage.sh $year 09
#./makepage_DD.sh $year 09
#./makeindex.sh $year 09

#./makepage.sh $year $month
#./makepage_DD.sh $year $month
#./makeindex.sh $year $month


# upload webpage
#cd /vol/hindcast3/waveverification
#./upload.sh

# secure data on local pc
#cd /disk4/waveverification/
#rsync -ua /vol/hindcast3/waveverification/data .



