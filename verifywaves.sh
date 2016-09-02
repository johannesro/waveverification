#!/bin/bash


#export PYTHONPATH="$PYTHONPATH:/home/johannesro/.fou_python/lib/python2.7/site-packages"
#export PYTHONPATH="$PYTHONPATH:/home/johannesro/python"
#export R_LIBS_SITE="/metno/vvfelles/pakker/R-Pakker"
#export R_LIBS="/disk1/local/R-packages"

#printf -v year %2.2i `date +%Y`
#printf -v month %2.2i `date +%m`

#cd /home/johannesro/waveverification

echo 'start wave verification script' >> wv.log

# gather data
./collectdata.py 

# make plots
./validate.py

# update webpage
cd /lustre/storeB/project/fou/hi/waveverification/website

./makepage.sh 
./makepage_DD.sh 
./makeindex.sh 


# upload webpage
cd /lustre/storeB/project/fou/hi/waveverification
./upload.sh

# secure data on local pc
#cd /disk4/waveverification/
#rsync -ua /vol/hindcast3/waveverification/data .



