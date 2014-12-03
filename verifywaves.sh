#!/bin/bash


export PYTHONPATH="$PYTHONPATH:/home/johannesro/.fou_python/lib/python2.7/site-packages"
export PYTHONPATH="$PYTHONPATH:/home/johannesro/python"
export R_LIBS_SITE="/metno/vvfelles/pakker/R-Pakker"
export R_LIBS="/disk1/local/R-packages"

printf -v year %2.2i `date +%Y`
printf -v month %2.2i `date +%m`

cd /home/johannesro/cronjobs/waveverification

echo 'start wave verification script' >> wv.log

# gather data
./collectdata.py 

# make plots
./validate.py

# update webpage
cd /disk4/waveverification/website 
./makepage.sh $year $month
./makepage_DD.sh $year $month
./makeindex.sh $year $month

# upload webpage
cd /disk4/waveverification
./upload.sh

