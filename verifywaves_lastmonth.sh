#!/bin/bash


export PYTHONPATH="$PYTHONPATH:/home/johannesro/.fou_python/lib/python2.7/site-packages"
export PYTHONPATH="$PYTHONPATH:/home/johannesro/python"
export R_LIBS_SITE="/metno/vvfelles/pakker/R-Pakker"
export R_LIBS="/disk1/local/R-packages"

# get current month
printf -v year %2.2i `date +%Y`
printf -v month %2.2i `date +%m`

# calculate which is the previous month:
lmonth=`date -d "${year}/${month}/01-1 days" +"%m"`
lyear=`date -d "${year}/${month}/01-1 days" +"%Y"`

cd /home/johannesro/cronjobs/waveverification

echo 'start wave verification script' >> wv.log

# gather data
./collectdata.py ${lyear} ${lmonth} 

# make plots
./validate.py ${lyear}${lmonth}

# update webpage
cd /disk4/waveverification/website 
./makepage.sh $lyear $lmonth
./makepage_DD.sh $lyear $lmonth
./makeindex.sh $lyear $lmonth

# upload webpage
cd /vol/hindcast3/waveverification
./upload.sh

