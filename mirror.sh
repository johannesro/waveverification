#!/bin/bash
# 
#Script to mirror data from StoreA to StoreB to have backup.
#Executed ones a month as a Cron-tab too secure the monthly 
#and quarterly Arc-MFC reports.
#
cp -r /lustre/storeA/project/fou/om/waveverification /lustre/storeB/project/fou/om/.
