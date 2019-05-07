#!/bin/bash
export PYTHONPATH=/nfs/production3/ma/home/atlas3-production/sw/lib/python3.7/site-packages
export LD_LIBRARY_PATH=$ATLAS_PROD/sw/ssl/lib
./hca2mtab.py $@
