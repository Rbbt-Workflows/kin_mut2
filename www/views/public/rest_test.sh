#!/bin/bash

jobname=$(curl -H"Expect:" -X POST http://kinmut2.bioinfo.cnio.es/KinMut2/predict --form mutations="P21802 S252W,O60674 D584E" --form jobname=REST-TEST --form _format=jobname -s)

while ! curl http://kinmut2.bioinfo.cnio.es/KinMut2/predict/$jobname/info?_format=json -s |grep '"status":"done"' >/dev/null; do
    sleep 1
done
    
curl http://kinmut2.bioinfo.cnio.es/KinMut2/predict/$jobname?_format=raw -s 
