#!/bin/bash


i=0
while(($i < $Nruns)); do

  while(($running_tasks < $max_tasks && $i < $Nruns)); do

    set -- ${_params[$i]}
    T=${1}
    MU=${2}
    n=${3}

    LastFold="T${T}_mu${MU}"

    f1="/${LastFold}/settings_num${n}.txt"

    f="${path}${f1}"

    ((i++))

  done
done
wait


