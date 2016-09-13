#!/bin/bash

nnodes=72720
thread=8

for n in $(seq 0 $((nnodes-1))); do
#for n in $(seq 2221 5000); do

  nf=$(printf %06d $n)

  while (($(jobs -p | wc -l) >= thread)); do
    sleep 1s
  done

  echo $nf

  ( cd input/${nf} && tar chf ../../input_tar/${nf}.tar * ) &

done
wait
