#!/bin/sh

for i in $(seq 1 10); do
    nohup julia optimize.jl $i >> logs/$i.log 2>&1 &
done
    