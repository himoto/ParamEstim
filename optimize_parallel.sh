#!/bin/sh

for i in $(seq 1 10); do
    nohup julia optimize.jl $i >> logs/$i.log 2>&1 &
done

# To terminate the process,
# $ pgrep -f optimize.jl | xargs kill -9

<<COMMENT_OUT
for i in $(seq 1 10); do
    julia optimize_continue.jl $i >> logs/$i.log 2>&1 &
done
COMMENT_OUT