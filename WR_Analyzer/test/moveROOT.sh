#! /bin/bash

for x in $(ls *.root); do
    a=$(date "+%M_%H_%m_%d_%y_")
    echo $a$x
    mv $x plots/rootfiles/$a$x
done
