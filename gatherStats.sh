#!/bin/bash

make

for value in {1..5}
do
    echo $value
    ./program
    sleep 0.5
done

