#!/bin/bash

make
./benchmark-blocked $1 1> results/blocked.out 2> results/blocked.error &
./benchmark-strassen $1 1> results/strassen.out 2> results/strassen.error &
./benchmark-naive $1 1> results/naive.out 2> results/naive.error &
