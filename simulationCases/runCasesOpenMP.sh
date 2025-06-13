#!/bin/bash

qcc -O2 -Wall -disable-dimensions -fopenmp -I$(PWD)/src-local JumpingBubbles.c -o JumpingBubbles -lm

export OMP_NUM_THREADS=4
./JumpingBubbles