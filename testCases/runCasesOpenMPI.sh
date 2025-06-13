#!/bin/bash

CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions -I$(PWD)/src-local JumpingBubbles.c -o JumpingBubbles -lm

./JumpingBubbles