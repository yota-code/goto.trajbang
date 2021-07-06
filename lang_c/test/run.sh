#!/usr/bin/env zsh

rm ./main.exe *.o

gcc -c ../src/trajbang3.c -I../include -o trajbang3.o
gcc main.c -I../include trajbang3.o -lm -o main.exe

./main.exe 34