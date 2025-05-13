#!/bin/bash 

List=("Dyna" ) 
NList=${#List[@]} 

for ((i=0; i<$NList; i++)) 
do 
	gfortran -O3 -o Run${List[i]} Get${List[i]}.f90
done 