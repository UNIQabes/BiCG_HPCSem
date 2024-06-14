#!/bin/zsh

gammaList=( 0.1 0.6 0.9 1 )

for gamma in $gammaList
do
	bin/bicg 1000 $gamma > report/N1000G${gamma}.csv
	bin/bicgstab 1000 $gamma > report/N1000G${gamma}_STAB.csv
done