#!/bin/sh

SRC=/home/andrew/GSSP/Radiative_transfer
DST=/scratch/$USER


if [ ! -e $DST ]
then
	mkdir $DST
fi

cp $SRC/* $DST
