#!/usr/bin/env bash



ARCHIVE=$1
TRAIN_DSPATH="${ARCHIVE}:train200"
TRAIN_DSPATH_H5="${ARCHIVE}/train200"
TEST_DSPATH="${ARCHIVE}:full"


for TRID in `h5ls ${TRAIN_DSPATH_H5} | grep tr | cut -d " " -f 1`
do
	for LAYER in 0 1 2 3
	do
		./mks.py --layer ${LAYER} ${TRAIN_DSPATH} ${TRID} ${TEST_DSPATH} >${TRID}.${LAYER}.repr
	done
done

