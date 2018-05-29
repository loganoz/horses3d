#!/bin/bash
# 
#	This script converts all *.hsol files found in the RESULTS folder
#	checking before if an existing *.tec is newer or older than the
#	*.hsol file.
#########################################################################
#
#
# Configuration:----------------------------------------
HORSES2PLT_FLAGS= #INSERT FLAGS FOR HORSES2PLT HERE!
MESH_NAME= #INSERT MESH FILE HERE!
FOLDER=./RESULTS
########################################################

########################################################
for entry in "$FOLDER"/*.hsol
do
	tecname=$(echo "${entry%.*}").tec
	if [ -f $tecname ]; then
		if [ $entry -nt $tecname ]; then
  			horses2plt $MESH_NAME $entry $HORSES2PLT_FLAGS
		else
			echo "$entry conversion skipped since a more recent *.tec file exists"
		fi
	else
  		horses2plt $MESH_NAME $entry $HORSES2PLT_FLAGS
	fi
done
########################################################
