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
SOLUTION_FILTER=
########################################################



########################################################
RED='\033[0;31m'
GREEN='\033[0;32m'
GRAY='\033[0m'
NC='\033[0m' # No Color
########################################################

i=0
entries="$FOLDER"/"$SOLUTION_FILTER"*.hsol
no_of_entries=0
for entry in $entries
do
	no_of_entries=$((no_of_entries+1))
done

for entry in $entries
do
	i=$((i+1))
	tecname=$(echo "${entry%.*}").tec
	if [ -f $tecname ]; then
		if [ $entry -nt $tecname ]; then
			printf "${GREEN}(%d/%d) %s...${NC}" "$i" "$no_of_entries" "Converting $entry"
			printf "${RED}"
  			horses2plt $MESH_NAME $entry $HORSES2PLT_FLAGS > /dev/null 
			printf "${NC}"
			if [ $? -eq 0 ]
			then
				printf "${GREEN}done!${NC}\n"
			fi
		else
			printf "${GRAY}(%d/%d) %s${NC}\n" "$i" "$no_of_entries" "$entry conversion skipped since a more recent *.tec file exists"
			
		fi
	else
		printf "${GREEN}(%d/%d) %s...${NC}" "$i" "$no_of_entries" "Converting $entry"
		printf "${RED}"
  		horses2plt $MESH_NAME $entry $HORSES2PLT_FLAGS > /dev/null
		printf "${NC}"
		if [ $? -eq 0 ]
		then
			printf "${GREEN}done!${NC}\n"
		fi
	fi
done
########################################################
