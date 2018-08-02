# Eric Erkela (COHERENT)

# This Bash script will generate the directory structure that we plan on using to store our 
# ORCA/Root files for our crystal characterization runs.

#!/bin/bash
echo ">> Crystal Serial Number? "
read SERIAL_NUM
DIR=$PWD

if [ ! -d $DIR/$SERIAL_NUM ]
then
	echo "making new crystal directory..."
	mkdir $SERIAL_NUM
fi

# generate position directories
cd $DIR/$SERIAL_NUM
if [ ! -d $DIR/$SERIAL_NUM/position ]
then 
	echo "making position directory..."
	mkdir position
else
	echo "Position directory already exists"
fi
cd position
for d in position_1 position_2 position_3 position_4 position_5
do
	if [ ! -d $DIR/$SERIAL_NUM/$d ]
	then
		echo "making $d directory..."
		mkdir $d
	else
		echo "$d directory already exists"
	fi
done


# Generate voltage directories
cd $DIR/$SERIAL_NUM
if [ ! -d $DIR/$SERIAL_NUM/voltage ]
then 
	echo "making voltage directory..."
	mkdir voltage
else
	echo "Voltage directory already exists"
fi
cd voltage
for d in 600_V 700_V 800_V 900_V 1000_V
do
	if [ ! -d $DIR/$SERIAL_NUM/$d ]
	then
		echo "Making $d Directory..."
		mkdir $d
	else
		echo "$d directory already exists"
	fi
done
