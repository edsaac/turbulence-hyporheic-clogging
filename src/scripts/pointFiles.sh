##  
## This script takes the postprocessing/boundaryProbes
## data and organizes it in three files for handling in
## python later
##
## ./pointFiles.sh <pathToFolder> <nameField>

INPUT_DIR=$1
OUTPUT_DIR=./organizedData
FIELD=$2

# Reset output directory
rm -r $OUTPUT_DIR
mkdir $OUTPUT_DIR

# Create a time file 
TIME_FILE=$OUTPUT_DIR/time.txt

TIMES=$(ls -v $INPUT_DIR)
LEN_TIMES=${#TIMES[@]}

echo $TIMES > $TIME_FILE
sed -i "s/ /\n/g" $TIME_FILE

# Create a coordinates file
COORD_FILE=$OUTPUT_DIR/xyz.txt

RAND_FILE=$INPUT_DIR/$(ls -v $INPUT_DIR | shuf -n 1)/points_$FIELD.xy
sed -e "s/\s\s*/ /g" $RAND_FILE | cut -d " " -f 1,2,3 > $COORD_FILE

# Read each time and generate a pressure array
# with each timestep as rows and coords points as cols
PRESSURE_FILE=$OUTPUT_DIR/$FIELD.txt

for T in $TIMES
do
  P_FILE=$INPUT_DIR/$T/points_$FIELD.xy
  sed -e "s/\s\s*/ /g" $P_FILE | cut --complement -d " " -f 1,2,3 | tr '\n' ' ' >> $PRESSURE_FILE
  echo "" >> $PRESSURE_FILE
done
sed -i "s/ $//g" $PRESSURE_FILE
