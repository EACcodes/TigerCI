#! /bin/sh

# TigerCI location 
export TIGER=$1


# other usual environmental variables for MOLCAS
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='methanol'
ref_energy=-115.06982215          # from the CASSCF calculation w/ MOLCAS
MRACPF_energy=-115.41990799              #from MOLCAS

# make/remove if necessary the workdir 
if [ -d $WorkDir ] 
then 
	'rm' -r $WorkDir
fi 
mkdir $WorkDir
cd $WorkDir


# copy over the input files
cp $ThisDir/$Project.molden $WorkDir
cp $ThisDir/*.input   $WorkDir

# set the input file correctly
sed -i".sed" "s:BLANK:$Project:"              $WorkDir/*.input
sed -i".sed" "s:SCRATCHPLACEHOLDER:$WorkDir:" $WorkDir/*.input
sed -i".sed" "s:THREADSPLACEHOLDER:$2:" $WorkDir/*.input

# Finally we run the actual calculation 
echo " "
echo "$Project-equ-MRACPF-test" 
echo "*********************************"
$TIGER/tiger_ci.exe $WorkDir/$Project.input > $WorkDir/$Project.output 2>&1 

# Test the result
echo "testing the final energy"
python $ThisDir/error-check.py $ref_energy $MRACPF_energy $WorkDir/$Project.output '1.0E-6'
if [ $? -ne 0  ]
then
    exit 2
fi

