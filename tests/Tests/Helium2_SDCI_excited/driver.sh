#! /bin/sh
#
# driver script for the Methane_exc test

# TigerCI location 
export TIGER=$1

# other usual environmental variables for TigerCI
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='Helium2_SDCI_excited'
export Basis='6-31G**'

benchmark_energy1=-5.7552182      # this value has been verified by an external code (MOLCAS)
benchmark_energy2=-3.8782561      # this value has NOT been verified by an external code

# make/remove if necessary the workdir 
if [ -d $WorkDir ] 
then 
	'rm' -r $WorkDir
fi 
mkdir $WorkDir
cp $Project.molden $WorkDir
cd $WorkDir

# here are the actual molcas input files
cat > $Project.input <<EOF
CARTESIAN
NUMBER OF ORBITALS
10
NUMBER OF BASIS FUNCTIONS
10
NUMBER OF INACTIVE ORBITALS
2
SPIN MULTIPLICITY
1
NUMBER OF ELECTRONS
4
NUMBER OF REFERENCES
1
REFERENCE OCCUPATIONS
22
SCRATCH DIRECTORY
$WorkDir/
INTEGRAL THRESHOLD
1.00000E-9
NONLOCAL
1
NUM ROOTS
2
NUM THREADS
$2
ORB FILE
$Project.molden
BASIS SET
$Basis
EOF

# Finally we run the actual calculation 
echo " "
echo "$Project" 
echo "*********************************"
$TIGER/tiger_ci.exe $WorkDir/$Project.input > $WorkDir/$Project.output 2>&1 

# Test the result
cd $ThisDir
echo "testing the final energy"
python error-check.py $benchmark_energy1 $benchmark_energy2 1.0E-6 $Project 
if [ -e workdir/FAILED ]
then
    exit 2
fi

