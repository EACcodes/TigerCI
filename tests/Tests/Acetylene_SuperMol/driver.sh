#! /bin/sh
#
# driver script for the Methane_SDCI test

# TigerCI location 
export TIGER=$1
export NUM_THREADS=$2
if [ $NUM_THREADS -gt 2 ]
then
    export NUM_THREADS=2
fi

# other usual environmental variables for TigerCI
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='Acetylene_SuperMol'
export Basis='cc-pVDZ'
benchmark_energy='-76.87830672'

# make/remove if necessary the workdir 
if [ -d $WorkDir ] 
then 
	'rm' -r $WorkDir
fi 
mkdir $WorkDir
cp $Project.molden $WorkDir
cd $WorkDir

# here are the actual TigerCI input files
cat > $Project.input <<EOF
NUMBER OF ORBITALS
40
NUMBER OF BASIS FUNCTIONS
40
NUMBER OF INACTIVE ORBITALS
6
NUMBER OF ACTIVE ORBITALS
2
SPIN MULTIPLICITY
1
NUMBER OF ELECTRONS
14
NUMBER OF REFERENCES
3
REFERENCE OCCUPATIONS
22222220
22222211
22222202
SCRATCH DIRECTORY
$WorkDir/
CD THRESHOLD
1.00000E-12
NONLOCAL
1
CD MEMORY
1
NUM THREADS
$NUM_THREADS
CYLINDERS
7: 1 3
8: 1 3
ORB FILE
$Project.molden
BASIS SET
$Basis
CARTESIAN
EOF

# Finally we run the actual calculation 
echo " "
echo "$Project" 
echo "*********************************"
$TIGER/tiger_ci.exe $WorkDir/$Project.input > $WorkDir/$Project.output 2>&1 

# Test the result
echo "testing the final energy"
energy=$(grep -i 'total energy' $WorkDir/$Project.output | awk '{print $10}')
difference=$(echo "$energy - $benchmark_energy" | bc | awk {'printf("%2.8E",$1)'} )
echo "final energy = $energy"
echo "benchmark    = $benchmark_energy"
echo "difference   = $difference"
result=$( echo $difference | awk ' function abs(x){if (x<0) return -x; return x} {if (abs($1) < 1.0E-6)  {print "good"}else{print "bad"}}')
if [ $result = 'good' ] 
then
	echo "Test is GOOD"
else
	echo "Test failed"
	exit 2
fi

