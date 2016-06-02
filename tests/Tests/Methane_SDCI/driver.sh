#! /bin/sh
#
# driver script for the Methane_SDCI test

# TigerCI location 
export TIGER=$1


# other usual environmental variables for MOLCAS
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='Methane_SDCI'
export Basis='6-31G**'
benchmark_energy='-40.38325482'


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
35
NUMBER OF BASIS FUNCTIONS
35
NUMBER OF INACTIVE ORBITALS
5
SPIN MULTIPLICITY
1
NUMBER OF ELECTRONS
10
NUMBER OF REFERENCES
1
REFERENCE OCCUPATIONS
22222
SCRATCH DIRECTORY
$WorkDir/
INTEGRAL THRESHOLD
1.00000E-9
NONLOCAL
1
NUM THREADS
$2
BASIS SET
$Basis
ORB FILE
$Project.molden
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
