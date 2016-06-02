#! /bin/sh 
#
# driver script for the Acetone_OUR_CD test
# runs LSDCI using our CD implementation

# TigerCI location
export TIGER=$1

# other usual environmental variables for TigerCI
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='Acetone_LSDCI'
export Basis='cc-pVDZ'
benchmark_energy='-114.186985403061911'

# make/remove if necessary the workdir 
if [ -d $WorkDir ] 
then 
	'rm' -r $WorkDir
fi 
mkdir $WorkDir
cp $Project.molden $WorkDir
cd $WorkDir

# here are the actual TigerCI input files
cat > $Project.input << EOF
NUMBER OF ORBITALS
40
NUMBER OF BASIS FUNCTIONS
40
NUMBER OF INACTIVE ORBITALS
8
SPIN MULTIPLICITY
1
NUMBER OF ELECTRONS
16
NUMBER OF REFERENCES
1
REFERENCE OCCUPATIONS
22222222
SCRATCH DIRECTORY
$WorkDir/
CD THRESHOLD
1.00000E-7
OCCUPATION THRESHOLD
0.8
VIRTUAL OCCUPATION THRESHOLD
0.8
WP DEFAULT RADIUS
0.8 
WP RADIUS MULTIPLIER
1.38
TOV VIRTUAL DEFAULT RADIUS
0.8
TOV VIRTUAL RADIUS MULTIPLIER
1.725
TOV OCCUPIED DEFAULT RADIUS
0.8
TOV OCCUPIED RADIUS MULTIPLIER
1.725
TOV CYLINDER RADIUS
0.8
CYLINDER RADIUS
0.8
CD MEMORY
1024
NUM THREADS
$2
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
