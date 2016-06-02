#! /bin/sh
#
# driver script for the Methane_SDCI test

# TigerCI location 
export TIGER=$1

# other usual environmental variables for TigerCI
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='Acetylene_SuperMol_MRACPF2'
export Basis='cc-pVDZ'
benchmark_energy='-76.9042985537'      # this value has NOT been verified with an external code

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
OCCUPATION THRESHOLD
0.5
VIRTUAL OCCUPATION THRESHOLD
0.5
WP DEFAULT RADIUS
0.5
WP RADIUS MULTIPLIER
1.2
TOV VIRTUAL DEFAULT RADIUS
0.5
TOV VIRTUAL RADIUS MULTIPLIER
1.2
TOV OCCUPIED DEFAULT RADIUS
0.5
TOV OCCUPIED RADIUS MULTIPLIER
1.2
TOV CYLINDER RADIUS
2.0
CYLINDER RADIUS
2.0
ACPF FLAG
2
CD MEMORY
512
FULLY INTEGRAL DIRECT MODE
NUM THREADS
$2
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

