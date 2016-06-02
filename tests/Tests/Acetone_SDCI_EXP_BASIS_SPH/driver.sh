#! /bin/sh 
#
# driver script for the Acetone_OUR_CD test
# runs LSDCI using our CD implementation

# TigerCI location
export TIGER=$1

# other usual environmental variables for TigerCI
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='Acetone_SDCI_EXP_BASIS'
export Basis='cc-pVTZ'
benchmark_energy='-114.317228454648685'
spherical_benchmark_energy='-114.299078614366792'

# make/remove if necessary the workdir 
if [ -d $WorkDir ] 
then 
	'rm' -r $WorkDir
fi 
mkdir $WorkDir
cp $Project.molden $WorkDir
cp $Project.spherical.molden $WorkDir
cd $WorkDir

# here are the actual TigerCI input files
cat > $Project.input << EOF
NUMBER OF ORBITALS
100
NUMBER OF BASIS FUNCTIONS
100
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
1.00000E-10
INTEGRAL THRESHOLD
1.00000E-9
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
NONLOCAL
1
ORB FILE
$Project.molden
BASIS SET
$Basis
CARTESIAN
EOF

cat > $Project.spherical.input << EOF
NUMBER OF ORBITALS
88
NUMBER OF BASIS FUNCTIONS
88
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
INTEGRAL THRESHOLD
1.00000E-9
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
$Project.spherical.molden
BASIS SET
$Basis
EOF

# Finally we run the actual calculation 
echo " "
echo "Acetone_SDCI_EXP_BASIS_SPH" 
echo "*********************************"
$TIGER/tiger_ci.exe $WorkDir/$Project.input > $WorkDir/$Project.output 2>&1
$TIGER/tiger_ci.exe $WorkDir/$Project.spherical.input > $WorkDir/$Project.spherical.output 2>&1

energy=$(grep -i 'total energy' $WorkDir/$Project.spherical.output | awk '{print $10}')
difference=$(echo "$energy - $spherical_benchmark_energy" | bc | awk {'printf("%2.8E",$1)'} )
echo "final energy = $energy"
echo "benchmark    = $spherical_benchmark_energy"
echo "difference   = $difference"
result=$( echo $difference | awk ' function abs(x){if (x<0) return -x; return x} {if (abs($1) < 1.0E-6)  {print "good"}else{print "bad"}}')
if [ $result = 'good' ] 
then
	echo "Test is GOOD"
else
	echo "Test failed"
#	exit 2
fi

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
