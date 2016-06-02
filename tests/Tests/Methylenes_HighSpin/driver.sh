#! /bin/sh 
#
# driver script for the Acetone_OUR_CD test
# runs LSDCI using our CD implementation

# TigerCI location
export TIGER=$1

# other usual environmental variables for TigerCI
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='Methylenes_HighSpin'
export Basis='cc-pVDZ'
export benchmark_energy_sing='-77.83628120'
export benchmark_energy_trip='-77.77532651'
export benchmark_energy_quint='-77.71560757'


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
48
NUMBER OF BASIS FUNCTIONS
48
NUMBER OF FROZEN ORBITALS
6
NUMBER OF INACTIVE ORBITALS
6
NUMBER OF ACTIVE ORBITALS
4
SPIN MULTIPLICITY
1
NUMBER OF ELECTRONS
16
NUMBER OF REFERENCES
19
REFERENCE OCCUPATIONS
2222222200
2222222020
2222222002
2222220220
2222220202
2222220022
2222220211
2222222011
2222221102
2222220112
2222221201
2222222101
2222221210
2222221120
2222220121
2222221021
2222222110
2222221012
2222221111
BASIS SET
$Basis
SCRATCH DIRECTORY
$WorkDir/
CD THRESHOLD
1E-15
AO INTEGRAL THRESHOLD
1.0E-14
NONLOCAL
1
REFERENCE CI FLAG
0
ORB FILE
$Project.molden
NUM THREADS
$2
END OF INPUT
EOF

# Finally we run the actual calculation 
echo " "
echo "$Project" 
echo "*********************************"
$TIGER/tiger_ci.exe $WorkDir/$Project.input > $WorkDir/$Project.output1 2>&1
cp $Project.output1 ..

# Test the result
echo "testing the final singlet energy"
energy=$(grep -i 'total energy' $WorkDir/$Project.output1 | awk '{print $10}')
difference=$(echo "$energy - $benchmark_energy_sing" | bc | awk {'printf("%2.8E",$1)'} )
echo "final singlet energy = $energy"
echo "benchmark    = $benchmark_energy_sing"
echo "difference   = $difference"
result=$( echo $difference | awk ' function abs(x){if (x<0) return -x; return x} {if (abs($1) < 1.0E-6)  {print "good"}else{print "bad"}}')
if [ $result = 'good' ] 
then
	echo "Test is GOOD"
else
	echo "Test failed"
	exit 2
fi

rm -rf $WorkDir/*
cp ../$Project.molden $WorkDir

# here are the actual TigerCI input files for singlet
cat > $Project.input << EOF
NUMBER OF ORBITALS
48
NUMBER OF BASIS FUNCTIONS
48
NUMBER OF FROZEN ORBITALS
6
NUMBER OF INACTIVE ORBITALS
6
NUMBER OF ACTIVE ORBITALS
4
SPIN MULTIPLICITY
3
NUMBER OF ELECTRONS
16
NUMBER OF REFERENCES
19
REFERENCE OCCUPATIONS
2222222200
2222222020
2222222002
2222220220
2222220202
2222220022
2222220211
2222222011
2222221102
2222220112
2222221201
2222222101
2222221210
2222221120
2222220121
2222221021
2222222110
2222221012
2222221111
BASIS SET
$Basis
SCRATCH DIRECTORY
$WorkDir/
CD THRESHOLD
1E-15
AO INTEGRAL THRESHOLD
1.0E-14
NONLOCAL
1
REFERENCE CI FLAG
0
ORB FILE
$Project.molden
NUM THREADS
$2
END OF INPUT
EOF

# Finally we run the actual calculation 
echo " "
$TIGER/tiger_ci.exe $WorkDir/$Project.input > $WorkDir/$Project.output3 2>&1
cp $Project.output3 ..

# Test the result
echo "testing the final triplet energy"
energy=$(grep -i 'total energy' $WorkDir/$Project.output3 | awk '{print $10}')
difference=$(echo "$energy - $benchmark_energy_trip" | bc | awk {'printf("%2.8E",$1)'} )
echo "final triplet energy = $energy"
echo "benchmark    = $benchmark_energy_trip"
echo "difference   = $difference"
result=$( echo $difference | awk ' function abs(x){if (x<0) return -x; return x} {if (abs($1) < 1.0E-6)  {print "good"}else{print "bad"}}')
if [ $result = 'good' ] 
then
	echo "Test is GOOD"
else
	echo "Test failed"
	exit 2
fi

rm -rf $WorkDir/*
cp ../$Project.molden $WorkDir

# here are the actual TigerCI input files
cat > $Project.input << EOF
NUMBER OF ORBITALS
48
NUMBER OF BASIS FUNCTIONS
48
NUMBER OF FROZEN ORBITALS
6
NUMBER OF INACTIVE ORBITALS
6
NUMBER OF ACTIVE ORBITALS
4
SPIN MULTIPLICITY
5
NUMBER OF ELECTRONS
16
NUMBER OF REFERENCES
19
REFERENCE OCCUPATIONS
2222222200
2222222020
2222222002
2222220220
2222220202
2222220022
2222220211
2222222011
2222221102
2222220112
2222221201
2222222101
2222221210
2222221120
2222220121
2222221021
2222222110
2222221012
2222221111
BASIS SET
$Basis
SCRATCH DIRECTORY
$WorkDir/
CD THRESHOLD
1E-15
AO INTEGRAL THRESHOLD
1.0E-14
NONLOCAL
1
REFERENCE CI FLAG
0
ORB FILE
$Project.molden
NUM THREADS
$2
END OF INPUT
EOF

# Finally we run the actual calculation 
echo " "
$TIGER/tiger_ci.exe $WorkDir/$Project.input > $WorkDir/$Project.output5 2>&1
cp $Project.output5 ..

# Test the result
echo "testing the final quintet energy"
energy=$(grep -i 'total energy' $WorkDir/$Project.output5 | awk '{print $10}')
difference=$(echo "$energy - $benchmark_energy_quint" | bc | awk {'printf("%2.8E",$1)'} )
echo "final quintet energy = $energy"
echo "benchmark    = $benchmark_energy_quint"
echo "difference   = $difference"
result=$( echo $difference | awk ' function abs(x){if (x<0) return -x; return x} {if (abs($1) < 1.0E-6)  {print "good"}else{print "bad"}}')
if [ $result = 'good' ] 
then
	echo "Test is GOOD"
else
	echo "Test failed"
	exit 2
fi
