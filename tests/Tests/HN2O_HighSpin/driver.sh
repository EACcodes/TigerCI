#! /bin/sh 
#
# driver script for the Acetone_OUR_CD test
# runs LSDCI using our CD implementation

# TigerCI location
export TIGER=$1

# other usual environmental variables for TigerCI
export ThisDir=`pwd`
export WorkDir="$ThisDir/workdir"
export Project='HN2O_HighSpin'
export Basis='cc-pVDZ'
export benchmark_energy_doub='-184.26682995'
export benchmark_energy_quart='-184.00569535'
export benchmark_energy_sext='-183.74176208'


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
47
NUMBER OF BASIS FUNCTIONS
47
NUMBER OF FROZEN ORBITALS
9
NUMBER OF INACTIVE ORBITALS
9
NUMBER OF ACTIVE ORBITALS
5
SPIN MULTIPLICITY
2
NUMBER OF ELECTRONS
23
NUMBER OF REFERENCES
51
REFERENCE OCCUPATIONS
22222222222100
22222222222010
22222222222001
22222222221200
22222222220210
22222222220201
22222222221020
22222222220120
22222222220021
22222222221002
22222222220102
22222222220012
22222222212200
22222222202210
22222222202201
22222222212020
22222222202120
22222222202021
22222222212002
22222222202102
22222222202012
22222222210220
22222222201220
22222222200221
22222222210202
22222222201202
22222222200212
22222222210022
22222222201022
22222222200122
22222222221110
22222222221101
22222222221011
22222222220111
22222222212110
22222222212101
22222222212011
22222222202111
22222222211210
22222222211201
22222222210211
22222222201211
22222222211120
22222222211021
22222222210121
22222222201121
22222222211102
22222222211012
22222222210112
22222222201112
22222222211111
BASIS SET
$Basis
SCRATCH DIRECTORY
$WorkDir/
CD THRESHOLD
1E-12
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
echo "testing the final doublet energy"
energy=$(grep -i 'total energy' $WorkDir/$Project.output1 | awk '{print $10}')
difference=$(echo "$energy - $benchmark_energy_doub" | bc | awk {'printf("%2.8E",$1)'} )
echo "final doublet energy = $energy"
echo "benchmark    = $benchmark_energy_doub"
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
47
NUMBER OF BASIS FUNCTIONS
47
NUMBER OF FROZEN ORBITALS
9
NUMBER OF INACTIVE ORBITALS
9
NUMBER OF ACTIVE ORBITALS
5
SPIN MULTIPLICITY
4
NUMBER OF ELECTRONS
23
NUMBER OF REFERENCES
51
REFERENCE OCCUPATIONS
22222222222100
22222222222010
22222222222001
22222222221200
22222222220210
22222222220201
22222222221020
22222222220120
22222222220021
22222222221002
22222222220102
22222222220012
22222222212200
22222222202210
22222222202201
22222222212020
22222222202120
22222222202021
22222222212002
22222222202102
22222222202012
22222222210220
22222222201220
22222222200221
22222222210202
22222222201202
22222222200212
22222222210022
22222222201022
22222222200122
22222222221110
22222222221101
22222222221011
22222222220111
22222222212110
22222222212101
22222222212011
22222222202111
22222222211210
22222222211201
22222222210211
22222222201211
22222222211120
22222222211021
22222222210121
22222222201121
22222222211102
22222222211012
22222222210112
22222222201112
22222222211111
BASIS SET
$Basis
SCRATCH DIRECTORY
$WorkDir/
CD THRESHOLD
1E-12
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
echo "testing the final quartet energy"
energy=$(grep -i 'total energy' $WorkDir/$Project.output3 | awk '{print $10}')
difference=$(echo "$energy - $benchmark_energy_quart" | bc | awk {'printf("%2.8E",$1)'} )
echo "final quartet energy = $energy"
echo "benchmark    = $benchmark_energy_quart"
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
47
NUMBER OF BASIS FUNCTIONS
47
NUMBER OF FROZEN ORBITALS
9
NUMBER OF INACTIVE ORBITALS
9
NUMBER OF ACTIVE ORBITALS
5
SPIN MULTIPLICITY
6
NUMBER OF ELECTRONS
23
NUMBER OF REFERENCES
51
REFERENCE OCCUPATIONS
22222222222100
22222222222010
22222222222001
22222222221200
22222222220210
22222222220201
22222222221020
22222222220120
22222222220021
22222222221002
22222222220102
22222222220012
22222222212200
22222222202210
22222222202201
22222222212020
22222222202120
22222222202021
22222222212002
22222222202102
22222222202012
22222222210220
22222222201220
22222222200221
22222222210202
22222222201202
22222222200212
22222222210022
22222222201022
22222222200122
22222222221110
22222222221101
22222222221011
22222222220111
22222222212110
22222222212101
22222222212011
22222222202111
22222222211210
22222222211201
22222222210211
22222222201211
22222222211120
22222222211021
22222222210121
22222222201121
22222222211102
22222222211012
22222222210112
22222222201112
22222222211111
BASIS SET
$Basis
SCRATCH DIRECTORY
$WorkDir/
CD THRESHOLD
1E-12
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
echo "testing the final sextet energy"
energy=$(grep -i 'total energy' $WorkDir/$Project.output5 | awk '{print $10}')
difference=$(echo "$energy - $benchmark_energy_sext" | bc | awk {'printf("%2.8E",$1)'} )
echo "final sextet energy = $energy"
echo "benchmark    = $benchmark_energy_sext"
echo "difference   = $difference"
result=$( echo $difference | awk ' function abs(x){if (x<0) return -x; return x} {if (abs($1) < 1.0E-6)  {print "good"}else{print "bad"}}')
if [ $result = 'good' ] 
then
	echo "Test is GOOD"
else
	echo "Test failed"
	exit 2
fi
