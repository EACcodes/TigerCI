#! /bin/bash
#PBS -q eac
#PBS -l nodes=1:ppn=8,walltime=5:00:00
#PBS -V
#PBS -N ACPF_1
#PBS -M dkrisilo@princeton.edu
#PBS -m abe

# Set the relevant variables ....
export ThisDir=$PBS_O_WORKDIR
export PID=$$
export WorkDir="/scratch/dkrisilo/${PID}"
export MOLCAS="$WorkDir/molcas74"
export Project='Methane'
export Basis1='6-31G**'
export Basis2='STO-3G'


echo "Starting CI run "
echo "-----------------------------"
echo "Variable check ...."
echo "This directory   = " $ThisDir
echo "Work directory   = " $WorkDir
echo "Molcas directory = " $MOLCAS
echo "Basis sets       = " $Basis1 $Basis2
echo "PID              = " $PID
echo "-----------------------------"
echo ""
echo "Time: " `date`
echo ""
 
mkdir -p $WorkDir 
cd $WorkDir 


# Here is the place to copy over anything we might need for later
cp ${ThisDir}/LOVO_prep.py $WorkDir

# Unpack MOLCAS
cp /home/dkrisilo/molcas.tgz $WorkDir
tar -xzf molcas.tgz


# Time to make the input files 
#--------------------------------
#   SEWARD FILE
#--------------------------------
cat > $Project.Seward.in << EOF
&SEWARD &END
Title
$Project with $Basis1 basis
Basis set
C.$Basis1....
C1	0.0000000000	0.0000000000	0.0000000000
End of Basis
Basis set
H.$Basis1....
H1	0.0000000	0.0000000	2.0579126
H2	1.9402185	0.0000000	-0.6859709
H3	-0.9701102	-1.6802790	-0.6859709
H4	-0.9701102	1.6802790	-0.6859709
End of Basis
CHOINPUT
THRCholesky
1.0D-15
NOPRESCREEN
ENDCHOINPUT
End of Input
EOF



#--------------------------------
#   SCF FILE
#--------------------------------
cat > $Project.Scf.in << EOF
&SCF &END
Title
$Project with $Basis1 basis
Occupied
5
Iterations
50
End of Input
EOF



# LOCALIZATION FILE
# --------------------------------
cat > $Project.Localization.in << EOF
&LOCALISATION &END
NORBITALS
5
End of Input
EOF


#   BREWIN_CI FILE
#--------------------------------
cat > $Project.Brewin_ci.in << EOF
&Brewin_CI &End
RESTART FLAG
0
NUMBER OF ORBITALS
35
NUMBER OF BASIS FUNCTIONS
35
CSFS ONLY FLAG
0
SPECIAL CI FLAG
1
REFERENCE CI FLAG
0
LOCAL CI FLAG
1
VIRTUAL TRUNCATION FLAG
1
DIRECT CI FLAG
933
NUMBER OF FROZEN ORBITALS
0
NUMBER OF INACTIVE ORBITALS
5
NUMBER OF ACTIVE ORBITALS
0
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
INTEGRAL BUFFER SIZE
500
ACPF FLAG
1
CD THRESHOLD
1.00000D-12
INTEGRAL THRESHOLD
1.00000E-9
OCCUPATION THRESHOLD
0.8
VIRTUAL OCCUPATION THRESHOLD
0.8
WP DEFAULT RADIUS
0.8 
WP RADIUS MULTIPLIER
1.2
TOV VIRTUAL DEFAULT RADIUS
0.8
TOV VIRTUAL RADIUS MULTIPLIER
1.5
TOV OCCUPIED DEFAULT RADIUS
0.8
TOV OCCUPIED RADIUS MULTIPLIER
1.5
TOV CYLINDER RADIUS
0.8
CYLINDER RADIUS
0.8
NONLOCAL
1
END OF BREWIN VARIABLES
END OF INPUT
EOF


# Start calculating ! 
echo "running seward"
molcas $WorkDir/$Project.Seward.in > $WorkDir/$Project.Seward.out

echo "running SCF"
molcas $WorkDir/$Project.Scf.in > $WorkDir/$Project.Scf.out

echo "running LOVO_prep"
python $WorkDir/LOVO_prep.py "$Basis2" > LOVO_prep.out

echo "running Localization"
molcas $WorkDir/$Project.Localization.in >$WorkDir/$Project.Localization.out
cp ${WorkDir}/*LocOrb PIPEKO

mv *.ChVec1 CHVEC1
mv *.ChRst CHORST
mv *.ChMap CHOMAP
mv *.ChRed CHRED


echo "running Brewin_ci"
molcas $WorkDir/$Project.Brewin_ci.in > $WorkDir/$Project.Brewin_ci.out



# Now we clean up and copy our outputs
cp ${WorkDir}/*out $ThisDir
'rm' -r  $WorkDir
echo "done "


echo ""
echo "Time: " `date`
echo ""



