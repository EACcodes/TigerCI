#! /usr/bin/python
# !
#
#  This script provides the prep work for a LOVO 
#  calculation. It asks MOLCAS to preform a second
#  seward calculation on an enlarged basis containing 
#  both the normal basis set and a secondary minimal   
#  basis. 
#
#  The reference for the LOVOs :
#  J. E. Subotnik, A. D. Dutoi, M. Head-Gordon, JCP, 123, 114108 (2005)
#
#  Important NOTE !!! 
#  While MOLCAS supports a variety of formats this script doesn't ! 
#  I'm basing this off Jeremy's molcas.tcsh script and his input files. 
#  If you prefer to use a different format (such as an xyz file) thats fine
#  ... but you'll need to modify the script , in particular add your own
#  method instead of processInput_Jeremy 
#
#  The script should be used before running tiger_ci ! 
#  The script also requires that the name of the minimal basis be passed as 
#  an argument when the script is called. For example
#
#  python LOVO_prep.py ${minimal_basis_set}  input_file_name > LOVO_prep.out
# 
#  The input file name is optional.
#
#  This script is compatible with python version 2.2 
#  though personally I prefer 2.6 :)   
#
#//////////////////////////////////////////////////////////////////////

import os
import sys
import re 
import struct
import time
import shutil


#************* Constants ********************#

min_block_length = 512
new_molcas_input = "LOVO.input"
molcas_output = "LOVO.molcas.out"
overlap_file = "LOVO_overlap"
basis_file_name = 'LOVO.basis.data' 


#********************************************#


#************* Class Defs *******************#
#   The following classes are useful for     #
#   data storage when reading inputs         #
#   ( yes using classes as glorified         #
#     data structures is a bit overkill)     #
#********************************************#

class System:
    # Holds all the information for our given system
    # atom types , geometry and basis data
        def __init__(self):
            self.atoms = [] 
    
class Atom:
    # Holds all the data on a single atom type
    # geomtery and basis data mostly
    def __init__(self):
        self.name = ''
        self.basis = ''
        self.geometry = []
        self.options = []    
    
#**************Main program *******************#

def main(minimal_basis,in_file_name):

	# Make room for our calculation
	workdir = setupWorkspace()
	
	# Build the new Seward input
	makeMOLCASInput(minimal_basis,in_file_name)
	
	# Run MOLCAS with the new input
	runMOLCAS(workdir)

	# Read the output file and grab the basis set data
	numBasis = readLovoOut() 

	# Now we get the overlap matrix
	readOverlap(numBasis)

	# Finally we can cleap up and return to the CI calculation
	os.system( "cp %s  %s .." % (overlap_file,basis_file_name) )
	os.chdir("..")
	return
	
    
    
    
def makeMOLCASInput(minimal_basis,in_file_name):
	#   This method reads the old Seward input
	#   then constructs a new input file based on
	#   the old one and the minimal basis set

	print("\tWriting the new MOLCAS input file .... ", end=' ') 
	
	#   we "sanitize" all inputs by moving everything into lower case
	minimal_basis = minimal_basis.lower()

	# find and open the seward input
	if  in_file_name == '' :
		sewardSearch = re.compile(r"Seward.in").search
		files = os.listdir(".") 
		seward_in_list = [ i for i in files if sewardSearch(i) ]
		if len(seward_in_list) >1 : err_func(2,"Found multiple seward input files ??? ")
		seward_in = ''.join(seward_in_list) 
	else:
		seward_in = in_file_name
	
	# get the text of the input file
	f = open(seward_in,'r') 

	# Now we need to process the input file and retrieve the geometry and
	# current basis data.... The following method is highly dependent on the
	# input format ( but it is the only part of this method which is input
	# format dependent ! )
	config = System()
	processInput_Jeremy(minimal_basis,f,config)
	f.close()


	# Now that the information on the system is stored in config
	# We need to produce an new output file containing the minimal
	# basis set

	out = open( new_molcas_input ,'w' )
	out.write("&GATEWAY \nTitle\noverlap for LOVO calculation\n")
	index = 1

	# First we print what we already have 
	for atomType in config.atoms:
		out.write("basis set \n")
		out.write( atomType.name + "." + atomType.basis + "....\n")
		print(type(atomType.geometry))
		for line in atomType.geometry :
			   out.write( reindex(line,index) + "\n" )
			   index = index + 1
		for line in atomType.options:
			out.write( line + "\n" )
		out.write("end of basis set \n")


	# Now we print the new minimal basis
	for atomType in config.atoms:
		out.write("basis set \n")
		out.write( atomType.name + "." + minimal_basis + "....\n")
		for line in atomType.geometry :
		    out.write( reindex(line,index) + "\n" )
		    index = index + 1 
		out.write("charge\n0.0\n")
		out.write("end of basis set\n")


	# Now we finish with a few optimizations ... 
	out.write("\n\n&SEWARD\nONEONLY\nMULTIPOLES\n0\n")
	out.close()
	print("Finished new MOLCAS input file")
	return
    
    
def processInput_Jeremy(minimal_basis,f,config):
	#   This method processes an old seward input format
	#   which is written in the "old" molcas input format
	#   which is currently used in Jeremy's script
	#
	#	In case you ever need to implement another function 
	#	for processing a different input format :
	#
	#  INPUTS:
	#	minimal_basis : a string containing the minimal_basis name
	#	f 	      : a file handle for the already opened input file
	#	config        : a system object containing an empty list 
	#
	#  OUTPUTS:
	#	config        : upon exiting config has a list containing an Atom object 
	#			for each atom type in the calculation. The Atom object contains
	#			the current basis set and geometry information. 
	#			
	#
	#//////////////////////////////////////////////////////////////////////////////////

	line = f.readline()
	record_flag = False
	counter = 0


	# The regex searches 
	basisSearch    = re.compile(r'basis' ).search
	endSearch      = re.compile(r'end').search
	atomBasisDef   = re.compile(r'[A-Za-z]+\.').search           # searches for lines like He.6-31G**
	while not line == '':


		# Sanitize the input
		if line[-1] == "\n" : 
			line = (line.rstrip()).lower()
		else:
			line = line.lower()


		if basisSearch(line) and not record_flag :
			# Start of basis record
			record_flag = True
			currentAtom = Atom()
			#print "Start ..." , line
	    
		if endSearch(line) and record_flag :
			# End of basis record
			record_flag = False
			config.atoms.append( currentAtom )
			#print "End ... " , line
	    
		if record_flag and not basisSearch(line) :

			# Ignore this line if it contains a comment or has 0 length(its blank)
			if len(line.split()) == 0 :
				pass
			elif line.split()[0] == "*":
				pass

			else:
                                # The following new bit of logic is designed to 
                                # allow additional options to be passed to MOLCAS inside the basis set definition

				toSplit = line[:]
				bin     = toSplit.split()
				if atomBasisDef(line):
					bin = line.split('.')
					currentAtom.name  = bin[0]
					currentAtom.basis = bin[1]

				elif  line.find(currentAtom.name) >= 0 and len(bin) >= 4 :
					if isFloat(bin[1]) and isFloat(bin[2]) and isFloat(bin[3]):
						# if the line conatins the atomic symbol 
						# and contains after the symbol 3 floating numbers 
						# it looks like H1 1.0 0.0 3.0 
						# and contains the geometry information 
						currentAtom.geometry.append(line)
					
					else:
						# its an extra option 
						currentAtom.options.append(line)	

				else:
					# its an extra option 
					currentAtom.options.append(line)
	
		# Get the next line
		line  = f.readline()
	
	return


def runMOLCAS(workdir) :
	# This function runs MOLCAS with our new input file 
	# First we move ourselves and our input into the LOVO directory

	try:
		os.system("cp %s LOVO" % new_molcas_input )
		os.chdir( os.getcwd() + "/LOVO" )

	except OSError as e : 
		err_func(3, "Having moving files / directory issues : " + e ) 

	# Now we construct a new script to run MOLCAS for us 
	f = open ( 'run.sh','w' ) 

	# Here we a previously set MOLCAS home directory
	MOLCAS = 0
	try :
		MOLCAS = os.environ.get("MOLCAS")
		print("\tUsing the molcas directory: "  + MOLCAS)
	except OSError :
		pass

	# if we have a MOLCAS home directory set we can use that or 
	# we use the default
	if MOLCAS:
		f.write("#! /bin/bash \nexport WorkDir=%s\nmolcas MOLCAS=%s MOLCAS_PRINT=VERBOSE %s  2>&1 >%s \n" \
				% ( workdir +"/LOVO" , MOLCAS , new_molcas_input , molcas_output ))
	else:
		f.write("#! /bin/bash \nexport WorkDir=%s\nmolcas MOLCAS_PRINT=VERBOSE %s  2>&1 >%s \n" \
				% ( workdir +"/LOVO" , new_molcas_input , molcas_output ))		
	f.write("export WorkDir=%s\n" % workdir )                  
	f.close()
	print("\tStarting MOLCAS run ... ", end=' ')
	try :
		os.system("chmod u+x run.sh")
		os.system("./run.sh")
	except Exception as e :
		err_func(6,"Issues when trying to run MOLCAS ... " + str(e) ) 
	print("... Finished MOLCAS run ")
	return 


def readLovoOut() :
	# Reads the output from MOLCAS and grabs the basis set data 
	f = open(molcas_output,'r')
	MolcasScrewedUp = True
	text = [ i.rstrip().lower() for i in f ]
	f.close() 

	# things to search for 
	numBasSearch = re.compile(r'number of basis').search 
	petiteSearch = re.compile(r'petite').search

	# Find the number of basis functions
	line = [ q for q in text if numBasSearch(q) ] 
	numBasis = (line[0]).split()[-1]
	numBasis = int(numBasis)

	# Now find the basis functions
	shellStruct = {}
	flag = 0 
	counter = 0 
	for line in text : 

		if petiteSearch(line) :
			flag = 1 
			continue 

		if flag and counter <  6 : 
			counter = counter + 1 
			continue 

		if flag and counter ==  6 : 
		  # We hit the actual basis set definitions ! 
			bin = line.split()
			if len(bin) != 4 : break
			#BasisSet.append(line)
			#functionList.append(( bin[1],bin[2],bin[3] ))
			label    = bin[1]
			center   = int(bin[3])
			ang      = bin[2][1]
			prinQNum = bin[2][0]
			if prinQNum == '*' : prinQNum = previous  		#MOLCAS is lazy ! 
			prinQNum = int(prinQNum)
			previous = prinQNum

			# Match the angular momentum
			if   ang == 's' : 
				l = 0 
			elif ang == 'p' :
				l = 1 
			elif ang == 'd' : 
				l = 2 
			elif ang == 'f' :
				l = 3 
			elif  ang == 'g' :
				l = 4 
			elif  ang == 'h' :
				l = 5 
			elif  isNum(ang):
				# some fun diffuse or polarized function
				# make a note of it and then move				
				pass			
			else :
				print("What in the world are you doing that you have angular momentum functions of l > 5 in your basis set ?? ") 

			if not center in shellStruct : shellStruct[center] = [] 
			shellStruct[center].append(( label , prinQNum , l  ))


	#print shellStruct
	numCenters = len(shellStruct)

	# Finding basis functions was the easy part ... now we process the basis set
	# we use the shell structure containing ( label , principle q number , ang momentum ) 
	# and we expect the output to be a list with each entry containing  ( num_minimal_basis , nShells , number_orbitals_per_l ) 
	# for each atom  ( in both the normal and minimal basis ) and the maximum length of the tuples

	output , maxLineLength = readNormalBasis(shellStruct)


	# Now we write this all to the output file ... but only the second half of the data 
	# involving the minimal basis set
	basis_file = open( basis_file_name ,'w')
	basis_file.write(str(numBasis) + "       " + str(maxLineLength) +  "    \n" ) 
	counter = 0
	for center in output:
		if counter < numCenters / 2 : 
			counter = counter + 1
			continue
		else:
			basis_file.write(str(center[0]) + " " ) 
			basis_file.write(str(center[1]) + " " )

		for j in center[2]:
			basis_file.write(str(j) + " " )

		basis_file.write("\n")

	basis_file.close()	
	return 	numBasis

def isNum( x ):
	try :
		int(x)
		return True
	except Exception :
		pass
	return False


def isFloat(x):
	try:
		float(x)
		return True
	except ValueError :
		return False
		

def readNormalBasis(shellStruct):

	# We are given the shell structure .... 
	# Now we need to add up all the functions on different centers
	numCenters = len(shellStruct)
	output = [ 0 for i in range(numCenters+1) ] 
	maxLineLength = 0 
	# loop over each center
	for center in shellStruct : 
		num_minimal_basis     = 0 
		nShells               = 0 
		number_orbitals_per_l = []
		# for each basis function on the center
		for basis in shellStruct[center] : 
			( label , prinQNum , l ) = basis
			num_minimal_basis = num_minimal_basis + 1 
			if prinQNum > nShells :
				nShells = prinQNum 
			try :
				number_orbitals_per_l[l] =  number_orbitals_per_l[l] + 1 
			except IndexError : 
				number_orbitals_per_l.append(1)
				assert ( len(number_orbitals_per_l) == (l+1) ) 					# since I'm assuming the l's increase by 1 

		output[center] = ( num_minimal_basis , nShells , number_orbitals_per_l ) 
		if maxLineLength < len(number_orbitals_per_l) + 2 :  maxLineLength = len(number_orbitals_per_l) + 2
	
        # Now we have finished our data processing ! Notice that output[0] = 0 is still there
	# lets ignore it and move on .... 
	output.pop(0)
	return output , maxLineLength


def readOverlap(num_basis) :
	# This method reads the OneInt file to get the overlap matrix

	print("\tStarting to process the one electron integrals ...")
	files = os.listdir(".")
	one_ints_list = [ i for i in files if re.search("OneInt",i) ] 
	if len(one_ints_list) >1 : err_func(7,"Found multiple one integral files ??? ")
	OneInt = ''.join(one_ints_list)
	OneInt = open ( OneInt , 'rb' ) 
	assert(OneInt)

	print("\t\tSearching for MLTPL  0 label")
	go = True
	offset = 0
	while go :
		label = read_bytes(OneInt,offset,'b',1,8)
		label = [ chr(i) for i in label if i <= 256 and i >= 0  ] 
		label = "".join(label)
		if re.search('MLTPL  0',label):
			#print label
			break
		offset = offset + 1 

	print("\t\tFound the MLTPL   0 label")
	print("\t\tCalculating the location of the overlap matrix")

	# move past the label
	offset = offset + 8 
	# The location on disk of the overlap matrix is stored 24 bytes after the label
	offset = offset + 24
	loc = read_bytes(OneInt,offset,'i',4,1)[0]
	loc = loc * min_block_length
	
	print("\t\tNow reading and writing the overlap matrix")
	Output = open( overlap_file , 'wb' )
	OneInt.seek(0)
	OneInt.seek(loc)
        
	buff_size = 100 * 1024 * 1024                               # max buff size of 100 Mb ... reasonable ?
	if buff_size > num_basis * (num_basis+1)/2 * 8 :
                # The buffer size is greater than the entire matrix !
		print("\t\tNo buffering for read/write is necessary ")
		elements = OneInt.read( num_basis * (num_basis+1)//2 * 8)
		Output.write(elements)
                
	else :
                #Time to start buffering the input and output
		print("\t\tBuffering LOVO overlap read/writes ") 
		full_reads = num_basis * (num_basis+1)/2 * 8 / buff_size
		clean_up = num_basis * (num_basis+1)/2 * 8 % buff_size
		for i in full_reads:
			elements = OneInt.read(buff_size)
			Output.write(elements)
		elements = OneInt.read(clean_up)
		Output.write(elements)
	OneInt.close()
	Output.close()
	return 


def read_bytes ( OneInt , offset , typ_code , size , buf_length ) : 
	# Reads from the OneInt file a series of 0's and 1's
	# returns the data as a list 
	#
	#  OneInt	the file object
	#  offset	the offset to begin reading from ( in bytes ) 
	#  typ_code	the code for the type of data I want ( d for double , etc... )  
	#  size 	the size of data to read in number of data types ( i.e. 4 bytes , 4 doubles , etc...)
	#  buf_length	the number of data needed ( i.e. 4 characters , 120 doubles , etc.. )  

	fmt = '<' + typ_code
	OneInt.seek(0)
	OneInt.seek(offset,0)
	result = [ struct.unpack(fmt , OneInt.read(size) ) for i in range(buf_length) ] 
	result = [ i[0] for i in result ]
	return result


		   
    
def reindex ( line , index ) :
	# Given a line with geometry info I need to "reindex" and give the atom a new index
	bin = line.split()
	name = split_name(bin.pop(0))
	everythingElse = [ i + " " for i in bin]
	everythingElse = ''.join(everythingElse)
	return name + str(index) + " " + everythingElse
    
    
def split_name ( name ) : 
	# a method to split the atomic symbol from the index 
	ctest = lambda x : (ord(x) < ord('0') or ord(x) > ord('9'))  		# test for a non-integer character
	name_list = [ c for c in name if ctest(c) ] 
	new_name = ''.join( name_list)
	return new_name

  
def setupWorkspace():

	#   print something nice to the output file
	print("\n\n\n\t\tStarting LOVO prep\n\n")
	print("\tBuilding the workspace for the LOVO calc")
	
	#   We move to the work directory
	workdir = os.environ.get("WorkDir")
	if not workdir : err_func(4 , "The work directory environmental variable was not set!")
	try:
		os.chdir(workdir)
		if  os.path.exists("LOVO"):
			# If the LOVO directory already exists remove it
			print("Warning old LOVO directory detected ... removing and continuing")
			shutil.rmtree("LOVO")
		os.mkdir( "LOVO" ) 	
	except OSError as err:
		err_func(5,"workspace construction error ...:" + str(err))
	#   We have now finished making the LOVO directory in the WorkDir
	return workdir
    
  
  
def err_func( index , message ):
	# a simple error message function
	print("\n\n***************************************")
	print("***************************************")
	print("LOVO Prep error " , index) 
	print(message)
	print("***************************************")
	print("***************************************\n\n")
	os._exit(2)
	return

    
if __name__ == "__main__":
	if len(sys.argv) < 2 : 
		err_func(1,"Basis set information was incorrectly passed to this script")
	a = time.time()
	minimal_basis = sys.argv[1]
	# If the name of the input file was passed use that
	in_file_name = ''
	if len(sys.argv) == 3 : 
		in_file_name = sys.argv[2]
		print("using input file = " + in_file_name)
	main(minimal_basis,in_file_name)
	print("\n\nFinished LOVO_prep")
	print("Total time = %s (sec) " % ( str(time.time() - a) ))


