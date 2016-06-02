/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "OrbitalParsers.h"
#include "StringUtils.h"

using namespace std;
using namespace arma;

static mat MOs;

MoldenOrbitalParser::MoldenOrbitalParser(bool spherical) :
		  nbas(0), norb(0), spherical_harmonics(spherical) {
	if (spherical_harmonics){
		numd = 5;
		numf = 7;
		static const size_t d[5] = {4,2,0,1,3};
		static const size_t f[7] = {6,4,2,0,1,3,5};
		d_permutation = vector<size_t>(d,d+numd);
		f_permutation = vector<size_t>(f,f+numf);
	}
	else{
		numd = 6;
		numf = 10;
		static const size_t d[6] = {0,3,4,1,5,2};
		static const size_t f[10] = {0,4,5,3,9,6,1,8,7,2};
		d_permutation = vector<size_t>(d,d+numd);
		f_permutation = vector<size_t>(f,f+numf);
	}
}

mat get_MOs(){
	return MOs;
}

void set_MOs(mat new_MOs){
	MOs = new_MOs;
}

void MoldenOrbitalParser::parse(std::string filename, BasisData *basis) {
	vector<string> text;
	unordered_map<string, pair<size_t, size_t >> sections;

	cout << endl << "Reading molecular orbitals from " << filename << endl;
	bool angs = false;
	readMoldenFile(filename, text, sections, angs);
	checkCenters(basis, text, sections, angs);
	checkSphericalHarmonics(sections);

	vector<AngMom> funcOrder = getFuncOrder(text, sections);
	vector<size_t> p = buildBasisPermut(funcOrder);
	nbas = p.size();
	checkNBas(basis, nbas);
	readMOs(text, sections, p);
	cout << "  Finished reading in " << norb << " molecular orbitals represented";
	cout << " in " << nbas << " atomic orbitals" << endl;
	cout << endl;
}

void MoldenOrbitalParser::checkNBas(BasisData *basis, const int moldenNBas){
	int erkaleNBas = basis->get_nBas();
	if (!(erkaleNBas==moldenNBas)){
		stringstream error;
		error << "Erkale reports " << erkaleNBas << " basis functions while the orbital file has " ;
		error << moldenNBas << " basis functions!" ;
		throw runtime_error(error.str());
	}
}

void MoldenOrbitalParser::readMOs(vector<string>& text, unordered_map<string,
		pair<size_t,size_t> >& sections, vector<size_t> p){
	if (sections.find("[mo]") == sections.end()) {
		throw runtime_error("No [mo] section in the molden file");
	}

	// find the start of each MO
	norb = 0;
	size_t start = sections["[mo]"].first;
	size_t end = sections["[mo]"].second;
	vector<size_t> orb_starts;
	for(size_t i=start; i<=end; i++){
		string line = trim(text.at(i));
		string label; stringstream ss(line);
		ss >> label;
		if (label.find("occup")==0){
			orb_starts.push_back(i+1);
			norb++;
		}
	}

	MOs = zeros(nbas, norb);
	vector<double> orbital(nbas);
	for(size_t orb=0; orb<norb; orb++){
		size_t start = orb_starts[orb];
		for (size_t j=0; j<nbas; j++){
			stringstream ss(text.at(j+start));
			double coeff;
			int addr;
			ss >> addr >> coeff;
			orbital.at(addr-1) = coeff ; // remember we are base 0!
		}

		permuteVector(orbital, p, 0);
		for (size_t j=0; j<nbas; j++){
			MOs(j,orb) = orbital.at(j);
		}
	}
}

template <class T>
void MoldenOrbitalParser::permuteVector(vector<T>& v, vector<size_t>& p, size_t start){
	vector<T> copy(p.size());

	if (start >= v.size())
		throw invalid_argument("Invalid permutation start location.");

	// make a copy
	for(size_t i=0; i<p.size(); i++){
		copy.at(i) = v.at(start+i);
	}

	for (size_t i=0; i<p.size(); i++){
		v.at(start+i) = copy.at(p.at(i));
	}
}

vector<size_t> MoldenOrbitalParser::buildBasisPermut(vector<AngMom> funcOrder) {
	vector<size_t> p(funcOrder.size());
	for (size_t i = 0; i < p.size(); i++) {
		p[i] = i;
	}

	size_t i = 0 ;
	while(i < p.size()){
		switch (funcOrder[i]){
		case AngMom::s:
			i++;
			break;
		case AngMom::p:
			i += 3 ;
			break;
		case AngMom::d:
			permuteVector(p, d_permutation, i);
			i += numd;
			break;
		case AngMom::f:
			permuteVector(p, f_permutation, i);
			i += numf;
			break;
		default:
			throw runtime_error("Only s, p, d, and f orbitals are supported.");
		}
	}
	cout << "  Built the basis function permutation" << endl;
	return p;
}

vector<AngMom> MoldenOrbitalParser::getFuncOrder(vector<string>& text,
		unordered_map<string, pair<size_t, size_t> >& sections) {
	if (sections.find("[gto]") == sections.end()) {
		throw runtime_error("No [gto] section in the molden file");
	}

	// Here we figure out the order of angular momentum functions when we list all of the basis functions
	// we need this information to reorder the d, f, and g functions
	// We make two assumptions:
	// (1) We are using the cartesian form of the spherical harmonics
	// (2) We have already checked the nuclei and I'm sure the nuclei are in the correct order
	vector<AngMom> funcOrder;
	size_t start = sections["[gto]"].first;
	size_t end = sections["[gto]"].second;
	for (size_t i = start; i <= end; i++) {
		stringstream ss(trim(text[i]));
		char first_char;
		ss >> first_char;
		switch (first_char) {
		case 's':
		funcOrder.push_back(AngMom::s);
		break;

		case 'p':
			for (size_t j = 0; j < 3; j++) {
				funcOrder.push_back(AngMom::p);
			}
			break;

		case 'd':
			for (size_t j = 0; j < numd; j++) {
				funcOrder.push_back(AngMom::d);
			}
			break;

		case 'f':
			for (size_t j = 0; j < numf; j++) {
				funcOrder.push_back(AngMom::f);
			}
			break;

		case 'g':
			throw runtime_error("G functions are not currently supported.");
			break;

		default:
			break;
		}
	}


	return funcOrder;
}

void MoldenOrbitalParser::checkSphericalHarmonics(unordered_map<string, pair<size_t, size_t> >& sections) {
	// In a molden file the sections  [5D], [5D7F] or [9G] indicate spherical harmonics
	// however, [5D10F] or (apparently) [7F] indicate some weird mixture of spherical harmonics
	// for some shells and cartesian functions for other shells
	//
	// We are only supporting pure spherical or pure cartesian orbitals

	string sph[3] = {"[5d]", "[5d7f]", "[9g]"};
	vector<string> spherical(sph,sph+3);
	string mix[2] = {"[5d10f]", "[7f]"};
	vector<string> mixed(mix,mix+2);
	bool sphere_found = false;

	for (string sym : spherical){
		if (sections.find(sym) != sections.end()){
			if (!spherical_harmonics){
				stringstream ss;
				ss << "Input file indicates cartesian orbitals, but orbital ";
				ss << "file contains spherical harmonics.";
				throw runtime_error(ss.str());
			}
			sphere_found = true;
		}
	}


	if (spherical_harmonics && !sphere_found){
		stringstream ss;
		ss << "Input file indicates spherical harmonics, but orbital ";
		ss << "file contains cartesian orbitals.";
		throw runtime_error(ss.str());
	}

	for (string sym : mixed) {
		if (sections.find(sym) != sections.end()) {
			stringstream ss;
			ss << "Found a symbol for mixed spherical and Cartesian orbitals ";
			ss << "in the molden file -> " << sym;
			ss << ". We only support pure Cartesian or spherical functions";
			if (sphere_found){
				cout << "WARNING! " << ss.str() << endl;
				cout << "According to official molden specification, [7F] indicates mixed ";
				cout << "orbitals. However, MOLCAS appears to use it for pure sphericals." << endl;
			}
			else
				throw runtime_error(ss.str());
		}
	}
}

void MoldenOrbitalParser::readMoldenFile(string filename, vector<string>& text,
		unordered_map<string, pair<size_t, size_t >> &sections, bool &angs) {
	ifstream infile(filename);
	if (!infile.is_open()) {
		cout << "Unable to open the orbital file: " << filename << endl;
		throw runtime_error("Failed to open the Molden file");
	}

	string line;
	string curr_section = "";
	size_t i = 0;
	while (getline(infile, line)) {
		line = trim(line);
		transform(line.begin(), line.end(), line.begin(), ::tolower); // everything in lowercase
		text.push_back(line);
		if (line[0] == '[') {
			if (line.find("[atoms]") != string::npos)
				angs = ((line.find("(au)") != string::npos) ? false : true);

			// MOLCAS adds text after the [] ... remove it
			line.erase(find(line.begin(), line.end(), ']') + 1, line.end());

			if (!(curr_section == "")) {
				sections[curr_section].second = i - 1;
			}
			curr_section = line;
			sections[curr_section] = pair<size_t, size_t>(i + 1, 0);
		}
		i++;
	}
	// set the last section to end at the end of the text
	sections[curr_section].second = text.size() - 1;

	cout << "  Molden file contains sections: ";
	for (auto s : sections) {
		cout << s.first << " ";
	}
	cout << endl;
}

void MoldenOrbitalParser::writeMoldenFile(std::string filename){
	cout << "writing file!" << endl;
	ifstream infile(filename);
	if (!infile.is_open()) {
		cout << "Unable to open the orbital file: " << filename << endl;
		throw runtime_error("Failed to open the Molden file");
	}
        
        
        string file = getFileNameFromPath(filename);
        
	string filestem = file.substr(0,file.find(".molden"));
	ofstream outfile(filestem + ".local.molden");
	if (!outfile.is_open()) {
		cout << "Unable to open the orbital output file: " << file + ".localized" << endl;
		throw runtime_error("Failed to open the Molden output file");
	}

	string line;
	//copy headers and stuff over to new output file
	while (getline(infile, line)){
		outfile << line << endl;
		line = trim(line);
		transform(line.begin(), line.end(), line.begin(), ::tolower);
		if (line.find("[mo]") != string::npos)
			break;
	}

	//add MOs to output file
	int i = -1;
	while (getline(infile, line)){
		string trimmed(line);
		trimmed = trim(trimmed);
		transform(trimmed.begin(), trimmed.end(), trimmed.begin(), ::tolower);

		//check if first character is numeric
		if (trimmed[0] < '0' || trimmed[0] > '9'){
			if (trimmed.find("sym") != string::npos)
				i++; //increment MO
			outfile << line << endl;
			continue;
		}

		size_t pos = trimmed.find(' ');
#ifdef __FreeBSD__
                int j;
                stringstream ss(trimmed.substr(0,pos));
                ss >> j;
                j--;
#else
		int j = stoi(trimmed.substr(0,pos)) - 1;
#endif
		outfile << trimmed.substr(0, pos) + " " << MOs(j,i) << endl;
	}

}

bool MoldenOrbitalParser::is_int(std::string & s){
    
    try{
        const size_t trial = std::stoi(s);
    } catch (std::exception& e){
        return false;
    }
    
    return true;
}


vector<atom_t> MoldenOrbitalParser::readCoordinates(string filename){
    
    std::vector<string> text;
    std::unordered_map<string, pair<size_t, size_t >> sections;
    bool angs = false;
    readMoldenFile(filename, text, sections, angs);
    
    vector<atom_t> atoms;
    
    double conv = (angs) ? ANGSTROMINBOHR : 1.0;
    
    size_t start = sections["[atoms]"].first;
    size_t end = sections["[atoms]"].second;
    for (size_t i = start; i <= end; i++) {
    	string line = trim(text[i]);
	if (line == "") continue;
	stringstream ss(line);
	string stmp;
	size_t itmp, num;
	double x, y, z;
	ss >> stmp >> itmp >> num >> x >> y >> z;
        
        // remove any numbers directly attached to the element symbol
        string lastChar = stmp.substr(stmp.length()-1,stmp.length());
        while(is_int(lastChar)){
            stmp = stmp.substr(0,stmp.length()-1); // strip off number
            lastChar = stmp.substr(stmp.length()-1,stmp.length()); // save last char
        }
        
	atom_t atom;
        atom.el = stmp;
        atom.num = num;
	atom.x = x * conv;
	atom.y = y * conv;
	atom.z = z * conv;
        atom.Q = 0.0;
        
	atoms.push_back(atom);
    }
    
    return atoms;
}

void MoldenOrbitalParser::checkCenters(BasisData *basis, vector<string>& text,
		unordered_map<string, pair<size_t, size_t> >& sections, bool angs) {
	if (sections.find("[atoms]") == sections.end()) {
		throw runtime_error("No [atoms] section in the molden file");
	}

	vector<tiger::atom_t> molden_centers;
	size_t start = sections["[atoms]"].first;
	size_t end = sections["[atoms]"].second;
	for (size_t i = start; i <= end; i++) {
		string line = trim(text[i]);
		if (line == "") continue;
		stringstream ss(line);
		string stmp;
		size_t itmp, num;
		double x, y, z;
		ss >> stmp >> itmp >> num >> x >> y >> z;
		tiger::atom_t atom;
		atom.Z = num;
		atom.coords.x = x;
		atom.coords.y = y;
		atom.coords.z = z;
		molden_centers.push_back(atom);
	}

	size_t erkaleNnuc = basis->get_nAtoms();
	size_t myNnuc = molden_centers.size();
	if (myNnuc != erkaleNnuc) {
		stringstream ss;
		ss << "The number of atoms in the orbital file and the xyz file  do not match!";
		ss << myNnuc << " " << erkaleNnuc;
		throw runtime_error(ss.str());
	}

	double bohr_to_angs = (angs ? 1.8897261 : 1.0);
	double thresh = 1.0E-3;
	vector<tiger::atom_t> basis_centers = basis->get_atoms();
	for (size_t i = 0; i < myNnuc; i++) {
		int erZ = basis_centers.at(i).Z;
		int moZ = molden_centers.at(i).Z;

		double erkale_x = basis_centers.at(i).coords.x / bohr_to_angs;
		double erkale_y = basis_centers.at(i).coords.y / bohr_to_angs;
		double erkale_z = basis_centers.at(i).coords.z / bohr_to_angs;
		double molden_x = molden_centers.at(i).coords.x;
		double molden_y = molden_centers.at(i).coords.y;
		double molden_z = molden_centers.at(i).coords.z;

		double diff_x = abs(erkale_x - molden_x);
		double diff_y = abs(erkale_y - molden_y);
		double diff_z = abs(erkale_z - molden_z);

		if ((erZ != moZ) || (diff_x > thresh) || (diff_y > thresh) || (diff_z > thresh))
			throw runtime_error("The nuclei are ordered differently in the xyz and orbital file");
	}
	cout << "  Nuclei are consistent" << endl;
}

// Double checks that all molecular orbitals are orthonormal
void MoldenOrbitalParser::check_norms(BasisData *basis){
	int nOrb = MOs.n_cols;
	mat overlap = basis->get_overlap();
	bool error = false;
	for (int i = 0; i < nOrb; i++){
		double norm = calc_norm(MOs.col(i), MOs.col(i), overlap);
		if (abs(norm - 1.0) > 5.0E-6){
			cout << "Orbital  " << i << " is not normalized, with a norm of: " << norm << endl;
			error = true;
		}
		for (int j = i + 1; j < nOrb; j++){
			double ortho = calc_norm(MOs.col(i), MOs.col(j), overlap);
			if (abs(ortho) > 5.0E-6){
				cout << "Orbitals " << i << " and " << j << " are not orthogonal, with overlap: " << ortho << endl;
				error = true;
			}
		}
	}

	if (error){
		stringstream err;
		err << "******************************************" << endl;
		err << "                                          " << endl;
		err << "SERIOUS WARNING FROM MOLECULAR_ORBITAL_MOD" << endl;
		err << "NORMALIZATION ISSUES !!!! " << endl;
		err << "                                          " << endl;
		err << "******************************************" << endl;
		throw runtime_error(err.str());
	}
	else
		cout << "Molecular orbitals are normalized and orthogonal!" << endl;
}

// Calculates the norm of a molecular orbital (AO basis is not orthogonal)
double MoldenOrbitalParser::calc_norm(vec vec1, vec vec2, mat overlap){
	return dot(vec1, (overlap * vec2));
}

// get MOs
void c_get_coefficients(int64_t& nBas, int64_t& nOrb, double *molecular_orbitals){    
	if ((MOs.n_rows != nBas) || (MOs.n_cols != nOrb))
		throw runtime_error("Number of molecular orbitals does not match " +
				string("number of basis functions and number of orbitals"));

	double * ptr = molecular_orbitals;

	//convert integral matrix to two-dimensional array
	for (int64_t i = 0; i < nOrb; i++){
		for (int64_t j = 0; j < nBas; j++){
			*ptr = MOs(j,i); //[j][i] for Fortran arrays
			ptr++;
		}
	}
}
