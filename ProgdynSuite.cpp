#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

#define avNum 6.0221415E23
#define RgasK 0.00198588
#define RgasJ 8.31447
#define conver1 4.184E26 //dividing by this converts amu angs^2 /s^2 to kcal/mol
#define constant_e 2.71828182846
#define constant_pi 3.14159265358979
#define c 29979245800
#define h 6.626075E-34

using namespace std;

template <typename T> string to_string(const T& n) {
	ostringstream stm;
	stm << n;
	return stm.str();
}

float Distance(unsigned int Atom1, unsigned int Atom2, vector<float>& A, vector<float>& B, vector<float>& C) {
	Atom1--; // Goes from a system of indices starting at 1 to a system starting at 0
	Atom2--; // Same thing
	if((A.size() < Atom1) && (A.size() < Atom2)) return 0.0; // prevents bounds errors
	else if(A.size() < Atom1) return sqrt((0-A[Atom2])*(0-A[Atom2]) + (0-B[Atom2])*(0-B[Atom2]) + (0-C[Atom2])*(0-C[Atom2]));
	else if(A.size() < Atom2) return sqrt((A[Atom1]-0)*(A[Atom1]-0) + (B[Atom1]-0)*(B[Atom1]-0) + (C[Atom1]-0)*(C[Atom1]-0));
	else return sqrt((A[Atom1]-A[Atom2])*(A[Atom1]-A[Atom2]) + (B[Atom1]-B[Atom2])*(B[Atom1]-B[Atom2]) + (C[Atom1]-C[Atom2])*(C[Atom1]-C[Atom2]));
}

class Atom {
private:
	float x;
	float y;
	float z;
	float orig_x_pos;
	float orig_y_pos;
	float orig_z_pos;
	float old_x;
	float old_y;
	float old_z;
	float older_x;
	float older_y;
	float older_z;
	float atomic_weight;
	float atomic_num;
	string atomic_symbol;
	float x_velocity;
	float y_velocity;
	float z_velocity;
	float new_x_force;
	float new_y_force;
	float new_z_force;
	float x_force;
	float y_force;
	float z_force;
	float old_x_force;
	float old_y_force;
	float old_z_force;
	float x_after_vel;
	float y_after_vel;
	float z_after_vel;
public:
	Atom(float new_x, float new_y, float new_z, float new_atomic_weight, string new_atomic_symbol) {
		x = new_x;
		y = new_y;
		z = new_z;
		atomic_weight = new_atomic_weight;
		atomic_symbol = new_atomic_symbol;
		return;
	}
	Atom(int new_atomic_num, float new_x, float new_y, float new_z) {
		atomic_num = new_atomic_num;
		if(atomic_num==1) {atomic_symbol="H";atomic_weight=1.00783;}
		else if(atomic_num==2) {atomic_symbol="He";atomic_weight=4.0026;}
		else if(atomic_num==3) {atomic_symbol="Li";atomic_weight=6.941;}
		else if(atomic_num==4) {atomic_symbol="Be";atomic_weight=9.012;}
		else if(atomic_num==5) {atomic_symbol="B";atomic_weight=10.811;}
		else if(atomic_num==6) {atomic_symbol="C";atomic_weight=12.;}
		else if(atomic_num==7) {atomic_symbol="N";atomic_weight=14.007;}
		else if(atomic_num==8) {atomic_symbol="O";atomic_weight=15.9994;}
		else if(atomic_num==9) {atomic_symbol="F";atomic_weight=18.9984;}
		else if(atomic_num==10) {atomic_symbol="Ne";atomic_weight=20.1797;}
		else if(atomic_num==11) {atomic_symbol="Na";atomic_weight=22.989;}
		else if(atomic_num==12) {atomic_symbol="Mg";atomic_weight=24.305;}
		else if(atomic_num==13) {atomic_symbol="Al";atomic_weight=26.98154;}
		else if(atomic_num==14) {atomic_symbol="Si";atomic_weight=28.0855;}
		else if(atomic_num==15) {atomic_symbol="P";atomic_weight=30.9738;}
		else if(atomic_num==16) {atomic_symbol="S";atomic_weight=32.066;}
		else if(atomic_num==17) {atomic_symbol="Cl";atomic_weight=35.4527;}
		else if(atomic_num==18) {atomic_symbol="Ar";atomic_weight=39.948;}
		else if(atomic_num==19) {atomic_symbol="K";atomic_weight=39.0983;}
		else if(atomic_num==20) {atomic_symbol="Ca";atomic_weight=40.078;}
		else if(atomic_num==21) {atomic_symbol="Sc";atomic_weight=44.96;}
		else if(atomic_num==22) {atomic_symbol="Ti";atomic_weight=47.867;}
		else if(atomic_num==23) {atomic_symbol="V";atomic_weight=50.94;}
		else if(atomic_num==24) {atomic_symbol="Cr";atomic_weight=51.9961;}
		else if(atomic_num==25) {atomic_symbol="Mn";atomic_weight=54.938;}
		else if(atomic_num==26) {atomic_symbol="Fe";atomic_weight=55.845;}
		else if(atomic_num==27) {atomic_symbol="Co";atomic_weight=58.933;}
		else if(atomic_num==28) {atomic_symbol="Ni";atomic_weight=58.693;}
		else if(atomic_num==29) {atomic_symbol="Cu";atomic_weight=63.546;}
		else if(atomic_num==30) {atomic_symbol="Zn";atomic_weight=65.38;}
		else if(atomic_num==31) {atomic_symbol="Ga";atomic_weight=69.723;}
		else if(atomic_num==32) {atomic_symbol="Ge";atomic_weight=72.64;}
		else if(atomic_num==33) {atomic_symbol="As";atomic_weight=74.9216;}
		else if(atomic_num==34) {atomic_symbol="Se";atomic_weight=78.96;}
		else if(atomic_num==35) {atomic_symbol="Br";atomic_weight=79.904;}
		else if(atomic_num==46) {atomic_symbol="Pd";atomic_weight=106.42;}
		else if(atomic_num==53) {atomic_symbol="I";atomic_weight=126.90447;}
		x = new_x;
		y = new_y;
		z = new_z;
		orig_x_pos = new_x;
		orig_y_pos = new_y;
		orig_z_pos = new_z;
		x_velocity = 0;
		y_velocity = 0;
		z_velocity = 0;
	}
	void set_x(float new_x) { x = new_x; }
	void set_y(float new_y) { y = new_y; }
	void set_z(float new_z) { z = new_z; }
	void set_atomic_weight(float new_atomic_weight) { atomic_weight = new_atomic_weight; }
	void set_old_x(float new_old_x) { old_x = new_old_x; }
	void set_old_y(float new_old_y) { old_y = new_old_y; }
	void set_old_z(float new_old_z) { old_z = new_old_z; }
	void set_orig_x(float new_orig_x_pos) { orig_x_pos = new_orig_x_pos; }
	void set_orig_y(float new_orig_y_pos) { orig_y_pos = new_orig_y_pos; }
	void set_orig_z(float new_orig_z_pos) { orig_z_pos = new_orig_z_pos; }
	void set_older_x(float new_older_x) { older_x = new_older_x; }
	void set_older_y(float new_older_y) { older_y = new_older_y; }
	void set_older_z(float new_older_z) { older_z = new_older_z; }
	void set_vel(float new_x_vel, float new_y_vel, float new_z_vel) {
		x_velocity = new_x_vel;
		y_velocity = new_y_vel;
		z_velocity = new_z_vel;
		return;
	}
	void set_x_vel(float new_x_vel) { x_velocity = new_x_vel; }
	void set_y_vel(float new_y_vel) { y_velocity = new_y_vel; }
	void set_z_vel(float new_z_vel) { z_velocity = new_z_vel; }
	void set_x_after_vel(float new_x_after_vel) { x_after_vel = new_x_after_vel; }
	void set_y_after_vel(float new_y_after_vel) { y_after_vel = new_y_after_vel; }
	void set_z_after_vel(float new_z_after_vel) { z_after_vel = new_z_after_vel; }
	void set_new_x_force(float new_new_x_force) { new_x_force = new_new_x_force; }
	void set_new_y_force(float new_new_y_force) { new_y_force = new_new_y_force; }
	void set_new_z_force(float new_new_z_force) { new_z_force = new_new_z_force; }
	void set_x_force(float new_x_force) { x_force = new_x_force; }
	void set_y_force(float new_y_force) { y_force = new_y_force; }
	void set_z_force(float new_z_force) { z_force = new_z_force; }
	void set_old_x_force(float new_old_x_force) { old_x_force = new_old_x_force; }
	void set_old_y_force(float new_old_y_force) { old_y_force = new_old_y_force; }
	void set_old_z_force(float new_old_z_force) { old_z_force = new_old_z_force; }
	float get_x() const { return x; }
	float get_y() const { return y; }
	float get_z() const { return z; }
	float get_orig_x() const { return orig_x_pos; }
	float get_orig_y() const { return orig_y_pos; }
	float get_orig_z() const { return orig_z_pos; }
	float get_old_x() const { return old_x; }
	float get_old_y() const { return old_y; }
	float get_old_z() const { return old_z; }
	float get_older_x() const { return older_x; }
	float get_older_y() const { return older_y; }
	float get_older_z() const { return older_z; }
	float get_atomic_weight() const { return atomic_weight; }
	float get_atomic_num() const { return atomic_num; }
	string get_atomic_symbol() const { return atomic_symbol; }
	float get_x_vel() const { return x_velocity; }
	float get_y_vel() const { return y_velocity; }
	float get_z_vel() const { return z_velocity; }
	float get_x_after_vel() const { return x_after_vel; }
	float get_y_after_vel() const { return y_after_vel; }
	float get_z_after_vel() const { return z_after_vel; }
	float get_new_x_force() const { return new_x_force; }
	float get_new_y_force() const { return new_y_force; }
	float get_new_z_force() const { return new_z_force; }
	float get_x_force() const { return x_force; }
	float get_y_force() const { return y_force; }
	float get_z_force() const { return z_force; }
	float get_old_x_force() const { return old_x_force; }
	float get_old_y_force() const { return old_y_force; }
	float get_old_z_force() const { return old_z_force; }
};

class ProgdynConf {
public:
	void set_method(string new_method) { method = new_method; }
	void set_meth2(string new_meth2) { meth2 = new_meth2; }
	void set_meth3(string new_meth3) { meth3 = new_meth3; }
	void set_meth4(string new_meth4) { meth4 = new_meth4; }
	void set_meth5(string new_meth5) { meth5 = new_meth5; }
	void set_meth6(string new_meth6) { meth6 = new_meth6; }
	void set_meth7(string new_meth7) { meth7 = new_meth7; }
	void set_charge(string new_charge) { charge = new_charge; }
	void set_multiplicity(string new_multiplicity) { multiplicity = new_multiplicity; }
	void set_memory(string new_memory) { memory = new_memory; }
	void set_processors(unsigned int new_processors) { processors = new_processors; }
	void set_checkpoint(string new_checkpoint) { checkpoint = new_checkpoint; }
	void set_initialDis(unsigned int new_initialDis) { initialDis = new_initialDis; }
	void set_diag(unsigned int new_diag) { diag = new_diag; }
	void set_timestep(float new_timestep) { timestep = new_timestep; }
	void set_scaling(float new_scaling) { scaling = new_scaling; }
	void set_temp(float new_temp) { temp = new_temp; }
	void set_searchdir(string new_searchdir) { searchdir = new_searchdir; }
	void set_classical(unsigned int new_classical) { classical = new_classical; }
	void set_numimag(unsigned int new_numimag) { numimag = new_numimag; }
	void set_geometry(string new_geometry) { geometry = new_geometry; }
	void set_highlevel(unsigned int new_highlevel) { highlevel = new_highlevel; }
	void set_boxon(bool new_boxon) { boxon = new_boxon; }
	void set_boxsize(float new_boxsize) { boxsize = new_boxsize; }
	void set_maxAtomMove(float new_maxAtomMove) { maxAtomMove = new_maxAtomMove; }
	void set_cannonball(float new_cannonball) { cannonball = new_cannonball; }
	void set_disMode(vector<int>& new_disMode) { disMode = new_disMode; }
	void set_controlPhase(map<int,string>& new_controlPhase) { controlPhase = new_controlPhase; }
	void set_rotationmode(unsigned int new_rotationmode) { rotationmode = new_rotationmode; }
	void set_linkatoms(unsigned int new_linkatoms) { linkatoms = new_linkatoms; }
	void set_fixedatom1(int new_fixedatom1) { fixedatom1 = new_fixedatom1; }
	void set_fixedatom2(int new_fixedatom2) { fixedatom2 = new_fixedatom2; }
	void set_fixedatom3(int new_fixedatom3) { fixedatom3 = new_fixedatom3; }
	void set_fixedatom4(int new_fixedatom4) { fixedatom4 = new_fixedatom4; }
	void set_DRP(bool new_DRP) { DRP = new_DRP; }
	void set_methodfilelines(unsigned int new_methodfilelines) { methodfilelines = new_methodfilelines; }
	void set_killcheck(bool new_killcheck) { killcheck = new_killcheck; }
	void set_damping(float new_damping) { damping = new_damping; }
	void set_etolerance(unsigned int new_etolerance) { etolerance = new_etolerance; }
	void set_nmrmethod(string new_nmrmethod) { nmrmethod = new_nmrmethod; }
	void set_nmrmethod2(string new_nmrmethod2) { nmrmethod2 = new_nmrmethod2; }
	void set_nmrmethod3(string new_nmrmethod3) { nmrmethod3 = new_nmrmethod3; }
	void set_nmrtype(unsigned int new_nmrtype) { nmrtype = new_nmrtype; }
	void set_nmrevery(unsigned int new_nmrevery) { nmrevery = new_nmrevery; }
	void set_nonstandard(bool new_nonstandard) { nonstandard = new_nonstandard; }
	void set_title1(string new_title1) { title1 = new_title1; }
	void set_title2(string new_title2) { title2 = new_title2; }
	void set_title3(string new_title3) { title3 = new_title3; }
	void set_title4(string new_title4) { title4 = new_title4; }
	void set_reversetraj(bool new_reversetraj) { reversetraj = new_reversetraj; }

	string get_method() const { return method; }
	string get_meth2() const { return meth2; }
	string get_meth3() const { return meth3; }
	string get_meth4() const { return meth4; }
	string get_meth5() const { return meth5; }
	string get_meth6() const { return meth6; }
	string get_meth7() const { return meth7; }
	string get_charge() const { return charge; }
	string get_multiplicity() const { return multiplicity; }
	string get_memory() const { return memory; }
	unsigned int get_processors() const { return processors; }
	string get_checkpoint() const { return checkpoint; }
	unsigned int get_initialDis() const { return initialDis; }
	unsigned int get_diag() const { return diag; }
	float get_timestep() const { return timestep; }
	float get_scaling() const { return scaling; }
	float get_temp() const { return temp; }
	string get_searchdir() const { return searchdir; }
	unsigned int get_classical() const { return classical; }
	unsigned int get_numimag() const { return numimag; }
	string get_geometry() const { return geometry; }
	unsigned int get_highlevel() const { return highlevel; }
	bool get_boxon() const { return boxon; }
	float get_boxsize() const { return boxsize; }
	float get_maxAtomMove() const { return maxAtomMove; }
	float get_cannonball() const { return cannonball; }
	vector<int> get_disMode() const { return disMode; }
	map<int,string> get_controlPhase() const { return controlPhase; }
	unsigned int get_rotationmode() const { return rotationmode; }
	unsigned int get_linkatoms() const { return linkatoms; }
	int get_fixedatom1() const { return fixedatom1; }
	int get_fixedatom2() const { return fixedatom2; }
	int get_fixedatom3() const { return fixedatom3; }
	int get_fixedatom4() const { return fixedatom4; }
	bool get_DRP() const { return DRP; }
	unsigned int get_methodfilelines() const { return methodfilelines; }
	bool get_killcheck() const { return killcheck; }
	float get_damping() const { return damping; }
	unsigned int get_etolerance() const { return etolerance; }
	string get_nmrmethod() const { return nmrmethod; }
	string get_nmrmethod2() const { return nmrmethod2; }
	string get_nmrmethod3() const { return nmrmethod3; }
	unsigned int get_nmrtype() const { return nmrtype; }
	unsigned int get_nmrevery() const { return nmrevery; }
	bool get_nonstandard() const { return nonstandard; }
	string get_title1() const { return title1; }
	string get_title2() const { return title2; }
	string get_title3() const { return title3; }
	string get_title4() const { return title4; }
	bool get_reversetraj() const { return reversetraj; }

	void read_from_file() {
		ifstream progdyn_conf;
		progdyn_conf.open("progdyn.conf");
		string line;
		while(getline(progdyn_conf, line)) {
			if(line == "") break;
			if(line[0] != '#') {
				stringstream ss;
				ss << line;
				string first_word;
				ss >> first_word;
				if(first_word=="method") ss >> method;
				else if(first_word=="method2") ss >> meth2;
				else if(first_word=="method3") ss >> meth3;
				else if(first_word=="method4") ss >> meth4;
				else if(first_word=="method5") ss >> meth5;
				else if(first_word=="method6") ss >> meth6;
				else if(first_word=="method7") ss >> meth7;
				else if(first_word=="charge") ss >> charge;
				else if(first_word=="multiplicity") ss >> multiplicity;
				else if(first_word=="memory") ss >> memory;
				else if(first_word=="processors") ss >> processors;
				else if(first_word=="checkpoint") ss >> checkpoint;
				else if(first_word=="initialdis") ss >> initialDis;
				else if(first_word=="diagnostics") ss >> diag;
				else if(first_word=="timestep") ss >> timestep;
				else if(first_word=="scaling") ss >> scaling;
				else if(first_word=="temperature") ss >> temp;
				else if(first_word=="searchdir") ss >> searchdir;
				else if(first_word=="classical") ss >> classical;
				else if(first_word=="numimag") ss >> numimag;
				else if(first_word=="geometry") ss >> geometry;
				else if(first_word=="highlevel") ss >> highlevel;
				else if(first_word=="boxon") ss >> boxon;
				else if(first_word=="boxsize") ss >> boxsize;
				else if(first_word=="maxAtomMove") ss >> maxAtomMove;
				else if(first_word=="cannonball") ss >> cannonball;
				else if(first_word=="displacements") {
					int pos_disMode;
					float disMode_value;
					ss >> pos_disMode;
					ss >> disMode_value;
					disMode.at(pos_disMode) = disMode_value;
				}
				else if(first_word=="controlphase") {
					int pos_controlPhase;
					string pos_or_neg;
					ss >> pos_controlPhase;
					ss >> pos_or_neg;
					controlPhase.insert(pair<int,string>(pos_controlPhase,pos_or_neg));
				}
				else if(first_word=="rotationmode") ss >> rotationmode;
				else if(first_word=="linkatoms") ss >> linkatoms;
				else if(first_word=="fixedatom1") ss >> fixedatom1;
				else if(first_word=="fixedatom2") ss >> fixedatom2;
				else if(first_word=="fixedatom3") ss >> fixedatom3;
				else if(first_word=="fixedatom4") ss >> fixedatom4;
				else if(first_word=="DRP") {
					ss >> DRP;
					if(DRP==1) classical = 2;
				}
				else if(first_word=="methodfile") ss >> methodfilelines;
				else if(first_word=="killcheck") ss >> killcheck;
				else if(first_word=="damping") ss >> damping;
				else if(first_word=="etolerance") ss >> etolerance;
				else if(first_word=="NMRmethod") ss >> nmrmethod;
				else if(first_word=="NMRmethod2") ss >> nmrmethod2;		
				else if(first_word=="NMRmethod3") ss >> nmrmethod2;
				else if(first_word=="NMRtype") ss >> nmrtype;
				else if(first_word=="NMRevery") ss >> nmrevery;
				else if(first_word=="nonstandard") ss >> nonstandard;
				else if(first_word=="title") {
					ss >> title1;
					ss >> title2;
					ss >> title3;
					ss >> title4;
				} 
				else if(first_word=="reversetraj") {
					string true_or_false;
					ss >> true_or_false;
					if(true_or_false == "true") reversetraj = true;
					else reversetraj = false;
				}// end if
			} // end if
		} // end while
		return;
	} // end function
	
	ProgdynConf() {
		method = "";
		meth2 = "";
		meth3 = "";
		meth4 = "";
		meth5 = "";
		meth6 = "";
		meth7 = "";
		charge = "0";
		multiplicity = "1";
		memory = "4gb";
		processors = 2;
		checkpoint = "g09.chk";
		initialDis = 2;
		diag = 0;
		timestep = 1E-15;
		scaling = 1.0;
		temp = 363.15;
		searchdir = "positive";
		classical = 1;
		numimag = 1;
		geometry = "linear";
		highlevel = 999;
		boxon = false;
		boxsize = 7.5;
		maxAtomMove = 0.1;
		cannonball = 0;
		rotationmode = 1;
		linkatoms = 0;
		fixedatom1 = -1;
		fixedatom2 = -1;
		fixedatom3 = -1;
		fixedatom4 = -1;
		DRP = false;
		methodfilelines = 0;
		killcheck = true;
		damping = 1.0;
		etolerance = 2;
		nmrmethod = "";
		nmrmethod2 = "";
		nmrmethod3 = "";
		nmrtype = 0;
		nmrevery = 9999999;
		nonstandard = false;
		title1 = "alame3";
		title2 = "XC";
		title3 = "UB3D2PCM";
		title4 = "363dis2";
		reversetraj = false;
	}
private:
	string method, meth2, meth3, meth4, meth5, meth6, meth7;
	string charge;
	string multiplicity;
	string memory;
	unsigned int processors;
	string checkpoint;
	unsigned int initialDis;
	unsigned int diag;
	float timestep;
	float scaling;
	float temp;
	string searchdir;
	unsigned int classical;
	unsigned int numimag;
	string geometry;
	unsigned int highlevel;
	bool boxon;
	float boxsize;
	float maxAtomMove;
	float cannonball;
	vector<int> disMode;
	map<int,string> controlPhase;
	unsigned int rotationmode;
	unsigned int linkatoms;
	int fixedatom1, fixedatom2, fixedatom3, fixedatom4;
	bool DRP;
	unsigned int methodfilelines;
	bool killcheck;
	float damping;
	unsigned int etolerance;
	string nmrmethod, nmrmethod2, nmrmethod3;
	unsigned int nmrtype;
	unsigned int nmrevery;
	bool nonstandard;
	string title1, title2, title3, title4;
	bool reversetraj;
};

template <class T1,class T2,class T3>
class Tuple {
private:
	T1 first_element;
	T2 second_element;
	T3 third_element;
public:
	inline bool operator < (const Tuple& rhs) const { 
		if(this->first_element < rhs.get_first()) return true;
		else if(this->first_element > rhs.get_first()) return false;
		if(this->second_element < rhs.get_second()) return true;
		else if(this->second_element > rhs.get_second()) return false;
		if(this->third_element < rhs.get_third()) return true;
		else return false;
	}
	T1 get_first() const { return first_element; }
	T2 get_second() const { return second_element; }
	T3 get_third() const { return third_element; }
	Tuple(T1 f_e, T2 s_e, T3 t_e) {
		first_element = f_e;
		second_element = s_e;
		third_element = t_e;
	}
};

void cat_append(string source_str, string destination_str) {
	ifstream source_file(source_str.c_str());
	ofstream destination_file;
	destination_file.open(destination_str.c_str(), fstream::app);
	string line;
    while(getline(source_file,line)) {
		destination_file << line << endl;
	}
	destination_file.close();
	return;
}

void cp(string source_str, string destination_str) {
	ifstream source_file(source_str.c_str());
	ofstream destination_file;
	destination_file.open(destination_str.c_str());
	string line;
	while(getline(source_file,line)) {
		destination_file << line << endl;
	}
	destination_file.close();
	return;
}

void mv(string source_str, string destination_str) {
	ifstream source_file(source_str.c_str());
	ofstream destination_file;
	destination_file.open(destination_str.c_str());
	string line;
	while(getline(source_file,line)) {
		destination_file << line << endl;
	}
	destination_file.close();
	remove(source_str.c_str());
	return;
}

void append_to_dynfollowfile(string line) {
	ofstream dynfollowfile;
	dynfollowfile.open("dynfollowfile", ios::app);
	dynfollowfile << line << endl;
	return;
}

void append_to_geoRecord(string line) {
	ofstream geoRecord;
	geoRecord.open("geoRecord", ios::app);
	geoRecord << line << endl;
	return;
}

void print_date_to_file(string filename) {
	ofstream file;
	file.open(filename.c_str());
	time_t rawtime;
	time(&rawtime);
	file << ctime(&rawtime);
	file.close();
	return;
}

string print_date_to_file() {
	stringstream s_stream;
	time_t rawtime;
	time(&rawtime);
	s_stream << ctime(&rawtime);
	return s_stream.str();
}

void print_to_NMRlistcc(string line) {
	ofstream NMRlistcc;
	NMRlistcc.open("NMRlistcc", ios::app);
	NMRlistcc << line << endl;
	return;
}

bool grep_Normal_termination(string file_path) {
	ifstream log_file(file_path.c_str());
	string line;
	while(getline(log_file,line)) {
		if(line.find("Normal termination") != string::npos) return true;
	}
	return false;
}

string str_find_XXXX_last_n_lines(string filename, unsigned int n) {
	ifstream file(filename.c_str());
	string lines[n];
	unsigned int size = 0;
	string line;
	while(getline(file,line)) lines[size++%n] = line;
	unsigned int start = size > n ? (size%n) : 0;
	unsigned int count = min(n,size);
	for(unsigned int i = 0; i < count; i++) {
		if(line.find("XXXX") != string::npos) return line; // cout << l[(start+i)%k] << endl ; // FIX THIS!!!
	}
	return false;
}

bool bool_find_XXXX_last_n_lines(string filename, unsigned int n) {
	ifstream file(filename.c_str());
	string lines[n];
	unsigned int size = 0;
	string line;
	while(getline(file,line)) lines[size++%n] = line;
	unsigned int start = size > n ? (size%n) : 0;
	unsigned int count = min(n,size);
	for(unsigned int i = 0; i < count; i++) {
		if(line.find("XXXX") != string::npos) return true; // cout << l[(start+i)%k] << endl ; // FIX THIS!!!
	}
	return false;
}

bool file_exists(string filename) {
	ifstream file(filename.c_str());
	if(file.good()) return true;
	else false;
}

void create_old() { //awk '/Input orientation/,/Distance matrix/ {print}' olddynrun | awk '/   0   / {print}' > old
	ifstream olddynrun("olddynrun");
	ofstream old;
	old.open("old");
	string line;
	while(getline(olddynrun,line)) {
		if(line.find("Input orientation") != string::npos) {
			while(getline(olddynrun,line)) {
				if(line.find("   0   ") != string::npos) old << line << endl;
				if(line.find("Distance matrix") != string::npos) break;
			}
		}
	}
	old.close();
	return;
}

void create_older() { //awk '/Input orientation/,/Distance matrix/ {print}' olderdynrun | awk '/   0   / {print}' > older
	ifstream olderdynrun("olderdynrun");
	ofstream older;
	older.open("older");
	string line;
	while(getline(olderdynrun,line)) {
		if(line.find("Input orientation") != string::npos) {
			while(getline(olderdynrun,line)) {
				if(line.find("   0   ") != string::npos) older << line << endl;
				if(line.find("Distance matrix") != string::npos) break;
			}
		}
	}
	older.close();
	return;
}

void process_x_log() { //awk '/Nuclear Magnetic Resonance/,/HF-SCF/ {if ($2=="C") print $1,$2,"Isotropic =",$3; if ($2=="H") print $1,$2,"Isotropic =",$3}' x.log >> NMRlistcc
	ifstream x_log("x.log");
	ofstream NMRlistcc;
	NMRlistcc.open("NMRlistcc");
	string line;
	while(getline(x_log,line)) {
		if(line.find("Nuclear Magnetic Resonance") != string::npos) {
			while(getline(x_log,line)) {
				stringstream s_stream;
				string pos_1, pos_2, pos_3;
				s_stream << line;
				s_stream >> pos_1;
				s_stream >> pos_2;
				s_stream >> pos_3;
				if(pos_2 == "C" || pos_2 == "H") {
					NMRlistcc << pos_1 << " " << pos_2 << " Isotropic = " << pos_3;
				}
				if(line.find("HF-SCF") != string::npos) break;
			}
		}
	}
	NMRlistcc.close();
}

void create_temp401(string JOB_NAME) {
	string line;
	ifstream log_file(JOB_NAME.c_str());
	ofstream temp401("temp401");
	while(getline(log_file,line)) {
		if(line.find("        1         2") != string::npos) {
			temp401 << line << endl;
			while(getline(log_file,line)) {
				temp401 << line << endl;
				if(line.find("Harmonic frequencies") != string::npos) break;
			}
		}
	}
	temp401.close();
	return;
}

void create_tempfreqs() {
	string line;
	ifstream temp401("temp401");
	ofstream tempfreqs("tempfreqs");
	while(getline(temp401,line)) {
		if(line.find("Frequencies --") != string::npos) {
			stringstream s_stream;
			string discard, pos_3, pos_4, pos_5, pos_6, pos_7;
			s_stream << line;
			s_stream >> discard; s_stream >> discard;
			s_stream >> pos_3 >> pos_4 >> pos_5 >> pos_6 >> pos_7;
			tempfreqs << pos_3 << endl << pos_4 << endl << pos_5 << endl << pos_6 << endl << pos_7 << endl;
		}
	}
	tempfreqs.close();
	return;
}

void create_tempredmass() {
	string line;
	ifstream temp401("temp401");
	ofstream tempredmass("tempredmass");
	while(getline(temp401,line)) {
		if(line.find("Reduced masses") != string::npos) {
			stringstream s_stream;
			string discard, pos_4, pos_5, pos_6, pos_7, pos_8;
			s_stream << line;
			s_stream >> discard; s_stream >> discard; s_stream >> discard;
			s_stream >> pos_4 >> pos_5 >> pos_6 >> pos_7 >> pos_8;
			tempredmass << pos_4 << endl << pos_5 << endl << pos_6 << endl << pos_7 << endl << pos_8 << endl;
		}
	}
	tempredmass.close();
	return;
}

void create_tempfrc() {
	string line;
	ifstream temp401("temp401");
	ofstream tempfrc("tempfrc");
	while(getline(temp401,line)) {
		if(line.find("Force constants") != string::npos) {
			stringstream s_stream;
			string discard, pos_4, pos_5, pos_6, pos_7, pos_8;
			s_stream << line;
			s_stream >> discard; s_stream >> discard; s_stream >> discard;
			s_stream >> pos_4 >> pos_5 >> pos_6 >> pos_7 >> pos_8;
			tempfrc << pos_4 << endl << pos_5 << endl << pos_6 << endl << pos_7 << endl << pos_8 << endl;
		}
	}
	tempfrc.close();
	return;
}

void create_tempmodes() {
	string line;
	ifstream temp401("temp401");
	ofstream tempmodes("tempmodes");
	while(getline(temp401,line)) {
		if(line.find("0") != string::npos) {
			stringstream s_stream;
			string pos_1_str;
			s_stream << line;
			s_stream >> pos_1_str;
			if(pos_1_str.size() == 1) {
				int pos_1_int;
				s_stream.str(string()); s_stream.clear();
				s_stream << pos_1_str;
				s_stream >> pos_1_int;
				if((!s_stream.fail()) && (pos_1_int < 4)) tempmodes << line << endl;
			}
		}
	}
	tempmodes.close();
	return;
}

void create_tempmasses(string JOB_NAME) {
	string line;
	ifstream log_file(JOB_NAME.c_str());
	ofstream tempmasses("tempmasses");
	while(getline(log_file,line)) {
		if(line.find("has atomic number") != string::npos) {
			tempmasses << line << endl;
		}
	}
	tempmasses.close();
	return;
}

void create_tempstangeos(string JOB_NAME) {
	string line;
	ifstream log_file(JOB_NAME.c_str());
	ofstream tempstangeos("tempstangeos");
	while(getline(log_file,line)) {
		if(line.find("Standard orientation:") != string::npos) {
			stringstream s_stream;
			string discard;
			int pos_3 = 1;
			s_stream << line;
			s_stream >> discard; s_stream >> discard;
			s_stream >> pos_3;
			if(pos_3 == 0) tempstangeos << line << endl;
			while(getline(log_file,line)) {
				stringstream s_stream;
				string discard;
				int pos_3 = 1;
				s_stream << line;
				s_stream >> discard; s_stream >> discard;
				s_stream >> pos_3;
				if(pos_3 == 0) tempstangeos << line << endl;
				if(line.find("tional const") != string::npos) break;
			}
		}
	}
	tempstangeos.close();
	return;
}

void create_tempinputgeos(string JOB_NAME) {
	string line;
	ifstream log_file(JOB_NAME.c_str());
	ofstream tempinputgeos("tempinputgeos");
	while(getline(log_file,line)) {
		if(line.find("Input orientation:") != string::npos) {
			stringstream s_stream;
			string discard;
			int pos_3 = 1;
			s_stream << line;
			s_stream >> discard; s_stream >> discard;
			s_stream >> pos_3;
			if(pos_3 == 0) tempinputgeos << line << endl;
			while(getline(log_file,line)) {
				stringstream s_stream;
				string discard;
				int pos_3 = 1;
				s_stream << line;
				s_stream >> discard; s_stream >> discard;
				s_stream >> pos_3;
				if(pos_3 == 0) tempinputgeos << line << endl;
				if(line.find("Stoichiometry") != string::npos) break;
			}
		}
	}
	tempinputgeos.close();
	return;
}

void create_temp_files(string JOB_NAME) {
	create_temp401(JOB_NAME);
	create_tempfreqs();
	create_tempredmass();
	create_tempfrc();
	create_tempmodes();
	create_tempmasses(JOB_NAME);
	create_tempstangeos(JOB_NAME);
	create_tempinputgeos(JOB_NAME);
	return;
}

void read_tempstangeos(vector<Atom>& vecAtoms) {
	int old_atom = 0;
	string line;
	ifstream tempstangeos;
	tempstangeos.open("tempstangeos");
	while(getline(tempstangeos, line)) {
		int atomic_num, atom;
		float x_pos, y_pos, z_pos;
		stringstream ss;
		ss << line;
		ss >> atom;
		ss >> atomic_num;
		ss >> x_pos;
		ss >> x_pos;
		ss >> y_pos;
		ss >> z_pos;
		if(atom < old_atom) break;
		else old_atom = atom;
		Atom new_atom(atomic_num, x_pos, y_pos, z_pos);
		vecAtoms.push_back(new_atom);
	}
	return;
}

void read_tempmasses(vector<Atom>& vecAtoms) {
	string line;
	unsigned int curr_atom = 0;
	ifstream tempmasses;
	tempmasses.open("tempmasses");
	while(getline(tempmasses, line)) {
		float weight = 0;
		string discard;
		stringstream ss;
		ss << line;
		for(int i = 0; i < 8; i++) {
			ss >> discard; // 1-8
		}
		ss >> weight; // 9
		if(weight > 0) {vecAtoms.at(curr_atom).set_atomic_weight(weight);}
		curr_atom++;
	}
	tempmasses.close();
	return;
}

void read_tempfreqs_redmass_frc(unsigned int numFreq, vector<float>& freq, vector<float>& redMass, vector<float>& frc, float scaling) {
	string line;
	ifstream tempfreqs, tempredmass, tempfrc;
	tempfreqs.open("tempfreqs");
	tempredmass.open("tempredmass");
	tempfrc.open("tempfrc");
	for(unsigned int i = 0; i < numFreq; i++) {
		float value;
		if(getline(tempfreqs, line)) {
			stringstream ss;
			value = 0;
			ss << line;
			ss >> value;
			value *= scaling;
			if(value < 0) freq.at(i) = 2;
			else freq.at(i) = value;
		}
		if(getline(tempredmass, line)) {
			stringstream ss;
			value = 1;
			ss << line;
			ss >> value;
			redMass.at(i) = value;
		}
		if(getline(tempfrc, line)) {
			stringstream ss;
			value = 0.0001;
			ss << line;
			ss >> value;
			if(value == 0) frc.at(i) = 0.0001;
			else frc.at(i) = value;
		}	
	}
	return;
}

void read_tempmodes(map< Tuple<int,int,int>, float >& map_tuples, unsigned int numFreq, unsigned int numAtoms) {
	string line, discard;
	ifstream tempmodes;
	tempmodes.open("tempmodes");
	for(unsigned int i = 0; i < numFreq; i += 5) {
		for(unsigned int j = 0; j < (3*numAtoms); j++) {
			if(getline(tempmodes,line)) {
				int pos_1 = 0, pos_2 = 0;
				float pos_4 = 0, pos_5 = 0, pos_6 = 0, pos_7 = 0, pos_8 = 0;
				stringstream ss;
				ss << line;
				ss >> pos_1;   // 1
				ss >> pos_2;   // 2
				ss >> discard; // 3
				ss >> pos_4;   // 4
				ss >> pos_5;   // 5
				ss >> pos_6;   // 6
				ss >> pos_7;   // 7
				ss >> pos_8;   // 8
				Tuple<int,int,int> tuple_1(i,pos_2,pos_1);
				Tuple<int,int,int> tuple_2(i+1,pos_2,pos_1);
				Tuple<int,int,int> tuple_3(i+2,pos_2,pos_1);
				Tuple<int,int,int> tuple_4(i+3,pos_2,pos_1);
				Tuple<int,int,int> tuple_5(i+4,pos_2,pos_1);
				map_tuples.insert(pair<Tuple<int,int,int>,float>(tuple_1,pos_4));
				map_tuples.insert(pair<Tuple<int,int,int>,float>(tuple_2,pos_5));
				map_tuples.insert(pair<Tuple<int,int,int>,float>(tuple_3,pos_6));
				map_tuples.insert(pair<Tuple<int,int,int>,float>(tuple_4,pos_7));
				map_tuples.insert(pair<Tuple<int,int,int>,float>(tuple_5,pos_8));
			}
		}
	}
	tempmodes.close();
	return;
}

void prog1stpoint(ProgdynConf& progdyn_conf, int isomernum, int runpointnum) {
	// input(ProgdynConf& progdyn_conf, int isomernum, int runpointnum, ifstream geoPlusVel)
	// output(ofstream diagnostics, ofstream g09.com)
	string method = progdyn_conf.get_method();
	string meth2 = progdyn_conf.get_meth2();
	string meth3 = progdyn_conf.get_meth3();
	string meth4 = progdyn_conf.get_meth4();
	string meth5 = progdyn_conf.get_meth5();
	string meth6 = progdyn_conf.get_meth6();
	string meth7 = progdyn_conf.get_meth7();
	string charge = progdyn_conf.get_charge();
	string multiplicity = progdyn_conf.get_multiplicity();
	string memory = progdyn_conf.get_memory();
	unsigned int processors =progdyn_conf.get_processors();
	string checkpoint = progdyn_conf.get_checkpoint();
	unsigned int diag = progdyn_conf.get_diag();
	float timestep = progdyn_conf.get_timestep();
	unsigned int highlevel = progdyn_conf.get_highlevel();
	unsigned int linkatoms = progdyn_conf.get_linkatoms();
	unsigned int methodfilelines = progdyn_conf.get_methodfilelines();
	bool killcheck = progdyn_conf.get_killcheck();
	string nmrmethod = progdyn_conf.get_nmrmethod();
	string nmrmethod2 = progdyn_conf.get_nmrmethod2();
	string nmrmethod3 = progdyn_conf.get_nmrmethod3();
	unsigned int nmrtype = progdyn_conf.get_nmrtype();
	unsigned int nmrevery = progdyn_conf.get_nmrevery();
	bool nonstandard = progdyn_conf.get_nonstandard();
	string title1 = progdyn_conf.get_title1();
	string title2 = progdyn_conf.get_title2();
	string title3 = progdyn_conf.get_title3();
	string title4 = progdyn_conf.get_title4();

	if(diag == 1) {
		ofstream diagnostics;
		diagnostics.open("diagnostics");
		diagnostics << "****************** starting prog1stpoint *****************" << endl;
		diagnostics << "method,charge,multiplicity,memory" << endl;
		diagnostics << method << "," << charge << "," << multiplicity << "," << memory << endl;
		diagnostics << "processors,checkpoint,title" << endl;
		diagnostics << processors << "," << checkpoint << "," << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		diagnostics.close();
	}

	vector<Atom> vecAtoms;
	unsigned int numAtoms = 0;
	string line;
	ifstream geoPlusVel;
	geoPlusVel.open("geoPlusVel");
	getline(geoPlusVel, line);
	stringstream ss;
	ss << line;
	ss >> numAtoms;
	for(unsigned int i = 0; i < numAtoms; i++) {
		float x, y, z;
		float atomic_weight;
		string atomic_symbol;
		getline(geoPlusVel, line);
		ss.str(string()); ss.clear();
		ss << line;
		ss >> atomic_symbol;
		ss >> x;
		ss >> y;
		ss >> z;
		ss >> atomic_weight;
		Atom new_atom(x, y, z, atomic_weight, atomic_symbol);
		vecAtoms.push_back(new_atom); 
	}
	for(unsigned int i = 0; i < numAtoms; i++) {
		float x_vel, y_vel, z_vel;
		getline(geoPlusVel, line);
		ss.str(string()); ss.clear();
		ss << line;
		ss >> x_vel;
		ss >> y_vel;
		ss >> z_vel;
		vecAtoms.at(i).set_vel(x_vel, y_vel, z_vel);
	}
	geoPlusVel.close();

	ofstream g09_com;
	g09_com.open("g09.com");
	g09_com << "%nproc=" << processors << endl;
	g09_com << "%mem=" << memory << endl;
	if(killcheck != 1) g09_com << "%chk=" << checkpoint << endl;
	if(nonstandard == 0) {
		g09_com << "# " << method << " force scf=(tight,nosym) " << endl;
		if(meth2 == "unrestricted") g09_com << "guess=mix" << endl;
		if(meth3.size() > 2) g09_com << meth3 << endl;
		if(meth4.size() > 2) g09_com << meth4 << endl;
	}
	else if(nonstandard == 1) {
		g09_com << "# " << endl;
		g09_com << "nonstd" << endl;
		ifstream nonstandard;
		nonstandard.open("nonstandard");
		while(getline(nonstandard, line)) {
			g09_com << line << endl;
		}
		nonstandard.close();
	}
	g09_com << endl;
	g09_com << title1 << " " << title2 << " " << title3 << " " << title4 << endl;
	g09_com << "runpoint  " << runpointnum << endl;
	g09_com << "runisomer  " << isomernum << endl;
	g09_com << endl;
	g09_com << charge << " " << multiplicity << endl;
	g09_com.precision(7);
	g09_com.setf(ios::fixed);
	for(unsigned int i = 0; i < vecAtoms.size(); i++) {
		g09_com << vecAtoms[i].get_atomic_symbol() << " " << vecAtoms[i].get_x() << " " << vecAtoms[i].get_y() << " " << vecAtoms[i].get_z() << endl;
		if((i > highlevel) && (i <= (highlevel+linkatoms))) g09_com << "M H" << endl;
		else if(i > (highlevel+linkatoms)) g09_com << "M" << endl;
	}
	g09_com << endl;
	if(meth5.size() > 2) g09_com << meth5 << endl;
	if(meth6.size() > 2) g09_com << meth6 << endl;
	if(methodfilelines >= 1) {
		ifstream methodfile;
		methodfile.open("methodfile");
		for(unsigned int i = 0; i < methodfilelines; i++) {
			getline(methodfile, line);
			g09_com << line << endl;
		}
		methodfile.close();
	}
	if((nmrtype > 0) && ((runpointnum % nmrevery) == 0)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod << " nmr=giao geom=check" << endl;
		if(nmrmethod == method) g09_com << "guess=tcheck" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;		
	}
	if((nmrtype > 1) && ((runpointnum % nmrevery) == 0)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod2 << " nmr=giao geom=check" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;		
	}
	if((nmrtype > 2) && ((runpointnum % nmrevery) == 0)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod3 << " nmr=giao geom=check" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;		
	}
	g09_com << endl << endl << endl;
	g09_com.close();
}

void prog2ndpoint(ProgdynConf& progdyn_conf, int isomernum, int runpointnum, string trajdirection, string g09_log_filepath) {
	// input(ProgdynConf& progdyn_conf, int isomernum, int runpointnum, string trajdirection, ifstream geoPlusVel, ifstream methodfile)
	// output(ofstream diagnostics, ofstream traj, ofstream Echeck, ofstream g09.com)
	string method = progdyn_conf.get_method();
	string meth2 = progdyn_conf.get_meth2();
	string meth3 = progdyn_conf.get_meth3();
	string meth4 = progdyn_conf.get_meth4();
	string meth5 = progdyn_conf.get_meth5();
	string meth6 = progdyn_conf.get_meth6();
	string meth7 = progdyn_conf.get_meth7();
	string charge = progdyn_conf.get_charge();
	string multiplicity = progdyn_conf.get_multiplicity();
	string memory = progdyn_conf.get_memory();
	unsigned int processors =progdyn_conf.get_processors();
	string checkpoint = progdyn_conf.get_checkpoint();
	unsigned int diag = progdyn_conf.get_diag();
	float timestep = progdyn_conf.get_timestep();
	unsigned int highlevel = progdyn_conf.get_highlevel();
	unsigned int linkatoms = progdyn_conf.get_linkatoms();
	int fixedatom1 = progdyn_conf.get_fixedatom1();
	int fixedatom2 = progdyn_conf.get_fixedatom2();
	int fixedatom3 = progdyn_conf.get_fixedatom3();
	int fixedatom4 = progdyn_conf.get_fixedatom4();
	bool DRP = progdyn_conf.get_DRP();
	unsigned int methodfilelines = progdyn_conf.get_methodfilelines();
	bool killcheck = progdyn_conf.get_killcheck();
	unsigned int etolerance = progdyn_conf.get_etolerance();
	string nmrmethod = progdyn_conf.get_nmrmethod();
	string nmrmethod2 = progdyn_conf.get_nmrmethod2();
	string nmrmethod3 = progdyn_conf.get_nmrmethod3();
	unsigned int nmrtype = progdyn_conf.get_nmrtype();
	unsigned int nmrevery = progdyn_conf.get_nmrevery();
	bool nonstandard = progdyn_conf.get_nonstandard();
	string title1 = progdyn_conf.get_title1();
	string title2 = progdyn_conf.get_title2();
	string title3 = progdyn_conf.get_title3();
	string title4 = progdyn_conf.get_title4();

	unsigned int numAtoms=0;
	float initialDis;
	float potentialE, newPotentialE;
	float newPotentialEK;
	int scfcount=0;
	float KEinitmodes, KEinittotal;
	float desiredModeEnK;
	float pddga, pddgb, pddgc, qm;
	string line;

	if(diag >= 1) {
		ofstream diagnostics;
		diagnostics.open("diagnostics");
		diagnostics << "****************** starting prog2ndpoint *****************" << endl;
		diagnostics << "method,charge,multiplicity,memory" << endl;
		diagnostics << method << "," << charge << "," << multiplicity << "," << memory << endl;
		diagnostics << "processors,checkpoint,title" << endl;
		diagnostics << processors << "," << checkpoint << "," << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		diagnostics.close();
	}

	ofstream g09_com;
	g09_com.open("g09.com");
	g09_com << "%nproc=" << processors << endl;
	g09_com << "%mem=" << memory << endl;
	if(killcheck != 1) g09_com << "%chk=" << checkpoint << endl;
	if(nonstandard == 0) {
		g09_com << "# " << method << " force scf=(tight,nosym) " << endl;
		if(meth2 == "unrestricted") g09_com << "guess=mix" << endl;
		if(meth3.size() > 2) g09_com << meth3 << endl;
		if(meth4.size() > 2) g09_com << meth4 << endl;
	}
	else if(nonstandard == 1) {
		g09_com << "# " << endl;
		g09_com << "nonstd" << endl;
		ifstream nonstandard;
		nonstandard.open("nonstandard");
		while(getline(nonstandard, line)) {
			g09_com << line << endl;
		}
		nonstandard.close();
	}
	g09_com << endl;
	g09_com << title1 << " " << title2 << " " << title3 << " " << title4 << endl;
	g09_com << "runpoint  " << runpointnum << endl;
	g09_com << "runisomer  " << isomernum << endl;
	g09_com << endl;
	g09_com << charge << " " << multiplicity << endl;

	vector<Atom> vecAtoms;
	ifstream geoPlusVel;
	geoPlusVel.open("geoPlusVel");
	getline(geoPlusVel, line);
	stringstream ss;
	ss << line;
	ss >> numAtoms;
	for(unsigned int i = 0; i < numAtoms; i++) {
		float x, y, z;
		float atomic_weight;
		string atomic_symbol;
		getline(geoPlusVel, line);
		ss.str(string()); ss.clear();
		ss << line;
		ss >> atomic_symbol;
		ss >> x;
		ss >> y;
		ss >> z;
		ss >> atomic_weight;
		Atom new_atom(x, y, z, atomic_weight, atomic_symbol);
		vecAtoms.push_back(new_atom); 
	}
	for(unsigned int i = 0; i < numAtoms; i++) {
		float x_vel, y_vel, z_vel;
		getline(geoPlusVel, line);
		ss.str(string()); ss.clear();
		ss << line;
		ss >> x_vel;
		ss >> y_vel;
		ss >> z_vel;
		vecAtoms.at(i).set_vel(x_vel, y_vel, z_vel);
	}
	for(unsigned int i = 0; i < numAtoms; i++) {
		if(trajdirection=="reverserestart") {
			vecAtoms.at(i).set_x_after_vel(vecAtoms.at(i).get_x() - vecAtoms.at(i).get_x_vel());
			vecAtoms.at(i).set_y_after_vel(vecAtoms.at(i).get_y() - vecAtoms.at(i).get_y_vel());
			vecAtoms.at(i).set_z_after_vel(vecAtoms.at(i).get_z() - vecAtoms.at(i).get_z_vel());
		}
		else {
			vecAtoms.at(i).set_x_after_vel(vecAtoms.at(i).get_x() + vecAtoms.at(i).get_x_vel());
			vecAtoms.at(i).set_y_after_vel(vecAtoms.at(i).get_y() + vecAtoms.at(i).get_y_vel());
			vecAtoms.at(i).set_z_after_vel(vecAtoms.at(i).get_z() + vecAtoms.at(i).get_z_vel());
		}
		if(diag > 1) {
			ofstream diagnostics;
			diagnostics.open("diagnostics", fstream::app);
			if(i==0) diagnostics << "geometry after adding velocities" << endl;
			diagnostics << vecAtoms.at(i).get_x_after_vel() << " " << vecAtoms.at(i).get_y_after_vel() << " " << vecAtoms.at(i).get_z_after_vel() << endl;
			diagnostics.close();
		}
	}
	while(getline(geoPlusVel, line)) {
		if(line == "") break;
		else {
			ss.str(string()); ss.clear();
			string discard, pos_4, pos_11;
			float pos_5, pos_9, pos_13;
			ss << line;
			ss >> discard; // 1
			ss >> discard; // 2
			ss >> discard; // 3
			ss >> pos_4;   // 4
			ss >> pos_5;   // 5
			ss >> discard; // 6
			ss >> discard; // 7
			ss >> discard; // 8
			ss >> pos_9;   // 9
			ss >> discard; // 10
			ss >> pos_11;  // 11
			ss >> discard; // 12
			ss >> pos_13;  // 13
			if(pos_4 == "desired=") desiredModeEnK = pos_5;
			else if(pos_4 == "modes=") {
				KEinitmodes = pos_5;
				KEinittotal = pos_9;
			}
			if(pos_11 == "potential") potentialE = pos_13;
		}
	}
	geoPlusVel.close();

	ofstream traj;
	traj.open("traj");
	traj << numAtoms << endl;
	traj << potentialE << " " << title1 << " " << title2 << " " << title3 << " " << title4 << " runpoint 1 runisomer " << isomernum << endl;
	for(unsigned int i = 0 ; i < numAtoms; i++) {
		traj << vecAtoms.at(i).get_atomic_symbol() << " " << to_string(vecAtoms.at(i).get_x()) << " " << to_string(vecAtoms.at(i).get_y()) << " " << to_string(vecAtoms.at(i).get_z()) << endl;
	}
	traj.close();

	ifstream g09_log;
	g09_log.open(g09_log_filepath.c_str());
	while(getline(g09_log, line)) {
		if(line.find("SCF DONE") != string::npos || line.find("EUMP2 =") != string::npos|| line.find("Energy=") != string::npos) {
			ss.str(string()); ss.clear();
			string pos_1, pos_3, pos_6;
			float pos_2, discard, pos_5;
			ss << line;
			ss >> pos_1;
			ss >> pos_2;
			ss >> pos_3;
			ss >> discard;
			ss >> pos_5;
			ss >> pos_6;
			if(pos_1 == "Energy=" && pos_3 == "NIter=") newPotentialE = pos_2;
			else if(pos_1 == "SCF" && scfcount == 0) newPotentialE = pos_5;
			else if(pos_1 == "E2") {
				vector<float> temp_vec;
				unsigned int beg_pos = 0;
				for(unsigned int i = 0; i < pos_6.size(); i++) {
					if(pos_6.at(i) == 'D' || i == pos_6.size() - 1) {
						stringstream str_float;
						str_float << pos_6.substr(beg_pos, i-beg_pos);
						float value;
						str_float >> value;
						temp_vec.push_back(value);
						beg_pos = i;
					}
				}
				newPotentialE = temp_vec.at(1)*pow(10,temp_vec.at(2));
			}
			newPotentialEK = (newPotentialE - potentialE) * 627.509;
			if(pos_1 == "SCF") {
				if(scfcount == 0) pddga = pos_5;
				else if(scfcount == 1) qm = pos_5;
				else if(scfcount == 2) {
					pddgb = pos_5;
					pddgc = pddga - pddgb;
					newPotentialE = qm + pddgc;
					newPotentialEK = (newPotentialE - potentialE) * 627.509;
				}
				scfcount++;
			}
		}
		ss.str(string()); ss.clear();
		int pos_1;
		string discard, pos_3_str;
		ss << line;
		ss >> pos_1;
		ss >> discard;
		ss >> pos_3_str;
		if(pos_3_str.size() > 9) {
			unsigned int atom_num = 1000;
			float pos_3, pos_4, pos_5;
			ss >> pos_4;
			ss >> pos_5;
			ss.str(string()); ss.clear();
			ss << pos_3_str;
			ss >> pos_3;
			if(line.find("      1    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("      2    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("      3    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("      4    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("      5    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("      6    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("      7    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("      8    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("      9    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     10    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     11    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     12    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     13    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     14    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     15    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     16    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     17    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     18    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     19    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     20    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     21    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     22    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     23    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     24    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     25    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     26    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     27    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     28    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     29    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("     30    ") != string::npos) atom_num = pos_1-1;
			if(atom_num != 1000 && atom_num < numAtoms) {
				vecAtoms.at(atom_num).set_x_force(pos_3);
				vecAtoms.at(atom_num).set_y_force(pos_4);
				vecAtoms.at(atom_num).set_z_force(pos_5);
				if(diag > 1) {
					ofstream diagnostics;
					diagnostics.open("diagnostics", fstream::app);
					if(atom_num == 0) diagnostics << "i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3]" << endl;
					diagnostics << atom_num << " " << vecAtoms.at(atom_num).get_atomic_weight() << " " << vecAtoms.at(atom_num).get_x_force();
					diagnostics << " " << vecAtoms.at(atom_num).get_y_force() << " " << vecAtoms.at(atom_num).get_z_force() << endl;
					diagnostics.close();
				}
			}
		}
	}
	g09_log.close();

	if(DRP == 0) {
		ofstream Echeck;
		Echeck.open("Echeck", fstream::app);
		Echeck << "trajectory # " << isomernum << endl;
		Echeck << "point 1 potential E= " << newPotentialEK << "   point 1 kinetic E= " << KEinitmodes << "  Total = " << newPotentialEK+KEinitmodes << endl;
		Echeck << "desired total energy= " << desiredModeEnK << endl;
		if(newPotentialEK+KEinitmodes > desiredModeEnK+etolerance) Echeck << "XXXX bad total Energy" << endl;
		else if(newPotentialEK+KEinitmodes < desiredModeEnK-etolerance) Echeck << "XXXX bad total Energy" << endl;
		Echeck.close();
	}
	for(unsigned int i = 0; i < numAtoms; i++) {
		if(DRP == 1) {
			vecAtoms.at(i).set_x_force(0);
			vecAtoms.at(i).set_y_force(0);
			vecAtoms.at(i).set_z_force(0);
		}
		else {
			//cout << i+1 << " 1 " << vecAtoms.at(i).get_x_force() << endl;
			vecAtoms.at(i).set_x_force(0.5*1E20*vecAtoms.at(i).get_x_force()*627.509*(4184/(0.529177*avNum))*(pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
			//cout << i+1 << " 1 " << vecAtoms.at(i).get_x_force() << endl;
			//cout << i+1 << " 2 " << vecAtoms.at(i).get_y_force() << endl;
			vecAtoms.at(i).set_y_force(0.5*1E20*vecAtoms.at(i).get_y_force()*627.509*(4184/(0.529177*avNum))*(pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
			//cout << i+1 << " 2 " << vecAtoms.at(i).get_y_force() << endl;
			//cout << i+1 << " 3 " << vecAtoms.at(i).get_z_force() << endl;
			vecAtoms.at(i).set_z_force(0.5*1E20*vecAtoms.at(i).get_z_force()*627.509*(4184/(0.529177*avNum))*(pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
			//cout << i+1 << " 3 " << vecAtoms.at(i).get_z_force() << endl;
		}
		vecAtoms.at(i).set_x_after_vel(vecAtoms.at(i).get_x_after_vel() + vecAtoms.at(i).get_x_force());
		vecAtoms.at(i).set_y_after_vel(vecAtoms.at(i).get_y_after_vel() + vecAtoms.at(i).get_y_force());
		vecAtoms.at(i).set_z_after_vel(vecAtoms.at(i).get_z_after_vel() + vecAtoms.at(i).get_z_force());
		if((int(i) == fixedatom1) || (int(i) == fixedatom2) || (int(i) == fixedatom3) || (int(i) == fixedatom4)) {
			cout << fixedatom1 << " " << fixedatom2 << " " << fixedatom3 << " " << fixedatom4 << endl;		
			vecAtoms.at(i).set_x_after_vel(vecAtoms.at(i).get_x());
			vecAtoms.at(i).set_y_after_vel(vecAtoms.at(i).get_y());
			vecAtoms.at(i).set_z_after_vel(vecAtoms.at(i).get_z());
		}
		if(diag > 1) {
			ofstream diagnostics;
			diagnostics.open("diagnostics", fstream::app);
			if(i == 0) diagnostics << "i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3]" << endl;
			diagnostics << i << " " << vecAtoms.at(i).get_atomic_weight() << " " << vecAtoms.at(i).get_x_force();
			diagnostics << " " << vecAtoms.at(i).get_y_force() << " " << vecAtoms.at(i).get_z_force() << endl;
			diagnostics.close();
		}
		g09_com.precision(7);
		g09_com.setf(ios::fixed);
		g09_com << vecAtoms[i].get_atomic_symbol() << " " << vecAtoms[i].get_x_after_vel() << " " << vecAtoms[i].get_y_after_vel() << " " << vecAtoms[i].get_z_after_vel() << endl;
		if((i > highlevel) && (i <= (highlevel+linkatoms))) g09_com << "M H" << endl;
		else if(i > (highlevel+linkatoms)) g09_com << "M" << endl;
	}
	g09_com << endl;
	if(meth5.size() > 2) g09_com << meth5 << endl;
	if(meth6.size() > 2) g09_com << meth6 << endl;
	if(methodfilelines >= 1) {
		ifstream methodfile;
		methodfile.open("methodfile");
		for(unsigned int i = 0; i < methodfilelines; i++) {
			getline(methodfile, line);
			g09_com << line << endl;
		}
		methodfile.close();
	}
	if((nmrtype > 0) && ((runpointnum % nmrevery) == 0)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod << " nmr=giao geom=check" << endl;
		if(nmrmethod == method) g09_com << "guess=tcheck" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;		
	}
	if((nmrtype > 1) && ((runpointnum % nmrevery) == 0)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod2 << " nmr=giao geom=check" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;		
	}
	if((nmrtype > 2) && ((runpointnum % nmrevery) == 0)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod2 << " nmr=giao geom=check" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;
	}
	g09_com << endl << endl << endl;
	g09_com.close();
	traj.open("traj", fstream::app);
	traj << numAtoms << endl;
	traj << newPotentialE << " " << title1 << " " << title2 << " " << title3 << " " << title4 << " runpoint " << runpointnum << " runisomer " << isomernum << endl;
	for(unsigned int i = 0; i < numAtoms; i++) {
		traj << vecAtoms.at(i).get_atomic_symbol() << " " << vecAtoms.at(i).get_x() << " " << vecAtoms.at(i).get_y() << " " << vecAtoms.at(i).get_z() << endl;
	}
	traj.close();
}

void proganal(int isomernum, int runpointnum, string trajdirection, bool& nogo, string log_filepath) {
	// input(int isomernum, int runpointnum, string trajdirection, bool& nogo, string log_filepath, ifstream Echeck)
	// ouput(ofstream dynfollowfile, bool nogo)
	vector<float> A,B,C;
	ofstream dynfollowfile;
	dynfollowfile.open("dynfollowfile", fstream::app);
	ifstream log_file;	
	log_file.open(log_filepath.c_str());
	string line;
	stringstream ss_line;
	int first_title = 1;
	while(getline(log_file, line)) {
		if(line.find(" alame3") != string::npos) {
			if(first_title == 1) {
				ss_line.str(string()); ss_line.clear();
				ss_line << line;
				vector<string> vec_strings;
				string temp_string;
				for(unsigned int i = 0; i < 8; i++) {
					if(i == 6) {
						ss_line >> runpointnum;
						vec_strings.push_back(to_string(runpointnum));
					}
					else {
						ss_line >> temp_string;
						vec_strings.push_back(temp_string);
					}
				}
				dynfollowfile << vec_strings[0] << vec_strings[1] << vec_strings[2] << vec_strings[3] << vec_strings[5] << vec_strings [6] << vec_strings[7] << " ";
			}
			first_title++;
		}
		if(line.find("Standard orientation") != string::npos) {
			while(getline(log_file, line)) {
				unsigned int atom_num;
				string discard;
				ss_line.str(string()); ss_line.clear();
				ss_line << line;
				ss_line >> atom_num; // 1
				if((atom_num > 0) && (atom_num < 30)) {
					float x,y,z;
					ss_line >> discard; // 2
					ss_line >> discard; // 3
					ss_line >> x;
					ss_line >> y;
					ss_line >> z;
					A.push_back(x);
					B.push_back(y);
					C.push_back(z);
				}
				if(line.find("Rotational constants") != string::npos) break;
			}
		}
		if(line.find("before annihilation") != string::npos) {
			string discard, string_value;
			float float_value;
			ss_line.str(string()); ss_line.clear();
			ss_line << line;
			ss_line >> string_value; // 1
			ss_line >> discard; // 2
			ss_line >> discard; // 3
			ss_line >> discard; // 4
			ss_line >> discard; // 5
			ss_line >> float_value; // 6
			dynfollowfile << string_value << " " << float_value << " ";
		}
	}
	
	float C1N6 = Distance(1,6,A,B,C);
	float C3C8 = Distance(3,8,A,B,C);
	float C1C8 = Distance(1,8,A,B,C);
	dynfollowfile << "C1N6 " <<  C1N6 <<  "  C3C8 " <<  C3C8 << "  C1C8 " << C1C8 << " ";
	if(runpointnum > 500) {
		dynfollowfile << "Too many points.  XXXXN" << endl;
		nogo = true;
	}
	if((C1N6 > 2.5) && (C3C8 < 1.8)) {
		if(trajdirection == "reverse") nogo = true;
		dynfollowfile << "2,3-Rearrangement XXXX12N" << endl;
	}
	if((C1N6 > 2.5) && (C1C8 < 1.8)) {
		if(trajdirection == "reverse") nogo = true;
		dynfollowfile << "1,2-Rearrangement XXXX12N" << endl;
	}
	if((C1N6 > 3.2) && (C3C8 > 3.4) && (C1C8 > 3.2)) {
		if(trajdirection == "reverse") nogo = true;
		dynfollowfile << "Dissociation XXXX12N" << endl;
	}
	if((C1N6 < 1.7) && (C3C8 > 2.8)) {
		if(trajdirection == "reverse") nogo = true;
		dynfollowfile << "reformed starting material XXXX12N" << endl;
	}
	dynfollowfile << print_date_to_file() << endl; // date '+%b:%d:%Y %T' >> dynfollowfile
	dynfollowfile << str_find_XXXX_last_n_lines("Echeck", 1) << endl; // tail -1 Echeck | grep XXXX >> dynfollowfile
	dynfollowfile.close();
}

void progdynb(ProgdynConf& progdyn_conf, int isomernum, int runpointnum, string g09_log_filepath) {
	// input(ProgdynConf& progdyn_conf, int isomernum, int runpointnum, string g09_log_filepath, ifstream geoPlusVel, ifstream old, ifstream older, ifstream oldAdjForces, ifstream maxMove, ifsream uptimelist) 
	// output(ofstream diagnostics, ofstream vellist, ofstream oldAdjForces, ofstream maxMove, ofstream ZMAT, ofstream dyn, ofstream traj, ofstream g09.com)
	string method = progdyn_conf.get_method();
	string meth2 = progdyn_conf.get_meth2();
	string meth3 = progdyn_conf.get_meth3();
	string meth4 = progdyn_conf.get_meth4();
	string meth5 = progdyn_conf.get_meth5();
	string meth6 = progdyn_conf.get_meth6();
	string meth7 = progdyn_conf.get_meth7();
	string charge = progdyn_conf.get_charge();
	string multiplicity = progdyn_conf.get_multiplicity();
	string memory = progdyn_conf.get_memory();
	unsigned int processors =progdyn_conf.get_processors();
	string checkpoint = progdyn_conf.get_checkpoint();
	unsigned int diag = progdyn_conf.get_diag();
	float timestep = progdyn_conf.get_timestep();
	float temp = progdyn_conf.get_temp();
	unsigned int highlevel = progdyn_conf.get_highlevel();
	unsigned int linkatoms = progdyn_conf.get_linkatoms();
	int fixedatom1 = progdyn_conf.get_fixedatom1();
	int fixedatom2 = progdyn_conf.get_fixedatom2();
	int fixedatom3 = progdyn_conf.get_fixedatom3();
	int fixedatom4 = progdyn_conf.get_fixedatom4();
	bool DRP = progdyn_conf.get_DRP();
	unsigned int methodfilelines = progdyn_conf.get_methodfilelines();
	bool killcheck = progdyn_conf.get_killcheck();
	float damping = progdyn_conf.get_damping();
	string nmrmethod = progdyn_conf.get_nmrmethod();
	string nmrmethod2 = progdyn_conf.get_nmrmethod2();
	string nmrmethod3 = progdyn_conf.get_nmrmethod3();
	unsigned int nmrtype = progdyn_conf.get_nmrtype();
	unsigned int nmrevery = progdyn_conf.get_nmrevery();
	bool boxon = progdyn_conf.get_boxon();
	float boxsize = progdyn_conf.get_boxsize();
	bool nonstandard = progdyn_conf.get_nonstandard();
	string title1 = progdyn_conf.get_title1();
	string title2 = progdyn_conf.get_title2();
	string title3 = progdyn_conf.get_title3();
	string title4 = progdyn_conf.get_title4();

	int keepevery=1;
	float maxAtomMove=0.0;
	float nmrrand=0, nmrdo, nmrcc=0, loadlimit=0, KEatomstotal=0, thermostat=0, thermostatmult=1.0, newPotentialE, maxForce, oscillTest;

	if(diag == 1) {
		ofstream diagnostics;
		diagnostics.open("diagnostics");
		diagnostics << "****************** starting progdynb *****************" << endl;
		diagnostics << "method,charge,multiplicity,memory" << endl;
		diagnostics << method << "," << charge << "," << multiplicity << "," << memory << endl;
		diagnostics << "processors,checkpoint,title" << endl;
		diagnostics << processors << "," << checkpoint << "," << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		diagnostics.close();
	}

	vector<Atom> vecAtoms;
	unsigned int numAtoms = 0;
	string line;
	ifstream geoPlusVel;
	geoPlusVel.open("geoPlusVel");
	getline(geoPlusVel, line);
	stringstream ss;
	ss << line;
	ss >> numAtoms;
	for(unsigned int i = 0; i < numAtoms; i++) {
		float atomic_weight;
		string atomic_symbol, discard;
		getline(geoPlusVel, line);
		ss.str(string()); ss.clear();
		ss << line;
		ss >> atomic_symbol; // 1
		ss >> discard;       // 2
		ss >> discard;       // 3
		ss >> discard;       // 4
		ss >> atomic_weight; // 5
		Atom new_atom(0, 0, 0, atomic_weight, atomic_symbol);
		vecAtoms.push_back(new_atom); 
	}
	geoPlusVel.close();

	ifstream old;
	old.open("old");
	for(unsigned int i = 0; i < numAtoms; i++) {
		string discard;
		float pos_4, pos_5, pos_6;
		getline(old, line);
		ss.str(string()); ss.clear();
		ss << line;
		ss >> discard; // 1
		ss >> discard; // 2
		ss >> discard; // 3
		ss >> pos_4;   // 4
		ss >> pos_5;   // 5
		ss >> pos_6;   // 6
		vecAtoms.at(i).set_old_x(pos_4);
		vecAtoms.at(i).set_old_y(pos_5);
		vecAtoms.at(i).set_old_z(pos_6);
	}
	old.close();

	ifstream older;
	older.open("older");
	for(unsigned int i = 0; i < numAtoms; i++) {
		string discard;
		float pos_4, pos_5, pos_6;
		getline(older, line);
		ss.str(string()); ss.clear();
		ss << line;
		ss >> discard; // 1
		ss >> discard; // 2
		ss >> discard; // 3
		ss >> pos_4;   // 4
		ss >> pos_5;   // 5
		ss >> pos_6;   // 6
		vecAtoms.at(i).set_older_x(pos_4);
		vecAtoms.at(i).set_older_y(pos_5);
		vecAtoms.at(i).set_older_z(pos_6);
	}
	older.close();

	if(DRP == 1) {
		ifstream oldAdjForces;
		oldAdjForces.open("oldAdjForces");
		for(unsigned int i = 0; i < numAtoms; i++) {
			float pos_1, pos_2, pos_3;
			getline(oldAdjForces, line);
			ss.str(string()); ss.clear();
			ss << line;
			ss >> pos_1; // 1
			ss >> pos_2; // 2
			ss >> pos_3; // 3
			vecAtoms.at(i).set_old_x_force(pos_1);
			vecAtoms.at(i).set_old_y_force(pos_2);
			vecAtoms.at(i).set_old_z_force(pos_3);
		}
		oldAdjForces.close();

		ifstream maxMove;
		maxMove.open("maxMove");
		float pos_1;
		getline(maxMove, line);
		ss.str(string()); ss.clear();
		ss << line;
		ss >> pos_1;
		if(pos_1 < maxAtomMove && pos_1 > 0) maxAtomMove = pos_1;
		if(maxAtomMove < 0.000001) maxAtomMove = 0.000001;
		maxMove.close();
	}

	if((nmrrand == 0) && (runpointnum % nmrevery == 0)) nmrdo = 1;
	else if((nmrrand == 1) && (rand() < (1/nmrevery))) nmrdo = 1;

	ifstream uptimelist;
	uptimelist.open("uptimelist");
	string discard, pos_10;
	getline(uptimelist, line);
	uptimelist.close();
	ss.str(string()); ss.clear();
	ss << line;
	ss >> discard >> discard >> discard >> discard >> discard >> discard >> discard >> discard >> discard >> pos_10;
	float temp_float;
	ss.str(string()); ss.clear();
	if(pos_10.size() >= 4) ss << pos_10.substr(1,3);
	else ss << 0.0;
	ss >> temp_float;
	float x = 1.0001*temp_float; if(x < 8) x = 8;

	if((nmrrand == 1) && (x > loadlimit)) nmrdo = 0;

	ofstream vellist;
	vellist.open("vellist");
	if(diag == 3) vellist << "runpoint " << runpointnum - 1 << " runisomer " << isomernum << endl;
	for(unsigned int i = 0; i < numAtoms; i++) {
		float atomVel = pow(pow((vecAtoms.at(i).get_old_x() - vecAtoms.at(i).get_older_x()),2) 
		+ pow(vecAtoms.at(i).get_old_y() - vecAtoms.at(i).get_older_y(),2) 
		+ pow(vecAtoms.at(i).get_old_z() - vecAtoms.at(i).get_older_z(),2),.5);
		KEatomstotal = KEatomstotal + 0.5*vecAtoms.at(i).get_atomic_weight()*pow(atomVel,2)/(pow(timestep,2)*conver1);
		if(diag == 3) vellist << atomVel << endl;
	}
	float apparentTemp = KEatomstotal*2/(3*RgasK*numAtoms);
	if(diag == 4) vellist << "KEatomstotal " << KEatomstotal << " apparent Temperature " << apparentTemp << endl;
	if(thermostat == 1) {
		if(diag < 4) vellist << "KEatomstotal " << KEatomstotal << " desired temperature " << temp << " apparent Temperature " << apparentTemp << endl;
		if(apparentTemp > temp) damping = thermostatmult;
		else if(apparentTemp < temp) damping = 1/thermostatmult;
	}
	vellist.close();

	ifstream g09_log;
	g09_log.open(g09_log_filepath.c_str());
	while(getline(g09_log, line)) {
		if(line.find("SCF DONE") != string::npos || line.find("EUMP2 =") != string::npos || line.find("Energy=") != string::npos || line.find("ONIOM:") != string::npos) {
			ss.str(string()); ss.clear();
			string pos_1, pos_3, pos_6;
			float pos_2, discard, pos_5;
			ss << line;
			ss >> pos_1;
			ss >> pos_2;
			ss >> pos_3;
			ss >> discard;
			ss >> pos_5;
			ss >> pos_6;
			if(pos_1 == "Energy=" && pos_3 == "NIter=") newPotentialE = pos_2;
			else if(pos_1 == "SCF") newPotentialE = pos_5;
			else if(pos_1 == "extrapolated") newPotentialE = pos_5;
			else if(pos_1 == "E2") {
				vector<float> temp_vec;
				unsigned int beg_pos = 0;
				for(unsigned int i = 0; i < pos_6.size(); i++) {
					if(pos_6.at(i) == 'D' || i == pos_6.size() - 1) {
						stringstream str_float;
						str_float << pos_6.substr(beg_pos, i-beg_pos);
						float value;
						str_float >> value;
						temp_vec.push_back(value);
						beg_pos = i;
					}
				}
				newPotentialE = temp_vec.at(1)*pow(10,temp_vec.at(2));
			}
		}
		ss.str(string()); ss.clear();
		int pos_1;
		string discard, pos_3_str;
		ss << line;
		ss >> pos_1;
		ss >> discard;
		ss >> pos_3_str;
		if(pos_3_str.size() > 9) {
			unsigned int atom_num = 1000;
			float pos_3, pos_4, pos_5;
			ss >> pos_4;
			ss >> pos_5;
			ss.str(string()); ss.clear();
			ss << pos_3_str;
			ss >> pos_3;
			if(line.find("        1    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("        2    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("        3    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("        4    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("        5    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("        6    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("        7    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("        8    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("        9    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       10    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       11    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       12    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       13    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       14    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       15    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       16    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       17    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       18    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       19    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       20    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       21    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       22    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       23    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       24    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       25    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       26    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       27    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       28    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       29    ") != string::npos) atom_num = pos_1-1;
			else if(line.find("       30    ") != string::npos) atom_num = pos_1-1;
			if(atom_num != 1000 && atom_num < numAtoms) {
				cout << atom_num << " " << pos_3 << " " << pos_4 << " " << pos_5 << endl;
				vecAtoms.at(atom_num).set_x_force(pos_3);
				vecAtoms.at(atom_num).set_y_force(pos_4);
				vecAtoms.at(atom_num).set_z_force(pos_5);
				if(diag > 1) {
					ofstream diagnostics;
					diagnostics.open("diagnostics", fstream::app);
					if(atom_num == 0) diagnostics << "i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3]" << endl;
					diagnostics << atom_num << " " << vecAtoms.at(atom_num).get_atomic_weight() << " " << vecAtoms.at(atom_num).get_x_force();
					diagnostics << " " << vecAtoms.at(atom_num).get_y_force() << " " << vecAtoms.at(atom_num).get_z_force() << endl;
					diagnostics.close();
				}
			}
		}
	}
	g09_log.close();

	if(DRP == 1) {
		maxForce = 0;
		oscillTest = 0;
		ofstream oldAdjForces;
		oldAdjForces.open("oldAdjForces");
		oldAdjForces.precision(8);
		oldAdjForces.setf(ios::fixed);
		for(unsigned int i = 0; i < numAtoms; i++) {
			vecAtoms.at(i).set_x_force(1E20*vecAtoms.at(i).get_x_force()*627.509*(4184/(0.529177*avNum)*pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
			vecAtoms.at(i).set_y_force(1E20*vecAtoms.at(i).get_y_force()*627.509*(4184/(0.529177*avNum)*pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
			vecAtoms.at(i).set_z_force(1E20*vecAtoms.at(i).get_z_force()*627.509*(4184/(0.529177*avNum)*pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
		oscillTest = vecAtoms.at(i).get_x_force()*vecAtoms.at(i).get_old_x_force() +  vecAtoms.at(i).get_x_force()*vecAtoms.at(i).get_old_x_force() + vecAtoms.at(i).get_x_force()*vecAtoms.at(i).get_old_x_force();
			if(vecAtoms.at(i).get_x_force() > maxForce) maxForce = vecAtoms.at(i).get_x_force();
			if(vecAtoms.at(i).get_y_force() > maxForce) maxForce = vecAtoms.at(i).get_y_force();
			if(vecAtoms.at(i).get_z_force() > maxForce) maxForce = vecAtoms.at(i).get_z_force();
			if((0 - vecAtoms.at(i).get_x_force()) > maxForce) maxForce = -vecAtoms.at(i).get_x_force();
			if((0 - vecAtoms.at(i).get_y_force()) > maxForce) maxForce = -vecAtoms.at(i).get_y_force();
			if((0 - vecAtoms.at(i).get_z_force()) > maxForce) maxForce = -vecAtoms.at(i).get_z_force();
			oldAdjForces << vecAtoms.at(i).get_x_force() << " " << vecAtoms.at(i).get_y_force() << " " << vecAtoms.at(i).get_z_force() << endl;
		}
		oldAdjForces << "oscillTest " << oscillTest << endl;

		ofstream maxMove;
		maxMove.open("maxMove");
		if(oscillTest < 0) {
			maxAtomMove = maxAtomMove * 0.5;
			maxMove << maxAtomMove << endl;
		}
		else if(oscillTest > 0) {
			maxAtomMove = maxAtomMove * 1.2;
			maxMove << maxAtomMove << endl;
		}
		oldAdjForces << "maxAtomMove " << maxAtomMove << endl;
		float forceMult = maxAtomMove/maxForce;
		for(unsigned int i = 0; i < numAtoms; i++) {
			vecAtoms.at(i).set_new_x_force(vecAtoms.at(i).get_old_x_force() + forceMult*vecAtoms.at(i).get_x_force());
			vecAtoms.at(i).set_new_y_force(vecAtoms.at(i).get_old_y_force() + forceMult*vecAtoms.at(i).get_y_force());
			vecAtoms.at(i).set_new_z_force(vecAtoms.at(i).get_old_z_force() + forceMult*vecAtoms.at(i).get_z_force());
		}
		oldAdjForces.close();
	}
	else if(DRP == 0) {
		ofstream diagnostics;
		diagnostics.open("diagnostics");
		for(unsigned int i = 0; i < numAtoms; i++) {
			vecAtoms.at(i).set_x_force(1E20*vecAtoms.at(i).get_x_force()*627.509*(4184/(0.529177*avNum)*pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
			vecAtoms.at(i).set_y_force(1E20*vecAtoms.at(i).get_y_force()*627.509*(4184/(0.529177*avNum)*pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
			vecAtoms.at(i).set_z_force(1E20*vecAtoms.at(i).get_z_force()*627.509*(4184/(0.529177*avNum)*pow(timestep,2))/(vecAtoms.at(i).get_atomic_weight()/(avNum*1000)));
			if(diag > 1 && i == 0) diagnostics << "i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3]" << endl;
			if(diag > 1) diagnostics << i << " " << vecAtoms.at(i).get_atomic_weight() << " " << vecAtoms.at(i).get_x_force() << " " << vecAtoms.at(i).get_y_force() << " " << vecAtoms.at(i).get_z_force() << endl;
			vecAtoms.at(i).set_new_x_force(vecAtoms.at(i).get_old_x() + damping*(vecAtoms.at(i).get_old_x() - vecAtoms.at(i).get_older_x()) + vecAtoms.at(i).get_x_force());
			vecAtoms.at(i).set_new_y_force(vecAtoms.at(i).get_old_y() + damping*(vecAtoms.at(i).get_old_y() - vecAtoms.at(i).get_older_y()) + vecAtoms.at(i).get_y_force());
			vecAtoms.at(i).set_new_z_force(vecAtoms.at(i).get_old_z() + damping*(vecAtoms.at(i).get_old_z() - vecAtoms.at(i).get_older_z()) + vecAtoms.at(i).get_z_force());
			if((i == fixedatom1) || (i == fixedatom2) || (i == fixedatom3) || (i = fixedatom4)) {
				vecAtoms.at(i).set_new_x_force(vecAtoms.at(i).get_old_x());
				vecAtoms.at(i).set_new_y_force(vecAtoms.at(i).get_old_y());
				vecAtoms.at(i).set_new_z_force(vecAtoms.at(i).get_old_z());
			}
			if(boxon == 1) {
				if(vecAtoms.at(i).get_new_x_force() > boxsize) {
					if(vecAtoms.at(i).get_old_x() > vecAtoms.at(i).get_older_x()) vecAtoms.at(i).set_new_x_force(vecAtoms.at(i).get_old_x() + damping*(vecAtoms.at(i).get_old_x() - vecAtoms.at(i).get_older_x()) + vecAtoms.at(i).get_x_force());
				}
				if(vecAtoms.at(i).get_new_y_force() > boxsize) {
					if(vecAtoms.at(i).get_old_y() > vecAtoms.at(i).get_older_y()) vecAtoms.at(i).set_new_y_force(vecAtoms.at(i).get_old_y() + damping*(vecAtoms.at(i).get_old_y() - vecAtoms.at(i).get_older_y()) + vecAtoms.at(i).get_y_force());
				}
				if(vecAtoms.at(i).get_new_z_force() > boxsize) {
					if(vecAtoms.at(i).get_old_z() > vecAtoms.at(i).get_older_z()) vecAtoms.at(i).set_new_z_force(vecAtoms.at(i).get_old_z() + damping*(vecAtoms.at(i).get_old_z() - vecAtoms.at(i).get_older_z()) + vecAtoms.at(i).get_z_force());
				}
				if(vecAtoms.at(i).get_new_x_force() > -1*boxsize) {
					if(vecAtoms.at(i).get_old_x() < vecAtoms.at(i).get_older_x()) vecAtoms.at(i).set_new_x_force(vecAtoms.at(i).get_old_x() + damping*(vecAtoms.at(i).get_old_x() - vecAtoms.at(i).get_older_x()) + vecAtoms.at(i).get_x_force());
				}
				if(vecAtoms.at(i).get_new_y_force() > -1*boxsize) {
					if(vecAtoms.at(i).get_old_y() < vecAtoms.at(i).get_older_y()) vecAtoms.at(i).set_new_y_force(vecAtoms.at(i).get_old_y() + damping*(vecAtoms.at(i).get_old_y() - vecAtoms.at(i).get_older_y()) + vecAtoms.at(i).get_y_force());
				}
				if(vecAtoms.at(i).get_new_z_force() > -1*boxsize) {
					if(vecAtoms.at(i).get_old_z() < vecAtoms.at(i).get_older_z()) vecAtoms.at(i).set_new_z_force(vecAtoms.at(i).get_old_z() + damping*(vecAtoms.at(i).get_old_z() - vecAtoms.at(i).get_older_z()) + vecAtoms.at(i).get_z_force());
				}
			}
		}
	}
	if((runpointnum % keepevery) == 0) cat_append("g09.log", "dyn"); // cat g09.log >> dyn

	ofstream g09_com;
	g09_com.open("g09.com");
	g09_com << "%nproc=" << processors << endl;
	g09_com << "%mem=" << memory << endl;
	if(killcheck != 1) g09_com << "%chk=" << checkpoint << endl;
	if(nonstandard == 0) {
		g09_com << "# " << method << " force scf=(tight,nosym) " << endl;
		if(meth2 == "unrestricted") g09_com << "guess=mix" << endl;
		if(meth3.size() > 2) g09_com << meth3 << endl;
		if(meth4.size() > 2) g09_com << meth4 << endl;
	}
	else if(nonstandard == 1) {
		g09_com << "# " << endl;
		g09_com << "nonstd" << endl;
		ifstream nonstandard;
		nonstandard.open("nonstandard");
		while(getline(nonstandard, line)) {
			g09_com << line << endl;
		}
		nonstandard.close();
	}
	g09_com << endl;
	g09_com << title1 << " " << title2 << " " << title3 << " " << title4 << endl;
	g09_com << "runpoint  " << runpointnum << endl;
	g09_com << "runisomer  " << isomernum << endl;
	g09_com << endl;
	g09_com << charge << " " << multiplicity << endl;
	g09_com.precision(7);
	g09_com.setf(ios::fixed);
	for(unsigned int i = 0; i < vecAtoms.size(); i++) {
		g09_com << vecAtoms[i].get_atomic_symbol() << " " << vecAtoms.at(i).get_new_x_force() << " " << vecAtoms.at(i).get_new_y_force() << " " << vecAtoms.at(i).get_new_z_force() << endl;
		if((i > highlevel) && (i <= (highlevel+linkatoms))) g09_com << "M H" << endl;
		else if(i > (highlevel+linkatoms)) g09_com << "M" << endl;
	}
	g09_com << endl;
	if(meth5.size() > 2) g09_com << meth5 << endl;
	if(meth6.size() > 2) g09_com << meth6 << endl;
	if(methodfilelines >= 1) {
		ifstream methodfile;
		methodfile.open("methodfile");
		for(unsigned int i = 0; i < methodfilelines; i++) {
			getline(methodfile, line);
			g09_com << line << endl;
		}
		methodfile.close();
	}
	if((nmrtype > 0) && (nmrdo == 1)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod << " nmr=giao geom=check" << endl;
		if(nmrmethod == method) g09_com << "guess=tcheck" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;		
	}
	else if((nmrtype > 1) && (nmrdo == 1)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod2 << " nmr=giao geom=check" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;		
	}
	else if((nmrtype > 2) && (nmrdo == 1)) {
		g09_com << "--link1--" << endl;
		g09_com << "%nproc=" << processors << endl;
		g09_com << "%mem=" << memory << endl;
		g09_com << "%chk=" << checkpoint << endl;
		g09_com << "# " << nmrmethod3 << " nmr=giao geom=check" << endl;
		if(meth7.size() > 2) g09_com << meth7 << endl;
		g09_com << endl;
		g09_com << title1 << "," << title2 << "," << title3 << "," << title4 << endl;
		g09_com << "runpoint ," << runpointnum << endl;
		g09_com << "runisomer ," << isomernum << endl;
		g09_com << endl;
		g09_com << charge << "," << multiplicity << endl;		
	}
	else if((nmrcc == 1) && (nmrdo == 1)) {
		ofstream ZMAT;
		ZMAT.open("ZMAT");
		ZMAT << "CCSD(T) NMR calculation" << endl;
		ZMAT.precision(7);
		ZMAT.setf(ios::fixed);
		for(unsigned int i = 0; i < numAtoms; i++) {
			ZMAT << vecAtoms.at(i).get_atomic_symbol() << " " << vecAtoms.at(i).get_new_x_force() << " " << vecAtoms.at(i).get_new_y_force() << " " << vecAtoms.at(i).get_new_z_force() << endl;	
		}
		ZMAT << endl;
		ZMAT << "*ACES2(CALC=CCSD[T],PROP=NMR,BASIS=dzp" << endl;
		ZMAT << "ABCDTYPE=AOBASIS,TREAT_PERT=SEQUENTIAL,CC_PROG=ECC" << endl;
		ZMAT << "COORD=CARTESIAN" << endl;
		ZMAT << "MEM_UNIT=GB,MEMORY=2)" << endl;
		ZMAT << endl;
		ZMAT.close();
	}
	g09_com << endl << endl << endl;
	g09_com.close();
}

void proggenHP(ProgdynConf& progdyn_conf, string g09_log_filepath) {
	// input(ProgdynConf& progdyn_conf, ifstream tempstangeos, ifstream tempmasses, ifstream tempfreqs, ifstream tempresmass, ifstream tempfrc, ifstream tempmodes, ifstream cannontraj, ifstream temp811, ifsteam tempinputgeos, string g09_log_filepath)
	// output(ofstream diagnostics, ofstream modesread, ofstream maxMove)
	vector<int> disMode; // FIX THIS
	for(unsigned int i = 0; i < 10000; i++) { disMode.push_back(-1); } // FIX THIS
	progdyn_conf.set_disMode(disMode); // FIX THIS
	progdyn_conf.read_from_file();
	string method = progdyn_conf.get_method();
	string charge = progdyn_conf.get_charge();
	string multiplicity = progdyn_conf.get_multiplicity();
	string memory = progdyn_conf.get_memory();
	unsigned int processors =progdyn_conf.get_processors();
	string checkpoint = progdyn_conf.get_checkpoint();
	unsigned int diag = progdyn_conf.get_diag();
	unsigned int initialDis = progdyn_conf.get_initialDis();
	float timestep = progdyn_conf.get_timestep();
	float scaling = progdyn_conf.get_scaling();
	float temp = progdyn_conf.get_temp();
	string searchdir = progdyn_conf.get_searchdir();
	unsigned int classical = progdyn_conf.get_classical();
	unsigned int numimag = progdyn_conf.get_numimag();
	string geometry = progdyn_conf.get_geometry();
	unsigned int highlevel = progdyn_conf.get_highlevel();
	bool boxon = progdyn_conf.get_boxon();
	float boxsize = progdyn_conf.get_boxsize();
	float maxAtomMove = progdyn_conf.get_maxAtomMove();
	float cannonball = progdyn_conf.get_cannonball();
	map<int,string> controlPhase = progdyn_conf.get_controlPhase();
	unsigned int rotationmode = progdyn_conf.get_rotationmode();
	bool DRP = progdyn_conf.get_DRP();
	string title1 = progdyn_conf.get_title1();
	string title2 = progdyn_conf.get_title2();
	string title3 = progdyn_conf.get_title3();
	string title4 = progdyn_conf.get_title4();

	int classicalSpacing=2;
	float zpeGauss=0;
	float zpeGaussK=0;
	float zpePlusE=0;
	float potentialE=0;

	vector<Atom> vecAtoms;

	read_tempstangeos(vecAtoms); // gets x y z position of each atom along with the atomic symbol and sets the atomic mass to default value

	read_tempmasses(vecAtoms); // gets the atomic mass of each atom (may be the same as the default atomic mass)

	if(diag>=1) {
		ofstream diagnostics_file;
		diagnostics_file.open("diagnostics", ios::app);
		diagnostics_file << "********************* starting proggen *********************" << endl
				 << "method,charge,multiplicity,memory" << endl
				 << method << "," << charge << "," << multiplicity << "," << memory << endl
				 << "processors,checkpoint,title1,title2,title3,title4,initialDis,timestep,scaling,temp" << endl
				 << processors << "," << checkpoint << "," << title1 << "," << title2 << "," << title3 << "," << title4 << ","
				 << initialDis << "," << timestep << "," << scaling << "," << temp << endl
				 << "classical,numimag,highlevel,boxon,boxsize,DRP,maxAtomMove,cannonball" << endl
				 << classical << "," << numimag << "," << highlevel << "," << boxon << "," << boxsize << "," << DRP << "," << maxAtomMove << "," << cannonball << endl;
		for(unsigned int i = 0; i < vecAtoms.size(); i++) {
			diagnostics_file << "vecAtoms[" << i << "].atomic_num,vecAtoms[" << i << "].atomic_symbol,vecAtoms[" << i << "].atomic_mass,vecAtoms[" << i << "].get_x_pos,vecAtoms["
					 << i << "].get_y_pos,vecAtoms[" << i << "].get_z_pos" << endl
					 << vecAtoms.at(i).get_atomic_num() << "," << vecAtoms.at(i).get_atomic_symbol() << "," << vecAtoms.at(i).get_atomic_weight() << ","
					 << vecAtoms.at(i).get_x() << "," << vecAtoms.at(i).get_y() << "," << vecAtoms.at(i).get_z() << endl;
		}
		diagnostics_file.close();
	}

	unsigned int numFreq = 3*vecAtoms.size()-6;
	if(geometry == "linear") numFreq = 3*vecAtoms.size()-5;
	vector<float> freq(numFreq), redMass(numFreq), frc(numFreq);

	read_tempfreqs_redmass_frc(numFreq, freq, redMass, frc, scaling);

	if(diag >= 1) {
		ofstream diagnostics_file;
		diagnostics_file.open("diagnostics", ios::app);
		for(unsigned int i = 0; i <= numFreq; i++) {
			diagnostics_file << "freq[" << i << "],redMass[" << i << "],frc[" << i << "]" << endl
					 << freq[i] << "," << redMass[i] << "," << frc[i] << endl;
		}
		diagnostics_file.close();
	}

	map< Tuple<int,int,int>, float > mode;
	if(classical != 2) {
		read_tempmodes(mode, numFreq, vecAtoms.size());
	}

	if(diag > 2) {
		ofstream modesread;
		modesread.open("modesread", ios::app);
		modesread << fixed << setprecision(5);
		for(map< Tuple<int,int,int>, float>::iterator it = mode.begin(); it != mode.end(); ++it) {
			if(it->first.get_first() < numFreq) {
				if(it->first.get_second() == 1) {
					if(it->first.get_third() == 1) modesread << it->second << " ";
					else if(it->first.get_third() == 2) modesread << it->second << " ";
					else if(it->first.get_third() == 3) modesread << it->second << endl;
				}
			}
		}	
	}

	vector< Tuple<float,float,float> > cannonArr;
	string line;
	if(cannonball > 0) {
		ifstream cannontraj;
		cannontraj.open("cannontraj");
		for(unsigned int i = 0; i < vecAtoms.size(); i++) {
			getline(cannontraj, line);
			stringstream ss;
			string discard;
			float x,y,z;
			ss << line;
			ss >> x;
			ss >> y;
			ss >> z;
			Tuple<float,float,float> new_tuple(x,y,z);
			cannonArr.push_back(new_tuple);
		}
		cannontraj.close();
	}

	srand(time(NULL));
	float tester = rand() % 1000;
	vector<float> randArrayA, randArrayB, randArrayC, randArrayD, randArrayE, randArrayR;
	ifstream temp811;
	temp811.open("temp811");
	for(unsigned int i = 0; i < tester; i++) {
		getline(temp811, line);
	}
	for(unsigned int i = 0; i < numFreq; i++) {
		stringstream ss1, ss2, ss3;
		float value;
		getline(temp811, line);
		ss1 << line;
		ss2 >> value;
		randArrayA.push_back(value);
		getline(temp811, line);
		ss2 << line;
		ss2 >> value;
		randArrayB.push_back(value);
		getline(temp811, line);
		ss3 << line;
		ss3 >> value;
		randArrayC.push_back(value);
	}
	for(unsigned int i = 0; i < 6; i++) {
		stringstream ss;
		float value;
		getline(temp811, line);
		ss << line;
		ss >> value;
		randArrayR.push_back(value);
	}

	unsigned int counter = 0;
	while(counter < numFreq) {
		if((initialDis==2) || (disMode[counter]==2)) {
			stringstream ss1, ss2;
			float value;
			getline(temp811, line);
			ss1 << line;
			ss1 >> value;
			float tempNum = 2*(value-0.5);
			float prob = pow(constant_e,(-1*(pow(tempNum,2))));
			getline(temp811, line);
			ss2 << line;
			ss2 >> value;
			if(value < prob) {
				randArrayD.push_back(tempNum);
				counter++;
			}
		}
		else if((initialDis!=2) && (disMode[counter]!=2)) counter++;
	}

	for(unsigned int i = 0; i < vecAtoms.size(); i++) {
		for(unsigned int j = 0; j < 3; j++) {
			stringstream ss;
			float value;
			getline(temp811, line);
			ss << line;
			ss >> value;
			if(value > 0.5) randArrayE.push_back(1);
			else randArrayE.push_back(-1);
		}
	}

	vector<float> zpeJ, zpeK, zpeRat, Q;
	vector<int> vibN;
	for(unsigned int i = 0; i < numFreq; i++) {
		if(classical == 1) zpeJ.push_back(0.5*h*c*classicalSpacing);
		else zpeJ.push_back(0.5*h*c*freq.at(i));
		zpeK.push_back(zpeJ.at(i)*avNum/4184);
		if(temp < 10) vibN.push_back(0.0);
		else {
			zpeRat.push_back(pow(constant_e, (-2*zpeK.at(i))/(RgasK*temp)));
			if(zpeRat.at(i) == 1) zpeRat.at(i) = 0.9999999999;
			Q.push_back(1/(1 - zpeRat.at(i)));
			float newRand = randArrayA.at(i);
			vibN.push_back(0);
			tester = 1/Q.at(i);
			for(unsigned int j = 0; j < (4000*zpeRat.at(i) + 2); j++) {
				if(newRand > tester) vibN.at(i)++;
				tester += (pow(zpeRat.at(i), j)/Q.at(i));
			}
		}
	}

	float desiredModeEnK = 0;
	vector<float> modeEn, modeEnK, maxShift, shift;
	for(unsigned int i = 0; i < numFreq; i++) {
		modeEn.push_back((zpeJ.at(i)*1E18)*(2*vibN.at(i) + 1));
		if(classical == 1) modeEn.at(i) = (zpeJ.at(i)*1E18)*2*vibN.at(i);
		modeEnK.push_back(zpeK.at(i)*(2*vibN.at(i) + 1));
		if(classical == 1) modeEnK.at(i) = zpeK.at(i)*2*vibN.at(i);
		desiredModeEnK += modeEnK.at(i);
		if(freq.at(i) < 10) modeEn.at(i) = (zpeJ.at(i)*1E18)*(2*vibN.at(i));
		maxShift.push_back(pow((2*modeEn.at(i)/frc.at(i)),0.5));
		shift.push_back(0);
		if(initialDis == 3) shift.at(i) = maxShift.at(i)*sin(randArrayC.at(i)*constant_pi*2);
		else if(initialDis == 2) shift.at(i) = maxShift.at(i)*randArrayD.at(i);
		else if(initialDis == 1) shift.at(i) = maxShift.at(i)*(2*(randArrayC.at(i) - 0.5));
		else if(initialDis == 0) shift.at(i) = 0;
		if(disMode.at(i) == 3) shift.at(i) = maxShift.at(i)*sin(randArrayC.at(i)*constant_pi*2);
		else if(disMode.at(i) == 2) shift.at(i) = maxShift.at(i)*randArrayD.at(i);
		else if(disMode.at(i) == 1) shift.at(i) = maxShift.at(i)*(2*(randArrayC.at(i) - 0.5));
		else if(disMode.at(i) == 0) shift.at(i) = 0;
		else if(disMode.at(i) == 10) shift.at(i) = 0; // kept for backward compatability
		if(freq.at(i) < 10) shift.at(i) = 0;
		if(numimag == 1) shift.at(1) = 0;
		else if(numimag == 2) shift.at(2) = 0;
	}

	if(diag > 1) {
		ofstream diagnostics;
		diagnostics.open("diagnostics", fstream::app);
		for(unsigned int i = 0; i < numFreq; i++) {
			if(i == 0) diagnostics << "zpeJ[i],zpeK[i],zpeRat[i],Q[i],vibN[i],modeEn[i],maxShift[i],shift[i]" << endl;
			diagnostics << zpeJ.at(i) << " " << zpeK.at(i) << " " << zpeRat.at(i) << " " << Q.at(i) << " " << vibN.at(i) << " " << modeEn.at(i) << " " << maxShift.at(i) << " " << shift.at(i) << endl;
		}
		diagnostics.close();
	}

	map< Tuple<int,int,int>, float > shiftMode;
	if(classical != 2) {
		for(unsigned int i = 0; i < numFreq; i++) {
			for(unsigned int j = 0; j < vecAtoms.size(); j++) {
				float x_value = shift.at(i), y_value = shift.at(i), z_value = shift.at(i);
				Tuple<int,int,int> x_tuple(i+1,j+1,1);
				Tuple<int,int,int> y_tuple(i+1,j+1,2);
				Tuple<int,int,int> z_tuple(i+1,j+1,3);
				map< Tuple<int,int,int>, float>::iterator it = mode.find(x_tuple); if(it != mode.end()) x_value += it->second;
				it = mode.find(y_tuple); if(it != mode.end()) y_value += it->second;
				it = mode.find(z_tuple); if(it != mode.end()) z_value += it->second;
				shiftMode.insert(pair<Tuple<int,int,int>,float>(x_tuple,x_value));
				shiftMode.insert(pair<Tuple<int,int,int>,float>(y_tuple,y_value));
				shiftMode.insert(pair<Tuple<int,int,int>,float>(z_tuple,z_value));
				vecAtoms.at(j).set_x(vecAtoms.at(j).get_x()+x_value);
				vecAtoms.at(j).set_y(vecAtoms.at(j).get_y()+y_value);
				vecAtoms.at(j).set_z(vecAtoms.at(j).get_z()+z_value);
			}
		}
	}

	vector<float> kinEn, vel;
	ofstream diagnostics;
	diagnostics.open("diagnostics", fstream::app);
	for(unsigned int i = 0; i < numFreq; i++) {
		kinEn.push_back(100000*(modeEn.at(i)-0.5*frc.at(i)*pow(shift.at(i),2)));
		vel.push_back(pow((2*kinEn.at(i)/(redMass.at(i)/avNum)),0.5));
		if(numimag > 1) numimag = 1;
		if(i+1 > numimag) {
			if(randArrayB.at(i) < 0.5) vel.at(i) *= -1;
		}
		if(i+1 == numimag) {
			if(searchdir == "negative") vel.at(i) *= -1;
		}
		if(diag > 1) {
			if(i == 0) diagnostics << "vel[i]" << endl;
			diagnostics << vel.at(i) << endl;
		}
	}
	diagnostics.close();

	for(map<int,string>::iterator it = controlPhase.begin(); it != controlPhase.end(); ++it) {
		if(it->first < numFreq) {
			if((it->second == "positive") && (vel.at(it->first) < 0)) vel.at(it->first) *= -1;
			else if((it->second == "negative") && (vel.at(it->first) > 0)) vel.at(it->first) *= -1;
		}
	}

	map< Tuple<int,int,int>, float > velMode;
	if(classical != 2) {
		for(unsigned int i = 0; i < numFreq; i++) {
			for(unsigned int j = 0; j < vecAtoms.size(); j++) {
				float x_value = 0, y_value = 0, z_value = 0;
				Tuple<int,int,int> x_tuple(i,j+1,1);
				Tuple<int,int,int> y_tuple(i,j+1,2);
				Tuple<int,int,int> z_tuple(i,j+1,3);
				map< Tuple<int,int,int>, float>::iterator it = mode.find(x_tuple);
				if(it != mode.end()) x_value = it->second*vel.at(i)*timestep;
				it = mode.find(y_tuple); if(it != mode.end()) y_value = it->second*vel.at(i)*timestep;
				it = mode.find(z_tuple); if(it != mode.end()) z_value = it->second*vel.at(i)*timestep;
				velMode.insert(pair<Tuple<int,int,int>,float>(x_tuple,x_value));
				velMode.insert(pair<Tuple<int,int,int>,float>(y_tuple,y_value));
				velMode.insert(pair<Tuple<int,int,int>,float>(z_tuple,z_value));
				vecAtoms.at(j).set_x_vel(vecAtoms.at(j).get_x_vel()+x_value);
				vecAtoms.at(j).set_y_vel(vecAtoms.at(j).get_y_vel()+y_value);
				vecAtoms.at(j).set_z_vel(vecAtoms.at(j).get_z_vel()+z_value);
			}
		}
	}

	float degFreedomEnK=0, degFreedomEnJ=0, cartEn=0, kinEnCart=0;
	if(classical == 2) {
		string oldline, discard;
		int atom;
		float pos_4, pos_5, pos_6;
		stringstream ss;
		ifstream tempinputgeos;
		tempinputgeos.open("tempinputgeos");
		while(getline(tempinputgeos, line)) {
			if(oldline == line) oldline = "";
			ss.str(string()); ss.clear();
			ss << line;
			ss >> atom;    // 1
			ss >> discard; // 2
			ss >> discard; // 3
			ss >> pos_4;   // 4
			ss >> pos_5;   // 5
			ss >> pos_6;   // 6
			vecAtoms.at(atom).set_x(pos_4);
			vecAtoms.at(atom).set_y(pos_5);
			vecAtoms.at(atom).set_z(pos_6);
			vecAtoms.at(atom).set_orig_x(pos_4);
			vecAtoms.at(atom).set_orig_y(pos_5);
			vecAtoms.at(atom).set_orig_z(pos_6);
		}
		degFreedomEnK = temp*RgasK;
		degFreedomEnJ = degFreedomEnK/(avNum/4184);
		cartEn = degFreedomEnJ*1E18;
		kinEnCart = 100000*cartEn;
		for(unsigned int i = 0; i < vecAtoms.size(); i++) {
			vecAtoms.at(i).set_x_vel(pow(randArrayE.at(i)*timestep*(2*kinEnCart/(vecAtoms.at(i).get_atomic_weight()/avNum)),0.5)); 
			vecAtoms.at(i).set_y_vel(pow(randArrayE.at(i+1)*timestep*(2*kinEnCart/(vecAtoms.at(i).get_atomic_weight()/avNum)),0.5));
			vecAtoms.at(i).set_z_vel(pow(randArrayE.at(i+2)*timestep*(2*kinEnCart/(vecAtoms.at(i).get_atomic_weight()/avNum)),0.5));
			if(DRP == 1) {
				vecAtoms.at(i).set_x_vel(0);
				vecAtoms.at(i).set_y_vel(0);
				vecAtoms.at(i).set_z_vel(0);
			}
		}
		tempinputgeos.close();
	}

	float KEinitmodes = 0;
	for(unsigned int i = 0; i < vecAtoms.size(); i++) {
		KEinitmodes += 0.5*vecAtoms.at(i).get_atomic_weight()*((pow(vecAtoms.at(i).get_x_vel(),2) + pow(vecAtoms.at(i).get_y_vel(),2) + pow(vecAtoms.at(i).get_z_vel(),2))/(pow(timestep,2)*conver1));
	}

	vector<float> rotateX1, rotateX2, rotateX3, rotateY1, rotateY2, rotateY3, rotateZ1, rotateZ2, rotateZ3;
	float eRotX = 0, eRotY = 0, eRotZ = 0, keRx, keRy, keRz, rotEdesired, scaleX, scaleY, scaleZ;
	int signX = 1, signY = 1, signZ = 1;
	if(rotationmode > 0) {
		for(unsigned int i = 0; i < vecAtoms.size(); i++) {
			rotateX1.push_back(0);
			rotateX2.push_back(-1*vecAtoms.at(i).get_orig_z());
			rotateX3.push_back(vecAtoms.at(i).get_orig_y());
			rotateY1.push_back(-1*vecAtoms.at(i).get_orig_z());
			rotateY2.push_back(0);
			rotateY3.push_back(vecAtoms.at(i).get_orig_x());
			rotateZ1.push_back(-1*vecAtoms.at(i).get_orig_y());
			rotateZ2.push_back(vecAtoms.at(i).get_orig_x());
			rotateZ3.push_back(0);
		}
		for(unsigned int i = 0; i < vecAtoms.size(); i++) {
			eRotX += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateX1.at(i)*rotateX1.at(i))/((timestep*timestep)*conver1);
			eRotX += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateX2.at(i)*rotateX2.at(i))/((timestep*timestep)*conver1);
			eRotX += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateX3.at(i)*rotateX3.at(i))/((timestep*timestep)*conver1);
			eRotY += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateY1.at(i)*rotateY1.at(i))/((timestep*timestep)*conver1);
			eRotY += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateY2.at(i)*rotateY2.at(i))/((timestep*timestep)*conver1);
			eRotY += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateY3.at(i)*rotateY3.at(i))/((timestep*timestep)*conver1);
			eRotZ += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateZ1.at(i)*rotateZ1.at(i))/((timestep*timestep)*conver1);
			eRotZ += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateZ2.at(i)*rotateZ2.at(i))/((timestep*timestep)*conver1);
			eRotZ += 0.5*vecAtoms.at(i).get_atomic_weight()*(rotateZ3.at(i)*rotateZ3.at(i))/((timestep*timestep)*conver1);
		}
		keRx = -0.5*0.001987*temp*log(1-randArrayR.at(0));
		keRy = -0.5*0.001987*temp*log(1-randArrayR.at(1));
		keRz = -0.5*0.001987*temp*log(1-randArrayR.at(2));
		if(eRotX < 1) keRx = 0;
		if(eRotY < 1) keRy = 0;
		if(eRotZ < 1) keRz = 0;
		rotEdesired = keRx + keRy + keRz;
		if(randArrayR.at(3) < 0.5) signX = -1;
		if(randArrayR.at(4) < 0.5) signY = -1;
		if(randArrayR.at(5) < 0.5) signZ = -1;
		if(eRotX < 1) eRotX = 1;
		if(eRotY < 1) eRotY = 1;
		if(eRotZ < 1) eRotZ = 1;
		scaleX = pow((keRx/eRotX),0.5);
		scaleY = pow((keRy/eRotY),0.5);
		scaleZ = pow((keRz/eRotZ),0.5);
		for(unsigned int i = 0; i < vecAtoms.size(); i++) {
			rotateX1.at(i) = rotateX1.at(i)*scaleX*signX;
			rotateX2.at(i) = rotateX2.at(i)*scaleX*signX;
			rotateX3.at(i) = rotateX3.at(i)*scaleX*signX;
			rotateY1.at(i) = rotateY1.at(i)*scaleX*signX;
			rotateY2.at(i) = rotateY2.at(i)*scaleX*signX;
			rotateY3.at(i) = rotateY3.at(i)*scaleX*signX;
			rotateZ1.at(i) = rotateZ1.at(i)*scaleX*signX;
			rotateZ2.at(i) = rotateZ2.at(i)*scaleX*signX;
			rotateZ3.at(i) = rotateZ3.at(i)*scaleX*signX;
		}
		for(unsigned int i = 0; i < vecAtoms.size(); i++) {
			vecAtoms.at(i).set_x_vel(vecAtoms.at(i).get_x_vel()+rotateX1.at(i)+rotateY1.at(i)+rotateZ1.at(i));
			vecAtoms.at(i).set_y_vel(vecAtoms.at(i).get_y_vel()+rotateX2.at(i)+rotateY2.at(i)+rotateZ2.at(i));
			vecAtoms.at(i).set_z_vel(vecAtoms.at(i).get_z_vel()+rotateX3.at(i)+rotateY3.at(i)+rotateZ3.at(i));
		}
	}

	if(cannonball > 0) {
		int multiplier = 1, tester = 0;
		float tolerance = 0.1;
		while(tester==0) {
			float KEinittotal=0;
			for(unsigned int i = 0; i < vecAtoms.size(); i++) {
				
			}
			if(KEinittotal>(KEinitmodes+cannonball+tolerance)) multiplier *= 0.98901364;
			if(KEinittotal<(KEinitmodes+cannonball-tolerance)) multiplier *= 1.01;
			if((KEinittotal<(KEinitmodes+cannonball+tolerance)) && (KEinittotal>(KEinitmodes+cannonball-tolerance))) tester = 1;
			for(unsigned int i = 0; i < vecAtoms.size(); i++) {
				
			}
		}	
	}

	ofstream geoPlusVel;
	geoPlusVel.open("geoPlusVel");
	geoPlusVel << vecAtoms.size() << endl;
	geoPlusVel << fixed;
	for(unsigned int i = 0; i < vecAtoms.size(); i++) {
		geoPlusVel << " " << vecAtoms.at(i).get_atomic_symbol() << "  " << setprecision(7) << vecAtoms.at(i).get_x() << " " << vecAtoms.at(i).get_y() << " " << vecAtoms.at(i).get_z() << "   " << setprecision(6) << vecAtoms.at(i).get_atomic_weight() << endl;
	}
	
	float KEinittotal = 0;
	geoPlusVel << setprecision(7);
	for(unsigned int i = 0; i < vecAtoms.size(); i++) {
		KEinittotal += 0.5*vecAtoms.at(i).get_atomic_weight()*(pow(vecAtoms.at(i).get_x_vel(),2)+pow(vecAtoms.at(i).get_y_vel(),2)+pow(vecAtoms.at(i).get_z_vel(),2))/(pow(timestep,2)*conver1);
		geoPlusVel << vecAtoms.at(i).get_x_vel() << " " << vecAtoms.at(i).get_y_vel() << " " << vecAtoms.at(i).get_z_vel() << endl;
	}

	geoPlusVel << setprecision(6);
	for(unsigned int i = 0; i < numFreq; i++) {
		if(initialDis==0) geoPlusVel << randArrayA.at(i) << "   " << randArrayB.at(i) << "    " << vibN.at(i) << "    " << scientific << vel.at(i) << "       " << fixed << shift.at(i) << "  " << disMode.at(i) << endl;
		else if(initialDis==1) geoPlusVel << randArrayA.at(i) << "   " << randArrayC.at(i) << "    " << vibN.at(i) << "    " << scientific << vel.at(i) << "       " << fixed << shift.at(i) << "  " << disMode.at(i) << endl;
		else if(initialDis==2) geoPlusVel << randArrayA.at(i) << "   " << randArrayD.at(i) << "    " << vibN.at(i) << "    " << scientific << vel.at(i) << "       " << fixed << shift.at(i) << "  " << disMode.at(i) << endl;
		else if(initialDis==4) geoPlusVel << randArrayA.at(i) << "   " << randArrayC.at(i) << "    " << vibN.at(i) << "    " << scientific << vel.at(i) << "       " << fixed << shift.at(i) << "  " << disMode.at(i) << "  " << sin(randArrayC.at(i)*constant_pi*2) << endl;
	}
	geoPlusVel.unsetf(ios_base::floatfield);
	geoPlusVel << "temp  " << temp << endl;
	geoPlusVel << "initialDis " << initialDis << endl;
	geoPlusVel << "classical " << classical << endl;
	geoPlusVel << scientific << setprecision(0) << "timestep " << timestep << endl;
	geoPlusVel << "numimag " << numimag << endl;
	geoPlusVel << fixed << setprecision(3) << "Total mode energy desired= " << desiredModeEnK << endl;
	geoPlusVel << "KE initial from modes= " << KEinitmodes << "   KE initital total= " << KEinittotal << "   Rotational Energy desired= " << rotEdesired << endl;
	if(cannonball>0) geoPlusVel << "cannonball " << cannonball << "   cannon Energy= " << KEinittotal-KEinitmodes << endl;
	if(boxon>0) geoPlusVel << "boxsize " << boxsize << endl;
	if(DRP>0) geoPlusVel << "DRP " << DRP << "   maxAtomMove " << maxAtomMove << endl;

	ifstream g09_log;
	g09_log.open(g09_log_filepath.c_str());
	while(getline(g09_log,line)) {
		if(line.find("Zero-point correction") != string::npos) {
			stringstream ss;
			string discard;
			ss << line;
			for(unsigned int i = 0; i < 2; i++) { ss >> discard; }
			ss >> zpeGauss;
		}
		if(line.find("zero-point Energies") != string::npos) {
			stringstream ss;
			string discard;
			ss << line;
			for(unsigned int i = 0; i < 6; i++) { ss >> discard; }
			ss >> zpePlusE;
		}
	}
	g09_log.close();

	zpeGaussK = zpeGauss*627.509;
	potentialE = zpePlusE - zpeGauss;
	geoPlusVel << setprecision(6) << "Gaussian zpe= " << zpeGauss << " or " << zpeGaussK << " kcal/mol  E + zpe= " << zpePlusE << "   potential E= " << potentialE << endl;
	geoPlusVel << endl;
	geoPlusVel.close();
}

void progcfour() {
	cout << "You need to define progcfour before it can be used" << endl;
	return;
}

int main() {
	const string LOG_FILE = "docslog";
	unsigned int iteration_number = 0;

// File System Setup
	const string LC_ALL = "C";
	const string JOB_NAME = "freqinHP";
	const string TEMPORARY_DIR = "$HOME/ProgdynSuite/RunOutputs/temp1";
	const string GAUSS_SCRDIR = TEMPORARY_DIR + "/temporary_files";
	const string PROG_SCRDIR = TEMPORARY_DIR + "/prog_files";
	const string PROG_HOME = TEMPORARY_DIR;
	const string G09_ROOT = "/blueapps/chem";
	const string RAND_DIR = TEMPORARY_DIR;
	const string COM_LOG_FILES = "FIX THIS";
	//. $G09_ROOT/gaussian/g09/bsd/g09.profile

	unsigned int runpointnumber = 1;
	unsigned int isomernumber = 1;
	bool nogo = false;
	bool goingwell = false;
	string skipstart = "";
	bool bypassproggen = false;
	remove("diagnostics");
	bool detour = false;
	bool tempdone = false;

	ProgdynConf progdyn_conf;
	progdyn_conf.read_from_file();
	
	chdir(TEMPORARY_DIR.c_str());
	
	while (true) { // start while::A
		goingwell = false;
		while (true) { // start while::B
			if (skipstart == "") cout << "skipping start and continuing from previous runs" << endl;
			else { // generate geoPlusVel and first input file
				if(runpointnumber == 1) {
					append_to_dynfollowfile("X did not complete first point so new isomer started");
				}
				else if(runpointnumber == 2) {
					append_to_dynfollowfile("X did not complete second point so new isomer started");
				}
				else if(runpointnumber == 3) {
					append_to_dynfollowfile("X did not complete third point so new isomer started");
				}
				if(bypassproggen) cout << "Taking starting conditions from pre-generated geoPlusVel" << endl;
				else {
					create_temp_files(JOB_NAME);
					proggenHP(progdyn_conf, JOB_NAME);
				}
				if(isomernumber) isomernumber++;
				else isomernumber = 1;
				runpointnumber = 1;
				remove("g09.com");
				prog1stpoint(progdyn_conf, isomernumber, runpointnumber);
				if(file_exists("g09.com")) {
					remove("tempfreqs");
					remove("tempredmass");
					remove("tempfrc");
					remove("tempmodes");
					remove("tempstangeos");
					remove("tempmasses");
					remove("temp401");
					remove("temp811");
					remove("tempinputgeos");
					append_to_geoRecord(to_string(isomernumber) + " ----trajectory isomer number----");
					cat_append("geoPlusVel", "geoRecord");
					goingwell = false;
					cp(TEMPORARY_DIR + "/g09.com", PROG_SCRDIR + "/g09.com");
					cp(PROG_SCRDIR + "/g09.com", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".com");
					system((G09_ROOT + "/bin/rung09 " + PROG_SCRDIR + "/g09.com > " + PROG_SCRDIR + "/g09.log").c_str());
					cp(PROG_SCRDIR + "/g09.log", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".log");
					iteration_number++;
					goingwell = grep_Normal_termination(PROG_SCRDIR + "/g09.log");
					if(goingwell) {
						cat_append(PROG_SCRDIR + "/g09.log", "dyn");
						cp(PROG_SCRDIR + "/g09.log", "olderdynrun");
					}
					else cp(PROG_SCRDIR + "/g09.log", TEMPORARY_DIR + "/g09.log");
					break;
				}
				else break;
				remove("g09.com");
				runpointnumber = 2;
				prog2ndpoint(progdyn_conf, isomernumber, runpointnumber, skipstart, TEMPORARY_DIR + "/g09.log");
				proganal(isomernumber, runpointnumber, skipstart, nogo, TEMPORARY_DIR + "/g09.log");
				tempdone = bool_find_XXXX_last_n_lines("dynfollowfile", 10); // tail -1 dynfollowfile | awk '/XXXX/ {print}' > $PROG_SCRDIR/tempdone
				if (tempdone) {
					remove("dyn");
					remove("traj");
					runpointnumber = 0;
					break;
				}
				if(file_exists("g09.com")) {
					goingwell = false;
					cp(TEMPORARY_DIR + "/g09.com", PROG_SCRDIR + "/g09.com");
					cp(PROG_SCRDIR + "/g09.com", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".com");
					system((G09_ROOT + "/bin/rung09 " + PROG_SCRDIR + "/g09.com > " + PROG_SCRDIR + "/g09.log").c_str());
					cp(PROG_SCRDIR + "/g09.log", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".log");
					iteration_number++;
					goingwell = grep_Normal_termination(PROG_SCRDIR + "/g09.log");
					if(goingwell) {
						cp(PROG_SCRDIR + "/g09.log", "olddynrun");
						cat_append(PROG_SCRDIR + "/g09.log", "dyn");
						proganal(isomernumber, runpointnumber, skipstart, nogo, TEMPORARY_DIR + "/g09.log");
						create_old(); //awk '/Input orientation/,/Distance matrix/ {print}' olddynrun | awk '/   0   / {print}' > old
						create_older(); //awk '/Input orientation/,/Distance matrix/ {print}' olderdynrun | awk '/   0   / {print}' > older
						runpointnumber = 3;
						progdynb(progdyn_conf, isomernumber, runpointnumber, PROG_SCRDIR +"/g09.log");
						remove("old");
						remove("older");
					}
					else {
						cp(PROG_SCRDIR + "/g09.log", TEMPORARY_DIR + "/g09.log");
						break;
					}
				}
				else break;
				skipstart = "forward";
			}
			if(skipstart == "reverserestart") {
				remove("g09.com");
				runpointnumber = 1;
				prog1stpoint(progdyn_conf, isomernumber, runpointnumber);
				if(file_exists("g09.com")) {
					goingwell = false;
					cp(TEMPORARY_DIR + "/g09.com", PROG_SCRDIR + "/g09.com");
					cp(PROG_SCRDIR + "/g09.com", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".com");
					system((G09_ROOT + "/bin/rung09 " + PROG_SCRDIR + "/g09.com > " + PROG_SCRDIR + "/g09.log").c_str());
					cp(PROG_SCRDIR + "/g09.log", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".log");
					iteration_number++;
					goingwell = grep_Normal_termination(PROG_SCRDIR + "/g09.log");
					if(goingwell) {
						cp(PROG_SCRDIR + "/g09.log", "olderdynrun");
					}
					else {
						cp(PROG_SCRDIR + "/g09.log", TEMPORARY_DIR + "/g09.log");
						break;
					}
				}
				else {
					break;
				}
				remove("g09.com");
				runpointnumber = 2;
				prog2ndpoint(progdyn_conf, isomernumber, runpointnumber, skipstart, TEMPORARY_DIR + "/g09.log");
				proganal(isomernumber, runpointnumber, skipstart, nogo, TEMPORARY_DIR + "/g09.log");
				tempdone = false;
				if(file_exists("g09.com")) {
					goingwell = false;
					cp(TEMPORARY_DIR + "/g09.com", PROG_SCRDIR + "/g09.com");
					cp(PROG_SCRDIR + "/g09.com", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".com");
					system((G09_ROOT + "/bin/rung09 " + PROG_SCRDIR + "/g09.com > " + PROG_SCRDIR + "/g09.log").c_str());
					cp(PROG_SCRDIR + "/g09.log", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".log");
					iteration_number++;
					goingwell = grep_Normal_termination(PROG_SCRDIR + "/g09.log");
					if(goingwell) {
						cp(PROG_SCRDIR + "/g09.log", "olddynrun");
						cat_append(PROG_SCRDIR + "/g09.log", "dyn");
						proganal(isomernumber, runpointnumber, skipstart, nogo, TEMPORARY_DIR + "/g09.log");
						create_old(); //awk '/Input orientation/,/Distance matrix/ {print}' olddynrun | awk '/   0   / {print}' > old
						create_older(); //awk '/Input orientation/,/Distance matrix/ {print}' olderdynrun | awk '/   0   / {print}' > older
						runpointnumber = 3;
						progdynb(progdyn_conf, isomernumber, runpointnumber, PROG_SCRDIR +"/g09.log");
						remove("old");
						remove("older");
					}
					else {
						cp(PROG_SCRDIR + "/g09.log", TEMPORARY_DIR + "/g09.log");
						break;
					}
				}
				else {
					break;
				}
				skipstart = "reverse";
			}
			while(true) { // start while::A
				runpointnumber++;
				goingwell = false;
				cp(TEMPORARY_DIR + "/g09.com", PROG_SCRDIR + "/g09.com");
				cp(PROG_SCRDIR + "/g09.com", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".com");
				system((G09_ROOT + "/bin/rung09 " + PROG_SCRDIR + "/g09.com > " + PROG_SCRDIR + "/g09.log").c_str());
				cp(PROG_SCRDIR + "/g09.log", COM_LOG_FILES + "/g09" + to_string(iteration_number) + ".log");
				iteration_number++;
				goingwell = grep_Normal_termination(PROG_SCRDIR + "/g09.log");
				if(goingwell) {
					proganal(isomernumber, runpointnumber, skipstart, nogo, TEMPORARY_DIR + "/g09.log"); //awk -f $TEMPORARY_DIR/proganal $PROG_SCRDIR/g09.log >> $TEMPORARY_DIR/dynfollowfile
					mv("olddynrun", "olderdynrun");
					create_old(); //awk '/Input orientation/,/Distance matrix/ {print}' $PROG_SCRDIR/g09.log | awk '/   0   / {print}' > old
					cp(PROG_SCRDIR + "/g09.log", "olddynrun");
					create_older(); //awk '/Input orientation/,/Distance matrix/ {print}' olderdynrun | awk '/   0   / {print}' > older
					progdynb(progdyn_conf, isomernumber, runpointnumber, PROG_SCRDIR +"/g09.log"); //awk -f $TEMPORARY_DIR/progdynb $PROG_SCRDIR/g09.log > g09.com
					remove("old");
					remove("older");
				}
				else {
					cp(PROG_SCRDIR + "/g09.log", TEMPORARY_DIR + "/g09.log");
					break;
				}
				if(file_exists("ZMAT")) {
					cp("ZMAT", PROG_SCRDIR.c_str());
					progcfour(); // $PROG_SCRDIR/progcfour $TEMPORARY_DIR $PROG_SCRDIR
					mv("ZMAT", "temp.ZMAT");
					print_to_NMRlistcc("generic one two three " + to_string(runpointnumber) + " runisomer " + to_string(isomernumber));
					process_x_log(); // awk '/Nuclear Magnetic Resonance/,/HF-SCF/ {if ($2=="C") print $1,$2,"Isotropic =",$3; if ($2=="H") print $1,$2,"Isotropic =",$3}' x.log >> NMRlistcc
				}
				if(detour) {
					detour = false;
					print_date_to_file(LOG_FILE); //date >> $LOG_FILE
					cat_append("run.com", LOG_FILE);
					cp("run.log", "temp.log");
					system((G09_ROOT + "/bin/rung09 " + TEMPORARY_DIR + "/run.com > " + TEMPORARY_DIR + "/run.log").c_str());
				}
				if(nogo) break;
				tempdone = false;
				tempdone = bool_find_XXXX_last_n_lines("dynfollowfile", 2); // tail -2 dynfollowfile | awk '/XXXX/ {print}' > $PROG_SCRDIR/tempdone
				if(tempdone) {
					if(progdyn_conf.get_reversetraj())
						if(skipstart == "reverse") {
							skipstart = "";
							remove("geoPlusVel");
							remove("olddynrun");
							remove("olderdynrun");
							mv("traj","traj"+to_string(isomernumber));
							mv("dyn","dyn"+to_string(isomernumber));
						}
						else if(skipstart == "forward" ) skipstart = "reverserestart";
					else {
						skipstart = "";
						remove("geoPlusVel");
						remove("olddynrun");
						remove("olderdynrun");
						mv("traj","traj"+to_string(isomernumber));
						mv("dyn","dyn"+to_string(isomernumber));
					}
					break;
				}
			} // end while::C
			if(nogo) break;
			if(goingwell) {
				cout << "starting a new point or a new direction" << endl;
			}
			else {
				break;
			}
		} // end while::B
		if(nogo) break;
		if(goingwell) cout << "starting a new point or a new direction2" << endl;
		else break;
	} // end while::A
	return 0;
}
