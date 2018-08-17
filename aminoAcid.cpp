#include "aminoAcid.h"

aminoAcid::aminoAcid(){
}

aminoAcid::aminoAcid(string csv) {
	
	deu = sod = car = 0;
	int pos = csv.find(",");
	/*
	name = csv.substr(0, pos);
	csv = csv.substr(pos + 1);
	pos = csv.find(",");

	adduct = csv.substr(0, pos);
	csv = csv.substr(pos + 1);
	pos = csv.find(",");

	fomular = csv.substr(0, pos);
	csv = csv.substr(pos + 1);
	pos = csv.find(",");
	*/
	derFomular = csv.substr(0, pos);
	csv = csv.substr(pos + 1);
	pos = csv.find(",");
	
	char * mptr = NULL;
	char * mptr2 = NULL;
	char * eptr = NULL;
	char * eptr2 = NULL;
	map<char, int> isotop;
	char estr[20] = "";
	char mstr[20] = "";
	strcpy_s(estr, derFomular.c_str());
	strcpy_s(mstr, derFomular.c_str());
	eptr = strtok_s(estr, "1234567890+", &eptr2);
	mptr = strtok_s(mstr, "ABCDEFGHIJKLMNOPQRSTUVWXYZ+", &mptr2);
	while (mptr != NULL) {
		if (eptr[0] == 'S')
			break;
		isotop[eptr[0]] = stoi(mptr);
		eptr = strtok_s(NULL, "1234567890+", &eptr2);
		mptr = strtok_s(NULL, "ABCDEFGHIJKLMNOPQRSTUVWXYZ+", &mptr2);
	}

	/*
	mw = atof(csv.substr(0, pos).c_str());
	csv = csv.substr(pos + 1);
	pos = csv.find(",");

	H = atof(csv.substr(0, pos).c_str());
	csv = csv.substr(pos + 1);
	pos = csv.find(",");

	Na = atof(csv.substr(0, pos).c_str());
	csv = csv.substr(pos + 1);
	pos = csv.find(",");
	*/

	for (int c = 0; c <= numCar; c++) {
		for (int d = 0; d <= numDeu; d++) {
			for (int n = 0; n <= numSod; n++) {
				if (n > isotop['N'] || d > isotop['H'] || c > isotop['C'])
					isotopes[d][n][c] = 0;
				else
					isotopes[d][n][c] = atof(csv.substr(0, pos).c_str());
				csv = csv.substr(pos + 1);
				pos = csv.find(",");
			}
		}
	}
}

aminoAcid::~aminoAcid(){
}

double aminoAcid::search(double m, double error) {
	double smallest = 1000;
	for (int d = 0; d <= numDeu; d++) {
		for (int n = 0; n <= numSod; n++) {
			for (int c = 0; c <= numCar; c++) {
				if (isotopes[d][n][c] == 0)
					continue;
				// need to weed out all of the impossible isotopes
				if (fabs(isotopes[d][n][c] - m) < smallest) {
					smallest = fabs(isotopes[d][n][c] - m);
					deu = d;
					sod = n;
					car = c;
				}
			}
		}
	}
	return smallest;
}

string aminoAcid::getName() {
	return name;
}

string aminoAcid::getDerFomular() {
	return derFomular;
}

double aminoAcid::getMass() {
	return isotopes[deu][sod][car];
}

double aminoAcid::getMass(int a, int b, int c) {
	return isotopes[a][b][c];
}

void aminoAcid::print() {
	cout << name << " || " << isotopes[deu][sod][car] << " || " << "D" << deu << "Na" << sod << "C" << car << endl;
}

string aminoAcid::printcsv() {
	return name + "," + adduct + "," + fomular + "," + derFomular + "," + to_string(mw) + "," + to_string(H) + "," + to_string(Na) + "," + to_string(isotopes[deu][sod][car]) + "," + to_string(deu) + "," + to_string(sod) + "," + to_string(car) + "\n";
}

void aminoAcid::retIsotopes(double myIsotopes[numDeu + 1][numSod + 1][numCar + 1]) {
	for (int a = 0; a <= numDeu; a++)
		for (int b = 0; b <= numSod; b++)
			for (int c = 0; c <= numCar; c++)
				myIsotopes[a][b][c] = isotopes[a][b][c];
}