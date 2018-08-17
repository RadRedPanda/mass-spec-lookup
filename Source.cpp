#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include "aminoAcid.h"

using namespace std;

#define refSheetPath ""
#define skipHeader 0
#define refSheetFileName "Database.csv"
//#define refSheetFileName "2016_10_4_ECF_mass list_Ye_Full.csv"
#define fileIn "newInput/PC9/" + to_string(sheet) + "_original exported out.csv"
#define fileOut "newInput/PC9/output/" + to_string(err).substr(2, 4) + "/" + to_string(sheet) + "_original exported out.csv"
#define errorLoop double err = 0.0001; err <= 0.0005; err += 0.0001

//vector<string> fileErr = {"0", "0001", "00025", "0005"};

struct data {
	double mass;
	double intensity;
	string extra = "";
};

struct isoForm {
	int r;
	int a;
	int b;
	int c;
	double m;
};

bool compareIso(isoForm i, isoForm j) {
	return i.m < j.m;
}

isoForm getClosest(isoForm val1, isoForm val2, double target, int &index)
{
	if (target - val1.m >= val2.m - target) {
		index++;
		return val2;
	}
	else
		return val1;
}

isoForm findClosest(vector<isoForm> arr, int n, double target, int &index)
{
	// Corner cases
	if (target <= arr[0].m) {
		index = 0;
		return arr[0];
	}
	if (target >= arr[n - 1].m) {
		index = n - 1;
		return arr[n - 1];
	}

	// Doing binary search
	int i = 0, j = n, mid = 0;
	while (i < j) {
		mid = (i + j) / 2;

		if (arr[mid].m == target) {
			index = mid;
			return arr[mid];
		}

		/* If target is less than array element,
		then search in left */
		if (target < arr[mid].m) {

			// If target is greater than previous
			// to mid, return closest of two
			if (mid > 0 && target > arr[mid - 1].m) {
				index = mid - 1;
				return getClosest(arr[mid - 1], arr[mid], target, index);
			}

			/* Repeat for left half */
			j = mid;
		}

		// If target is greater than mid
		else {
			if (mid < n - 1 && target < arr[mid + 1].m) {
				index = mid;
				return getClosest(arr[mid], arr[mid + 1], target, index);
			}
			// update i
			i = mid + 1;
		}
	}

	// Only single element left after search
	return arr[mid];
}

//	reads the reference sheet and puts the data in a vector of amino acids
vector<aminoAcid> readReferenceSheet(string refPath, int sh) {
	vector<aminoAcid> refDat;
	ifstream referenceFile;
	referenceFile.open(refPath);
	if (referenceFile.is_open()) {
		string line;
		for(int i=0; i<sh; i++)
			getline(referenceFile, line);
		while (!referenceFile.eof()) {
			aminoAcid temp(line);
			refDat.push_back(temp);
			getline(referenceFile, line);
		}
	}
	return refDat;
}

vector<isoForm> sortDat(vector<aminoAcid> refDat) {
	vector<isoForm> sortedIso;
	for (int row = 0; row < (int)refDat.size(); row++)
		for (int a = 0; a < 7; a++)
			for (int b = 0; b < 5; b++)
				for (int c = 0; c < 21; c++) {
					isoForm temp = { row, a, b, c, refDat[row].getMass(a, b, c) };
					if (temp.m > 0) {
						sortedIso.push_back(temp);
					}
				}
	return sortedIso;
}

int main() {
	string referencePath = refSheetPath;
	referencePath += refSheetFileName;
	cout << "Reference Path: " << referencePath << endl;
	cout << "Reading Reference Sheet...";

	vector<aminoAcid> referenceData = readReferenceSheet(referencePath, 2);
	cout << "DONE" << endl;

	cout << "Parsing Data...";
	vector<isoForm> sortedIso = sortDat(referenceData);
	cout << "DONE" << endl;

	cout << "Sorting reference sheet...";
	sort(sortedIso.begin(), sortedIso.end(), compareIso);
	cout << "DONE" << endl;

	//for (int a = 0; a < fileErr.size(); a++) {
		for (errorLoop) {
			cout << "Error: " << err << endl;
			for (int sheet = 0; sheet <= 13; sheet++) {
				string filePath = fileIn;
				cout << "Reading " << filePath << "...";

				ifstream file;
				file.open(filePath);
				vector<string> list;
				vector<struct data> values;
				if (file.is_open()) {
					string line;
					for (int i = 0; i < skipHeader; i++)
						getline(file, line);
					while (!file.eof()) {
						getline(file, line);
						int pos = line.find(",");
						if (pos != -1) {
							int pos2 = line.substr(pos + 1).find(",");
							struct data temp;
							temp.mass = atof(line.substr(0, pos).c_str());
							temp.intensity = atof(line.substr(pos + 1, pos + 1 + pos2).c_str());
							if (pos2 != -1) {
								temp.extra = line.substr(pos + 1 + pos2 + 1);
							}
							values.push_back(temp);
							list.push_back(line);
						}
					}
				}
				file.close();

				cout << "DONE" << endl;

				vector<aminoAcid> refTot;
				vector<string> csvOutput;
				//csvOutput.push_back("mass,first,second,third\n");

				cout << "Starting comparisons...";
				for (int i = 0; i < (int)values.size(); i++) {
					int index = -1, secIndex, thiIndex;
					isoForm closest = findClosest(sortedIso, sortedIso.size(), values[i].mass, index);
					if (fabs(closest.m - values[i].mass) <= err) {
						string output = to_string(values[i].mass) + "," + values[i].extra + ",";
						output += referenceData[closest.r].getDerFomular() + "+D" + to_string(closest.a) + "+13C" + to_string(closest.c) + "+15N" + to_string(closest.b) + ";" + to_string(closest.m - values[i].mass) + ":" + to_string(closest.m) + ", ";

						//edges
						if (index == 0) {
							secIndex = 1;
							thiIndex = 2;
						}
						else if (index == sortedIso.size() - 1) {
							secIndex = index - 1;
							thiIndex = index - 2;
						}
						//4 cases
						else if (fabs(sortedIso[index - 1].m - values[i].mass) > fabs(sortedIso[index + 1].m - values[i].mass)) {
							secIndex = index + 1;
							if (fabs(sortedIso[index - 1].m - values[i].mass) > fabs(sortedIso[index + 2].m - values[i].mass))
								thiIndex = index + 2;
							else
								thiIndex = index - 1;
						}
						else {
							secIndex = index - 1;
							if (fabs(sortedIso[index - 2].m - values[i].mass) > fabs(sortedIso[index + 1].m - values[i].mass))
								thiIndex = index + 1;
							else
								thiIndex = index - 2;
						}
						if (fabs(sortedIso[secIndex].m - values[i].mass) <= err)
							output += referenceData[sortedIso[secIndex].r].getDerFomular() + "+D" + to_string(sortedIso[secIndex].a) + "+13C" + to_string(sortedIso[secIndex].c) + "+15N" + to_string(sortedIso[secIndex].b) + ";" + to_string(sortedIso[secIndex].m - values[i].mass) + ":" + to_string(sortedIso[secIndex].m);
						output += ",";
						if (fabs(sortedIso[thiIndex].m - values[i].mass) <= err)
							output += referenceData[sortedIso[thiIndex].r].getDerFomular() + "+D" + to_string(sortedIso[thiIndex].a) + "+13C" + to_string(sortedIso[thiIndex].c) + "+15N" + to_string(sortedIso[thiIndex].b) + ";" + to_string(sortedIso[thiIndex].m - values[i].mass) + ":" + to_string(sortedIso[thiIndex].m);
						output += "\n";
						csvOutput.push_back(output);
					}
					else {
						string output = to_string(values[i].mass) + ",q," + values[i].extra + "\n";
						csvOutput.push_back(output);
					}
				}
				cout << "DONE" << endl;

				ofstream threeBestFile;
				string tbfPath = fileOut;

				cout << "Uploading to " << tbfPath << "...";

				threeBestFile.open(tbfPath);
				for (int i = 0; i < (int)csvOutput.size(); i++)
					threeBestFile << csvOutput[i];
				threeBestFile.close();

				cout << "DONE" << endl;

			}
		}
	//}

	system("PAUSE");
	return 0;
}