#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <map>

#define numDeu 15
#define numSod 15
#define numCar 15

#pragma once
using namespace std;

class aminoAcid
{
public:
	aminoAcid();
	aminoAcid(string csv);
	~aminoAcid();
	double search(double m, double error);
	string getName();
	string getDerFomular();
	double getMass();
	double getMass(int a, int b, int c);
	void print();
	string printcsv();
	void retIsotopes(double myIsotopes[numDeu + 1][numSod + 1][numCar + 1]);
private:
	int deu = 0, sod = 0, car = 0;
	string name;
	string adduct;
	string fomular;
	string derFomular;
	double mw;
	double H;
	double Na;
	double isotopes[numDeu + 1][numSod + 1][numCar + 1];
};

