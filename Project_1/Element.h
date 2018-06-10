#pragma once
#include <vector>
using namespace std;

class Element
{
public:
	vector<int> ID;
	vector<bool> pow;
	double NpoDx[4][4]; //Pochodne funkcji kszta³tu po x
	double NpoDy[4][4]; //Pochodne funkcji kszta³tu po y
	double hLocal[4][4];
	double cLocal[4][4];
	double pLocal[4];
	Element()
	{
		ID.resize(4);
		pow.resize(4);
	}
};