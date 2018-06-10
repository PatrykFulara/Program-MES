#pragma once
#include <math.h>
#include<stdio.h>
#include<iostream>
#include "GlobalData.h"
#include "Node.h"
#include "Element.h"

using namespace std;

class Grid
{
public:
	vector<Node> nodes; //Wektor wêz³ów
	vector<Element> elements; //Wektor elementów
	Grid(GlobalData *globalData);

private:

	double **hGlobal;
	double **cGlobal;
	double *pGlobal;
	double N1(double ksi, double eta);
	double N2(double ksi, double eta);
	double N3(double ksi, double eta);
	double N4(double ksi, double eta);
	double Matrix();
	int _numberOfNodes;
	int _numberOfElements;
	double _height;
	double _width;
	double _condK;
	double _cp;
	double _densityV;
	double _ambientTemp;
	double _alfa;
	double _simTime;
	double _stepTime;
	double _initTemp;
	int _heightsNodes; 
	int _widthNodes;
	void setNodes();
	void setElements();
	void Solve();
	bool tempVec(int n, double ** A, double * B, double *& X);
	double max(double *vec, int n);
	double min(double *vec, int n);
};