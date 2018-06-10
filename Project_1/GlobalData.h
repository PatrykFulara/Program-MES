#pragma once
#include <fstream>
using namespace std;

class GlobalData
{
public:
	double height, width, condK, cp, densityV, ambientTemp, alfa, simTime, stepTime, initTemp;
	int heightsNodes, widthNodes, numberOfNodes, numberOfElements;
	GlobalData();

private:
	fstream data;
	void countTotalNumberOfNodesAndElements();
};