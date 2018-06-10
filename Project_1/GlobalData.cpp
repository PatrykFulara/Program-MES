#include "GlobalData.h"
#include <iostream>

GlobalData::GlobalData()
{
	data.open("data.txt", ios::in);
	if (data.good())
	{
		data >> height; //D³ugoœæ [m]
		data >> width; //Szerokoœæ [m]
		data >> heightsNodes; //Iloœæ wêz³ów w pionie
		data >> widthNodes; //Iloœæ wêz³ów w poziomie
		data >> condK; //Wspó³czynnik przewodzenia
		data >> cp; //Ciep³o w³aœciwe
		data >> densityV; //Gêstoœæ 
		data >> ambientTemp; //Temperatura otoczenia
		data >> alfa; //Wspó³czynnik przenikania
		data >> simTime; //Granica czasu [s]
		data >> stepTime; //Krok czasowy [s]
		data >> initTemp; //Temperatura pocz¹tkowa [C]
		data.clear();
		data.close();
	}

	else
		cout << "Error reading file: data.txt" << endl;
	
	countTotalNumberOfNodesAndElements();
}

void GlobalData::countTotalNumberOfNodesAndElements() 
{
	numberOfNodes = heightsNodes * widthNodes;
	numberOfElements = (heightsNodes - 1) * (widthNodes - 1);
}
