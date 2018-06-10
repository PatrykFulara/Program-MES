#include "GlobalData.h"
#include <iostream>

GlobalData::GlobalData()
{
	data.open("data.txt", ios::in);
	if (data.good())
	{
		data >> height; //D�ugo�� [m]
		data >> width; //Szeroko�� [m]
		data >> heightsNodes; //Ilo�� w�z��w w pionie
		data >> widthNodes; //Ilo�� w�z��w w poziomie
		data >> condK; //Wsp�czynnik przewodzenia
		data >> cp; //Ciep�o w�a�ciwe
		data >> densityV; //G�sto�� 
		data >> ambientTemp; //Temperatura otoczenia
		data >> alfa; //Wsp�czynnik przenikania
		data >> simTime; //Granica czasu [s]
		data >> stepTime; //Krok czasowy [s]
		data >> initTemp; //Temperatura pocz�tkowa [C]
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
