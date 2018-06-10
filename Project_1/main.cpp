#include <iostream>
#include "GlobalData.h"
#include "Grid.h"

int main()
{
	GlobalData *globalData = new GlobalData();
	Grid grid(globalData);
	system("PAUSE");
}