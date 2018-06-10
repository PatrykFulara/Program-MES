#pragma once

class Node
{
public:
	double x, y, temperature;
	double localH[2][4];
	bool status; //Czy le¿y na krawêdzi
};
