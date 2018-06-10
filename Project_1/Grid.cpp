#include "Grid.h"

Grid::Grid(GlobalData *globalData)
{
	_numberOfNodes = globalData->numberOfNodes;
	_numberOfElements = globalData->numberOfElements;
	_height = globalData->height;
	_width = globalData->width;
	_heightsNodes = globalData->heightsNodes;
	_widthNodes = globalData->widthNodes;
	_condK = globalData->condK;
	_cp = globalData->cp;
	_densityV = globalData->densityV;
	_ambientTemp = globalData->ambientTemp;
	_alfa = globalData->alfa;
	_simTime = globalData->simTime;
	_stepTime = globalData->stepTime;
	_initTemp = globalData->initTemp;
	nodes.resize(_numberOfNodes);	//Ustalenie wielkoœci wektora wêz³ów
	elements.resize(_numberOfElements); //Ustalenie wielkoœci wektora elementów
	setNodes();
	setElements();

	//Inicjalizacja i zerowanie macierzy H,C i wektora P
	hGlobal = new double*[_numberOfNodes];
	cGlobal = new double*[_numberOfNodes];
	pGlobal = new double[_numberOfNodes];
	for (int i = 0; i < _numberOfNodes; i++) {
		hGlobal[i] = new double[_numberOfNodes];
		cGlobal[i] = new double[_numberOfNodes];
		for (int j = 0; j < _numberOfNodes; j++) {
			hGlobal[i][j] = 0.0;
			cGlobal[i][j] = 0.0;

		}
	}
	for (int i = 0; i < _numberOfNodes; i++) {
		pGlobal[i] = 0.0;
	}

	Matrix();
	Solve();
}

//Ustalenie wartoœci X,Y,status dla wêz³ów
void Grid::setNodes()
{
	for (int i = 0, column = 0; i < _numberOfNodes; i++) 
	{
		if (i % _heightsNodes == 0 && i != 0)
			column++;

		if (i % _heightsNodes == 0){ 
			nodes[i].y = 0;
			nodes[i].status = true;

		} else { 
			nodes[i].y = (_height / ((double)_heightsNodes - 1)) * (i % _heightsNodes); 
		}
		if (i == 0) {
			nodes[i].x = 0;
		} else {
			nodes[i].x = (_width / ((double)_widthNodes - 1)) * column;
		}

		if (column == 0 || column == (_widthNodes - 1) || (i % _heightsNodes == _heightsNodes - 1) )
			nodes[i].status = true;
	}
}

//Ustalenie ID wêz³ów dla ka¿dego elementu i powierzchni brzegowych
void Grid::setElements()
{
	for (int i = 0; i < _numberOfElements; i++)
	{
		//Wyzerowanie wartoœci
		for (int k = 0; k < 4; k++) {
			elements[i].pLocal[k] = 0.0;
			for (int j = 0; j < 4; j++) {
				elements[i].hLocal[k][j] = 0.0;
				elements[i].cLocal[k][j] = 0.0;
				
				elements[i].NpoDx[k][j] = 0.0;
				elements[i].NpoDy[k][j] = 0.0;
			}
		}

		if (i == 0)
			elements[i].ID[0] = 0;

		else
			if (i % (_heightsNodes - 1) != 0)
				elements[i].ID[0] = elements[i - 1].ID[3];
			else
				elements[i].ID[0] = elements[i - 1].ID[3] + 1;
		
		elements[i].ID[1] = elements[i].ID[0] + _heightsNodes;
		elements[i].ID[2] = elements[i].ID[1] + 1;
		elements[i].ID[3] = elements[i].ID[0] + 1;
		
	}
	
	//Ustalenie powierzchni brzegowych
	for (int i = 0; i < _numberOfElements; i++)
	{
		elements[i].pow[0] = nodes[elements[i].ID[0]].status && nodes[elements[i].ID[3]].status ? 1 : 0;		
		elements[i].pow[1] = nodes[elements[i].ID[0]].status && nodes[elements[i].ID[1]].status ? 1 : 0;		
		elements[i].pow[2] = nodes[elements[i].ID[1]].status && nodes[elements[i].ID[2]].status ? 1 : 0;
		elements[i].pow[3] = nodes[elements[i].ID[2]].status && nodes[elements[i].ID[3]].status ? 1 : 0;
		
	}
}

//Funkcje kszta³tu
double Grid::N1(double ksi, double eta) {
	return 0.25*(1 - ksi)*(1 - eta);
}
double Grid::N2(double ksi, double eta) {
	return 0.25*(1 + ksi)*(1 - eta);
}
double Grid::N3(double ksi, double eta) {
	return 0.25*(1 + ksi)*(1 + eta);
}
double Grid::N4(double ksi, double eta) {
	return 0.25*(1 - ksi)*(1 + eta);
}


double Grid::Matrix() {

	//Wspó³rzedne punktów ca³kowania Gaussa
	double pierw[4][2] = {
		{ -0.5773502691896258, -0.5773502691896258 },
		{ 0.5773502691896258, -0.5773502691896258 },
		{ 0.5773502691896258,  0.5773502691896258 },
		{ -0.5773502691896258,  0.5773502691896258 },
	};

	//Obliczanie funkcji kszta³tu i pochodnych
	double funkcjeKsztaltu[4][4] = {
		{ this->N1(pierw[0][0],pierw[0][1]), this->N2(pierw[0][0],pierw[0][1]), this->N3(pierw[0][0],pierw[0][1]), this->N4(pierw[0][0],pierw[0][1]) },
		{ this->N1(pierw[1][0],pierw[1][1]), this->N2(pierw[1][0],pierw[1][1]), this->N3(pierw[1][0],pierw[1][1]), this->N4(pierw[1][0],pierw[1][1]) },
		{ this->N1(pierw[2][0],pierw[2][1]), this->N2(pierw[2][0],pierw[2][1]), this->N3(pierw[2][0],pierw[2][1]), this->N4(pierw[2][0],pierw[2][1]) },
		{ this->N1(pierw[3][0],pierw[3][1]), this->N2(pierw[3][0],pierw[3][1]), this->N3(pierw[3][0],pierw[3][1]), this->N4(pierw[3][0],pierw[3][1]) },
	};

	double funkcjeKsztaltuPoKsi[4][4] = {
		{ -0.25*(1 - pierw[0][1]), 0.25*(1 - pierw[0][1]), 0.25*(1 + pierw[0][1]),	-0.25*(1 + pierw[0][1]) },
		{ -0.25*(1 - pierw[1][1]), 0.25*(1 - pierw[1][1]), 0.25*(1 + pierw[1][1]),	-0.25*(1 + pierw[1][1]) },
		{ -0.25*(1 - pierw[2][1]), 0.25*(1 - pierw[2][1]), 0.25*(1 + pierw[2][1]),	-0.25*(1 + pierw[2][1]) },
		{ -0.25*(1 - pierw[3][1]), 0.25*(1 - pierw[3][1]), 0.25*(1 + pierw[3][1]),	-0.25*(1 + pierw[3][1]) }
	};

	double funkcjeKsztaltuPoEta[4][4] = {
		{ -0.25*(1 - pierw[0][0]), -0.25*(1 + pierw[0][0]), 0.25*(1 + pierw[0][0]),	0.25*(1 - pierw[0][0]) },
		{ -0.25*(1 - pierw[1][0]), -0.25*(1 + pierw[1][0]), 0.25*(1 + pierw[1][0]),	0.25*(1 - pierw[1][0]) },
		{ -0.25*(1 - pierw[2][0]), -0.25*(1 + pierw[2][0]), 0.25*(1 + pierw[2][0]),	0.25*(1 - pierw[2][0]) },
		{ -0.25*(1 - pierw[3][0]), -0.25*(1 + pierw[3][0]), 0.25*(1 + pierw[3][0]),	0.25*(1 - pierw[3][0]) }
	};
	double XpoKsi = 0.0, XpoEta = 0.0, YpoKsi = 0.0, YpoEta = 0.0, detJ_inverse = 0.0, detJ = 0.0;


	double pom[2][2];
	double jacobian1D = 0;
	double npc[2][2];
	double tmpLocal[4][4],
		tmpPlocal[4][4],
		pomDlaH[4][4];
	double tmpNP[2][4];
	bool warunek;
	for (vector<Element>::iterator elem = this->elements.begin(); elem != this->elements.end(); elem++) {
					
			XpoKsi =(0.5773502691896258 / 4.0)*(this->nodes[elem->ID[0]].x - (this->nodes[elem->ID[1]].x) + (this->nodes[elem->ID[2]].x) - this->nodes[elem->ID[3]].x) +
					0.25*(-this->nodes[elem->ID[0]].x + (this->nodes[elem->ID[1]].x) + (this->nodes[elem->ID[2]].x) - this->nodes[elem->ID[3]].x);
			
			XpoEta =
				(0.5773502691896258 / 4.0)*(this->nodes[elem->ID[0]].x - (this->nodes[elem->ID[1]].x) + (this->nodes[elem->ID[2]].x) - this->nodes[elem->ID[3]].x) +
						 0.25*(-this->nodes[elem->ID[0]].x - (this->nodes[elem->ID[1]].x) + (this->nodes[elem->ID[2]].x) + this->nodes[elem->ID[3]].x);

			YpoKsi =
				(0.5773502691896258 / 4.0)*(this->nodes[elem->ID[0]].y - this->nodes[elem->ID[1]].y + (this->nodes[elem->ID[2]].y ) - (this->nodes[elem->ID[3]].y )) +
				0.25*(-this->nodes[elem->ID[0]].y + this->nodes[elem->ID[1]].y + (this->nodes[elem->ID[2]].y) - (this->nodes[elem->ID[3]].y));

			YpoEta =
				(0.5773502691896258 / 4.0)*(this->nodes[elem->ID[0]].y - this->nodes[elem->ID[1]].y + (this->nodes[elem->ID[2]].y) - (this->nodes[elem->ID[3]].y)) +
				0.25*(-this->nodes[elem->ID[0]].y - this->nodes[elem->ID[1]].y + (this->nodes[elem->ID[2]].y) + (this->nodes[elem->ID[3]].y ));
			
			//Wyznacznik [J]
			detJ = ((XpoKsi*YpoEta) - (XpoEta*YpoKsi));
			detJ_inverse = 1 / detJ;

			//Pomocnicza macierz wyniku mnozenia 1/det * [J]^-1
			pom[0][0] = detJ_inverse*YpoEta;
			pom[0][1] = -detJ_inverse*YpoKsi;
			pom[1][0] = -detJ_inverse*XpoEta;
			pom[1][1] = detJ_inverse*XpoKsi;

			//Mnozenie przez pochodne po funkcjach ksztaltu dla kazdego wêz³a
				for (int j = 0; j < 4; j++) {
					for (int k = 0; k < 4; k++) {
						elem->NpoDx[j][k] = pom[0][0] * funkcjeKsztaltuPoKsi[j][k] + pom[0][1] * funkcjeKsztaltuPoEta[j][k];
						elem->NpoDy[j][k] = pom[1][0] * funkcjeKsztaltuPoKsi[j][k] + pom[1][1] * funkcjeKsztaltuPoEta[j][k];
					}
				}
			
			
				for (int k = 0; k < 4; k++) {
					//Tworzenie lokalnej macierzy H
					for (int i = 0; i < 4; i++) {
						for (int j = 0; j < 4; j++) {
							elem->hLocal[i][j] += detJ*this->_condK*((elem->NpoDx[k][j] * elem->NpoDx[k][i]) + (elem->NpoDy[k][j] * elem->NpoDy[k][i]));						
						}			
					}

					//Tworzenie lokalnej macierzy C	
					for (int i = 0; i < 4; i++) {
						for (int j = 0; j < 4; j++) {
							elem->cLocal[i][j] += detJ*(funkcjeKsztaltu[k][j] * funkcjeKsztaltu[k][i]);
						}	
					}

				}
			
			//Przemna¿anie elementów macierzy C przez ciep³o w³aœciwe i gêstoœæ
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						
						elem->cLocal[i][j] *= (this->_cp*this->_densityV);
				
					}
				
				}

			//Zerowanie
				jacobian1D = 0;
				warunek = false;
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						if (i < 2)
							tmpNP[i][j];
						tmpLocal[i][j] = tmpPlocal[i][j] = pomDlaH[i][j] = 0.0;
					}
				}
			
			//Jacobian 1D dla ka¿dej powierzchni
				for (int i = 0; i < 4; i++) {
					if (elem->pow[i] == 1) {
						warunek = true;
						switch (i) {
						case 0:
							jacobian1D = (this->nodes[elem->ID[3]].y - this->nodes[elem->ID[0]].y) / 2;
							npc[0][0] = -1.0;
							npc[0][1] = -0.5773502691896258;
							npc[1][0] = -1.0;
							npc[1][1] = 0.5773502691896258;
							break;
						case 1:
							jacobian1D = (this->nodes[elem->ID[1]].x - this->nodes[elem->ID[0]].x) / 2;
							npc[0][0] = -0.5773502691896258;
							npc[0][1] = -1.0;
							npc[1][0] = 0.5773502691896258;
							npc[1][1] = -1.0;
							break;
						case 2:
							jacobian1D = (this->nodes[elem->ID[2]].y - this->nodes[elem->ID[1]].y) / 2;
							npc[0][0] = 1.0;
							npc[0][1] = -0.5773502691896258;
							npc[1][0] = 1.0;
							npc[1][1] = 0.5773502691896258;
							break;
						case 3:
							jacobian1D = (this->nodes[elem->ID[2]].x - this->nodes[elem->ID[3]].x) / 2;
							npc[0][0] = -0.5773502691896258;
							npc[0][1] = 1.0;
							npc[1][0] = 0.5773502691896258;
							npc[1][1] = 1.0;
							break;
						}
						for (int n = 0; n < 2; n++) {
							tmpNP[n][0] = this->N1(npc[n][0], npc[n][1]);
							tmpNP[n][1] = this->N2(npc[n][0], npc[n][1]);
							tmpNP[n][2] = this->N3(npc[n][0], npc[n][1]);
							tmpNP[n][3] = this->N4(npc[n][0], npc[n][1]);	
						}
					
						for (int j = 0; j < 4; j++) {
							tmpPlocal[i][j] += jacobian1D*(tmpNP[0][j] + tmpNP[1][j]);
						}

						for (int k = 0; k < 4; k++) {
							for (int j = 0; j < 4; j++) {
								pomDlaH[j][k] = pomDlaH[j][k] + ((tmpNP[0][k] * tmpNP[0][j])*jacobian1D);
								pomDlaH[j][k] = pomDlaH[j][k] + ((tmpNP[1][k] * tmpNP[1][j])*jacobian1D);
							}
						}
					}
				}
			


			//Tworzenie wektora P lokalnego
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						elem->pLocal[i] += tmpPlocal[j][i];
					}
				}
			
			//Uwzglêdnienie wspó³czynnika przenikania
				if (warunek) {
					for (int i = 0; i < 4; i++) {
						for (int j = 0; j < 4; j++) {
							elem->hLocal[i][j] += this->_alfa*(pomDlaH[i][j]);					
						}
					}
				}
			
			//Przemna¿anie lokalnego wektora P przez wspó³czynnik przenikania i temperature otoczenia
				for (int i = 0; i < 4; i++) {
					elem->pLocal[i] = elem->pLocal[i] * this->_alfa * this->_ambientTemp;
				}
	}
	 
	//Uzupe³nianie macierzy globalnej H i C
	for (int k = 0; k < this->elements.size();k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				this->hGlobal[this->elements[k].ID[i]][this->elements[k].ID[j]] += this->elements[k].hLocal[i][j];
				this->cGlobal[this->elements[k].ID[i]][this->elements[k].ID[j]] += this->elements[k].cLocal[i][j];
			}
		}
	}

	//Uzupe³nianie globalnego wektora P
	for (int k = 0; k < this->elements.size(); k++) {
		for (int j = 0; j < 4; j++) {
			this->pGlobal[this->elements[k].ID[j]] += this->elements[k].pLocal[j];
		}
	}

	//printf("\nP Globalne\n");
	/*for (int i = 0; i < 16; i++) {
		printf("%.3f  ", this->pGlobal[i]);
	}
	cout << endl;
	printf("\nH Globalne\n");
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			printf("%.3f  ", this->hGlobal[i][j]);
		}
		cout << endl;
	}*/


	printf("\n\n");
	return 0.0;
}


void Grid::Solve() {
	double **h = new double*[this->_numberOfNodes];
	double **cDivided = new double*[this->_numberOfNodes];
	double *p = new double[this->_numberOfNodes];
	double *t0 = new double[this->_numberOfNodes];
	double *ct0 = new double[this->_numberOfNodes];
	double *x = new double[this->_numberOfNodes];
	for (int i = 0; i < this->_numberOfNodes; i++) {
		h[i] = new double[this->_numberOfNodes];
		t0[i] = this->_initTemp;
		cDivided[i] = new double[this->_numberOfNodes];
		for (int j = 0; j < this->_numberOfNodes;j++) {
			h[i][j] = this->hGlobal[i][j] + (this->cGlobal[i][j] / this->_stepTime);
			cDivided[i][j] = this->cGlobal[i][j] / this->_stepTime;		
		}	
	}

	double minTemp = 0;
	double maxTemp = 0;
	cout << "Time[s]\tMinTemp[s]\tMaxTemp[s]\n";
	for (int i = this->_stepTime; i <= this->_simTime; i+=this->_stepTime) {
		for (int l = 0; l < this->_numberOfNodes; l++) {
			ct0[l] = p[l] = x[l] = 0.0;
		}

		for (int j = 0; j < this->_numberOfNodes; j++) {

			for (int k = 0; k < this->_numberOfNodes; k++) {
					p[j] += cDivided[j][k] * t0[k];
				}
				p[j] += this->pGlobal[j];
		}

		if (tempVec(this->_numberOfNodes, h, p, x)) {
			minTemp = min(x, this->_numberOfNodes);
			maxTemp = max(x, this->_numberOfNodes);
			
			for (int l = 0; l < this->_numberOfNodes; l++) {
				t0[l] = x[l];
			}
		}

		cout << i << "\t" << minTemp << "\t\t" << maxTemp << endl;
	}

}

const double eps = 1e-12;
//Obliczanie wektora temperatury (rozwi¹zanie uk³adu równañ)
bool Grid::tempVec(int n, double ** A, double * B, double *& X)
{
	int    i, j, k;
	double s;

	double **tab = new double*[n];
	for (int z = 0; z < n; z++) {
		tab[z] = new double[n];
		for (int x = 0; x < n; x++) {
			tab[z][x] = A[z][x];
		}
	}

	for (k = 0; k < n - 1; k++) {
		if (fabs(tab[k][k]) < eps) return false;

		for (i = k + 1; i < n; i++)
			tab[i][k] /= tab[k][k];

		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				tab[i][j] -= tab[i][k] * tab[k][j];
	}

	X[0] = B[0];

	for (i = 1; i < n; i++) {
		s = 0;

		for (j = 0; j < i; j++) s += tab[i][j] * X[j];

		X[i] = B[i] - s;
	}

	if (fabs(tab[n - 1][n - 1]) < eps) return false;

	X[n - 1] /= tab[n - 1][n - 1];

	for (i = n - 2; i >= 0; i--) {
		s = 0;

		for (j = i + 1; j < n; j++) s += tab[i][j] * X[j];

		if (fabs(tab[i][i]) < eps) return false;

		X[i] = (X[i] - s) / tab[i][i];
	}
	for (int i = 0; i < n; i++) {
		delete tab[i];
	}
	delete tab;
	return true;
}

double Grid::max(double *vec, int n) {
	double max = vec[0];
	for (int i = 1; i < n-1; i++) {
		if (vec[i] > max) {
			max = vec[i];
		}
	}
	return max;
}

double Grid::min(double *vec, int n) {
	double min = vec[0];
	for (int i = 1; i < n - 1; i++) {
		if (vec[i] < min) {
			min = vec[i];
		}
	}
	return min;
}