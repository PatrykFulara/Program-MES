#include <cmath>
const double eps = 1e-12;
bool ludist(int n, double ** A)
{
	int i, j, k;

	for (k = 0; k < n - 1; k++) {
		if (fabs(A[k][k]) < eps) return false;

		for (i = k + 1; i < n; i++)
			A[i][k] /= A[k][k];

		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				A[i][j] -= A[i][k] * A[k][j];
	}

	return true;
}

// Funkcja wyznacza wektor X na podstawie A i B
//---------------------------------------------
bool lusolve(int n, double ** A, double * B, double * X)
{
	int    i, j;
	double s;

	X[0] = B[0];

	for (i = 1; i < n; i++) {
		s = 0;

		for (j = 0; j < i; j++) s += A[i][j] * X[j];

		X[i] = B[i] - s;
	}

	if (fabs(A[n - 1][n - 1]) < eps) return false;

	X[n - 1] /= A[n - 1][n - 1];

	for (i = n - 2; i >= 0; i--) {
		s = 0;

		for (j = i + 1; j < n; j++) s += A[i][j] * X[j];

		if (fabs(A[i][i]) < eps) return false;

		X[i] = (X[i] - s) / A[i][i];
	}

	return true;
}