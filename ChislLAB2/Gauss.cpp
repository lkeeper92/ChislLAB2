#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <limits>
#include "Gauss.h"

using namespace std;
bool is_zero(double x) {
	return std::abs(x) < std::numeric_limits<double>::epsilon();
}

bool divison(vector<vector<double>>& matrix, int step)
{
	double temp = matrix[step][step];
	if (is_zero(temp)) return 0;
	for (int j = step; j < matrix[step].size(); j++)
		matrix[step][j] = matrix[step][j] / temp;
	return 1;
}

double module(vector<double> ar)
{
	double mod = 0;
	for (double i : ar)
		mod += i * i;
	return sqrt(mod);
}

double max_abs(vector<double> ar)
{
	double max = 0;
	for (double i : ar)
		if (max < abs(i)) max = abs(i);
	return max;

}

void minus_(vector<vector<double>>& matrix, int step)
{
	for (int i = step + 1; i < matrix.size(); i++)
	{
		double temp = matrix[i][step];
		for (int j = step; j < matrix[i].size(); j++)
		{
			matrix[i][j] = matrix[i][j] - matrix[step][j] * temp;
		}
	}
}

void lineswap(vector<vector <double>>& matrix, int step)
{
	double max = 0;
	int position = 0;
	for (int i = step; i < matrix.size(); i++)
	{
		if (max < abs(matrix[i][step])) {
			max = abs(matrix[i][step]);
			position = i;
		}
	}
	if (step != position) {
		vector<double> temp = matrix[step];
		matrix[step] = matrix[position];
		matrix[position] = temp;
	}
}
void output(vector<vector<double>> const& matrix)
{
	if (matrix.size() == 0) {
		cout << "Matrix is empty ! ";
		exit(0);
	}
	for (int i = 0; i < matrix.size(); i++)
	{
		int j;
		cout << "[ ";
		for (j = 0; j < matrix[i].size() - 1; j++)
			cout << setw(12) << matrix[i][j];
		cout << setw(3) << " | " << setw(6) << matrix[i][j] << " ]";
		cout << endl;
	}
	cout << endl;
}

void output(vector<double> const& matrix)
{
	if (matrix.size() == 0) {
		cout << "Matrix is empty ! ";
		exit(0);
	}
	cout << "{ ";
	for (int i = 0; i < matrix.size(); i++)
	{
		cout << setprecision(12) << setw(13) << matrix[i] << "  ";
	}
	cout << "} ";
	cout << endl;
	cout << resetiosflags(ios::fixed);
}
bool triangular_view(vector<vector<double>>& matrix, int const& numberOfEquations)
{
	for (int i = 0; i < numberOfEquations - 1; i++)
	{
		lineswap(matrix, i);
		if (!divison(matrix, i)) return 0;
		minus_(matrix, i);
	}
	if (!divison(matrix, numberOfEquations - 1)) return 0;
	return 1;
}
void findX(vector<vector<double>> matrix, vector <double>& result)
{
	int numberOfEquations = matrix.size();
	if (result.empty()) result.resize(numberOfEquations, 1);
	if (!triangular_view(matrix, numberOfEquations)) {
		cerr << "Matrix cant been transform in triangular view !!!\n";
		exit(1);
	}
	for (int i = numberOfEquations - 1; i >= 0; i--)
	{
		double temp = matrix[i][numberOfEquations];
		for (int c = numberOfEquations - 1; c > i; c--)
			temp -= matrix[i][c] * result[c];
		result[i] = temp;
	}
}