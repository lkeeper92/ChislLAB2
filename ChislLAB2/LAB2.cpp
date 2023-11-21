#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <limits>
#include "Gauss.h"

using namespace std;
const vector<double> operator-(const vector<double> x1, const vector<double> x2)
{
	vector<double> result;
	if (x1.size() > x2.size())
		for (int i = 0; i < x2.size(); i++)
			result.push_back(x1[i] - x2[i]);
	for (int i = 0; i < x1.size(); i++)
		result.push_back(x1[i] - x2[i]);
	return result;
}
const vector<double> operator/(const vector<double> x1, const vector<double> x2)
{
	vector<double> result;
	if (x1.size() == x2.size())
	{
		for (int i = 0; i < x2.size(); i++)
			result.push_back(x1[i] / x2[i]);

		return result;
	}
	cerr << "ATTENTION !!!";
	exit(0);
}

// ‘ункци€, вычисл€юща€ значение f1
double f1(double x1, double x2) {
	return 2 * pow(x1, 2) - x1 * x2 - 5 * x1 + 1;
}

// ‘ункци€, вычисл€юща€ значение f2
double f2(double x1, double x2) {
	return x1 + 3 * log10(x1) - pow(x2, 2);
}

// ‘ункци€, вычисл€юща€ значение производной f1 по x1
double f1x1(double x1, double x2) {
	return 4 * x1 - x2 - 5;
}

// ‘ункци€, вычисл€юща€ значение производной f1 по x2
double f1x2(double x1, double x2) {
	return  -x1;
}

// ‘ункци€, вычисл€юща€ значение производной f2 по x1
double f2x1(double x1, double x2) {
	return 1 + 3 * (1 / (x1 * log(10)));
}

// ‘ункци€, вычисл€юща€ значение производной f2 по x2
double f2x2(double x1, double x2) {
	return -2 * x2;
}

vector <double> solve(double const x1, double const x2, double& e, double& e2, int const& NIT) {
	// »нициализаци€ переменных дл€ хранени€ предыдущих значений
	vector <double> prev;
	prev.push_back(x1);
	prev.push_back(x2);
	// »нициализаци€ переменных дл€ хранени€ текущих значений
	vector <double> curr;
	curr.resize(2, 0);
	curr = prev;
	cout << "Start: (" << x1 << ", " << x2 << ")" << endl;
	cout << "Accuracy: " << e << endl;
	vector <double> dx;
	int k = 0;
	for (k; k <= NIT; k++)
	{
		prev = curr;
		vector<double> Fk;
		Fk.push_back(f1(prev[0], prev[1]));
		Fk.push_back(f2(prev[0], prev[1]));
		vector<vector<double>> Jk;
		Jk.resize(2);
		Jk[0].push_back(f1x1(prev[0], prev[1]));
		Jk[0].push_back(f1x2(prev[0], prev[1]));
		Jk[1].push_back(f2x1(prev[0], prev[1]));
		Jk[1].push_back(f2x2(prev[0], prev[1]));
		// matrix for gauss
		vector<vector<double>> matrix = Jk;
		matrix[0].push_back(-Fk[0]);
		matrix[1].push_back(-Fk[1]);
		findX(matrix, dx);

		curr[0] = prev[0] + dx[0];
		curr[1] = prev[1] + dx[1];
		double b2 = max_abs((curr - prev) / curr);
		if (module(curr) < 1) double b2 = max_abs(curr - prev);
		vector<double> b_1;
		b_1.push_back(f1(curr[0], curr[1]));
		b_1.push_back(f2(curr[0], curr[1]));
		double b1 = max_abs(b_1);
		cout << setw(5) << "b1: " << b1 << setw(5) << "b2: " << b2 << endl;
		if (b1 <= e && b2 <= e2) {
			cout << "\n ------------------------------------------------------------------------------\n";
			cout << "Done !\n";
			break;
		}
	}
	cout << "\n=========================================================================\n";
	cout << "Solution: (" << curr[0] << ", " << curr[1] << ")" << endl;
	cout << "Iter: " << k << endl;

	return  curr;
}

vector<vector<double>> Jkm(double const& M, vector<double> curr)
{
	vector<vector<double>> Jk;
	Jk.resize(2);
	Jk[0].push_back((f1(curr[0] + M, curr[1]) - f1(curr[0], curr[1])) / M);
	Jk[0].push_back((f1(curr[0], M + curr[1]) - f1(curr[0], curr[1])) / M);
	Jk[1].push_back((f2(curr[0] + M, curr[1]) - f2(curr[0], curr[1])) / M);
	Jk[1].push_back((f2(curr[0], M + curr[1]) - f2(curr[0], curr[1])) / M);
	return Jk;
}

vector <double> solve1(double const x1, double const x2, double& e, double& e2, int const& NIT) {
	// »нициализаци€ переменных дл€ хранени€ предыдущих значений
	vector <double> prev;
	prev.push_back(x1);
	prev.push_back(x2);
	// »нициализаци€ переменных дл€ хранени€ текущих значений
	vector <double> curr;
	curr = prev;
	curr.resize(2, 0);
	cout << "Start: (" << x1 << ", " << x2 << ")" << endl;
	cout << "Accuracy: " << e << endl;
	vector <double> dx;
	cout << "Enter M:";
	double M;
	cin >> M;
	int k = 0;
	for (k; k <= NIT; k++)
	{
		prev = curr;
		vector<double> Fk;
		Fk.push_back(f1(prev[0], prev[1]));
		Fk.push_back(f2(prev[0], prev[1]));
		vector<vector<double>> Jk;
		Jk = Jkm(M, curr);
		// matrix for gauss
		vector<vector<double>> matrix = Jk;
		matrix[0].push_back(-Fk[0]);
		matrix[1].push_back(-Fk[1]);
		findX(matrix, dx);

		curr[0] = prev[0] + dx[0];
		curr[1] = prev[1] + dx[1];
		double b2 = max_abs((curr - prev) / curr);
		if (module(curr) < 1) double b2 = max_abs(curr - prev);
		vector<double> b_1;
		b_1.push_back(f1(curr[0], curr[1]));
		b_1.push_back(f2(curr[0], curr[1]));
		double b1 = max_abs(b_1);
		cout << setw(5) << "b1: " << b1 << setw(5) << "b2: " << b2 << endl;
		if (b1 <= e && b2 <= e2) {
			cout << "\n ------------------------------------------------------------------------------\n";
			cout << "Done !\n";
			break;
		}
	}
	// ¬ывод результата
	cout << "\n=========================================================================\n";
	cout << "Solution: (" << curr[0] << ", " << curr[1] << ")" << endl;
	cout << "Iter: " << k << endl;
	return  curr;
}


int main()
{
	int NIT;
	std::cout << "Enter NIT: ";
	cin >> NIT;
	double x1_1 = 3;
	double x2_1 = 2;
	double x1_2 = 3;
	double x2_2 = -2;

	double e = 1e-9;
	double e2 = e;
	// –ешение системы дл€ каждого начального приближени€
	solve(x1_1, x2_1, e, e2, NIT);
	solve1(x1_1, x2_1, e, e2, NIT);

	return 0;
}

