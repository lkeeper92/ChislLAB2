#include "Newton.h"

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

// Функция, вычисляющая значение f1
double f1(double x1, double x2) {
	return 2 * pow(x1, 2) - x1 * x2 - 5 * x1 + 1;
}

double f2(double x1, double x2) {
	return x1 + 3 * log10(x1) - pow(x2, 2);
}

double f1x1(double x1, double x2) {
	return 4 * x1 - x2 - 5;
}

double f1x2(double x1, double x2) {
	return  -x1;
}

double f2x1(double x1, double x2) {
	return 1 + 3 * (1 / (x1 * log(10)));
}

double f2x2(double x1, double x2) {
	return -2 * x2;
}
vector <double> solve(double const x1, double const x2, double& e, double& e2, int const& NIT) {
	// Инициализация переменных для хранения предыдущих значений
	vector <double> prev;
	prev.push_back(x1);
	prev.push_back(x2);
	// Инициализация переменных для хранения текущих значений
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
		double b2 = 0;
		for (double i : curr)
		{
			if (abs(i) < 1) b2 = maxAbs(curr - prev);
			else b2 = maxAbs((curr - prev) / curr);
		}
		vector<double> b_1;
		b_1.push_back(f1(curr[0], curr[1]));
		b_1.push_back(f2(curr[0], curr[1]));
		double b1 = maxAbs(b_1);
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
	Jk[0].push_back((f1(curr[0] + M*curr[0], curr[1]) - f1(curr[0], curr[1])) / M*curr[0]);
	Jk[0].push_back((f1(curr[0], M*curr[1] + curr[1]) - f1(curr[0], curr[1])) / M*curr[1]);
	Jk[1].push_back((f2(curr[0] + M*curr[0], curr[1]) - f2(curr[0], curr[1])) / M*curr[0]);
	Jk[1].push_back((f2(curr[0], M*curr[1]	 + curr[1]) - f2(curr[0], curr[1])) / M*curr[1]);
	return Jk;
}
	
vector <double> solve1(double const x1, double const x2, double& e, double& e2, int const& NIT, double M) {
	// Инициализация переменных для хранения предыдущих значений
	vector <double> prev;
	prev.push_back(x1);
	prev.push_back(x2);
	// Инициализация переменных для хранения текущих значений
	vector <double> curr;
	curr = prev;
	curr.resize(2, 0);
	cout <<endl<< "Start: (" << x1 << ", " << x2 << ")" << endl;
	cout << "Accuracy: " << e << endl;
	vector <double> dx;
	cout << " \t For M = " << M << endl;
	double b1 = 0, b2 = 0, k=0;
	do
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
		for (double i : curr)
		{
			if (abs(i) < 1) b2 = maxAbs(curr - prev);
			else b2 = maxAbs((curr - prev) / curr);
		}
		vector<double> b_1;
		b_1.push_back(f1(curr[0], curr[1]));
		b_1.push_back(f2(curr[0], curr[1]));
		b1 = maxAbs(b_1);

		cout << setw(5) << "b1: " << b1 << setw(5) << "b2: " << b2 << endl;
		k++;
		if ((k - NIT)==1)
			cout << "Iter more than NIT" << endl;
	} while (b1 > e && b2 > e2);

				cout << "\n ------------------------------------------------------------------------------\n";
				cout << "Done !\n";
				cout << "Solution: (" << curr[0] << ", " << curr[1] << ")" << endl;
				cout << "Iter: " << k << endl;
				cout << "\n=========================================================================\n";
	return  curr;
}