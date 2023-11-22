#pragma once
#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <limits>
#include "Gauss.h"

using namespace std;
const vector<double> operator-(const vector<double> x1, const vector<double> x2);
const vector<double> operator/(const vector<double> x1, const vector<double> x2);

double f1(double x1, double x2);

double f2(double x1, double x2);

double f1x1(double x1, double x2);

double f1x2(double x1, double x2);

double f2x1(double x1, double x2);

double f2x2(double x1, double x2);
vector <double> solve(double const x1, double const x2, double& e, double& e2, int const& NIT);

vector<vector<double>> Jkm(double const& M, vector<double> curr);

vector <double> solve1(double const x1, double const x2, double& e, double& e2, int const& NIT, double M);
