#pragma once
#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <limits>

using namespace std;
bool is_zero(double x);

bool divison(vector<vector<double>>& matrix, int step);

double module(vector<double> ar);
double max_abs(vector<double> ar);

void minus_(vector<vector<double>>& matrix, int step);

void lineswap(vector<vector <double>>& matrix, int step);
void output(vector<vector<double>> const& matrix);

void output(vector<double> const& matrix);
bool triangular_view(vector<vector<double>>& matrix, int const& numberOfEquations);
void findX(vector<vector<double>> matrix, vector <double>& result);