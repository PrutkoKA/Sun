#pragma once
#include <functional>
#include "DataBase.h"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

using namespace Eigen;
using Eigen::SparseMatrix;

typedef double Scalar;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

class Loop : public DataBase {
public:

	vector< double > GetCoordinates();
	vector< double > GetSection();
	void SetTime(double time_);
	void IncrementTime(double dt);
	double GetTime();
	void Refresh(string col_name, function<double(vector < double >, int)> foo, vector < double > parameter);
	void Refresh(string col_name, DataBase db);
	void RefreshX(string col_name);
	void RefreshR(string col_name);
	void RefreshN(string col_name);
	void CalculateResolution(double X, double F, string function_name, string x_name);
	void CalculateConcentration(double X, string x_name);
	vector <double>& GetConcentrationRef() { return n; };
	void CalculateConcentrationWave();
	SparseMatrix< double > ILU_0(SparseMatrix< double >& SM);
	Vector SolveFromLU(SparseMatrix<double>& LU, Vector& b);
	vector < double > Redistribute(string col_name, vector < double >& conc);
	double FindArea(double* conc, double* x_, int method);
	double FindCubicFunctionValue(double a, double b, double c, double d, double x);
	void SetNewPoint(double S_left, double& sec_area, double& conc_cur, int searchMethod, vector < double >& new_x, double prev_x, int& id, double slope = 0., double* conc = NULL, double* x_ = NULL);
	double dxSearch(double S_left, double conc_, int method, double slope = 0., double prev_x = 0., double* conc = NULL, double* x_ = NULL);
	vector< double > RefineMesh();
	void FindMesh(double dt, string x_name);
	void SmoothN(double coef);
	const vector < double > smooth_least_square (const vector < double > x, const vector < double > f, unsigned int poly_degree, unsigned int half_count_points);
	void AddCollumnR();
	void AddCollumnN();

	vector < double > GetConcentration();
	vector < double > GetResolution();
	double GetAlpha();

	double TotalSum = 0.;

private:
	Eigen::SparseMatrix<double> S;
	vector < double > x;
	vector < double > R;
	vector < double > n;
	vector < double > n_w;	//wave (~)
	vector < double > n_w_old;
	vector < double > n_t;  //tilda? cap? (^)
	double alpha = 3.;
	const string COORD_COL = "coordinate";
	const string SECT_COL = "section";
	const string TIME_COL = "time";
	double time;
};

int gauss(vector < vector < double > > &a, vector < double >& ans);