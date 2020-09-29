// Solver
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>

//include (ExternalProject)
#include "yaml-cpp/yaml.h"
//
//target_link_libraries(playground
//	yaml-cpp
//)

#include "Loop.h"

using namespace std;

// set < string > solvers;
struct sol_struct
{
	string gas;
	double Cp;
	double cfl;
	int max_iter_num;
	double tolerance;
	double gamma;

	vector < double > RK_stage_coeffs;
	vector < int > res_smooth_flag;
	double eps_impl_res_smooth;
};

class Solver
{
public:
	string solver_name;

	vector < vector < double > > cv;
	vector < double > p;
	vector < double > dt;
	vector < double > dummy;

	vector < vector < double > > cvold;
	vector < vector < double > > diss;
	vector < vector < double > > rhs;

	vector < double > a;			///< cross sections of tube
	vector < double > x;			///< coordinates of grid (with dummy points)
	vector < double > vol;			///< volumes of cells
	int imax;						///< number of grid points
	int ib2;						///< last physical grid point
	int iter;
	double drho1;

	bool NoOutput;

	// int direction;

	Solver(sol_struct& sol_init_);

	// int ReadConfigFile(string file_name);	///< Reading Config file
	void ReadBoundaries(string file_name);		///< Reading Boundary file
	void InitFlow(double rho, double mass, double e, double p2);
	void RefreshBoundaries();
	void SetGasType(bool gas);		///< Sets gas type. true - Ideal, false - not Ideal
	void SetCp(double Cp_);			///< Setting Cp for Ideal gas
	double GetCp();					///< Getting Cp (ideal or not?)
	bool IsIdeal();					///< Checking if the gas is ideal or not

	void SetCFL(double cfl_);		///< Setting Courant number
	void SetMaxIterNum(int max_iter_num_); 	///< Setting maximal number of iterations
	void SetTolerance(double tolerance_);	///< Setting general tolerance for established convergence
	void CalculateVolumes();		///< Recalculate columes in grid

	double GetGamma();
	// double GetCp();
	double GetInflowT();
	double GetInflowP();
	double GetOutflowP();
	int GetInflowId();
	int GetOutflowId();
	double GetTolerance();
	vector < double > GetSection();
	int GetMaxIterNum();

	void SetGrid(Loop& loop);
	void InverseGrid();
	void ShowGrid();

	double Solve();
	virtual const vector < int >& GetDissFlag() = 0;
	virtual const vector < double >& GetDissBlend() = 0;
	virtual void Dissipation(double beta) = 0;
	virtual void LRState() = 0;
	virtual void Fluxes() = 0;
	void TimeSteps();
	void SourceTerm();
	void ImplResidualSmooth();
	double Convergence();
	void PrintResult();
	void SetOutputFile(string output_file_);
	void ShowOutput();
	void HideOutput();

private:
	Loop grid;
	string output_file;

	bool Ideal;						///< Type of gas, ideal or not.
	double Cp;						///< Cp

	double cfl;						///< Courant number
	int max_iter_num;				///< Maximal number of iterations
	double tolerance;				///< General tolerance for established convergence
	double gamma;					///< adiabatic number

	vector < double > RK_stage_coeffs; 	///< Runge-Kutta stage coefficients
	vector < int > res_smooth_flag;	///< residual smoothing (1=yes)
	double eps_impl_res_smooth;		///< coefficient of implicit residual smoothing
	int ncells;						///< number of cells

	int inflow_id;
	double p_b_in;
	double T_b_in;

	int outflow_id;
	double p_b_out;

	int RK_stages_num;					///< number of Runge-Kutta stages

	// vector < vector < double > > cv;
	// vector < double > p;
	// vector < double > dt;
	// vector < double > dummy;

	// vector < vector < double > > cvold;
	// vector < vector < double > > diss;

	// map < 
};

Solver* CreateReadConfigFile(string file_name);

double Sign(double value);
