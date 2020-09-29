// Solver
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include "yaml-cpp/yaml.h"

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
	map < string, int > vars { {"RhoA", 0}, {"RhoUA", 1}, {"RhoEA", 2} };
	map < string, int > vars_o { {"Rho", 0}, {"U", 1}, {"p", 2}, {"H", 3} };
	enum Vars_o {
		RHO,
		U,
		P,
		H
	};

	enum Vars {
		RHO_A,
		RHO_U_A,
		RHO_E_A
	};

	struct equation
	{
		string eq_name;
		int dt_var;
		vector < vector < int > > cur_dt;
		vector < vector < int > > cur_dx;

		equation(string eq_name_, string dt_term_, string dx_term_, map < string, int > vars_, map < string, int > vars_o_);
		vector < vector < int > > dt_term(string dt_term_, map < string, int > vars_o);
		vector < vector < int > > dx_term(string dx_term_, map < string, int > vars_o);
	};

	string solver_name;

	int eq_num, var_num;
	int cur_Q;
	vector < equation > equations;

	vector < vector < double > > cv;
	vector < vector < double > > fv;	// field var
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

	Solver(sol_struct& sol_init_);

	void SetEquation(string eq_name, string dt_term_, string dx_term_, map < string, int > vars_, map < string, int > vars_o_);

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
	virtual void LRState(string var_) = 0;
	virtual void Fluxes() = 0;
	virtual void RHS(int i) = 0;
	void RhoUPH();
	void TimeSteps();
	void SourceTerm();
	void ImplResidualSmooth();
	double Convergence();
	void PrintResult();
	void SetOutputFile(string output_file_);
	void ShowOutput();
	void HideOutput();

	double CalcH(double p_, double Rho_, double U_);

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

};

Solver* CreateReadConfigFile(string file_name);

double Sign(double value);
