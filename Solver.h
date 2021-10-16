// Solver
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include "yaml-cpp/yaml.h"
#include "adept.h"
//#include "adept/adept.h"
//#include "adept/adept_source.h"

#include "Loop.h"

#include "mass_matrix.h"

#include "Eigen/Dense"

#include <functional>

#include "EigenThings.h"

using namespace std::placeholders;

using namespace Eigen;

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
	bool time_expl;

	vector < double > RK_stage_coeffs;
	vector < vector < double > > RK_alpha;
	vector < vector < double > > RK_beta;
	vector < int > res_smooth_flag;
	double eps_impl_res_smooth;

	double TSRK_teta;
	vector < double > TSRK_d;
	vector < double > TSRK_eta;
	vector < vector < double > > TSRK_q;

	vector < vector < double > > ak;
	vector < double > bk;
	vector < double > ck;
	vector < vector < double > > alpha;

	bool lts; // local time step
	bool steadiness; // steadiness

	int time_stepping;	// time stepping
};

class Solver
{
public:
	//map < string, int > vars { {"RhoA", 0}, {"RhoUA", 1}, {"RhoEA", 2} };
	map < string, int > vars{ {"RhoA", 0}, {"RhoUA", 1}, {"RhoEA", 2} /*, {"nA", 3}*/ };
	//map < string, int > vars_o { {"Rho", 0}, {"U", 1}, {"p", 2}, {"H", 3} };
	map < string, int > vars_o{ {"Rho", 0}, {"U", 1}, {"p", 2}, {"H", 3}/*, {"n", 4}*/ };
	enum Vars_o {
		RHO,
		U,
		P,
		H/*,
		n*/
	};

	/*enum Vars {
		RHO_A,
		RHO_U_A,
		RHO_E_A
	};*/

	double omega = 1.;
	double beta = 2;
	double prev_drho = 1.;

	//Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner > gmres;
	/*Eigen::GMRES<SparseMatrix<double>, Eigen::DiagonalPreconditioner < double > > gmres;*/
	ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;

	//Eigen::BiCGSTAB<MatrixReplacement, Eigen::DiagonalPreconditioner <SparseMatrix<double> > > gmres;
	Eigen::SparseMatrix<double> S;
	//S.reserve(VectorXi::Constant((ib2 - 1)* eq_num, eq_num* eq_num));

	string b_type;
	bool lts;			// local time step
	bool steadiness;	// steadiness

	int time_stepping;

	int RHO_A = 0,
		RHO_U_A = 1,
		RHO_E_A = 2/*,
		N_A = 3*/;

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

	using DiagonalFunc = std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int,
		std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >, std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >) >;

	using LUFunc = std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >;

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
	vector < vector < double > > cvn;
	vector < vector < double > > cvnm1;
	vector < vector < double > > diss;
	vector < vector < double > > rhs;

	mass_matrix M_matrix;

	vector < vector < double > > Q_star;

	vector < double > a;			///< cross sections of tube
	vector < double > x;			///< coordinates of grid (with dummy points)
	vector < double > y;			///< coordinates of cell centers (for rhs)
	vector < double > vol;			///< volumes of cells
	int imax;						///< number of grid points
	int ib2;						///< last physical grid point
	int iter;
	double drho1;
	double dconc1;

	bool NoOutput;

	bool time_expl;
	double Global_Time;

	/*vector < vector < vector < double > > > L_SGS;
	vector < vector < vector < double > > > U_SGS;
	vector < vector < vector < double > > > D_SGS;*/
	vector < MatrixXd > L_SGS;
	vector < MatrixXd > U_SGS;
	vector < MatrixXd > D_SGS;

	Solver(sol_struct& sol_init_);

	void SetEquation(string eq_name, string dt_term_, string dx_term_, map < string, int > vars_, map < string, int > vars_o_);

	void AdjustMesh(double* rho_, double* mass_, double* e_, double* p_, double x_);

	void ReadBoundaries(string file_name);		///< Reading Boundary file
	void InitFlow(double rho, double mass, double e, double p2);
	void InitFlow(double* rho_, double* mass_, double* e_, double* p_, double x_);
	void InitFlowAG2(double* rho_, double* mass_, double* e_, double* p_, double x_);
	void InitFlowAG(double* rho_, double* mass_, double* e_, double* p_, double x_);
	void ResetDummy();
	void RefreshBoundaries();
	void SetGasType(bool gas);		///< Sets gas type. true - Ideal, false - not Ideal
	void SetCp(double Cp_);			///< Setting Cp for Ideal gas
	double GetCp();					///< Getting Cp (ideal or not?)
	bool IsIdeal();					///< Checking if the gas is ideal or not

	void SetCFL(double cfl_);		///< Setting Courant number
	double GetCFL();		///< Setting Courant number
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

	void calculate_mass_matrix() { M_matrix.calculate_mass_matrix(++x.begin(), --x.end(), x.size() - 2); };
	void fill_inverse_mass_matrix() { M_matrix.fill_inverse_mass_matrix(); }
	void print_inversed_mass_matrix() { M_matrix.print_inverse_matrix(); }

	void CalculateTimeSource(vector < vector < double > >& cvn_, vector < vector < double > >& cvnm1_, double physDt_);
	void CalculateUnsteadyRHS(double physDt_);

	double Solve(double physDt = 0.);
	//double SolveExplicit(double physDt = 0.);

	void UpdateP(int i);

	void StagesMemoryAlloc_Init(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage);
	void RHSProcessing(vector < vector < vector < double > > >& rhsstage, int rks, double physDt, vector < vector < double > >& rhsold);
	void New_Steady_UnsteadyDualTime_CV(vector < vector < vector < double > > >& cvstage, int rks);
	void New_UnsteadyTwoStep_CV(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, int rks);
	void CVProcessing(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, int rks, DiagonalFunc D_Func = NULL, LUFunc L_Func = NULL, LUFunc U_Func = NULL);
	void FinalTSRK_CV_Calculation(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage);
	void FinalImplicitRHS_CV_Calculation(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage);

	//using DiagonalFunc = MatrixXd(*)(std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < std::vector < double > >&, double, double, double, int);
	//typedef MatrixXd(* DiagonalFunc)(void*, std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < std::vector < double > >&, double, double, double, int);

	//using DiagonalFunc = std::function < MatrixXd (std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >;
	//using LUFunc2 = std::function < MatrixXd(std::vector < MatrixXd >&, int, std::vector < std::vector < double > >&) >;

	void GetWn(vector < vector < double > >& Wn_);
	MatrixXd Diagonal(vector < MatrixXd >& LeftFluxJac, vector <MatrixXd >& RightFluxJac, vector < MatrixXd >& SourceJac, vector < vector < double > >& cv_, int i, LUFunc L_Func, LUFunc U_Func);
	MatrixXd Lower(std::vector < MatrixXd >& LeftFluxJac, std::vector < std::vector < double > >& cv_, int i);
	MatrixXd Lower2(std::vector < MatrixXd >& LeftFluxJac, std::vector < std::vector < double > >& cv_, int i);
	MatrixXd Upper(std::vector < MatrixXd >& RightFluxJac, std::vector < std::vector < double > >& cv_, int i);
	MatrixXd Upper2(std::vector < MatrixXd >& RightFluxJac, std::vector < std::vector < double > >& cv_, int i);
	void LUSGS(DiagonalFunc D_Func, LUFunc L_Func, LUFunc U_Func, int rks, vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, double physDt_ = 1.);
	void GMRES(DiagonalFunc D_Func, LUFunc L_Func, LUFunc U_Func, int rks, vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, double physDt_ = 1.);
	SparseMatrix< double > ILU_0(SparseMatrix< double >& SM);

	//void LUSGS(int rks, MatrixXd& Li, MatrixXd& Ui, MatrixXd& Di, vector < vector < vector < double > > >& rhsstage, vector < vector < double > >& cvstage);
	double SolveExplImpl(double physDt = 0.);
	double ForwardEuler(double physDt);
	MatrixXd ToEigen(vector < vector < double > >& M);
	VectorXd ToEigen(vector < double >& V);
	vector < double > ToVector(VectorXd& EV);
	vector < vector < double > > ToMatrix(const MatrixXd& EM);
	VectorXd GetEigenVector(vector < vector < double > >& Vec, int i);
	vector < double > MatVecMult(vector < vector < double > >& Mat, vector < double >& Vec);
	vector < double > VecSum(double a, vector < double >& Vec1, double b, vector < double >& Vec2);
	void AreaMult(vector < vector < double > >& Mi, int k, int rks);
	MatrixXd MakeD(MatrixXd& Li, MatrixXd& Ui, MatrixXd& Di, double dt, int i, int rks);
	virtual void GetFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac, bool POS_NEG, bool simple = false) = 0;
	MatrixXd GetJacobian(vector < double >& jac);
	virtual void ComputeRHSandJacobian(bool NO_JAC = false) = 0;
	virtual const vector < int >& GetDissFlag() = 0;
	virtual const vector < double >& GetDissBlend() = 0;
	virtual void Dissipation(double beta) = 0;
	virtual void LRState() = 0;
	virtual void LRState(string var_) = 0;
	virtual void LRState(vector < vector < double > >& cv_, vector < vector < double > >& ls_, vector < vector < double > >& rs_) = 0;
	virtual void Fluxes() = 0;
	virtual void RHS(int i) = 0;
	void RhoUPH();
	void RhoUPH(vector < vector < double > > & cv_, vector < vector < double > > & fv_);
	void TimeSteps(bool local_time = true, double dt_ = 0.);
	double SpectralRadius(vector< vector < double > >& cv_, int i);
	void SourceTerm();
	void SourceTerm(int i);
	void ImplResidualSmooth();
	void ImplResidualSmooth(vector < vector < vector < double > > >& rhs);
	double Convergence();
	void PrintResult();
	void SetOutputFile(string output_file_);
	void ShowOutput();
	void HideOutput();

	void InitFirsStages(vector < vector < vector < double > > >& cvstage, vector < vector < vector < double > > >& rhsstage);

	vector < vector < double > > MakeCV(vector < vector < double > >& fv_);

	double CalcH(double p_, double Rho_, double U_);

	//void FillJacobian(vector < vector < double > >& M_SGS, vector < double >& jac, double s);
	void FillJacobian(MatrixXd& M_SGS, vector < double >& jac, double s);

	Loop grid;
private:
	//Loop grid;
	string output_file;

	bool Ideal;						///< Type of gas, ideal or not.
	double Cp;						///< Cp

	double cfl;						///< Courant number
	int max_iter_num;				///< Maximal number of iterations
	double tolerance;				///< General tolerance for established convergence
	double gamma;					///< adiabatic number

	vector < double > RK_stage_coeffs; 	///< Runge-Kutta stage coefficients
	vector < vector < double > > RK_alpha;
	vector < vector < double > > RK_beta;
	vector < int > res_smooth_flag;	///< residual smoothing (1=yes)
	double eps_impl_res_smooth;		///< coefficient of implicit residual smoothing
	int ncells;						///< number of cells

	double TSRK_teta;
	vector < double > TSRK_d;
	vector < double > TSRK_eta;
	vector < vector < double > > TSRK_q;

	int inflow_id;
	double p_b_in;
	double T_b_in;

	int outflow_id;
	double p_b_out;

	int RK_stages_num;					///< number of Runge-Kutta stages

	vector < vector < double > > ak;
	vector < double > bk;
	vector < double > ck;
	vector < vector < double > > alpha;

};

Solver* CreateReadConfigFile(string file_name);

double Sign(double value);

//adept::adouble aSign(adept::adouble value);

int count_(vector < int >& vec, int val);