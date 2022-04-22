// Solver
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <string>
#include <set>
#include "yaml-cpp/yaml.h"
#include "adept.h"

#include "Loop.h"

#include "mass_matrix.h"

#include "Eigen/Dense"

#include <functional>

#include "EigenThings.h"

using namespace std::placeholders;

using namespace Eigen;

using namespace std;

struct sol_struct
{
	string gas;
	double Cp;
	double cfl;
	int max_iter_num;
	double tolerance;
	double gamma;
	double free_g;
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
	bool remesh; // remeshing for unsteady problems

	int time_stepping;	// time stepping

	string solver_name;
};

enum class operation
{
	plus,
	minus,
	mult,
	div
};

struct eq_term
{
	enum class var_type
	{
		conservative,
		field,
		function,
		gamma,
		gammam,
		area,
		coord,
		cp,
		g,
		k_planc,
		not_defined
	};
	operation op;
	string name;
	double degree;
	double coef;

	var_type v_type = var_type::not_defined;
	int var_id = -1;

	eq_term(const operation op_, const string name_, const double degree_, const double coef_ = 1.) : op(op_), name(name_), degree(degree_), coef(coef_) {};
	eq_term(const string& term_s);
};

class Solver
{
public:
	const double k_planc = 6.626e-34;
	int g_RHO = -1, g_U = -1, g_P = -1, g_H = -1, g_E = -1, g_A = -1, g_X = -1;
	int g_RHO_A = -1, g_RHO_U_A = -1, g_RHO_E_A = -1;

	map<string, int> func_name_ids;

	enum Vars {
		RHO_A,
		RHO_U_A,
		RHO_E_A,
		CONS_VAR_COUNT
	};
	map < string, int > vars;
	map<int, string> c_var_name;

	enum Vars_o {
		RHO,
		U,
		E,
		P,
		H,
		A,
		//N,
		TEMP,
		DA,
		DP,
		DTEMP,
		//FT,
		//RAD,
		//RADFUNC,
		FIELD_VAR_COUNT
	};
	map < string, int > vars_o;
	map<int, string> var_name;

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
	bool remesh; // remeshing for unsteady problems

	int time_stepping;

	
	map<operation, string> op_name = { {operation::plus, "+"}, {operation::minus, "-"}, {operation::mult, "*"}, {operation::div, "/"} };
	map<string, operation> name_to_op = { {"+", operation::plus}, {"-", operation::minus}, {"*", operation::mult}, {"/", operation::div} };
	
	struct equation
	{
		enum class term_name
		{
			dt,
			dx,
			source
		};

		string eq_name;
		pair<string, vector<eq_term>> cur_dt;
		pair<string, vector<eq_term>> cur_dx;
		pair<string, vector<eq_term>> cur_source;

		equation(string eq_name_, vector<string> dt_term_, vector<string> dx_term_, vector<string> source_term, map < string, int > vars_, map < string, int > vars_o_);
		pair<string, vector<eq_term>> get_equation(const string& eq_name, const vector<string>& eq_terms_s);
		pair<string, vector<eq_term>> get_equation(const string& eq_name, const vector<eq_term>& eq_terms);
	};

	using custom_f = std::function<void(const std::vector<std::vector<double>>&, const vector<adept::adouble>&, const std::map < std::string, int >&, const std::vector<std::string>&, int i, adept::adouble&)>;

	class custom_func
	{
	public:
		string func_name;
		const map<string, int>& var_names;
		vector<string> param_names;
		custom_f func;
		
		custom_func(const string& func_name_,
			custom_f function_,
			const map<string, int>& var_names_,
			const vector<string>& param_names_) :
			func_name(func_name_),
			func(function_),
			var_names(var_names_),
			param_names(param_names_) {};

		template <typename T>
		T get_function_value(const std::vector<std::vector<double>>& params, const vector<adept::adouble>& params_a, int i);
	};

	using DiagonalFunc = std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int,
		std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >, std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >) >;

	using LUFunc = std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >;

	string solver_name;

	int eq_num, var_num;
	int cur_Q;
	set<string> delayed_fv_equations;
	vector < equation > equations;
	map<string, vector<eq_term>> fv_equation;
	
	map<string, custom_func/*<double>*/> functions;

	vector < vector < double > > cv;
	vector < vector < double > > fv;	// field var
	vector < double > p;
	vector < double > dt;
	vector < double > dummy;

	vector < vector < double > > cvold;
	vector < vector < double > > cvn_old;
	vector < vector < double > > cvnm1_old;
	vector < vector < double > > cvn;
	vector < vector < double > > cvnm1;
	vector < vector < double > > diss;
	vector < vector < double > > rhs;

	mass_matrix M_matrix;

	vector < vector < double > > Q_star;

	vector<string> RemeshFuncs;
	vector<double> MaxOfRemeshFuncs;
	double RemeshTau;
	double MaxX;
	double MaxF;
	double MaxFn;
	string RemeshVar;

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

	vector < MatrixXd > L_SGS;
	vector < MatrixXd > U_SGS;
	vector < MatrixXd > D_SGS;

	Solver(sol_struct& sol_init_);

	void AddDelayedFvEquation(const set<string>& equations_to_delay);
	void SetEquation(string eq_name, const vector<string>& dt_term_, const vector<string>& dx_term_, const vector<string>& source_term, map < string, int > vars_, map < string, int > vars_o_);
	void set_fv_equation(const string& eq_name, const vector<string>& eq_terms_s);
	void set_fv_equation(const string& eq_name, const vector<eq_term>& eq_terms);
	template<typename T>
	T make_equation(const int eq, const equation::term_name term_name, const vector<T*>& f_vars);
	template<typename T>
	T make_fv_equation(const string& eq_name, const int point, const vector<T*>& field_var = vector<T*>(), const T* cons_var = nullptr);

	void set_function(const string& func_name, 
		custom_f function,
		const map<string, int>& var_names, const vector<string>& param_names);

	enum class filling_type
	{
		common,		// fill fv as usual
		only_vars,	// fill only fv without differentials
		only_diff	// fill only differentials
	};
	template<typename T>
	void fill_fv_equations(const filling_type& f_type, vector<T*>& fv_a, int i, bool compute_differential = true, const T* x_ = nullptr);

	void ReadBoundaries(string file_name);		///< Reading Boundary file
	void InitFlow(double rho, double mass, double e, double p2);
	void InitFlow(double* rho_, double* mass_, double* e_, double* p_, double x_);
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

	void UpdateP(int i);

	void StagesMemoryAlloc_Init(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage);
	void RHSProcessing(vector < vector < vector < double > > >& rhsstage, int rks, double physDt, vector < vector < double > >& rhsold);
	void New_Steady_UnsteadyDualTime_CV(vector < vector < vector < double > > >& cvstage, int rks);
	void New_UnsteadyTwoStep_CV(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, int rks);
	void CVProcessing(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, int rks, DiagonalFunc D_Func = NULL, LUFunc L_Func = NULL, LUFunc U_Func = NULL);
	void FinalTSRK_CV_Calculation(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage);
	void FinalImplicitRHS_CV_Calculation(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage);

	void GetWn(vector < vector < double > >& Wn_);
	MatrixXd Diagonal(vector < MatrixXd >& LeftFluxJac, vector <MatrixXd >& RightFluxJac, vector < MatrixXd >& SourceJac, vector < vector < double > >& cv_, int i, LUFunc L_Func, LUFunc U_Func);
	MatrixXd Lower(std::vector < MatrixXd >& LeftFluxJac, std::vector < std::vector < double > >& cv_, int i);
	MatrixXd Lower2(std::vector < MatrixXd >& LeftFluxJac, std::vector < std::vector < double > >& cv_, int i);
	MatrixXd Upper(std::vector < MatrixXd >& RightFluxJac, std::vector < std::vector < double > >& cv_, int i);
	MatrixXd Upper2(std::vector < MatrixXd >& RightFluxJac, std::vector < std::vector < double > >& cv_, int i);
	void LUSGS(DiagonalFunc D_Func, LUFunc L_Func, LUFunc U_Func, int rks, vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, double physDt_ = 1.);
	void GMRES(DiagonalFunc D_Func, LUFunc L_Func, LUFunc U_Func, int rks, vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, double physDt_ = 1.);
	SparseMatrix< double > ILU_0(SparseMatrix< double >& SM);

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
	void TimeSteps(bool local_time = true, double dt_ = 0.);
	double SpectralRadius(vector< vector < double > >& cv_, int i);
	void SourceTerm(int i);
	void ImplResidualSmooth();
	void ImplResidualSmooth(vector < vector < vector < double > > >& rhs);
	double Convergence();
	void PrintResult();
	void SetOutputFile(string output_file_);
	void ShowOutput();
	void HideOutput();

	void InitFirsStages(vector < vector < vector < double > > >& cvstage, vector < vector < vector < double > > >& rhsstage);

	double CalcH(double p_, double Rho_, double U_);

	//void FillJacobian(vector < vector < double > >& M_SGS, vector < double >& jac, double s);
	void FillJacobian(MatrixXd& M_SGS, vector < double >& jac, double s);

	vector<adept::adouble> construct_side_flux_array(const vector<adept::adouble*>& vars, const int i);

	virtual void deactivate_adept_stack() = 0;

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
	double free_g;

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
	double T_b_in, T_b_out;

	int outflow_id;
	double p_b_out;

	int RK_stages_num;					///< number of Runge-Kutta stages

	vector < vector < double > > ak;
	vector < double > bk;
	vector < double > ck;
	vector < vector < double > > alpha;

	template<typename T>
	T make_equation_general(vector<eq_term>& eq_terms, int point, const vector<T*>& field_var, const T* cons_var);
	template<typename T>
	T get_var_value(const string& var_name_, const int point, eq_term::var_type& v_type, int& v_id, const vector<T*>& field_var = vector<T*>(), const T* cons_var = nullptr, bool fv_eq = true);
	template<typename T>
	T calculate_term_value(eq_term& term, int point, const vector<T*>& field_var = vector<T*>(), const T* cons_var = nullptr, bool fv_eq = true);
	template<typename T>
	void compute_differential_var(int i, const vector<int>& skipped, vector<T*>& fv);
	template<typename T>
	void fill_fv_underneath(const filling_type& f_type, int i, vector<T*>& fv_new, vector<int>& skipped, bool compute_differential = true, const T* cons_var = nullptr);
	template<typename T>
	void apply_operation(T& term, operation& op, const T& value);
};

Solver* CreateReadConfigFile(string file_name);

double Sign(double value);

//adept::adouble aSign(adept::adouble value);

int count_(vector < int >& vec, int val);