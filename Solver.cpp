// Solver

#include "Solver.h"
#include "CDS.h"
#include "CUSP.h"
#include "VL.h"
#include "HLLE.h"
//#include "Eigen/Dense"
//
//using namespace Eigen;

//double_type aSign(double_type value);

using namespace YAML;

set < string > solvers { "cds", "cusp", "vanleer", "hlle" };

Solver::Solver(sol_struct& sol_init_) :
	Ideal((sol_init_.gas == "Ideal" ? true : false)),
	Cp(sol_init_.Cp),
	cfl(sol_init_.cfl),
	max_iter_num(sol_init_.max_iter_num),
	tolerance(sol_init_.tolerance),
	gamma(sol_init_.gamma),
	free_g(sol_init_.free_g),
	RK_stage_coeffs(sol_init_.RK_stage_coeffs),
	RK_alpha(sol_init_.RK_alpha),
	RK_beta(sol_init_.RK_beta),
	res_smooth_flag(sol_init_.res_smooth_flag),
	eps_impl_res_smooth(sol_init_.eps_impl_res_smooth),
	ak(sol_init_.ak),
	bk(sol_init_.bk),
	ck(sol_init_.ck),
	alpha(sol_init_.alpha),
	time_expl(sol_init_.time_expl),
	TSRK_teta(sol_init_.TSRK_teta),
	TSRK_d(sol_init_.TSRK_d),
	TSRK_eta(sol_init_.TSRK_eta),
	TSRK_q(sol_init_.TSRK_q),
	lts(sol_init_.lts),
	remesh(sol_init_.remesh),
	steadiness(sol_init_.steadiness),
	time_stepping(sol_init_.time_stepping),
	solver_name(sol_init_.solver_name)
{
	if (false)
	{
		make_fv_equation<double_type>(var_name[0], 0);	// To resolve double_type template in HLLE. Whaaat?? It even may be not executed.
		make_equation<double>(0, Solver::equation::term_name::dt, vector<double*>());
		vector<double_type*> b;
		fill_fv_equations<double_type>(filling_type::common, b, 0);
	}

	eq_num = vars.size();
	var_num = vars_o.size();
	iter = 0;

	NoOutput = false;

	if (Ideal) {
		cout << "\tgas: " << "Ideal" << endl;

	}
	else {
		cout << "\tgas: " << "not Ideal" << endl;
	}

	cout << "\tCp: " << Cp << endl;
	cout << "\tcfl: " << cfl << endl;
	cout << "\tmax_iter_num: " << max_iter_num << endl;
	cout << "\ttolerance: " << tolerance << endl;
	cout << "\tgamma: " << gamma << endl;
	cout << "\teps_impl_res_smooth: " << eps_impl_res_smooth << endl;

	cout << "\tres_smooth_flag: ";
	for (auto flag_ : res_smooth_flag)
		cout << flag_ << "   ";

	cout << endl;

	cout << "\tRK_stage_coeffs: ";
	for (auto coef_ : RK_stage_coeffs)
		cout << coef_ << "   ";

	cout << endl;
	cout << (time_expl ? "Explicit solver\n" : "Implicit solver\n");
	cout << (steadiness ? "Steady problem\n" : "Unsteady problem\n");

	RK_stages_num = RK_stage_coeffs.size();
}

void Solver::AddDelayedFvEquation(const set<string>& equations_to_delay)
{
	delayed_fv_equations = equations_to_delay;
}

void Solver::SetEquation(string eq_name, const vector<string>& dt_term_, const vector<string>& dx_term_, const vector<string>& source_term, map < string, int > vars_, map < string, int > vars_o_)
{
	equations.push_back( equation (eq_name, dt_term_, dx_term_, source_term, vars_, vars_o_) );
}

void Solver::set_fv_equation(const string& eq_name, const vector<string>& eq_terms_s)
{
	vector<eq_term> eq_terms;

	eq_terms.reserve(eq_terms_s.size());

	for (auto& term_s : eq_terms_s)
		if (term_s.size() == 0)
			continue;
		else
			eq_terms.push_back(term_s);

	set_fv_equation(eq_name, eq_terms);
}

void Solver::set_function(const string& func_name, 
	custom_f function,
	const map<string, int>& var_names, const vector<string>& param_names)
{
	functions.insert(pair<string, custom_func>(func_name, custom_func(func_name, function, var_names, param_names)));
	set_fv_equation(
		func_name,
		{ func_name, "" }
	);
}

void Solver::set_fv_equation(const string& eq_name, const vector<eq_term>& eq_terms)
{
	fv_equation.insert(pair<string, vector<eq_term>>(eq_name, eq_terms));
}

pair<string, vector<eq_term>> Solver::equation::get_equation(const string& eq_name, const vector<string>& eq_terms_s)
{
	vector<eq_term> eq_terms;

	eq_terms.reserve(eq_terms_s.size());

	for (auto& term_s : eq_terms_s)
		if (term_s.size() == 0)
			continue;
		else
			eq_terms.push_back(term_s);

	pair<string, vector<eq_term>> result(eq_name, eq_terms);
	return result;
}

pair<string, vector<eq_term>> Solver::equation::get_equation(const string& eq_name, const vector<eq_term>& eq_terms)
{
	return pair<string, vector<eq_term>>(eq_name, eq_terms);
}

template<typename T>
T Solver::make_equation(const int eq, const equation::term_name term_name, const vector<T*>& f_vars)
{
	vector<eq_term> &eq_terms = term_name == equation::term_name::dt ? equations[eq].cur_dt.second :
								term_name == equation::term_name::dx ? equations[eq].cur_dx.second :
																	   equations[eq].cur_source.second;

	if (eq_terms.empty())
		return 0.;

	return make_equation_general<T>(eq_terms, -1, f_vars, nullptr); 	// fv anc cv will not be queried. Conservative variables will not be queried, but set, if it exists -> arrays are provided
}

template<typename T>
void Solver::apply_operation(T& term, operation& op, const T& value)
{
	switch (op)
	{
	case operation::plus:
		term += value;
		break;
	case operation::minus:
		term -= value;
		break;
	case operation::mult:
		term *= value;
		break;
	case operation::div:
		term /= value;
		break;
	default:
		break;
	}
}

template<typename T>
T Solver::get_var_value (const string& var_name_, const int point, eq_term::var_type& v_type, int& v_id, const vector<T*>& field_var, const T* cons_var, bool fv_eq)
{
	bool arrays_are_provided = cons_var || !fv_eq;
	if (v_type == eq_term::var_type::not_defined)
	{
		if (var_name_ == var_name[g_X])
		{
			v_type = eq_term::var_type::coord;
			return x[point];
		}
		if (var_name_ == var_name[g_A])
		{
			v_type = eq_term::var_type::area;
			return (arrays_are_provided ? 1. : a[point]);
		}
		if (var_name_ == "GAMMAM")
		{
			v_type = eq_term::var_type::gammam;
			return gamma - 1.;
		}
		if (var_name_ == "GAMMA")
		{
			v_type = eq_term::var_type::gamma;
			return gamma;
		}
		if (var_name_ == "CP")
		{
			v_type = eq_term::var_type::cp;
			return Cp;
		}
		if (var_name_ == "g")
		{
			v_type = eq_term::var_type::g;
			return free_g;
		}
		if (var_name_ == "k")
		{
			v_type = eq_term::var_type::k_planc;
			return k_planc;
		}
		if (vars.find(var_name_) != vars.end())
		{
#ifdef USE_OMP
#pragma omp critical
#endif
			{
				v_id = vars[var_name_];
				v_type = eq_term::var_type::conservative;
			}
			return (arrays_are_provided ? cons_var[v_id] : cv[v_id][point]);
		}
		else if (functions.find(var_name_) != functions.end())
		{
			v_type = eq_term::var_type::function;
			v_id = func_name_ids.size();
			func_name_ids[var_name_] = v_id;
			custom_func& call = functions.at(var_name_);

			vector<double_type> field_var_val;
			if (arrays_are_provided)
			{
				field_var_val.resize(field_var.size());
				for (int i = 0; i < field_var.size(); ++i)
					field_var_val[i] = *field_var[i];
			}

			return call.get_function_value<T>(fv, field_var_val, point, Global_Time);
		}
		else
		{
#ifdef USE_OMP
#pragma omp critical
#endif
			{
				v_id = vars_o[var_name_];
				v_type = eq_term::var_type::field;
			}
			return (arrays_are_provided ? *field_var[v_id] : fv[v_id][point]);
		}
	}
	else
	{
		switch (v_type)
		{
		case eq_term::var_type::coord:
			return x[point];
		case eq_term::var_type::area:
			return (arrays_are_provided ? 1. : a[point]);
		case eq_term::var_type::gammam:
			return gamma - 1.;
		case eq_term::var_type::gamma:
			return gamma;
		case eq_term::var_type::cp:
			return Cp;
		case eq_term::var_type::g:
			return free_g;
		case eq_term::var_type::k_planc:
			return k_planc;
		case eq_term::var_type::conservative:
			return arrays_are_provided ? cons_var[v_id] : cv[v_id][point];
		case eq_term::var_type::field:
			return arrays_are_provided ? *field_var[v_id] : fv[v_id][point];
		case eq_term::var_type::function:
		{
			custom_func& call = functions.at(var_name_);

			vector<double_type> field_var_val;
			if (arrays_are_provided)
			{
				field_var_val.resize(field_var.size());
				for (int i = 0; i < field_var.size(); ++i)
					field_var_val[i] = *field_var[i];
			}

			return call.get_function_value<T>(fv, field_var_val, point, Global_Time);
		}
		}
	}
};

template<typename T>
T Solver::calculate_term_value(eq_term& term, int point, const vector<T*>& field_var, const T* cons_var, bool fv_eq)
{
	operation op = term.op;
	string& var_name_ = term.name;
	auto& v_type = term.v_type;
	auto& v_id = term.var_id;
	T value = get_var_value(var_name_, point, v_type, v_id, field_var, cons_var, fv_eq);

	double degree = term.degree;
	double coef = term.coef;
	if (fabs(degree - 1.) > 1e-5)
		value = pow(value, degree);
		
	value *= coef;
	return value;
};

template<typename T>
T Solver::make_fv_equation(const string& eq_name, const int point, const vector<T*>& field_var, const T* cons_var)
{
	if (fv_equation.empty())
		return 0.;

	auto& eq_terms = fv_equation.at(eq_name);
	return make_equation_general<T>(eq_terms, point, field_var, cons_var);
}

template<typename T>
T Solver::make_equation_general(vector<eq_term>& eq_terms, int point, const vector<T*>& field_var, const T* cons_var)
{
	bool fv_eq = point != -1;			// If point isn't -1 - then it is a field value equation, otherwise it is a conservative equation
	operation op = eq_terms[0].op;

	T term = calculate_term_value<T>(eq_terms[0], point, field_var, cons_var, fv_eq);
	if (op != operation::plus)
		term = -term;

	auto eq_terms_it = eq_terms.begin();
	++eq_terms_it;
	while (eq_terms_it != eq_terms.end())
	{
		op = (*eq_terms_it).op;
		T value = calculate_term_value<T>(*eq_terms_it, point, field_var, cons_var, fv_eq);
		apply_operation(term, op, value);
		++eq_terms_it;
	}
	return term;
}

template<typename T>
void Solver::fill_fv_equations(const filling_type &f_type, vector<T*>& fv_a, int i, bool compute_differential, const T* x_)
{
	vector<int> skipped;

	fill_fv_underneath<T>(f_type, i, fv_a, skipped, compute_differential, x_);
	if (compute_differential && f_type != filling_type::only_vars)
	{
		compute_differential_var<T>(i, skipped, fv_a);
		for (auto it = delayed_fv_equations.begin(); it != delayed_fv_equations.end() && f_type == filling_type::common; ++it)
		{
			string v_name = *it;
			*fv_a[vars_o[v_name]] = make_fv_equation<T>(v_name, i, fv_a, x_);
		}
	}
}

template<typename T>
void Solver::compute_differential_var(int i, const vector<int>& skipped, vector<T*>& fv)
{
	int left_id = max(0, i - 1);
	int right_id = min(imax - 1, i + 1);
	double alpha1(0.);
	double alpha2(0.);
	vector<double> fv_left, fv_right;
	vector<double*> fv_left_ref, fv_right_ref;

	if (!skipped.empty())
	{
		fv_left.resize(var_num);
		fv_right.resize(var_num);
		fv_left_ref.resize(var_num);
		fv_right_ref.resize(var_num);
		for (int j = 0; j < var_num; ++j)
		{
			fv_left_ref[j] = &fv_left[j];
			fv_right_ref[j] = &fv_right[j];
		}
		if (fabs(x[right_id] - x[left_id]) > 1e-20)
		{
			alpha1 = (x[right_id] - x[i]) / (x[right_id] - x[left_id] + 1e-20) / (x[i] - x[left_id] + 1e-20);
			alpha2 = (x[i] - x[left_id]) / (x[right_id] - x[left_id] + 1e-20) / (x[right_id] - x[i] + 1e-20);
		
			fill_fv_equations<double>(filling_type::only_vars, fv_left_ref, left_id, false);		// Extra calculations if several d vars...
			fill_fv_equations<double>(filling_type::only_vars, fv_right_ref, right_id, false);
		}
	}

	for (const int& var : skipped)
	{
		string name = var_name[var];
		name.erase(name.begin());
		int prev_var = vars_o[name];

		*fv[var] = (*fv[prev_var] - fv_left[prev_var]) * alpha1 + (fv_right[prev_var] - *fv[prev_var]) * alpha2;
	}
}

template<typename T>
void Solver::fill_fv_underneath(const filling_type& f_type, int i, vector<T*>& fv_new, vector<int>& skipped, bool compute_differential, const T* cons_var)
{
	for (int var = 0; var < var_num; ++var)
	{
		string& v_name = var_name[var];
		bool is_diff = v_name[0] == 'd';
		if (is_diff && f_type != filling_type::only_vars)
		{
			if (compute_differential)
			{
				skipped.push_back(var);
			}
			continue;
		}
		if (!is_diff && f_type != filling_type::only_diff)
		{
			if (f_type == filling_type::only_vars && delayed_fv_equations.contains(v_name))
				continue;

			if (var == g_A)
				*fv_new[var] = a[i];
			else if (var == g_X)
				*fv_new[var] = x[i];
			else
				*fv_new[var] = make_fv_equation<T>(v_name, i, fv_new, cons_var);
		}
	}
}

void Solver::ReadBoundaries(string file_name)
{
	cout << "Opening file " << "\"" << file_name.c_str() << "\"..." << endl;
	Node config = LoadFile(file_name.c_str());

	if (config["type"].as< string >() == "dirichlet")
	{
		b_type = "dirichlet";
		if (config["first"]) {
			if (config["first"]["inflow"]) {
				inflow_id = 0;
				p_b_in = config["first"]["inflow"]["p01"].as< double >();
				T_b_in = config["first"]["inflow"]["t01"].as< double >();

				cout << "Inflow total pressure: " << p_b_in << endl;
				cout << "Inflow temperature: " << T_b_in << endl;

			}
			else if (config["first"]["outflow"]) {
				outflow_id = 0;
				p_b_out = config["first"]["outflow"]["p2"].as< double >();
				cout << "Outflow pressure: " << p_b_out << endl;

			}
		}

		if (config["last"]) {
			if (config["last"]["inflow"]) {
				inflow_id = imax;
				p_b_in = config["last"]["inflow"]["p01"].as< double >();
				T_b_in = config["last"]["inflow"]["t01"].as< double >();

				cout << "Inflow total pressure: " << p_b_in << endl;
				cout << "Inflow temperature: " << T_b_in << endl;

			}
			else if (config["last"]["outflow"]) {
				outflow_id = imax;
				p_b_out = config["last"]["outflow"]["p2"].as< double >();
				cout << "Outflow pressure: " << p_b_out << endl;

			}
		}
	}
	if (config["type"].as< string >() == "transmissive")
	{
		b_type = "transmissive";
		inflow_id = 0;
		outflow_id = imax;
	}
	if (config["type"].as< string >() == "walls")
	{
		b_type = "walls";
		if (config["first"]) {
			inflow_id = 0;
			T_b_in = config["first"]["T"].as< double >();
			cout << "First wall temperature: " << T_b_in << endl;
		}
		if (config["last"]) {
			outflow_id = imax;
			T_b_out = config["last"]["T"].as< double >();
			cout << "Last wall temperature: " << T_b_out << endl;
		}

				
	}
}

void Solver::InitFlow(double rho, double mass, double e, double p2)
{
	cout << "Initializing flow" << endl;

	cv.resize(eq_num);
	cvold.resize(eq_num);
	diss.resize(eq_num);
	rhs.resize(eq_num);
	Q_star.resize(eq_num);
	for (int i = 0; i < eq_num; ++i) {
		cv[i].resize(imax, 0.);
		cvold[i].resize(imax, 0.);
		diss[i].resize(imax, 0.);
		rhs[i].resize(imax, 0.);
		Q_star[i].resize(imax, 0.);
	}
	p.resize(imax, 0.);
	dt.resize(imax, 0.);
	dummy.resize((imax + 1) * (var_num + 1), 0.);

	for (int i = 0; i < imax; ++i)
	{
		cv[RHO_A][i] = rho * a[i];		//	$\rho S$
		cv[RHO_U_A][i] = mass;					//	$\rho u S$
		cv[RHO_E_A][i] = rho * e * a[i];	//	$\rho E S$
		p[i] = p2;
	}

	//if (bk.size() > 0) {
	if (!time_expl) {
		L_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
			/*vector < vector < double > > (eq_num, 
				vector < double > (eq_num, 0)));*/

		U_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
			/*vector < vector < double > > (eq_num, 
				vector < double > (eq_num, 0)));*/

		D_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
			/*vector < vector < double > > (eq_num,
				vector < double > (eq_num, 0)));*/
	}
}

void Solver::InitFlow(double* rho_, double* mass_, double* e_, double* p_, double x_)
{
	cout << "Initializing flow" << endl;

	cv.resize(eq_num);
	cvold.resize(eq_num);
	diss.resize(eq_num);
	rhs.resize(eq_num);
	Q_star.resize(eq_num);
	for (int i = 0; i < eq_num; ++i) {
		cv[i].resize(imax, 0.);
		cvold[i].resize(imax, 0.);
		diss[i].resize(imax, 0.);
		rhs[i].resize(imax, 0.);
		Q_star[i].resize(imax, 0.);
	}
	p.resize(imax, 0.);
	dt.resize(imax, 0.);
	dummy.resize((imax + 1) * (var_num + 1), 0.);

	int id;

	for (int i = 0; i < imax; ++i)
	{
		id = 0;
		if (x[i] > x_) { id = 1; }
		cv[RHO_A][i] = rho_[id] * a[i];		//	$\rho S$
		cv[RHO_U_A][i] = mass_[id];					//	$\rho u S$
		cv[RHO_E_A][i] = rho_[id] * e_[id] * a[i];	//	$\rho E S$
		p[i] = p_[id];
	}

	//if (bk.size() > 0) {
	if (!time_expl) {
		L_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
		/*vector < vector < double > > (eq_num,
			vector < double > (eq_num, 0)));*/

		U_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
		/*vector < vector < double > > (eq_num,
			vector < double > (eq_num, 0)));*/

		D_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
		/*vector < vector < double > > (eq_num,
			vector < double > (eq_num, 0)));*/
	}
}

void Solver::ResetDummy()
{
	dummy.assign(dummy.size(), 0.);
}

void Solver::RefreshBoundaries()
{
	double rgas;
	double u;
	double cs2;
	double c02;
	double rinv;
	double dis;
	double cb;
	double cc02;
	double tb;
	double pb;
	double rhob;
	double ub;
	double rho;
	double cs;

	int in_id = GetInflowId();
	in_id = (in_id > 0 ? in_id - 1 : in_id);
	int in_id_p1 = (in_id > 0 ? in_id - 1 : in_id + 1);

	int out_id = GetOutflowId();
	out_id = (out_id > 0 ? out_id - 1 : out_id);
	int out_id_p1 = (out_id > 0 ? out_id - 1 : out_id + 1);

	double gam1 = gamma - 1.;
	double gap1 = gamma + 1.;

	if (b_type == "dirichlet") {
		// Inflow boundaries
		rgas = gam1 * Cp / gamma;
		u = fv[g_U][in_id_p1]; // cv[1][in_id_p1] / cv[0][in_id_p1];		///< $u = \frac{(\rho u S)_{[1]}} {(\rho  S)_{[1]}}$
		cs2 = gamma * fv[g_P][in_id_p1] / fv[g_RHO][in_id_p1]; // p[in_id_p1] * a[in_id_p1] / cv[0][in_id_p1];	// Speed of sound $c_s^2 = \frac{\gamma_ p_{[1]} S_{[1]}} { \rho  S_{[1]}}$
		c02 = cs2 + 0.5 * gam1 * u * u;			// Stagnation speed of sound$c_0^2 = c_s^2 + 0.5 (\gamma_ - 1) u^2$
		rinv = abs(u) - 2. * sqrt(cs2) / gam1;		// Riemann invariant $R^-$, $r_{inv} = u - 2 \frac{\sqrt{c_s^2}}{\gamma_ - 1} $ 		// fix is here
		dis = gap1 * c02 / (gam1 * rinv * rinv) - 0.5 * gam1;	// $dis = \frac{(\gamma_ + 1) c_0^2}{(\gamma_ - 1) r_{inv}^2}$
		if (dis < 0)
			dis = 1e-20;

		cb = -rinv * (gam1 / gap1) * (1 + sqrt(dis));	// Boundary speed of sound $c_b = -r_{inv} \frac{\gamma_ - 1}{\gamma_ + 1} (1 + \sqrt{dis})$
		cc02 = min(cb * cb / c02, 1.);		// It is like cb shoud be not greater than 1. // stagnation speed of sound, $cc_0^2 = min\left(\frac{cb^2}{c_0^2}, 2\right)$
		tb = T_b_in * cc02;			// $T_b = T_{01} cc_0^2$
		pb = p_b_in * pow((tb / T_b_in), gamma / gam1);	// $p_b = p_{01} \left(\frac{T_b}{T_{01}}\right)^\frac{\gamma_}{\gamma_ - 1}$
		rhob = pb / (rgas * tb);		//  $\rho_b = \frac{p_b}{r_{gas} T_b}$
		ub = Sign(u) * sqrt(2. * Cp * (T_b_in - tb));		// $u_b = \sqrt{2 c_{p, gas} \left(T_{01} - T_b\right)}$ 		// fix is here

		cv[g_RHO_A][in_id] = rhob * a[in_id_p1];		//	$\rho S = \rho_b S_{[1]}$
		cv[g_RHO_U_A][in_id] = rhob * a[in_id_p1] * ub;					//	$\rho u S = \rho_b S_{[1]} u_b$
		cv[g_RHO_E_A][in_id] = (pb / gam1 + 0.5 * rhob * ub * ub) * a[in_id_p1];	//	$\rho E S = \left(\frac{p_b}{\gamma_ - 1} + 0.5 \rho_b u_b^2\right) S_{[1]}$
		p[in_id] = pb;

		// Outflow boundaries

		

		rho = fv[g_RHO][out_id_p1]; // cv[0][out_id_p1] / a[out_id_p1];
		u = fv[g_U][out_id_p1]; // cv[1][out_id_p1] / cv[0][out_id_p1];
		cs = sqrt(gamma * p[out_id_p1] / rho);

		if (abs(u) >= cs) {		// supersonic flow 		// fix is here
			pb = p[out_id_p1];
			rhob = rho;
			ub = u;
		}
		else {				//subsonic flow
			pb = p_b_out;
			rhob = rho + Sign(u) * (p_b_out - p[out_id_p1]) / (cs * cs);		// $\rho_b = \rho + \frac{(p_2 - p_{[ib2 - 1]})}{c_s^2}$ 		// fix is here
			ub = u - Sign(u) * (p_b_out - p[out_id_p1]) / (cs * rho);		// $u_b = u - \frac{(p_2 - p_{[ib2 - 1]})}{c_s \rho}$ 		// fix is here
		}

		cv[g_RHO_A][out_id] = rhob * a[out_id_p1];
		cv[g_RHO_U_A][out_id] = rhob * ub * a[out_id_p1];
		cv[g_RHO_E_A][out_id] = (pb / gam1 + 0.5 * rhob * ub * ub) * a[out_id_p1];
		p[out_id] = pb;
	}
	if (b_type == "transmissive")
	{
		for (int eq = 0; eq < eq_num; ++eq)
			cv[eq][in_id] = cv[eq][in_id_p1];

		p[in_id] = (gamma - 1.) * (cv[g_RHO_E_A][in_id] / a[in_id] - 0.5 * cv[g_RHO_U_A][in_id] * cv[g_RHO_U_A][in_id] / (cv[g_RHO_A][in_id] * a[in_id]));

		for (int eq = 0; eq < eq_num; ++eq)
			cv[eq][out_id] = cv[eq][out_id_p1];

		p[out_id] = (gamma - 1.) * (cv[g_RHO_E_A][out_id] / a[out_id] - 0.5 * cv[g_RHO_U_A][out_id] * cv[g_RHO_U_A][out_id] / (cv[g_RHO_A][out_id] * a[out_id]));
	}
	if (b_type == "walls")
	{
		constexpr double m_p_re?? = 5.9786e26;	// Reciprocal protan mass
		cv[g_RHO_A][in_id] = cv[g_RHO_A][in_id_p1];
		cv[g_RHO_U_A][in_id] = 0.;
		cv[g_RHO_U_A][in_id_p1] = 0.;

		p[in_id] = 2. * m_p_re?? * cv[g_RHO_A][in_id] / a[in_id] * k_planc * T_b_in;
		cv[g_RHO_E_A][in_id] = p[in_id] * a[in_id] * gam1;
		cv[g_RHO_E_A][in_id_p1] = cv[g_RHO_E_A][in_id];

		cv[g_RHO_A][out_id] = cv[g_RHO_A][out_id_p1];
		cv[g_RHO_U_A][out_id] = 0.;
		cv[g_RHO_U_A][out_id_p1] = 0.;

		p[out_id] = 2. * m_p_re?? * cv[g_RHO_A][out_id] / a[out_id] * k_planc * T_b_out;
		cv[g_RHO_E_A][out_id] = p[out_id] * a[out_id] * gam1;
		cv[g_RHO_E_A][out_id_p1] = cv[g_RHO_E_A][out_id];
	}
	//cv[N_A][0] = cv[N_A][1];
	//cv[N_A][imax - 1] = cv[N_A][ib2 - 1];
}

void Solver::SetGasType(bool gas)
{
	Ideal = gas;
}

void Solver::SetCp(double Cp_)
{
	Cp = Cp_;
}

bool Solver::IsIdeal()
{
	return Ideal;
}

void Solver::SetCFL(double cfl_)
{
	cfl = cfl_;
}

double Solver::GetCFL()
{
	return cfl;
}

void Solver::SetMaxIterNum(int max_iter_num_)
{
	max_iter_num = max_iter_num_;
}

void Solver::SetTolerance(double tolerance_)
{
	tolerance = tolerance_;
}

void Solver::CalculateVolumes()
{
	vol.resize(imax, 0.);

	for (int i = 1; i < ib2; ++i)
	{
		vol[i] = (x[i + 1] - x[i]) / 4. * (a[i] + a[i + 1]) + 
				 (x[i] - x[i - 1]) / 4. * (a[i] + a[i - 1]);
	}
	vol[0] = vol[1];
	vol[imax - 1] = vol[ib2 - 1];
}

double Solver::GetGamma()
{
	return gamma;
}

double Solver::GetCp()
{
	return Cp;
}

double Solver::GetInflowT()
{
	return T_b_in;
}

double Solver::GetInflowP()
{
	return p_b_in;
}

double Solver::GetOutflowP()
{
	return p_b_out;
}

int Solver::GetInflowId()
{
	return inflow_id;
}

int Solver::GetOutflowId()
{
	return outflow_id;
}

double Solver::GetTolerance()
{
	return tolerance;
}

vector < double > Solver::GetSection()
{
	return grid.GetSection();
}

int Solver::GetMaxIterNum()
{
	return max_iter_num;
}

void Solver::SetGrid(Loop& loop)
{
	grid.SetData(loop.GetTable(), loop.GetColumnNames());

	x = grid.GetCoordinates();
	a = grid.GetSection();
	y.resize(x.size());
	for (int i = 0; i < y.size(); ++i)
	{
		if (i == 0 || i == y.size() - 1)
			y[i] = x[i];
		else
			y[i] = 0.5 * (x[i] + 0.5 * (x[i - 1] + x[i + 1]));
	}

	imax = x.size();
	ib2 = imax - 1;
	ncells = imax - 3;

	CalculateVolumes();

	grid.AddColumn("volume", vol);
}

void Solver::InverseGrid()
{
	grid.InverseData("coordinate");

	x = grid.GetCoordinates();
	a = grid.GetSection();
	vol = grid.GetValues("volume");
}

void Solver::ShowGrid()
{
	grid.ShowData();
}

void Solver::CalculateTimeSource(vector < vector < double > >& cvn_, vector < vector < double > >& cvnm1_, double physDt_)
{
	for (int i = 1; i < ib2; ++i)
	{
		for (int eq = 0; eq < eq_num; ++eq) {
			Q_star[eq][i] = (2. / physDt_ * vol[i] * cvn_[eq][i] - 1. / (2. * physDt_) * vol[i] * cvnm1_[eq][i]);
		}
	}
}

void Solver::CalculateUnsteadyRHS(double physDt_)
{
	for (int i = 1; i < ib2; ++i)
	{
		for (int eq = 0; eq < eq_num; ++eq) {
			//rhs[eq][i] += 3. / (2. * physDt_) * vol[i] * cv_[eq][i] - Q_star_[eq][i];
			rhs[eq][i] = (rhs[eq][i] + (3. * vol[i] * cv[eq][i] - 4. * vol[i] * cvn[eq][i] + vol[i] * cvnm1[eq][i]) / (2. * physDt_));	// Good one
			//Some variants for little physDt (not working)
			/*//rhsstage[rks][eq][i] = (rhs[eq][i] - 3. / (2. * physDt) * vol[i] * beta * cv[eq][i]) / (1 + 3. / (2. * physDt) * RK_beta[rks][rks] * cfl * dt[i] * beta);
			//rhsstage[rks][eq][i] = (rhs[eq][i] + (3. * vol[i] * cvstage[rks][eq][i] - 4. * vol[i] * cvn[eq][i] + vol[i] * cvnm1[eq][i] - 3. * vol[i] * beta * cvstage[rks][eq][i]) / (2. * physDt)) /
			//	(1 + 3. / (2. * physDt) * RK_beta[rks][rks] * cfl * dt[i] * beta);
			//rhsstage[rks][eq][i] = (rhs[eq][i] + (3. * vol[i] * cvstage[0][eq][i] - 4. * vol[i] * cvn[eq][i] + vol[i] * cvnm1[eq][i]) / (2. * physDt)) /
			//	(1 + 3. / (2. * physDt) * RK_beta[rks][rks] * cfl * dt[i] /** beta* /);*/
		}
	}
}

double Solver::Solve(double physDt)
{
	if (!steadiness && !time_expl && time_stepping == 1) {
		cout << "There is no implicit unsteady two step Runge-Kutta scheme" << endl;
		return -404.;
	}
	return SolveExplImpl(physDt);
}

double Solver::ForwardEuler(double physDt)
{
	double fac, adtv;
	double Time_ = 0.;
	double Euler_CFL = 0.25;
	double DT_;

	while (Time_ < physDt) 
	{
		iter++;

		// Store previous solution; set dissipation = 0
		//cvold = cv;

		// Calculate time step
		TimeSteps(false);

		if (dt[0] * Euler_CFL > physDt - Time_) {
			DT_ = (physDt - Time_) / Euler_CFL;
		}
		else {
			DT_ = dt[0];
		}

		RhoUPH();
		LRState("");
		ComputeRHSandJacobian();

		fac = Euler_CFL;
		for (int i = 1; i < ib2; ++i)
		{
			adtv = fac * DT_ / vol[i];
			for (int eq = 0; eq < eq_num; ++eq) {
				rhs[eq][i] += adtv * rhs[eq][i];		// Calculating of minus part in Eq. 6.7
			}
		}

		// implicit residual smoothing
		//if (res_smooth_flag[rks] > 0 && eps_impl_res_smooth > 0.)
			ImplResidualSmooth();

			// update (conservative variables and pressure)
		for (int i = 1; i < ib2; ++i)
		{
			for (int eq = 0; eq < eq_num; ++eq) {
				cv[eq][i] -= rhs[eq][i];		// Calculating of conservative variables (Eq. 6.7)
			}
			UpdateP(i);
		}

		// boundary conditions

		RhoUPH();
		RefreshBoundaries();

		Time_ += Euler_CFL * DT_;

		cout << Time_ << endl;
	}

	return 0.; // Convergence();
}

void Solver::UpdateP(int i)
{
	p[i] = make_fv_equation<double>(var_name[g_P], i);
	p[i] = max(p[i], 1e-20);
}

void Solver::StagesMemoryAlloc_Init(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage)
{
	int add_st = 0;

	// Store previous solution; set dissipation = 0
	cvold = cv;

	if (time_expl == true) {
		
		if (!steadiness && time_stepping == 1) {
			RK_stages_num = TSRK_d.size() - 1;
			add_st = 2;
		}
		
	}
	else {
		if (bk.size() > 0) {
			RK_stages_num = bk.size();
		}
	}

	cvstage.resize(RK_stages_num + add_st,
		vector < vector < double > >(eq_num,
			vector < double >(imax, 0.)));

	rhsstage.resize(RK_stages_num + add_st,
		vector < vector < double > >(eq_num,
			vector < double >(imax, 0.)));

	if (time_expl == true) {
		if (!steadiness && time_stepping == 1) {
			InitFirsStages(cvstage, rhsstage);
		}
		else {
			cvstage[0] = cvold;
		}
	}
}

double Solver::SolveExplImpl(double physDt)
{
	double fac, adtv, rrho, rhou, rhoe;
	vector < int > diss_flag_;
	vector < double > diss_blend_;
	vector < vector < vector < double > > > cvstage;
	vector < vector < vector < double > > > rhsstage;
	vector < vector < double > > rhsold;
	bool SkipJacComputation;

	auto Diag = std::bind(&Solver::Diagonal, this, _1, _2, _3, _4, _5, _6, _7);
	auto Low = std::bind(&Solver::Lower2, this, _1, _2, _3);
	auto Up = std::bind(&Solver::Upper2, this, _1, _2, _3);

	iter++;
	if (solver_name == "cds") {
		diss_flag_ = GetDissFlag();
		diss_blend_ = GetDissBlend();
	}
	for (int i = 0; i < eq_num && solver_name == "cds"; ++i) {
		diss[i].assign(diss[i].size(), 0.);
	}

	// Calculate time step
	{if (!steadiness && time_stepping == 1)
		if (Global_Time == 0.) {
			TimeSteps(false);
			physDt = dt[0];
		}
		else {
			TimeSteps(false, physDt);
		}
	else
		//TimeSteps(lts, physDt/1.); }		// local = false!					***********   DO SOMETHING WITH THIS !!!   *****************
		TimeSteps(lts); }		// local = false!

	if (!steadiness && physDt < dt[0]) {			// **   IMPORTANT   **
		TimeSteps(lts, physDt / 1.);
	}

	StagesMemoryAlloc_Init(rhsstage, cvstage);
	RhoUPH();

	SkipJacComputation = false;
	// loop over the R. - K. stages:	Blazek, Section 6.1.2, p. 170, "Hybrid multistage schemes"
	for (int rks = 0; rks < RK_stages_num; ++rks)
	{
		if (solver_name == "cds" && diss_flag_[rks] == 1) {	// artificial dissipation
			Dissipation(diss_blend_[rks]);
		}
		else {
			LRState("");
		}
		if (time_expl == false && rks == 1)		{
			SkipJacComputation = true;		}

		ComputeRHSandJacobian(SkipJacComputation);					// convective flux residual
		RHSProcessing(rhsstage, rks, physDt, rhsold);		// RHS preparations, unsteady modification if needed
		CVProcessing(rhsstage, cvstage, rks, Diag, Low, Up);
	}

	if (time_expl == false) {
		LRState("");
		ComputeRHSandJacobian(true);
		RHSProcessing(rhsstage, RK_stages_num, physDt, rhsold);		// RHS preparations, unsteady modification if needed
		FinalImplicitRHS_CV_Calculation(rhsstage, cvstage);
		Global_Time += physDt;		// Is it true?
	} else if (!steadiness && time_stepping == 1) {	// If unsteady and explicit two step RK
		FinalTSRK_CV_Calculation(rhsstage, cvstage);
		Global_Time += physDt * cfl;
		return 0.;
	}
	else
		Global_Time += physDt;

	//Remesh();
	return Convergence();
}

void Solver::Remesh()
{
	double tau = RemeshTau;		// Manage this properly
	if (remesh && true) {
		vector <double> old_coords = grid.GetCoordinates();
		// New adaptive grid
		RhoUPH();

		for (int i = 0; i < 1; ++i) {
			for (int var = 0; var < eq_num; ++var)
				grid.SetRow(c_var_name[var], cv[var]);

			vector< vector <double> > functions;
			vector< double > Fs;

			functions.push_back(cvn[vars[RemeshVar]]);
			Fs.push_back(MaxFn);
			if (!RemeshFuncs.empty())
			{
				auto get_var_column = [this](string& col_name) -> vector<double>
				{
					if (vars.contains(col_name))
						return grid.GetValues(col_name);
					if (vars_o.contains(col_name))
						return fv[vars_o[col_name]];
				};

				int old_size = functions.size();
				functions.resize(old_size + RemeshFuncs.size());
				Fs.insert(Fs.end(), MaxOfRemeshFuncs.begin(), MaxOfRemeshFuncs.end());
				for (int i = 0; i < RemeshFuncs.size(); ++i)
					functions[old_size + i] = get_var_column(RemeshFuncs[i]);
			}

			grid.CalculateResolution(MaxX, MaxF, RemeshVar, "coordinate", functions, Fs);
			grid.CalculateConcentration(MaxX, "coordinate");

			vector <string> ignore;
			ignore.push_back("old_coords");
			x = grid.RefineMesh(/*dt[0], tau*/1., 10., 1., ignore);

			for (int var = 0; var < eq_num; ++var)
				cv[var] = grid.GetValues(c_var_name[var]);

			RefreshBoundaries();						// Refresh boundary conditions
			RhoUPH();

			CalculateVolumes();
			grid.SetRow("volume", vol);
		}
		// Now we will reevaluate cvn, cvnm1, cvold
		vector < vector < vector < double > >* > cvs{ &cvn, &cvnm1, &cvold };
		vector < vector < vector < double > > > cvs_old{ cvn_old, cvnm1_old, cvold };
		vector <vector <double> > old_coords_v;
		old_coords_v.push_back(grid.GetValues("old_coords"));
		old_coords_v.push_back(grid.GetValues("old_coords"));
		old_coords_v.push_back(old_coords);
		int i = 0;

		for (auto it : cvs)
		{
			grid.SetRow("coordinate", old_coords_v[i]);
			for (int var = 0; var < eq_num; ++var)
				grid.SetRow(c_var_name[var], cvs_old[i][var]);
			++i;

			vector < string > ignore;
			vector < vector < double > > new_tab;
			ignore.push_back("old_coords");
			ignore.push_back(grid.TYPE_COL);
			new_tab = grid.NewTable("coordinate", x, ignore, false);
			grid.SetData(new_tab);

			for (int var = 0; var < eq_num; ++var)
				(*it)[var] = grid.GetValues(c_var_name[var]);
		}

		grid.SetRow("coordinate", x);
		for (int var = 0; var < eq_num; ++var)
			grid.SetRow(c_var_name[var], cv[var]);

		calculate_mass_matrix();
		fill_inverse_mass_matrix();
	}
}

void Solver::RHSProcessing(vector < vector < vector < double > > >& rhsstage, int rks, double physDt, vector < vector < double > >& rhsold)
{
	double adtv;

	// If not steady-state and dual-time stepping: creating ensteady residual
	if (!steadiness) {
		CalculateUnsteadyRHS(physDt);
		if (time_expl == true) {
			if (time_stepping == 0) {
				rhsstage[rks] = rhs;
			}
			else {
				rhsstage[rks + 1] = rhs;
			}
		}
	}
	else {
		if (time_expl == true) {
			rhsstage[rks] = rhs;
		}
	}
	if (time_expl == false) {
		if (rks == 0) {
			rhsold = rhs;
		}
		else {
			rhsstage[rks - 1] = rhs;
			rhs = rhsold;
		}
	}
	// If steady or (unsteady and dual-time stepping)    Preparation of RHS for future CV calculation
	for (int i = 1; i < ib2 && (steadiness || time_stepping == 0) && time_expl == true; ++i)
	{
		for (int eq = 0; eq < eq_num; ++eq) {
			rhs[eq][i] = 0.;
			for (int stage = 0; stage <= rks; ++stage) {
				if (RK_beta[rks][stage] == 0.) continue;

				adtv = RK_beta[rks][stage] * cfl * dt[i] / vol[i];
				rhs[eq][i] += adtv * rhsstage[stage][eq][i];
			}
		}
	}

	// implicit residual smoothing
	if (res_smooth_flag[rks] > 0 && eps_impl_res_smooth > 0. && steadiness && time_expl == true)
		ImplResidualSmooth();
}

void Solver::New_Steady_UnsteadyDualTime_CV(vector < vector < vector < double > > >& cvstage, int rks)
{
	auto scalar_prod = [](const vector < double >::const_iterator& vec1, const vector < double >::const_iterator& vec2, int size) -> double
	{
		vector <double> prods(size);
		std::transform(vec1, vec1 + size, vec2, prods.begin(), std::multiplies<double>());
		return std::accumulate(prods.begin(), prods.end(), 0.);
	};

	for (int i = 1; i < ib2; ++i)
	{
		for (int eq = 0; eq < eq_num; ++eq) {
			cv[eq][i] = 0.;
			for (int stage = 0; stage <= rks; ++stage) {
				if (RK_alpha[rks][stage] == 0.) continue;

				cv[eq][i] += RK_alpha[rks][stage] * cvstage[stage][eq][i];
			}
			if (steadiness)
				cv[eq][i] -= rhs[eq][i];		// Calculating of conservative variables (Eq. 6.7)
			else
				cv[eq][i] -= scalar_prod(rhs[eq].begin() + 1, M_matrix.inversed_mass_matrix[i - 1].begin(), ib2 - 1);
		}
		UpdateP(i);
	}
}

void Solver::New_UnsteadyTwoStep_CV(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, int rks)
{
	double r_ = 1.;		// Don't understand this one

	for (int i = 1; i < ib2; ++i)
	{
		for (int eq = 0; eq < eq_num; ++eq) {
			cv[eq][i] = 0.;
			double q_ = 0.;
			for (int st_ = 0; st_ < rks + 2; ++st_) {
				q_ += TSRK_q[rks + 1][st_];
			}
			if (TSRK_d[rks + 1] != 0.)
				cv[eq][i] += TSRK_d[rks + 1] * (cvstage[0][eq][i] - cvstage[1][eq][i]);

			cv[eq][i] += (1. - q_) * cvstage[1][eq][i];

			for (int st_ = 0; st_ < rks + 2; ++st_) {
				if (TSRK_q[rks + 1][st_] == 0.)
					continue;

				cv[eq][i] += TSRK_q[rks + 1][st_] * (cvstage[st_][eq][i] - cfl * dt[i] / r_ * rhsstage[st_][eq][i]);
			}
		}
		UpdateP(i);
	}
}

void Solver::CVProcessing(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, int rks, DiagonalFunc D_Func, LUFunc L_Func, LUFunc U_Func)
{
	if (time_expl == true) {
		if (steadiness || time_stepping == 0) {			// If steady or (unsteady and dual-time stepping)
			New_Steady_UnsteadyDualTime_CV(cvstage, rks);
		} 
		else {
			New_UnsteadyTwoStep_CV(rhsstage, cvstage, rks);		// If unsteady and two step RK
		}
	} else {
		if (1 == 1) {
			LUSGS(D_Func, L_Func, U_Func, rks, rhsstage, cvstage);
		}
		else {
			//LUSGS(D_Func, L_Func, U_Func, rks, rhsstage, cvstage);
			GMRES(D_Func, L_Func, U_Func, rks, rhsstage, cvstage);
		}
	}

	RefreshBoundaries();						// Refresh boundary conditions
	RhoUPH();

	if (time_expl == true) {
		if (steadiness || time_stepping == 0) {
			if (rks < RK_stages_num - 1 + 0) {	// Storing new CV
				cvstage[rks + 1 + 0] = cv;
			}
		}
		else {
			cvstage[rks + 1 + 1] = cv;
		}
	}
}

void Solver::FinalTSRK_CV_Calculation(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage)
{
	double r_ = 1.;		// Don't understand this one

	RhoUPH();
	LRState("");
	ComputeRHSandJacobian();
	rhsstage[RK_stages_num + 1] = rhs;

	// update (conservative variables and pressure)
	for (int i = 1; i < ib2; ++i)
	{
		for (int eq = 0; eq < eq_num; ++eq) {
			cv[eq][i] = 0.;
			double eta_ = 0.;
			for (int st_ = 0; st_ < RK_stages_num + 2; ++st_) {
				eta_ += TSRK_eta[st_];
			}
			cv[eq][i] += TSRK_teta * (cvstage[0][eq][i] - cvstage[1][eq][i]);

			cv[eq][i] += (1. - eta_) * cvstage[1][eq][i];

			for (int st_ = 0; st_ < RK_stages_num + 2; ++st_) {
				if (TSRK_eta[st_] == 0.)
					continue;

				cv[eq][i] += TSRK_eta[st_] * (cvstage[st_][eq][i] - cfl * dt[i] / r_ * rhsstage[st_][eq][i]);
			}
		}
		UpdateP(i);
	}

	// boundary conditions

	RefreshBoundaries();
	RhoUPH();
}

void Solver::FinalImplicitRHS_CV_Calculation(vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage)
{
	cv = cvold;
	for (int i = 1; i < ib2; ++i) {
		for (int eq = 0; eq < eq_num; ++eq) {
			rhs[eq][i] = 0.;
			for (int rks = 0; rks < RK_stages_num; ++rks)
				rhs[eq][i] += rhsstage[rks][eq][i] * dt[i] * cfl / vol[i] * ak[RK_stages_num][rks];
		}
	}

	for (int i = 1; i < ib2; ++i) {
		for (int eq = 0; eq < eq_num; ++eq)
			for (int rks = 0; rks < RK_stages_num; ++rks)
				cv[eq][i] -= alpha[RK_stages_num][rks] * (cvold[eq][i] - cvstage[rks][eq][i] * a[i]);
	}

	for (int i = 1; i < ib2; ++i) {
		for (int eq = 0; eq < eq_num; ++eq)
			cv[eq][i] -= rhs[eq][i];

		UpdateP(i);
	}

	RefreshBoundaries();
	RhoUPH();
}

void Solver::GetWn(vector < vector < double > >& Wn_)
{
	for (int eq = 0; eq < eq_num; ++eq)
	{
		for (int i = 0; i < imax; ++i)
		{
			Wn_[eq][i] = cvold[eq][i] / a[i];
		}
	}
}

MatrixXd Solver::Diagonal(std::vector < MatrixXd >& LeftFluxJac, std::vector <MatrixXd >& RightFluxJac, std::vector < MatrixXd >& SourceJac, std::vector < std::vector < double > > & cv_, int i, LUFunc L_Func, LUFunc U_Func)
{
	MatrixXd Di;
	MatrixXd Li, Ui;

	int Lind = i + 1;
	int Uind = i - 1;

	Li = Lower(LeftFluxJac, cv_, Lind);		// :_Func ??
	Ui = Upper(RightFluxJac, cv_, Uind);	// U_Func ??
	Di = SourceJac[i];

	Di = (Li + Ui - Di * vol[i]);

	return Di;
}

MatrixXd Solver::Lower(std::vector < MatrixXd >& LeftFluxJac, std::vector < std::vector < double > >& cv_, int i)
{
	MatrixXd Li;

	Li.setZero(eq_num, eq_num);
	if (i > 1) {
		Li = 0.5 * ((LeftFluxJac[i - 1] - LeftFluxJac[i - 2]) * 0.5 * (a[i] + a[i - 1]) + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cv_, i - 1) * omega);
	}

	return Li;
}

MatrixXd Solver::Lower2(std::vector < MatrixXd >& LeftFluxJac, std::vector < std::vector < double > >& cv_, int i)
{
	MatrixXd Li;

	Li.setZero(eq_num, eq_num);
	if (i > 1) {
		Li = LeftFluxJac[i - 1] * (a[i - 1] + a[i]) / 2.;
		Li -= LeftFluxJac[i - 2] * (a[i - 2] + a[i - 1]) / 2.;
	}

	return Li;
}

MatrixXd Solver::Upper(std::vector < MatrixXd >& RightFluxJac, std::vector < std::vector < double > >& cv_, int i)
{
	MatrixXd Ui;

	Ui.setZero(eq_num, eq_num);
	if (i < ib2 - 1) {
		Ui = 0.5 * (-(RightFluxJac[i + 1] - RightFluxJac[i]) * 0.5 * (a[i] + a[i + 1]) - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cv_, i + 1) * omega);
	}

	return Ui;
}

MatrixXd Solver::Upper2(std::vector < MatrixXd >& RightFluxJac, std::vector < std::vector < double > >& cv_, int i)
{
	MatrixXd Ui;

	Ui.setZero(eq_num, eq_num);
	if (i < ib2 - 1) {
		Ui = -RightFluxJac[i] * (a[i + 1] + a[i]) / 2.;
		Ui += RightFluxJac[i + 1] * (a[i + 2] + a[i + 1]) / 2.;
	}

	return Ui;
}

SparseMatrix< double > Solver::ILU_0(SparseMatrix< double >& SM)
{
	SparseMatrix< double > LU;
	int k_start, k_finish, j_start, j_finish;

	//LU.reserve(VectorXi::Constant((ib2 - 1) * eq_num, eq_num * eq_num));
	LU = SM;

	// There is more to optimize
	for (int i = 1; i < (ib2 - 1) * eq_num; ++i) {
		k_start = max(0, i / eq_num - 1) * eq_num;
		k_finish = i;
		for (int k = k_start; k < k_finish; ++k) {
			LU.coeffRef(i, k) = LU.coeff(i, k) / LU.coeff(k, k);
			j_start = k + 1;
			j_finish = min(ib2 - 1, i / eq_num - 1 + 3) * eq_num;
			for (int j = j_start; j < j_finish; ++j) {
				LU.coeffRef(i, j) = LU.coeff(i, j) - LU.coeff(i, k) * LU.coeff(k, j);		// <- Here for instance. We don't need compute if LU(k, j) is zero (No, it's OL)
			}
		}
	}
	LU.isCompressed();

	return LU;
}

void Solver::GMRES(DiagonalFunc D_Func, LUFunc L_Func, LUFunc U_Func, int rks, vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, double physDt_)
{
	if (time_expl == false /*&& rks == 0 && false*/) {
		MatrixXd M(eq_num, eq_num);
		VectorXd R1(eq_num);
		vector < vector < double > > Wn(eq_num, vector < double >(imax, 0.));

		//Eigen::SparseMatrix<double> S((ib2 - 1) * eq_num, (ib2 - 1) * eq_num);
		MatrixReplacement A;
	
		//S.reserve(VectorXi::Constant((ib2 - 1) * eq_num, eq_num * eq_num));
		S.setZero();
		for (int i = 1; i < ib2; ++i) {
			if (i > 1) {
				M = ak[rks][rks] * U_Func(U_SGS, cv, i);
				for (int j = 0; j < eq_num; ++j) {
					for (int k = 0; k < eq_num; ++k) {
						S.insert((i - 1 - 1) * eq_num + j, (i - 1) * eq_num + k) = M(j, k);
					}
				}
			}

			M = ak[rks][rks] * D_Func(L_SGS, U_SGS, D_SGS, cv, i, L_Func, U_Func);
			for (int j = 0; j < eq_num; ++j) {
				M(j, j) += vol[i] / (dt[i] * ck[rks] * cfl);
			}
			for (int j = 0; j < eq_num; ++j) {
				for (int k = 0; k < eq_num; ++k) {
					S.insert((i - 1) * eq_num + j, (i - 1) * eq_num + k) = M(j, k);
				}
			}

			if (i < ib2 - 1) {
				M = ak[rks][rks] * L_Func(L_SGS, cv, i);
				for (int j = 0; j < eq_num; ++j) {
					for (int k = 0; k < eq_num; ++k) {
						S.insert((i - 1 + 1) * eq_num + j, (i - 1) * eq_num + k) = M(j, k);
					}
				}
			}
		}
		S.makeCompressed();

		A.attachMyMatrix(S);
		Eigen::VectorXd b((ib2 - 1) * eq_num), x;
		Eigen::VectorXd approx((ib2 - 1) * eq_num);

		GetWn(Wn);
		for (int i = 1; i < ib2; ++i) {
			R1 = ak[rks][rks] * GetEigenVector(rhs, i);
			for (int l = 0; l < rks; ++l)
			{
				R1 += ak[rks][l] * GetEigenVector(rhsstage[l], i);
				R1 += alpha[rks][l] * (GetEigenVector(Wn, i) - GetEigenVector(cvstage[l], i)) * a[i];
			}
			for (int j = 0; j < eq_num; ++j) {
				b((i - 1) * eq_num + j) = -R1[j];
				//approx((i - 1) * eq_num + j) = -0*(cvstage[rks][j][i] * a[i] - cvold[j][i]);
			}
		}

		//Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
		//gmres.compute(A);
		//gmres.setTolerance(1e-10);

		//Eigen::GMRES<SparseMatrix<double>, Eigen::DiagonalPreconditioner < double > > gmres;

		//Eigen::GMRES<SparseMatrix<double>, Eigen::MyPreconditioner2 < double > > gmres;
		//Eigen::GMRES<SparseMatrix<double>, Eigen::MyPreconditioner3 < double > > gmres;		// It works though - Like Identity I suppose
		//gmres.preconditioner().setInvDiag(approx);

		Eigen::SparseMatrix<double> LU = ILU_0(S);
		//Eigen::GMRES<SparseMatrix<double>, Eigen::ILU_0_Preconditioner < double > > gmres;
		Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::ILU_0_Preconditioner < double > > gmres;
		gmres.preconditioner().setLUmatrix(LU, eq_num);

		gmres.compute(S);
		x = gmres.solve(b);
		//solver.compute(S);
		//x = solver.solve(b);

		for (int i = 1; i < ib2; ++i) {
			for (int eq = 0; eq < eq_num; ++eq) {
				cv[eq][i] = cvold[eq][i] + x((i - 1) * eq_num + eq);
				cvstage[rks][eq][i] = cv[eq][i] / a[i];
			}
		}
	}
}

void Solver::LUSGS(DiagonalFunc D_Func, LUFunc L_Func, LUFunc U_Func, int rks, vector < vector < vector < double > > >& rhsstage, vector < vector < vector < double > > >& cvstage, double physDt_)
{
	double dsi;
	double adtv;

	vector < MatrixXd > invD(imax, MatrixXd(eq_num, eq_num));
	vector < VectorXd > DdWijk(imax, VectorXd(eq_num));
	VectorXd DdWijkn(eq_num);
	vector < VectorXd > Wijk(imax, VectorXd(eq_num));
	VectorXd LdW1(eq_num);
	VectorXd R1(eq_num);
	VectorXd UdWn(eq_num);
	VectorXd dW(eq_num);
	VectorXd Wi(eq_num);

	MatrixXd Li(eq_num, eq_num);
	MatrixXd Ui(eq_num, eq_num);
	MatrixXd Di(eq_num, eq_num);

	vector < vector < double > >	Wn(eq_num, vector < double >(imax, 0.));

	auto scalar_prod = [](const vector < double >::const_iterator& vec1, const vector < double >::const_iterator& vec2, int size) -> double
	{
		vector <double> prods(size);
		std::transform(vec1, vec1 + size, vec2, prods.begin(), std::multiplies<double>());
		return std::accumulate(prods.begin(), prods.end(), 0.);
	};

	GetWn(Wn);

	Wijk[0](0) = Wn[0][0];
	Wijk[0](1) = Wn[1][0];
	Wijk[0](2) = Wn[2][0];

	// First sweep
	for (int i = 1; i < ib2; ++i)
	{
		Di = ak[rks][rks] * D_Func(L_SGS, U_SGS, D_SGS, cvold, i, L_Func, U_Func);
		for (int j = 0; j < eq_num; ++j) {
			/*for (int k = 0; k < eq_num; ++k) {
				Di(j, k) += vol[i] / (dt[i] * ck[rks] * cfl);
			}*/
			Di(j, j) += vol[i] / (dt[i] * ck[rks] * cfl);
			if (!steadiness) {
				Di(j, j) += 3. * vol[i] / (2. * physDt_);
			}
		}
		Li = ak[rks][rks] * L_Func(L_SGS, cvold, i);
		
		dW = (Wijk[i - 1] - GetEigenVector(Wn, i - 1)) * a[i - 1];
		LdW1 = Li * dW;
		//if (steadiness)
			R1 = ak[rks][rks] * GetEigenVector(rhs, i);
		/*else
		{
			for (int eq = 0; eq < eq_num; ++eq)
				R1(eq) = ak[rks][rks] * scalar_prod(rhs[eq].begin() + 1, M_matrix.inversed_mass_matrix[i - 1].begin(), ib2 - 1);
		}*/

		for (int l = 0; l < rks; ++l)
		{
			R1 += ak[rks][l] * GetEigenVector(rhsstage[l], i);
			R1 += alpha[rks][l] * (GetEigenVector(Wn, i) - GetEigenVector(cvstage[l], i)) * a[i];
		}

		DdWijk[i] = - R1 - LdW1;

		invD[i] = Di.inverse();
		/*for (int j = 0; j < eq_num; ++j) {
			for (int k = 0; k < eq_num; ++k) {
				invD[i](j, k) = 1. / Di(j, k);
			}
		}*/
		Wijk[i] = (invD[i] * DdWijk[i] + GetEigenVector(cvold, i)) / a[i];
	}

	// Second sweep
	cvstage[rks] = Wn;
	for (int i = ib2 - 1; i >= 1; --i)
	{
		Ui = ak[rks][rks] * U_Func(U_SGS, cvold, i);

		dW = (GetEigenVector(cvstage[rks], i + 1) - GetEigenVector(Wn, i + 1))* a[i + 1];
		UdWn = Ui * dW;

		DdWijkn = DdWijk[i] - UdWn;

		DdWijk[i] = invD[i] * DdWijkn;
		Wi = (GetEigenVector(cvold, i) + DdWijk[i]) / a[i];
		cvstage[rks][0][i] = Wi(0);
		cvstage[rks][1][i] = Wi(1);
		cvstage[rks][2][i] = Wi(2);
	}

	for (int i = 0; i < imax; ++i) {
		for (int eq = 0; eq < eq_num; ++eq) {
			cv[eq][i] = cvstage[rks][eq][i] * a[i];
		}
	}

	//for (int i = 1; i < ib2; ++i)
	//{
	//	VectorXd rhsadd;

	//	Li = L_Func2(L_SGS, cvold, i);		// cvold will not used
	//	dW = (GetEigenVector(cvstage[rks], i - 1) - GetEigenVector(Wn, i - 1))* a[i - 1];
	//	rhsadd = Li * dW;

	//	Ui = U_Func2(U_SGS, cvold, i);		// cvold will not used
	//	dW = (GetEigenVector(cvstage[rks], i + 1) - GetEigenVector(Wn, i + 1)) * a[i + 1];
	//	rhsadd += Ui * dW;

	//	Di = D_Func(L_SGS, U_SGS, D_SGS, cvold, i, L_Func2, U_Func2);
	//	dW = (GetEigenVector(cvstage[rks], i) - GetEigenVector(Wn, i)) * a[i];
	//	rhsadd += Di * dW;

	//	rhsstage[rks][0][i] = rhs[0][i] + rhsadd(0) * 1.;
	//	rhsstage[rks][1][i] = rhs[1][i] + rhsadd(1) * 1.;
	//	rhsstage[rks][2][i] = rhs[2][i] + rhsadd(2) * 1.;
	//}
}

MatrixXd Solver::ToEigen(vector < vector < double > >& M)
{
	MatrixXd EM(M.size(), M[0].size());
	
	for (int i = 0; i < M.size(); ++i)
		for (int j = 0; j < M[0].size(); ++j)
			EM(i,j) = M[i][j];

	return EM;
}

VectorXd Solver::ToEigen(vector < double >& V)
{
	VectorXd EV(V.size());

	for (int i = 0; i < V.size(); ++i)
		EV(i) = V[i];

	return EV;
}

vector < double > Solver::ToVector(VectorXd& EV)
{
	vector < double > V (EV.size());

	for (int i = 0; i < EV.size(); ++i) {
		V[i] = EV(i);
	}

	return V;
}

vector < vector < double > > Solver::ToMatrix(const MatrixXd& EM)
{
	vector < vector < double > > M(EM.rows(), vector < double > (EM.cols(), 0));

	for (int i = 0; i < EM.rows(); ++i) {
		for (int j = 0; j < EM.cols(); ++j) {
			M[i][j] = EM(i, j);
		}
	}

	return M;
}

VectorXd Solver::GetEigenVector(vector < vector < double > >& Vec, int i)
{
	VectorXd EVec(Vec.size());

	for (int j = 0; j < Vec.size(); ++j)
		EVec(j) = Vec[j][i];

	return EVec;
}

vector < double > Solver::MatVecMult(vector < vector < double > >& Mat, vector < double >& Vec)
{
	vector < double > ResVec(Mat.size(), 0.);

	for (int i = 0; i < Mat.size(); ++i) {
		for (int j = 0; j < Vec.size(); ++j) {
			ResVec[i] += Mat[i][j] * Vec[j];
		}
	}

	return ResVec;
}

vector < double > Solver::VecSum(double a, vector < double >& Vec1, double b, vector < double >& Vec2)
{
	vector < double > ResVec(Vec1.size(), 0.);

	for (int i = 0; i < Vec1.size(); ++i) {
		ResVec[i] = a * Vec1[i] + b * Vec2[i];
	}

	return ResVec;
}

void Solver::AreaMult(vector < vector < double > >& M, int i, int rks)
{
	double si = 0.5 * (a[i] + a[i + 1]);

	for (int j = 0; j < eq_num; j++)
		for (int k = 0; k < eq_num; k++)
			M[j][k] *= si * ak[rks][rks];
}

MatrixXd Solver::MakeD(MatrixXd& Li, MatrixXd& Ui, MatrixXd& Di, double dt, int i, int rks)
{
	MatrixXd D(eq_num, eq_num);

	D = Li + Ui - Di * vol[i] * ak[rks][rks];
	for (int j = 0; j < eq_num; ++j) {
		//for (int k = 0; k < eq_num; ++k) {
			//D(j, k) = 0.;
			//if (j == k) {
				D(j, j) += vol[i] / dt;
			//}
			//D(j, k) += Li(j, k) + Ui(j, k) - Di(j, k) * vol[i] * ak[rks][rks];
		//}
	}
	
	return D;
}

void Solver::RhoUPH()
{
	double gamma_ = GetGamma();

	if (fv.size() == 0) {
		fv.resize(var_num);
		for (int i = 0; i < var_num; ++i) {
			fv[i].resize(imax, 0.);
		}
	}

	vector<double*> fv_ref(var_num);
	vector<filling_type> f_types = { filling_type::only_vars, filling_type::only_diff };
	for (auto f_type : f_types)
	{
#ifdef USE_OMP
#pragma omp parallel for num_threads(MAX_THREAD_NUM) /*private(fv_ref)*/ firstprivate(f_type, fv_ref)
#endif
		for (int i = 0; i < imax; ++i)
		{
			for (int v = 0; v < var_num; ++v)
				fv_ref[v] = &fv[v][i];

			fill_fv_equations<double>(f_type, fv_ref, i);
		}
	}
	if (!delayed_fv_equations.empty())
		for (auto eq_name : delayed_fv_equations)
		{
			int var_id = vars_o.at(eq_name);
#ifdef USE_OMP
#pragma omp parallel for num_threads(MAX_THREAD_NUM) firstprivate(var_id, eq_name)
#endif
			for (int i = 0; i < imax; ++i)
				fv[var_id][i] = make_fv_equation<double>(eq_name, i);
		}
}

void Solver::TimeSteps(bool local_time, double dt_)		// Blazek, Section 6.1.4, p. 173, "Determination of the maximum time step"
{
	double rho;
	double u;
	double p;
	double cs;
	double dx;
	double sprad;
	double dt_min = 1e20;

	if (dt_ > 0.) {
		for (int i = 0; i < imax; ++i)
		{
			dt[i] = dt_;
		}
		return;
	}
	for (int i = 1; i < ib2; ++i)
	{
		rho = make_fv_equation<double>(var_name[g_RHO], i);
		u = make_fv_equation<double>(var_name[g_U], i);
		p = make_fv_equation<double>(var_name[g_P], i);
		cs = sqrt(gamma * p / rho);
		dx = 0.5 * (x[i + 1] - x[i - 1]);
		sprad = cs * sqrt(dx * dx + pow(a[i], 2)) + abs(u) * a[i];
		dt[i] = vol[i] / sprad;
		if (dt_min > dt[i])
			dt_min = dt[i];

	}
	dt[0] = dt[1];
	dt[imax - 1] = dt[ib2 - 1];
	for (int i = 0; i < imax && !local_time; ++i)
	{
		dt[i] = dt_min;
	}
}

double Solver::SpectralRadius(vector< vector < double > >& cv_, int i)
{
	double rho;
	double u;
	double cs;
	double rhoe;

	vector<double> cons_var(eq_num);
	for (int eq = 0; eq < eq_num; ++eq)
		cons_var[eq] = cv_[eq][i];

	double a = make_fv_equation<double>(var_name[g_A], i);
	rho = make_fv_equation<double>(var_name[g_RHO], i, vector<double*>(), cons_var.data()) / a;		// Here we need to devide by a, as ususaly cons vars are without area, but here it is with! And in this equation the area is = 1, but it shouldn't here. Sorry.
	u = make_fv_equation<double>(var_name[g_U], i, vector<double*>(), cons_var.data());
	rhoe = make_fv_equation<double>(var_name[g_E], i, vector<double*>(), cons_var.data()) * rho;
	cs = sqrt(gamma * (gamma - 1.) * (rhoe - 0.5 * rho * u * u) / rho);

	return (cs  + abs(u)) * a;
}

void Solver::SourceTerm(int i)
{
	double da = 0.5 * (a[i + 1] - a[i - 1]);
	rhs[RHO_U_A][i] = rhs[RHO_U_A][i] - p[i] * da;
}

void Solver::ImplResidualSmooth()
{
	double eps2, t;
	int i;
	double* d = dummy.data();

	eps2 = 2. * eps_impl_res_smooth + 1.;
	d[0] = 0.;
	for (int eq = 0; eq < eq_num; ++eq) 
		rhs[eq][0] = 0.;

	// elimination step

	for (int i = 1; i < ib2; ++i)
	{
		t = 1. / (eps2 - eps_impl_res_smooth * d[i - 1]);
		d[i] = t * eps_impl_res_smooth;
		for (int eq = 0; eq < eq_num; ++eq)
			rhs[eq][i] = t * (rhs[eq][i] + eps_impl_res_smooth * rhs[eq][i - 1]);
			
	}

	// backward substitution

	i = ib2 - 1;
	for (int ii = 2; ii < ib2; ++ii)
	{
		i = i - 1;

		for (int eq = 0; eq < eq_num; ++eq)
			rhs[eq][i] = rhs[eq][i] + d[i] * rhs[eq][i + 1];

	}
}

void Solver::ImplResidualSmooth(vector < vector < vector < double > > >& rhs_)
{
	double eps2, t;
	int i;
	double* d = dummy.data();

	eps2 = 2. * eps_impl_res_smooth + 1.;
	d[0] = 0.;
	for (int rks = 0; rks < RK_stages_num; ++rks)
		for (int eq = 0; eq < eq_num; ++eq)
			rhs_[rks][eq][0] = 0.;

	// elimination step

	for (int i = 1; i < ib2; ++i)
	{
		t = 1. / (eps2 - eps_impl_res_smooth * d[i - 1]);
		d[i] = t * eps_impl_res_smooth;

		for (int rks = 0; rks < RK_stages_num; ++rks)
			for (int eq = 0; eq < eq_num; ++eq)
				rhs_[rks][eq][i] = t * (rhs_[rks][eq][i] + eps_impl_res_smooth * rhs_[rks][eq][i - 1]);

	}

	// backward substitution

	i = ib2 - 1;
	for (int ii = 2; ii < ib2; ++ii)
	{
		i = i - 1;

		for (int rks = 0; rks < RK_stages_num; ++rks)
			for (int eq = 0; eq < eq_num; ++eq)
				rhs_[rks][eq][i] = rhs_[rks][eq][i] + d[i] * rhs_[rks][eq][i + 1];

	}
}

double Solver::Convergence()
{
	int idrho, nsup;
	double dr, drmax, rho, u, c, avms, drho, dconc, dn;

	// get the change of density and the mass flow

	drho = 0.;
	dconc = 0.;
	drmax = -1e20;
	avms = 0.;
	nsup = 0;
	double M_max = 0.;

	for (int i = 1; i < ib2; ++i)
	{
		dr = cv[g_RHO_A][i] - cvold[g_RHO_A][i];
		//dn = cv[N_A][i] - cvold[N_A][i];

		avms = avms + cv[g_RHO_U_A][i];
		drho = drho + dr * dr;
		//dconc = dconc + dn * dn;
		if (abs(dr) >= drmax) {
			drmax = abs(dr);
			idrho = i;
		}
		rho = fv[g_RHO][i];		// cv[RHO_A][i] / a[i];
		u = fv[g_U][i];			// cv[RHO_U_A][i] / cv[RHO_A][i];
		c = sqrt(gamma* p[i] / rho);
		if (u > c) nsup++;
		if (u / c > M_max) M_max = u / c;
	}
	avms = avms / double(ncells + 1);

	// print convergence history

	if (iter == 1) {
		drho1 = sqrt(drho / double(ncells + 1));
		//dconc1 = sqrt(dconc / double(ncells + 1));
	}

	drho = sqrt(drho / double(ncells + 1)) / drho1;
	//dconc = sqrt(dconc / double(ncells + 1)) / dconc1;

	double factor_;

	factor_ = prev_drho / drho;
	factor_ = min(2., factor_);
	//factor_ = min(1.1, factor_);
	//factor_ = max(0.9, factor_);
	factor_ = max(0.1, factor_);
	//cfl = cfl * factor_;

	prev_drho = drho;
	if (isnan(drho)) {
		cout << "nan" << endl;
		return -100.;
	}

	if (!NoOutput) {
		cout << "convergence" << "\t";
		cout << iter << "\t";
		cout << log10(drho) << "\t";
		//cout << log10(dconc) << "\t";
		cout << drmax << "\t";
		cout << idrho << "\t";
		cout << avms << "\t";
		cout << cv[1][1] - cv[1][ib2 - 1] << "\t";
		cout << nsup << "\t";
		cout << M_max << "\n";
	}

	return /*max(*/drho/*, dconc * 1.)*/;
}

void Solver::PrintResult()
{
	double rho, u, temp, c, mach;
	bool AG = false;

	FILE* file;

	double rgas = (gamma - 1.) * Cp / gamma;

	if (eq_num == 4)
		AG = true;

	file = fopen(output_file.c_str(), "w");
	fprintf(file, "x\tA\trho\tu\tp\tT\tM\tmass_flow");
	//if (AG) {
	//	fprintf(file, "\tn");
	//	//fv[n] = grid.GetResolution();
	//}
	fprintf(file, "\n");
	for (int i = 1; i < ib2; ++i)
	{
		rho = fv[g_RHO][i]; // cv[0][i] / a[i];
		u = fv[g_U][i]; // cv[1][i] / cv[0][i];
		temp = make_fv_equation<double>("T", i); // p[i] / (rgas * rho);
		c = sqrt(gamma * p[i] / rho);
		mach = /*abs*/(u) / c;
		fprintf(file, "%14.10le\t%le\t%le\t%le\t%le\t%le\t%le\t%le", 
			x[i], a[i], rho, u, p[i], temp, make_fv_equation<double>("dT", i), cv[g_RHO_U_A][i]);

		//if (AG) {
		//	fprintf(file, "\t%lf", fv[n][i]/* / grid.GetResolution()[i]*/);
		//}
		fprintf(file, "\n");
	}
}

void Solver::SetOutputFile(string output_file_)
{
	output_file = output_file_;
}

void Solver::ShowOutput()
{
	NoOutput = false;
}

void Solver::HideOutput()
{
	NoOutput = true;
}

void Solver::InitFirsStages(vector < vector < vector < double > > >& cvstage, vector < vector < vector < double > > >& rhsstage)
{
	cv = cvnm1;
	RhoUPH();
	LRState("");
	ComputeRHSandJacobian();
	cvstage[0] = cvnm1;
	rhsstage[0] = rhs;

	cv = cvn;
	cvstage[1] = cvn;
}

double Solver::CalcH(double p_, double Rho_, double U_)
{
	return gamma / (gamma - 1.) * p_ / Rho_ + 0.5 * U_ * U_;
}

Solver::equation::equation(string eq_name_, vector<string> dt_term_, vector<string> dx_term_, vector<string> source_term, map < string, int > vars_, map < string, int > vars_o_)
{
	eq_name = eq_name_;
	cur_dt = get_equation(eq_name_, dt_term_);
	cur_dx = get_equation(eq_name_, dx_term_);
	cur_source = get_equation(eq_name, source_term);
}

template <typename T>
T Solver::custom_func::get_function_value(const std::vector<std::vector<double>>& params, const vector<double_type>& params_a, int i, double time)
{
	double_type result(0.);
	func(params, params_a, var_names, param_names, i, time, result);
	return result;
}

#ifdef USE_ADEPT
template<>
double Solver::custom_func::get_function_value<double>(const std::vector<std::vector<double>>& params, const vector<double_type>& params_a, int i, double time)
{
	double_type result(0.);
	func(params, params_a, var_names, param_names, i, time, result);
	return result.value();
}
#endif

eq_term::eq_term(const string& term_s)
{
	auto get_var_name = [](const string& term_s) -> string
	{
		string var_name = "";
		bool var_is_started = false;
		for (unsigned int i = 0; i < term_s.size(); ++i)
		{
			if (isalpha(term_s[i]))
			{
				if (!(term_s[i] == 'e' && i < term_s.size() - 1 && !isalpha(term_s[i + 1])) || var_is_started)
				{
					var_is_started = true;
					var_name += term_s[i];
				}
			}
			else
			{
				if (var_is_started)
					break;
			}
		}
		return var_name;
	};
	unsigned op_pos = 0;
	name = get_var_name(term_s);
	char sym = term_s[0];
	switch (sym)
	{
	case '-':
		op = operation::minus;
		op_pos = 1;
		break;
	case '*':
		op = operation::mult;
		op_pos = 1;
		break;
	case '/':
		op = operation::div;
		op_pos = 1;
		break;
	case '+':
		op = operation::plus;
		op_pos = 1;
		break;
	default:
		op = operation::plus;
		break;
	}
	if (!sscanf(term_s.c_str() + op_pos, "%le", &coef))
		coef = 1.;
	auto deg_pos = term_s.find('^');
	if (deg_pos != string::npos)
	{
		if (!sscanf(term_s.c_str() + deg_pos + 1, "%le", &degree))
			degree = 1.;
	}
	else
		degree = 1.;
}

//void Solver::FillJacobian(vector < vector < double > >& M_SGS, vector < double >& jac, double s)
void Solver::FillJacobian(MatrixXd& M_SGS, vector < double >& jac, double s)
{
	for (int i = 0; i < eq_num; ++i)		// i - is a Conservative Variable number
	{
		for (int j = 0; j < eq_num; ++j)	// j - is a Flux number
		{
			//M_SGS[j][i] = jac[i * eq_num + j];// *s;
			M_SGS(j, i) = jac[i * eq_num + j];// *s;
		}
	}
	/*for (int eq = 0; eq < eq_num; ++eq) {
		M_SGS[eq].assign(jac.begin() + eq * eq_num, jac.begin() + (eq + 1) * eq_num);
	}*/
}

MatrixXd Solver::GetJacobian(vector < double >& jac)
{
	MatrixXd M_SGS(eq_num, eq_num);

	for (int i = 0; i < eq_num; ++i)
	{
		for (int j = 0; j < eq_num; ++j)
		{
			//M_SGS[j][i] = jac[i * eq_num + j];// *s;
			M_SGS(j, i) = jac[i * eq_num + j];// *s;
		}
	}

	return M_SGS;
}

Solver* CreateReadConfigFile(string file_name)
{
	string solver_s;

	cout << "Opening file " << "\"" << file_name.c_str() << "\"..." << endl;
	Node config = LoadFile(file_name.c_str());

	if (config["solver"]) {
		solver_s = config["solver"].as< string >();
		if ( solver_s == "cds") {
			// cout << "Solver is '" << solver_s << "'" << endl;
		} else if (solver_s == "cusp") {

		}
		else if (solver_s == "vanleer") {

		}
		else if (solver_s == "hlle" || solver_s == "hllc") {

		}
		else {
			cout << "There is no such a solver '" << solver_s << "'" << endl;
			cout << "Try these ones:" << endl;
			for (auto sol = solvers.begin(); sol != solvers.end(); ++sol) {
				cout << "\t" << *sol << endl;
			}

			return NULL;
		}
	} else {
		cout << "There is no a word 'solver'" << endl;

		return NULL;
	}

	sol_struct sol_init;
	cds_struct cds_init;

	sol_init.solver_name = solver_s;
	if (config["time_treat"].as< string >() == "implicit") {
		string IRKC_name = "IRKC_" + config["IRK"].as< string >();

		//sol_init.ak.push_back(vector < double >());
		sol_init.ak = config[IRKC_name]["ak"].as< vector < vector < double > > >();
		sol_init.bk = config[IRKC_name]["bk"].as< vector < double > >();
		sol_init.ck = config[IRKC_name]["ck"].as< vector < double > >();
		sol_init.alpha = config[IRKC_name]["alpha"].as< vector < vector < double > > >();
	}

	sol_init.gas = config["gas"].as< string >();
	sol_init.Cp = config["Cp"].as< double >();
	sol_init.cfl = config["cfl"].as< double >();
	sol_init.max_iter_num = config["max_iter_num"].as< int >();
	sol_init.tolerance = config["tolerance"].as< double >();
	sol_init.gamma = config["gamma"].as< double >();
	sol_init.free_g = config["g"].as< double >();

	// Explicit/Implicit time Scheme
	if (config["time_treat"]) {
		if (config["time_treat"].as< string >() == "implicit") {
			sol_init.time_expl = false;
		}
		else {
			sol_init.time_expl = true;
		}
	}
	else {
		sol_init.time_expl = true;
	}

	sol_init.RK_stage_coeffs = config["RK_stage_coeffs"].as< vector < double > >();

	if (config["steadiness"].as< string >() == "steady") {
		sol_init.steadiness = true;
	}
	else {
		sol_init.steadiness = false;
		if (config["time_stepping"].as< string >() == "dual_time_stepping") {
			sol_init.time_stepping = 0;
		}
		if (config["time_stepping"].as< string >() == "two_step_RK") {
			sol_init.time_stepping = 1;
		}
	}
	if (config["time_treat"].as< string >() == "explicit") {
		if (sol_init.steadiness == true || sol_init.time_stepping == 0) {
			string ERKC_name = "ERKC_" + config["ERK"].as< string >();

			sol_init.RK_alpha = config[ERKC_name]["alpha"].as< vector < vector < double > > >();
			sol_init.RK_beta = config[ERKC_name]["beta"].as< vector < vector < double > > >();
		}
		else if (sol_init.time_stepping == 1) {
			string TSRK_name = "TSRK_" + config["TSRK"].as< string >();

				sol_init.TSRK_teta = config[TSRK_name]["TSRK_teta"].as< double >();
				sol_init.TSRK_d = config[TSRK_name]["TSRK_d"].as< vector < double > >();
				sol_init.TSRK_eta = config[TSRK_name]["TSRK_eta"].as< vector < double > >();
				sol_init.TSRK_q = config[TSRK_name]["TSRK_q"].as< vector < vector < double > > >();
		}
	}

	sol_init.res_smooth_flag = config["res_smooth_flag"].as< vector < int > >();
	sol_init.eps_impl_res_smooth = config["eps_impl_res_smooth"].as< double >();

	sol_init.lts = config["lts"].as< bool >();

	sol_init.remesh = config["remesh"].as< bool >();

	if (solver_s == "cds") {
		cds_init.diss_blend = config["cds"]["diss_blend"].as< vector < double > >();
		cds_init.diss_flag = config["cds"]["diss_flag"].as< vector < int > >();
		cds_init.vis2 = config["cds"]["vis2"].as< double >();
		cds_init.vis4 = config["cds"]["vis4"].as< double >();

		return new CDS(sol_init, cds_init);
	}
	if (solver_s == "cusp") {
		return new CUSP(sol_init);
		
	}
	if (solver_s == "vanleer") {
		return new VL(sol_init);

	}
	if (solver_s == "hlle" || solver_s == "hllc") {
		return new HLLE(sol_init);

	}
}

double Sign(double value)
{
	if (value > 0.) return 1.;
	if (value < 0.) return -1.;
	return 0.;
}

int count_(vector < int >& vec, int val)
{
	return count(vec.begin(), vec.end(), val);
}

vector<double_type> Solver::construct_side_flux_array(const vector<double_type*>& vars, const int i)
{
	vector<double_type> flux(eq_num);
	unsigned int eq = 0;
	for (const auto &equation : equations)
	{
		double_type term = 1.;
		flux[eq] += make_equation<double_type>(eq, Solver::equation::term_name::dx, vars) * term;
		++eq;
	}
	return flux;
}