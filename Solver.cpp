// Solver

#include "Solver.h"
#include "CDS.h"
#include "CUSP.h"

using namespace YAML;

set < string > solvers { "cds", "cusp" };

Solver::Solver(sol_struct& sol_init_) :
					Ideal((sol_init_.gas == "Ideal" ? true : false)),
					Cp(sol_init_.Cp),
					cfl(sol_init_.cfl),
					max_iter_num(sol_init_.max_iter_num),
					tolerance(sol_init_.tolerance),
					gamma(sol_init_.gamma),
					RK_stage_coeffs(sol_init_.RK_stage_coeffs),
					res_smooth_flag(sol_init_.res_smooth_flag),
					eps_impl_res_smooth(sol_init_.eps_impl_res_smooth)
{
	eq_num = vars.size();
	var_num = vars_o.size();
	iter = 0;

	NoOutput = false;

	if (Ideal) {
		cout << "\tgas: " << "Ideal" << endl;

	} else {
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

	RK_stages_num = RK_stage_coeffs.size();

	SetEquation("mass", "RhoA", "RhoU", vars, vars_o);
	SetEquation("impulse", "RhoUA", "RhoUU+p", vars, vars_o);
	SetEquation("energy", "RhoEA", "RhoHU", vars, vars_o);
}

void Solver::SetEquation(string eq_name, string dt_term_, string dx_term_, map < string, int > vars_, map < string, int > vars_o_)
{
	equations.push_back( equation (eq_name, dt_term_, dx_term_, vars_, vars_o_) );
}

void Solver::ReadBoundaries(string file_name)
{
	cout << "Opening file " << "\"" << file_name.c_str() << "\"..." << endl;
	Node config = LoadFile(file_name.c_str());

	if (config["first"]) {
		if (config["first"]["inflow"]) {
			inflow_id = 0;
			p_b_in = config["first"]["inflow"]["p01"].as< double >();
			T_b_in = config["first"]["inflow"]["t01"].as< double >();

			cout << "Inflow total pressure: " << p_b_in << endl;
			cout << "Inflow temperature: " << T_b_in << endl;

		} else if (config["first"]["outflow"]) {
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

		} else if (config["last"]["outflow"]) {
			outflow_id = imax;
			p_b_out = config["last"]["outflow"]["p2"].as< double >();
			cout << "Outflow pressure: " << p_b_out << endl;

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
	for (int i = 0; i < eq_num; ++i) {
		cv[i].resize(imax, 0.);
		cvold[i].resize(imax, 0.);
		diss[i].resize(imax, 0.);
		rhs[i].resize(imax, 0.);
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
	double gam1 = gamma - 1.;
	double gap1 = gamma + 1.;

	// Inflow boundaries
	rgas = gam1 * Cp / gamma;
	u = fv[U][in_id_p1]; // cv[1][in_id_p1] / cv[0][in_id_p1];		///< $u = \frac{(\rho u S)_{[1]}} {(\rho  S)_{[1]}}$
	cs2 = gamma * fv[P][in_id_p1] / fv[RHO][in_id_p1]; // p[in_id_p1] * a[in_id_p1] / cv[0][in_id_p1];	// Speed of sound $c_s^2 = \frac{\gamma_ p_{[1]} S_{[1]}} { \rho  S_{[1]}}$
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

	cv[RHO_A][in_id] = rhob * a[in_id_p1];		//	$\rho S = \rho_b S_{[1]}$
	cv[RHO_U_A][in_id] = rhob * a[in_id_p1] * ub;					//	$\rho u S = \rho_b S_{[1]} u_b$
	cv[RHO_E_A][in_id] = (pb / gam1 + 0.5 * rhob * ub * ub) * a[in_id_p1];	//	$\rho E S = \left(\frac{p_b}{\gamma_ - 1} + 0.5 \rho_b u_b^2\right) S_{[1]}$
	p[in_id] = pb;

	// Outflow boundaries

	int out_id = GetOutflowId();
	out_id = (out_id > 0 ? out_id - 1 : out_id);
	int out_id_p1 = (out_id > 0 ? out_id - 1 : out_id + 1);

	rho = fv[RHO][out_id_p1]; // cv[0][out_id_p1] / a[out_id_p1];
	u = fv[U][out_id_p1]; // cv[1][out_id_p1] / cv[0][out_id_p1];
	cs = sqrt(gamma * p[out_id_p1] / rho);

	if (abs(u) >= cs) {		// supersonic flow 		// fix is here
		pb = p[out_id_p1];
		rhob = rho;
		ub = u;
	}
	else {				//subsonic flow
		pb = p_b_out;
		rhob = rho + Sign(u)*(p_b_out - p[out_id_p1]) / (cs * cs);		// $\rho_b = \rho + \frac{(p_2 - p_{[ib2 - 1]})}{c_s^2}$ 		// fix is here
		ub = u - Sign(u)*(p_b_out - p[out_id_p1]) / (cs * rho);		// $u_b = u - \frac{(p_2 - p_{[ib2 - 1]})}{c_s \rho}$ 		// fix is here
	}

	cv[RHO_A][out_id] = rhob * a[out_id_p1];
	cv[RHO_U_A][out_id] = rhob * ub * a[out_id_p1];
	cv[RHO_E_A][out_id] = (pb / gam1 + 0.5 * rhob * ub * ub) * a[out_id_p1];
	p[out_id] = pb;

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

double Solver::Solve()
{
	double fac, adtv, rrho, rhou, rhoe;
	vector < int > diss_flag_;
	vector < double > diss_blend_;

	iter++;

	if (solver_name == "cds") {
		diss_flag_ = GetDissFlag();
		diss_blend_ = GetDissBlend();
	}

	// Store previous solution; set dissipation = 0
	cvold = cv;
	for (int i = 0; i < eq_num && solver_name == "cds"; ++i)
		diss[i].assign(diss[i].size(), 0.);

	// Calculate time step
	TimeSteps();

	// loop over the R. - K. stages:
	for (int rks = 0; rks < RK_stages_num; ++rks)
	{
		  if (solver_name == "cds" && diss_flag_[rks] == 1) {	// artificial dissipation
		  	Dissipation(diss_blend_[rks]);
		  } else {
		  	//LRState();
			RhoUPH();
			LRState("");
		  }

		 // convective flux
		 Fluxes();
		 //RHS(0);

		// source term
		SourceTerm();

		// residual * time step
		fac = RK_stage_coeffs[rks] * cfl;
		for (int i = 1; i < ib2; ++i)
		{
			adtv = fac * dt[i] / vol[i];
			for (int eq = 0; eq < eq_num; ++eq)
				rhs[eq][i] = adtv * rhs[eq][i];
		}

		// implicit residual smoothing
		if (res_smooth_flag[rks] > 0 && eps_impl_res_smooth > 0.)
			ImplResidualSmooth();

		// update (conservative variables and pressure)
		for (int i = 1; i < ib2; ++i)
		{
			for (int eq = 0; eq < eq_num; ++eq)
				cv[eq][i] = cvold[eq][i] - rhs[eq][i];

			rrho = a[i] / cv[RHO_A][i];
			rhou = cv[RHO_U_A][i] / a[i];
			rhoe = cv[RHO_E_A][i] / a[i];
			p[i] = (gamma - 1.) * (rhoe - 0.5 * rhou * rhou *rrho);
		}

		// boundary conditions

		RhoUPH();
		RefreshBoundaries();
	}

	return Convergence();
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
	for (int i = 0; i < imax; ++i)
	{
		fv[RHO][i] = cv[0][i] / a[i];
		fv[U][i] = cv[1][i] / cv[0][i];
		fv[P][i] = p[i];
		fv[H][i] = gamma_ / (gamma_ - 1.) * fv[P][i] / fv[RHO][i] + 0.5 * pow(fv[U][i], 2);
	}
}

void Solver::TimeSteps()
{
	double rho;
	double u;
	double cs;
	double dx;
	double sprad;

	for (int i = 1; i < ib2; ++i)
	{
		rho = cv[RHO_A][i] / a[i];
		u = cv[RHO_U_A][i] / cv[RHO_A][i];
		cs = sqrt(gamma * p[i] / rho);
		dx = 0.5 * (x[i + 1] - x[i - 1]);
		sprad = cs * sqrt(dx * dx + pow(a[i], 2)) + abs(u) * a[i];
		dt[i] = vol[i] / sprad;
	}
	dt[0] = dt[1];
	dt[imax - 1] = dt[ib2 - 1];
}

void Solver::SourceTerm()
{
	double da;

	for (int i = 1; i < ib2; ++i)
	{
		da = 0.5 * (a[i + 1] - a[i - 1]);
		rhs[RHO_U_A][i] = rhs[RHO_U_A][i] - p[i] * da;
	}
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

double Solver::Convergence()
{
	int idrho, nsup;
	double dr, drmax, rho, u, c, avms, drho;

	// get the change of density and the mass flow

	drho = 0.;
	drmax = -1e20;
	avms = 0.;
	nsup = 0;

	for (int i = 1; i < ib2; ++i)
	{
		dr = cv[RHO_A][i] - cvold[RHO_A][i];
		avms = avms + cv[RHO_U_A][i];
		drho = drho + dr * dr;
		if (abs(dr) >= drmax) {
			drmax = abs(dr);
			idrho = i;
		}
		rho = fv[RHO][i];		// cv[RHO_A][i] / a[i];
		u = fv[U][i];			// cv[RHO_U_A][i] / cv[RHO_A][i];
		c = sqrt(gamma* p[i] / rho);
		if (u > c) nsup++;
	}
	avms = avms / double(ncells + 1);

	// print convergence history

	if (iter == 1) drho1 = sqrt(drho / double(ncells + 1));

	drho = sqrt(drho / double(ncells + 1)) / drho1;

	if (!NoOutput) {
		cout << "convergence" << "\t";
		cout << iter << "\t";
		cout << log10(drho) << "\t";
		cout << drmax << "\t";
		cout << idrho << "\t";
		cout << avms << "\t";
		cout << cv[1][1] - cv[1][ib2 - 1] << "\t";
		cout << nsup << "\n";
	}

	return drho;
}

void Solver::PrintResult()
{
	double rho, u, temp, c, mach;

	FILE* file;

	double rgas = (gamma - 1.) * Cp / gamma;

	file = fopen(output_file.c_str(), "w");
	fprintf(file, "x\tA\trho\tu\tp\tT\tM\tmass_flow\n");
	for (int i = 1; i < ib2; ++i)
	{
		rho = fv[RHO][i]; // cv[0][i] / a[i];
		u = fv[U][i]; // cv[1][i] / cv[0][i];
		temp = p[i] / (rgas * rho);
		c = sqrt(gamma * p[i] / rho);
		mach = abs(u) / c;
		fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
			x[i], a[i], rho, u, p[i], temp, mach, cv[RHO_U_A][i]);
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

double Solver::CalcH(double p_, double Rho_, double U_)
{
	return gamma / (gamma - 1.) * p_ / Rho_ + 0.5 * U_ * U_;
}

Solver::equation::equation(string eq_name_, string dt_term_, string dx_term_, map < string, int > vars_, map < string, int > vars_o_)
{
	eq_name = eq_name_;
	dt_var = vars_[dt_term_];
	cur_dt = dt_term(dt_term_, vars_o_);
	cur_dx = dx_term(dx_term_, vars_o_);
}
vector < vector < int > > Solver::equation::dt_term(string dt_term_, map < string, int > vars_o)
{
	vector < vector < int > > cur_dt;
	int cur_pos = 0;
	int id = 0;

	if (dt_term_.find_first_of("E") != string::npos) {
		dt_term_[dt_term_.find_first_of("E")] = 'H';
	}

	cur_dt.push_back(vector < int >());

	while (cur_pos < dt_term_.size()) {
		if (dt_term_.find_first_of("(") == cur_pos ||
			dt_term_.find_first_of(")") == cur_pos ||
			dt_term_.find_first_of("A") == cur_pos) {
			cur_pos++;
			if (cur_pos >= dt_term_.size()) {
				break;
			}
		}
		if (dt_term_.find_first_of("+") == cur_pos) {
			cur_dt.push_back(vector < int >());
			id++;
			cur_pos++;
		}
		for (auto it = vars_o.begin(); it != vars_o.end(); ++it)
		{
			if (dt_term_.find_first_of(it->first, cur_pos) == cur_pos) {
				cur_dt[id].push_back(it->second);
				cur_pos += (it->first).size();
				break;
			}
		}
	}

	return cur_dt;
}

vector < vector < int > > Solver::equation::dx_term(string dx_term_, map < string, int > vars_o)
{
	vector < vector < int > > cur_dx;
	int cur_pos = 0;
	int id = 0;

	cur_dx.push_back(vector < int >());

	while (cur_pos < dx_term_.size()) {
		if (dx_term_.find_first_of("(") == cur_pos ||
			dx_term_.find_first_of(")") == cur_pos ||
			dx_term_.find_first_of("A") == cur_pos) {
			cur_pos++;
			if (cur_pos >= dx_term_.size()) {
				break;
			}
		}
		if (dx_term_.find_first_of("+") == cur_pos) {
			cur_dx.push_back(vector < int >());
			id++;
			cur_pos++;

		}
		for (auto it = vars_o.begin(); it != vars_o.end(); ++it)
		{
			if (dx_term_.find_first_of(it->first, cur_pos) == cur_pos) {
				cur_dx[id].push_back(it->second);
				cur_pos += (it->first).size();
				break;
			}
		}
	}

	return cur_dx;
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

		} else {
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

	sol_init.gas = config["gas"].as< string >();
	sol_init.Cp = config["Cp"].as< double >();
	sol_init.cfl = config["cfl"].as< double >();
	sol_init.max_iter_num = config["max_iter_num"].as< int >();
	sol_init.tolerance = config["tolerance"].as< double >();
	sol_init.gamma = config["gamma"].as< double >();

	sol_init.RK_stage_coeffs = config["RK_stage_coeffs"].as< vector < double > >();;
	sol_init.res_smooth_flag = config["res_smooth_flag"].as< vector < int > >();
	sol_init.eps_impl_res_smooth = config["eps_impl_res_smooth"].as< double >();

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
}

double Sign(double value)
{
	if (value > 0.) return 1.;
	if (value < 0.) return -1.;
	return 0.;
}
