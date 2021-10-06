// Solver

#include "Solver.h"
#include "CDS.h"
#include "CUSP.h"
#include "VL.h"
#include "HLLE.h"
//#include "Eigen/Dense"
//
//using namespace Eigen;

//adept::adouble aSign(adept::adouble value);

using namespace YAML;

set < string > solvers { "cds", "cusp", "vanleer", "hlle" };

Solver::Solver(sol_struct& sol_init_) :
					Ideal((sol_init_.gas == "Ideal" ? true : false)),
					Cp(sol_init_.Cp),
					cfl(sol_init_.cfl),
					max_iter_num(sol_init_.max_iter_num),
					tolerance(sol_init_.tolerance),
					gamma(sol_init_.gamma),
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
					steadiness(sol_init_.steadiness), 
					time_stepping(sol_init_.time_stepping)
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

	SetEquation("mass", "RhoA", "RhoUA", vars, vars_o);
	SetEquation("impulse", "RhoUA", "(RhoUU+p)A", vars, vars_o);
	SetEquation("energy", "RhoEA", "RhoHUA", vars, vars_o);
	//SetEquation("grid", "nA", "nA", vars, vars_o);
}

void Solver::SetEquation(string eq_name, string dt_term_, string dx_term_, map < string, int > vars_, map < string, int > vars_o_)
{
	equations.push_back( equation (eq_name, dt_term_, dx_term_, vars_, vars_o_) );
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

void Solver::InitFlowAG2(double* rho_, double* mass_, double* e_, double* p_, double x_)
{
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
}

void Solver::InitFlowAG(double* rho_, double* mass_, double* e_, double* p_, double x_)
{
	cout << "Initializing flow" << endl;

	vector < double > Res, Conc;

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

	grid.AddColumn("rho", cv[RHO_A]);
	grid.AddColumn("rhoU", cv[RHO_U_A]);
	grid.AddColumn("rhoE", cv[RHO_E_A]);
	grid.CalculateResolution(1., 1., "rho", "coordinate");
	grid.CalculateConcentration(1., "coordinate");
	Res = grid.GetResolution();
	Conc = grid.GetConcentration();
	double alpha_ = grid.GetAlpha();
	//for (int i = 1; i < ib2; ++i)
	//{	
	//	cv[N_A][i] = Conc[i] * a[i]; // (Conc[i] - alpha_ * (alpha_ + 1.) * (Conc[i + 1] - 2. * Conc[i] + Conc[i - 1])) / Res[i];
	//}
	//cv[N_A][0] = cv[N_A][1];
	//cv[N_A][imax - 1] = cv[N_A][ib2 - 1];

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
			rhob = rho + Sign(u) * (p_b_out - p[out_id_p1]) / (cs * cs);		// $\rho_b = \rho + \frac{(p_2 - p_{[ib2 - 1]})}{c_s^2}$ 		// fix is here
			ub = u - Sign(u) * (p_b_out - p[out_id_p1]) / (cs * rho);		// $u_b = u - \frac{(p_2 - p_{[ib2 - 1]})}{c_s \rho}$ 		// fix is here
		}

		cv[RHO_A][out_id] = rhob * a[out_id_p1];
		cv[RHO_U_A][out_id] = rhob * ub * a[out_id_p1];
		cv[RHO_E_A][out_id] = (pb / gam1 + 0.5 * rhob * ub * ub) * a[out_id_p1];
		p[out_id] = pb;
	}
	if (b_type == "transmissive")
	{
		double sigma = 1.;
		cv[RHO_A][in_id] = cv[RHO_A][in_id_p1] * sigma * 1;
		cv[RHO_U_A][in_id] = cv[RHO_U_A][in_id_p1] * sigma * 1;
		cv[RHO_E_A][in_id] = cv[RHO_E_A][in_id_p1] * sigma * 1;

		p[in_id] = (gamma - 1.) * (cv[RHO_E_A][in_id] / a[in_id] - 0.5 * cv[RHO_U_A][in_id] * cv[RHO_U_A][in_id] / (cv[RHO_A][in_id] * a[in_id]));


		cv[RHO_A][out_id] = cv[RHO_A][out_id_p1] * sigma;
		cv[RHO_U_A][out_id] = cv[RHO_U_A][out_id_p1] * sigma;
		cv[RHO_E_A][out_id] = cv[RHO_E_A][out_id_p1] * sigma;

		p[out_id] = (gamma - 1.) * (cv[RHO_E_A][out_id] / a[out_id] - 0.5 * cv[RHO_U_A][out_id] * cv[RHO_U_A][out_id] / (cv[RHO_A][out_id] * a[out_id]));
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
	double rrho, rhou, rhoe;

	rrho = a[i] / cv[RHO_A][i];
	rhou = cv[RHO_U_A][i] / a[i];
	rhoe = cv[RHO_E_A][i] / a[i];
	p[i] = (gamma - 1.) * (rhoe - 0.5 * rhou * rhou * rrho);
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

//double Solver::SolveExplicit(double physDt)
//{
//	double fac, adtv, rrho, rhou, rhoe;
//	vector < int > diss_flag_;
//	vector < double > diss_blend_;
//	vector < vector < vector < double > > > cvstage;
//	vector < vector < vector < double > > > rhsstage;
//	vector < vector < double > > rhsold;
//	bool SkipJacComputation;
//
//	iter++;
//
//	if (solver_name == "cds") {
//		diss_flag_ = GetDissFlag();
//		diss_blend_ = GetDissBlend();
//	}
//	for (int i = 0; i < eq_num && solver_name == "cds"; ++i) {
//		diss[i].assign(diss[i].size(), 0.);
//	}
//
//	// Calculate time step
//	{if (!steadiness && time_stepping == 1)
//		if (Global_Time == 0.) {
//			TimeSteps(false);
//			physDt = dt[0];
//		}
//		else {
//			TimeSteps(false, physDt);
//		}
//	else
//		TimeSteps(lts); }		// local = false!
//
//	StagesMemoryAlloc_Init(rhsstage, cvstage);
//	RhoUPH();
//
//	SkipJacComputation = false;
//	// loop over the R. - K. stages:	Blazek, Section 6.1.2, p. 170, "Hybrid multistage schemes"
//	for (int rks = 0; rks < RK_stages_num; ++rks)
//	{
//		if (solver_name == "cds" && diss_flag_[rks] == 1) {	// artificial dissipation
//			Dissipation(diss_blend_[rks]);
//		} else {
//			LRState("");
//		}
//		if (time_expl == false && rks == 1)
//		{
//			SkipJacComputation = true;
//		}
//
//		ComputeRHSandJacobian();					// convective flux residual
//		RHSProcessing(rhsstage, rks, physDt, rhsold);		// RHS preparations, unsteady modification if needed
//		CVProcessing(rhsstage, cvstage, rks);		// update (conservative variables and pressure)
//	}
//
//	// If unsteady and two step RK
//	if (!steadiness && time_stepping == 1) {
//		FinalTSRK_CV_Calculation(rhsstage, cvstage);
//		Global_Time += physDt * cfl;
//		return 0.;
//	}
//
//	//FILE* file;
//	//file = fopen("Output/time_dep.txt", "a");
//	//for (int i = 1; i < ib2; ++i)
//	//{
//	//	fprintf(file, "%lf\t", cv[RHO_U_A][i]);
//	//	// fprintf(file, "%lf\t", p[i]);
//	//}
//	//fprintf(file, "\n");
//	//fclose(file);
//
//	return Convergence();
//}

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
	} else if (!steadiness && time_stepping == 1) {	// If unsteady and explicit two step RK
		FinalTSRK_CV_Calculation(rhsstage, cvstage);
		Global_Time += physDt * cfl;
		return 0.;
	}

	if (false) {

		double rho_[2] = { 1, 0.125 };
		double u_[2] = { 0.75, 0 };
		double mass_[2];
		double p_[2] = { 1., 0.1 };
		double e_[2];

		for (int i = 0; i < 2; ++i) {
			e_[i] = p_[i] / (gamma - 1.) / rho_[i] + 0.5 * u_[i] * u_[i];
			mass_[i] = rho_[0] * u_[i];
		}

		// cusp_s->InverseGrid();
		InitFlowAG2(rho_, mass_, e_, p_, 0.3);

		grid.SetRow("rho", cv[RHO_A]);
		grid.SetRow("rhoU", cv[RHO_U_A]);
		grid.SetRow("rhoE", cv[RHO_E_A]);

		grid.CalculateResolution(1., 1., "rho", "coordinate");
		grid.CalculateConcentration(1., "coordinate");

		x = grid.RefineMesh();

		cv[RHO_A] = grid.GetValues("rho");
		cv[RHO_U_A] = grid.GetValues("rhoU");
		cv[RHO_E_A] = grid.GetValues("rhoE");

		//InitFlowAG2(rho_, mass_, e_, p_, 0.3);

		RefreshBoundaries();						// Refresh boundary conditions
		RhoUPH();
	}

	if (false) {
		for (int i = 0; i < 100; ++i) {
			grid.SetRow("rho", cv[RHO_A]);
			grid.SetRow("rhoU", cv[RHO_U_A]);
			grid.SetRow("rhoE", cv[RHO_E_A]);

			grid.CalculateResolution(1., 1., "rho", "coordinate");
			grid.CalculateConcentration(1., "coordinate");

			x = grid.RefineMesh();

			cv[RHO_A] = grid.GetValues("rho");
			cv[RHO_U_A] = grid.GetValues("rhoU");
			cv[RHO_E_A] = grid.GetValues("rhoE");

			RefreshBoundaries();						// Refresh boundary conditions
			RhoUPH();

			CalculateVolumes();
			grid.SetRow("volume", vol);
		}
	}

	//if (true)
	//{		// New piece of work
	//	vector < double > Res, Conc;

	//	grid.SetRow("rho", cv[RHO_A]);
	//	grid.SetRow("rhoU", cv[RHO_U_A]);
	//	grid.SetRow("rhoE", cv[RHO_E_A]);

	//	Res = grid.GetResolution();

	//	vector < double > rem_cvN_A = cv[N_A];
	//	double coef = 0.5;
	//	coef = min(max(0., coef), 1.);
	//	double side_coef = (1. - coef) / 1.;

	//	for (int j = 1; j < ib2 - 1; ++j)
	//	{
	//		if (j > 1 && j < ib2 - 2) {
	//			cv[N_A][j] = (rem_cvN_A[j - 1] / (x[j] - x[j - 1]) + rem_cvN_A[j + 1] / (x[j + 1] - x[j])) / (1. / (x[j] - x[j - 1]) + 1. / (x[j + 1] - x[j])) * side_coef + rem_cvN_A[j] * coef;
	//		}
	//		if (j == 1) {
	//			cv[N_A][j] = (rem_cvN_A[j + 1] / (x[j + 1] - x[j])) / ( 1. / (x[j + 1] - x[j])) * side_coef * 1. + rem_cvN_A[j] * coef;
	//		}
	//		if (j == ib2 - 2) {
	//			cv[N_A][j] = (rem_cvN_A[j - 1] / (x[j] - x[j - 1])) / (1. / (x[j] - x[j - 1])) * side_coef * 1. + rem_cvN_A[j] * coef;
	//		}
	//		//cv[N_A][j] = (pow(rem_cvN_A[j - 1], 1) + pow(rem_cvN_A[j + 1], 1)) * side_coef + rem_cvN_A[j] * coef;
	//	}
	//	cv[N_A][0] = cv[N_A][0];
	//	cv[N_A][ib2 - 1] = cv[N_A][ib2 - 1];

	//	vector < double > new_x = grid.GetValues("coordinate");
	//	for (int i = 1; i < ib2; ++i) {
	//		double max_dx, dx;

	//		if (i > 1 || i < ib2 - 1) {
	//			max_dx = min(x[i] - x[i - 1], x[i + 1] - x[i]);
	//		}
	//		else {
	//			if (i == 1) {
	//				max_dx = x[i + 1] - x[i];
	//			}
	//			else {
	//				max_dx = x[i] - x[i - 1];
	//			}
	//		}

	//		dx = 0.5 * (x[i] - x[i - 1]) * ((cv[N_A][i] - cvold[N_A][i]) - (cv[N_A][i - 1] - cvold[N_A][i - 1])) / cv[N_A][i] *Sign(pow(Res[i] - Res[i - 1], 1)) * 1.;
	//		dx += 0.5 * (x[i + 1] - x[i]) * ((cv[N_A][i] - cvold[N_A][i]) - (cv[N_A][i + 1] - cvold[N_A][i + 1])) / cv[N_A][i]  *Sign(pow(Res[i + 1] - Res[i], 1))* 1.;
	//		//dx = 0.5 * (x[i] - x[i - 1]) * ((cv[N_A][i] - cvold[N_A][i]) - 0. * (cv[N_A][i - 1] - cvold[N_A][i - 1])) / cv[N_A][i] * Res[i];
	//		dx = min(abs(dx), max_dx / grid.GetAlpha()) * Sign(dx);
	//		new_x[i] += dx;
	//		//new_x[i] -= 0.5 * (x[i] - x[i - 1]) * (cv[N_A][i] - cvold[N_A][i]) / cv[N_A][i] * Res[i];
	//	}
	//	vector < string > ignore;
	//	vector < vector < double > > new_tab;

	//	ignore.push_back(grid.TYPE_COL);
	//	new_tab = grid.NewTable("coordinate", new_x, ignore, false);

	//	grid.SetData(new_tab);
	//	// Restruction of grid;
	//	//grid.FindMesh(0.1, "coordinate");
	//	cv[RHO_A] = grid.GetValues("rho");
	//	cv[RHO_U_A] = grid.GetValues("rhoU");
	//	cv[RHO_E_A] = grid.GetValues("rhoE");

	//	grid.CalculateResolution(1., 1., "rho", "coordinate");
	//	grid.CalculateConcentration(1., "coordinate");
	//	Res = grid.GetResolution();
	//	Conc = grid.GetConcentration();


	//	double alpha_ = grid.GetAlpha();
	//	for (int i = 1; i < ib2; ++i)
	//	{
	//		cv[N_A][i] = (Conc[i]) * a[i]; // -alpha_ * (alpha_ + 1.) * (Conc[i + 1] - 2. * Conc[i] + Conc[i - 1])) / Res[i];
	//	}
	//	cv[N_A][0] = cv[N_A][1];
	//	cv[N_A][imax - 1] = cv[N_A][ib2 - 1];
	//	x = new_x;

	//	RefreshBoundaries();						// Refresh boundary conditions
	//	RhoUPH();
	//}

	//FILE* file;

	//file = fopen("Output/time_dep.txt", "a");
	//for (int i = 1; i < ib2; ++i)
	//{
	//	fprintf(file, "%lf\t", cv[RHO_U_A][i]);
	//	// fprintf(file, "%lf\t", p[i]);
	//}
	//fprintf(file, "\n");
	//fclose(file);

	return Convergence();
}

void Solver::AdjustMesh(double* rho_, double* mass_, double* e_, double* p_, double x_)
{
	InitFlowAG2(rho_, mass_, e_, p_, x_);

	grid.SetRow("rho", cv[RHO_A]);
	grid.SetRow("rhoU", cv[RHO_U_A]);
	grid.SetRow("rhoE", cv[RHO_E_A]);

	grid.CalculateResolution(1., 1., "rho", "coordinate");
	grid.CalculateConcentration(1., "coordinate");

	x = grid.RefineMesh();

	cv[RHO_A] = grid.GetValues("rho");
	cv[RHO_U_A] = grid.GetValues("rhoU");
	cv[RHO_E_A] = grid.GetValues("rhoE");

	InitFlowAG2(rho_, mass_, e_, p_, 0.3);

	RefreshBoundaries();						// Refresh boundary conditions
	RhoUPH();
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
	for (int i = 1; i < ib2; ++i)
	{
		for (int eq = 0; eq < eq_num; ++eq) {
			cv[eq][i] = 0.;
			for (int stage = 0; stage <= rks; ++stage) {
				if (RK_alpha[rks][stage] == 0.) continue;

				cv[eq][i] += RK_alpha[rks][stage] * cvstage[stage][eq][i];
			}
			auto find_lower_index = [&] (int x_index) -> int {
				int y_index = x_index;
				if (x[x_index] == y[x_index])
					return x_index;
				while (x[x_index] < y[y_index] && y_index > 0)
					--y_index;
				while (y_index < y.size() - 1 && y[y_index + 1] <= x[x_index])
					++y_index;
				return y_index;
			};
			int y_ind = find_lower_index(i);
			if (x[i] == y[y_ind])
				cv[eq][i] -= rhs[eq][i];		// Calculating of conservative variables (Eq. 6.7)
			else
				cv[eq][i] -= (y[y_ind + 1] - x[i]) / (y[y_ind + 1] - y[y_ind]) * rhs[eq][y_ind]
						   + (x[i] - y[y_ind])     / (y[y_ind + 1] - y[y_ind]) * rhs[eq][y_ind + 1];
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

//void Solver::LUSGS(int rks, MatrixXd& Li, MatrixXd& Ui, MatrixXd& Di, vector < vector < vector < double > > >& rhsstage, vector < vector < double > >& cvstage)
//{
//	double dsi;
//	double adtv;
//
//	vector < MatrixXd > invD(imax, MatrixXd(eq_num, eq_num));
//	vector < VectorXd > DdWijk(imax, VectorXd(eq_num));
//	VectorXd DdWijkn(eq_num);
//	vector < VectorXd > Wijk(imax, VectorXd(eq_num));
//	VectorXd LdW1(eq_num);
//	VectorXd UdWn(eq_num);
//	VectorXd dW(eq_num);
//	VectorXd Wi(eq_num);
//
//	vector < double > jac(eq_num * eq_num);
//	vector < double > flux(eq_num);
//
//	vector < vector < double > >	ls(var_num, vector < double >(imax, 0.)),
//									rs(var_num, vector < double >(imax, 0.));
//	vector < vector < double > >	ls_(var_num, vector < double >(imax, 0.)),
//									rs_(var_num, vector < double >(imax, 0.));
//	vector < vector < double > >	cvl(eq_num, vector < double >(imax, 0.)),
//									cvr(eq_num, vector < double >(imax, 0.));
//	vector < vector < double > >	cvl_(eq_num, vector < double >(imax, 0.)),
//									cvr_(eq_num, vector < double >(imax, 0.));
//	vector < vector < double > >	cvold_(eq_num, vector < double >(imax, 0.));
//
//	for (int eq = 0; eq < eq_num; ++eq)
//	{
//		for (int i = 0; i < imax; ++i)
//		{
//			cvold_[eq][i] = cvold[eq][i] / a[i];
//		}
//	}
//	cvstage = cvold_;
//
//	Wijk[0](0) = cvold_[0][0];
//	Wijk[0](1) = cvold_[1][0];
//	Wijk[0](2) = cvold_[2][0];
//
//	//LRState(cvold, ls, rs);
//	//cvl = MakeCV(ls);
//	//cvr = MakeCV(rs);
//
//	// implicit residual smoothing
//	/*if (res_smooth_flag[rks] > 0 && eps_impl_res_smooth > 0.)
//		ImplResidualSmooth();*/
//
//		// First sweep
//	for (int i = 1; i < ib2; ++i)
//	{
//		//LRState(cvstage, ls_, rs_);
//		//cvl_ = MakeCV(ls_);
//		//cvr = MakeCV(rs);
//
//		dsi = 0.5 * (a[i] + a[i + 1]);
//
//		adtv = ck[rks] * cfl * dt[i];
//		//Li.hasNaN()
//		Li.setZero(eq_num, eq_num);
//		Ui.setZero(eq_num, eq_num);
//		Li += 0.5*(L_SGS[i] + L_SGS[i].cwiseAbs()) * dsi;
//		Ui += 0.5*(U_SGS[i] - U_SGS[i].cwiseAbs()) * dsi;
//		dsi = -0.5 * (a[i] + a[i - 1]);
//		Li += 0.5*(L_SGS[i - 1] + L_SGS[i].cwiseAbs()) * dsi;
//		Ui += 0.5*(U_SGS[i - 1] - U_SGS[i - 1].cwiseAbs()) * dsi;
//		Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//		Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//
//		Li = 0.5 * ( (L_SGS[i] - L_SGS[i - 1]) * 0.5 * (a[i] + a[i + 1]) + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//		Ui = 0.5 * (-(U_SGS[i] - U_SGS[i - 1]) * 0.5 * (a[i] + a[i - 1]) - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//
//		/*GetFluxAndJacobian(i, flux, cvold_, jac, true, true);
//		Li = GetJacobian(jac);
//		Ui = Li;
//		dsi = 0.5 * (a[i] + a[i + 1]);
//		Li = 0.5 * (Li + Li.cwiseAbs()) * dsi;
//		Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//		dsi = -0.5 * (a[i] + a[i - 1]);
//		Ui = 0.5 * (Ui - Ui.cwiseAbs()) * dsi;
//		Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);*/
//
//		//GetFluxAndJacobian(i, flux, cvl, jac, true);
//		//GetFluxAndJacobian(i, flux, cvold_, jac, true);
//		//Li = GetJacobian(jac);
//		//Li = (Li + Li.cwiseAbs()) / 2.;
//		Li *= /*dsi **/ ak[rks][rks];
//		//ToMatrix(Li);
//
//		/*GetFluxAndJacobian(i, flux, cvr, jac, false);
//		Ui = GetJacobian(jac);
//		Ui = (Ui - Ui.cwiseAbs()) / 2.;
//		Ui *= dsi * ak[rks][rks];
//		ToMatrix(Ui);*/
//
//		dsi = -0.5 * (a[i - 1] + a[i]);		//  ** MINUS ?? **
//
//		/*GetFluxAndJacobian(i - 1, flux, cvl, jac, true);
//		Li += dsi * ak[rks][rks] * ((GetJacobian(jac) + GetJacobian(jac).cwiseAbs()) / 2.);*/
//
//		//Ui = U_SGS[i];
//		//GetFluxAndJacobian(i, flux, cvold_, jac, false);
//		//Ui = GetJacobian(jac);
//		//Ui *= dsi * ak[rks][rks];
//
//		/*GetFluxAndJacobian(i - 1, flux, cvr, jac, false);
//		Ui += dsi * ak[rks][rks] * ((GetJacobian(jac) - GetJacobian(jac).cwiseAbs()) / 2.);*/
//
//		//GetFluxAndJacobian(i, flux, cvold_, jac, false);
//		//Ui = GetJacobian(jac);
//		//Ui = (Ui - Ui.cwiseAbs()) / 2.;
//		Ui *= /*dsi **/ ak[rks][rks];
//		//ToMatrix(Ui);
//
//		Di = D_SGS[i];
//		Di = MakeD(Li, Ui, Di, adtv, i, rks);
//
//		invD[i] = Di.inverse();
//		//ToMatrix(invD[i]);
//		//ToMatrix(Di * invD[i]);
//
//		//dW = Wijk[i - 1] - GetEigenVector(cvold, i - 1);
//		//dW = (GetEigenVector(cvl_, i - 1) - GetEigenVector(cvl, i - 1));
//
//		Li.setZero(eq_num, eq_num);
//		if (i > 0) {
//			dsi = 0.5 * (a[i - 1] + a[i]);
//			Li += 0.5*(L_SGS[i - 1] + L_SGS[i - 1].cwiseAbs()) * dsi;
//		}
//
//		if (i > 1) {
//			dsi = -0.5 * (a[i - 2] + a[i - 1]);
//			Li += 0.5*(L_SGS[i - 2] + L_SGS[i - 2].cwiseAbs()) * dsi;
//			Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i - 1) * omega);
//		}
//		/*Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i - 1) * omega);*/
//		//GetFluxAndJacobian(i - 1, flux, cvold_, jac, true);
//
//		//GetFluxAndJacobian(i - 1, flux, cvold_, jac, true, true);
//		Li = GetJacobian(jac);
//		dsi = 0.5 * (a[i] + a[i - 1]);
//		Li = 0.5 * (Li + Li.cwiseAbs()) * dsi;
//		Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i - 1) * omega);
//
//		Li.setZero(eq_num, eq_num);
//		if (i > 1)
//			Li = 0.5 * ((L_SGS[i - 1] - L_SGS[i - 2]) * 0.5 * (a[i] + a[i - 1]) + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i - 1) * omega);
//		
//		//GetFluxAndJacobian(i - 1, flux, cvl, jac, true);
//		//GetFluxAndJacobian(i - 1, flux, cvold_, jac, true);
//		//Li = GetJacobian(jac);
//		//Li = (Li + Li.cwiseAbs()) / 2.;
//		//Li *= dsi * ak[rks][rks] * (GetEigenVector(cvl_, i - 1) - GetEigenVector(cvl, i - 1)); // dW; // (GetEigenVector(cvl_, i - 1) - GetEigenVector(cvl, i - 1));
//		//if (i == 1) Li = Li.setZero();
//		Li *= /*dsi **/ ak[rks][rks] * (Wijk[i - 1] - GetEigenVector(cvold_, i - 1)) * a[i - 1]; // dW; // (GetEigenVector(cvl_, i - 1) - GetEigenVector(cvl, i - 1));
//		//if (i > 1) {
//		//	dsi = -0.5 * (a[i - 2] + a[i - 1]);		//  ** MINUS ?? **
//
//		//	GetFluxAndJacobian(i - 2, flux, cvl, jac, true);
//		//	Li += dsi * ak[rks][rks] * ((GetJacobian(jac) + GetJacobian(jac).cwiseAbs()) / 2.) * (GetEigenVector(cvl_, i - 2) - GetEigenVector(cvl, i - 2)); // dW; // (GetEigenVector(cvl_, i - 2) - GetEigenVector(cvl, i - 2));
//		//}
//		//else {
//		//	Li *= 0.;
//		//}
//
//		//dW = Wijk[i - 1] - GetEigenVector(cvold, i - 1);
//		//ToVector(dW);
//
//		//LdW1 = Li * dW;
//		LdW1 = Li;
//		DdWijk[i] = -GetEigenVector(rhs, i) * ak[rks][rks] - LdW1;
//		for (int l = 0; l < rks; ++l)
//		{
//			DdWijk[i] -= ak[rks][l] * GetEigenVector(rhsstage[l], i);
//		}
//
//		Wijk[i] = (invD[i] * DdWijk[i] + GetEigenVector(cvold, i)) / a[i];
//		cvstage[0][i] = Wijk[i](0);
//		cvstage[1][i] = Wijk[i](1);
//		cvstage[2][i] = Wijk[i](2);
//	}
//
//	// Second sweep
//	cvstage = cvold_;
//	for (int i = ib2 - 1; i >= 1; --i)
//	{
//		//LRState(cvold, ls_, rs_);
//		//LRState(cvstage, ls_, rs_);
//		//cvl_ = MakeCV(ls_);
//		cvr_ = MakeCV(rs_);
//
//		dsi = -0.5 * (a[i] + a[i + 1]);
//
//		//dW = GetEigenVector(cvstage, i + 1) - GetEigenVector(cvold, i + 1);
//		//dW = (GetEigenVector(cvr_, i) - GetEigenVector(cvr, i));
//
//		Ui = 0.5*(U_SGS[i] - U_SGS[i].cwiseAbs()) * dsi;
//		if (i < ib2 - 1) {
//			dsi = 0.5 * (a[i + 1] + a[i + 2]);
//			Ui += 0.5*(U_SGS[i + 1] - U_SGS[i + 1].cwiseAbs()) * dsi;
//			Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i + 1) * omega);
//		}
//		/*Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i + 1) * omega);*/
//		////GetFluxAndJacobian(i + 1, flux, cvold_, jac, false);
//		//GetFluxAndJacobian(i, flux, cvr, jac, false);
//		//Ui = GetJacobian(jac);
//		//Ui = (Ui - Ui.cwiseAbs()) / 2.;
//		//Ui *= dsi * ak[rks][rks] * (GetEigenVector(cvr_, i) - GetEigenVector(cvr, i)); // dW;
//
//		//GetFluxAndJacobian(i + 1, flux, cvold_, jac, true, true);
//		Ui = GetJacobian(jac);
//		dsi = -0.5 * (a[i] + a[i + 1]);
//		Ui = 0.5 * (Ui - Ui.cwiseAbs()) * dsi;
//		Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i + 1) * omega);
//
//		Ui.setZero(eq_num, eq_num);
//		if (i < ib2 - 1)
//			Ui = 0.5 * (-(U_SGS[i + 1] - U_SGS[i]) * 0.5 * (a[i] + a[i + 1]) - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i + 1) * omega);
//
//		//GetFluxAndJacobian(i + 1, flux, cvr_, jac, false);
//		//GetFluxAndJacobian(i + 1, flux, cvold_, jac, false);
//		//Ui = GetJacobian(jac);
//		//Ui = (Ui - Ui.cwiseAbs()) / 2.;
//		//if (i == ib2 - 1) Ui = Ui.setZero();
//		Ui *= /*dsi **/ ak[rks][rks] * (GetEigenVector(cvstage, i + 1) - GetEigenVector(cvold_, i + 1)) * a[i + 1]; // dW;
//
//		//if (i < ib2 - 1) {
//		//	dsi = 0.5 * (a[i + 1] + a[i + 2]);
//
//		//	GetFluxAndJacobian(i + 1, flux, cvr, jac, false);
//		//	Ui += dsi * ak[rks][rks] * ((GetJacobian(jac) - GetJacobian(jac).cwiseAbs()) / 2.) * (GetEigenVector(cvr_, i + 1) - GetEigenVector(cvr, i + 1)); // dW;
//		//}
//		//else {
//		//	Ui *= 0.;
//		//}
//		//dW = GetEigenVector(cvstage, i + 1) - GetEigenVector(cvold, i + 1);
//		//ToVector(dW);
//
//		//UdWn = Ui * dW;
//		UdWn = Ui;
//		DdWijkn = DdWijk[i] - UdWn;
//
//		dW = invD[i] * DdWijkn;
//		DdWijk[i] = dW;
//		Wi = (GetEigenVector(cvold, i) + dW) / a[i];
//		cvstage[0][i] = Wi(0);
//		cvstage[1][i] = Wi(1);
//		cvstage[2][i] = Wi(2);
//	}
//
//	/*for (int eq = 0; eq < eq_num; ++eq)
//	{
//		for (int i = 0; i < imax; ++i)
//		{
//			cvold_[eq][i] = cvstage[eq][i] * a[i];
//		}
//	}*/
//
//	//LRState(cvold_, ls_, rs_);		// not old ofcourse, it is new
//	//LRState(cvold, ls, rs);
//	
//	/*vector < vector < double > >	cvl_(eq_num, vector < double >(imax, 0.)),
//									cvr_(eq_num, vector < double >(imax, 0.));*/
//
//	/*cvl_ = MakeCV(ls_);
//	cvr_ = MakeCV(rs_);
//	cvl = MakeCV(ls);
//	cvr = MakeCV(rs);*/
//
//	for (int i = 1; i < ib2; ++i)
//	{
//		VectorXd rhsadd;
//
//		//GetFluxAndJacobian(i, flux, cvl_, jac, true, true);
//		//Li = (GetJacobian(jac) + GetJacobian(jac).cwiseAbs()) / 2. * (GetEigenVector(cvl_, i) - GetEigenVector(cvl, i)) * (a[i + 1] + a[i]) / 2.;
//		Li = /*0.5*(L_SGS[i - 1] + L_SGS[i - 1].cwiseAbs())*/ L_SGS[i - 1] * (a[i - 1] + a[i]) / 2.;// *(GetEigenVector(cvstage, i - 1) - GetEigenVector(cvold_, i - 1))* a[i - 1];
//		if (i > 1) {
//			Li -= /*0.5*(L_SGS[i - 2] + L_SGS[i - 2].cwiseAbs())*/ L_SGS[i - 2] * (a[i - 2] + a[i - 1]) / 2.;// *(GetEigenVector(cvstage, i - 1) - GetEigenVector(cvold_, i - 1))* a[i - 1];
//			//Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i - 1) * omega);
//		}
//		else {
//			Li.setZero(eq_num, eq_num);
//		}
//		/*Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i - 1) * omega);*/
//
//		/*GetFluxAndJacobian(i - 1, flux, cvold_, jac, true, true);
//		Li = GetJacobian(jac);
//		dsi = 0.5 * (a[i] + a[i - 1]);
//		Li = 0.5 * (Li + Li.cwiseAbs()) * dsi;*/
//		//Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i - 1) * omega);
//
//		Li *= (GetEigenVector(cvstage, i - 1) - GetEigenVector(cvold_, i - 1)) * a[i - 1];
//		//if (i == 1) Li = Li.setZero();
//
//		//GetFluxAndJacobian(i, flux, cvstage, jac, true);
//		//Li = (GetJacobian(jac) + GetJacobian(jac).cwiseAbs()) / 2. * /*(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))*/DdWijk[i] * (a[i + 1] + a[i]) / 2.;
//		rhsadd = Li*1.;
//
//		//GetFluxAndJacobian(i, flux, cvr_, jac, false, true);
//		//Ui = (GetJacobian(jac) - GetJacobian(jac).cwiseAbs()) / 2. * (GetEigenVector(cvr_, i) - GetEigenVector(cvr, i)) * (a[i + 1] + a[i]) / 2.;
//
//		Ui = -/*0.5*(U_SGS[i] - U_SGS[i].cwiseAbs())*/ U_SGS[i] * (a[i + 1] + a[i]) / 2.;// *(GetEigenVector(cvstage, i + 1) - GetEigenVector(cvold_, i + 1))* a[i + 1];
//		if (i < ib2 - 1) {
//			Ui += /*0.5*(U_SGS[i + 1] - U_SGS[i + 1].cwiseAbs())*/ U_SGS[i + 1] * (a[i + 2] + a[i + 1]) / 2.;// *(GetEigenVector(cvstage, i + 1) - GetEigenVector(cvold_, i + 1))* a[i + 1];
//			//Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i + 1) * omega);
//		}
//		else {
//			Ui.setZero(eq_num, eq_num);
//		}
//
//		/*GetFluxAndJacobian(i + 1, flux, cvold_, jac, true, true);
//		Ui = GetJacobian(jac);
//		dsi = -0.5 * (a[i] + a[i + 1]);
//		Ui = 0.5 * (Ui - Ui.cwiseAbs()) * dsi;*/
//		//Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i + 1) * omega);
//
//		/*Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i + 1) * omega);*/
//		Ui *= (GetEigenVector(cvstage, i + 1) - GetEigenVector(cvold_, i + 1)) * a[i + 1];
//		//if (i == ib2 - 1) Li = Ui.setZero();
//
//		//GetFluxAndJacobian(i, flux, cvstage, jac, false);
//		//Ui = -(GetJacobian(jac) - GetJacobian(jac).cwiseAbs()) / 2. * /*(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))*/DdWijk[i] * (a[i - 1] + a[i]) / 2.;
//		rhsadd += Ui*1;
//
//		//GetFluxAndJacobian(i - 1, flux, cvl_, jac, true, true);
//		//Li = (GetJacobian(jac) + GetJacobian(jac).cwiseAbs()) / 2. * (GetEigenVector(cvl_, i - 1) - GetEigenVector(cvl, i - 1)) * (a[i - 1] + a[i]) / 2.;
//
//		Li = /*0.5*(L_SGS[i] + L_SGS[i].cwiseAbs())*/ L_SGS[i] * (a[i + 1] + a[i]) / 2.;// *(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))* a[i];
//		Li -= /*0.5*(L_SGS[i - 1] + L_SGS[i - 1].cwiseAbs())*/ L_SGS[i - 1] * (a[i - 1] + a[i]) / 2.;// *(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))* a[i];
//		//Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//		Li *= (GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i)) * a[i];
//
//		//GetFluxAndJacobian(i, flux, cvold_, jac, true, true);
//		//Li = GetJacobian(jac);
//		//Ui = Li;
//		//dsi = 0.5 * (a[i] + a[i + 1]);
//		//Li = 0.5 * (Li + Li.cwiseAbs()) * dsi;
//		////Li = 0.5 * (Li + MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//		//Li *= (GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i)) * a[i];
//
//		//GetFluxAndJacobian(i - 1, flux, cvstage, jac, true);
//		//Li = (GetJacobian(jac) + GetJacobian(jac).cwiseAbs()) / 2. * /*(GetEigenVector(cvstage, i - 1) - GetEigenVector(cvold_, i - 1))*/DdWijk[i] * (a[i - 1] + a[i]) / 2.;
//		rhsadd += Li*1;
//
//		//GetFluxAndJacobian(i - 1, flux, cvr_, jac, false, true);
//		//Ui = (GetJacobian(jac) - GetJacobian(jac).cwiseAbs()) / 2. * (GetEigenVector(cvr_, i - 1) - GetEigenVector(cvr, i - 1)) * (a[i - 1] + a[i]) / 2.;
//
//		Ui = /*0.5*(U_SGS[i] - U_SGS[i].cwiseAbs())*/ U_SGS[i]  * (a[i + 1] + a[i]) / 2.;// *(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))* a[i];
//		Ui -= /*0.5*(U_SGS[i - 1] - U_SGS[i - 1].cwiseAbs())*/ U_SGS[i - 1] * (a[i - 1] + a[i]) / 2.;// *(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))* a[i];
//		//Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//		Ui *= (GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i)) * a[i];
//
//		//dsi = -0.5 * (a[i] + a[i - 1]);
//		//Ui = 0.5 * (Ui - Ui.cwiseAbs()) * dsi;
//		////Ui = 0.5 * (Ui - MatrixXd::Identity(eq_num, eq_num) * SpectralRadius(cvold, i) * omega);
//		//Ui *= (GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i)) * a[i];
//
//		//GetFluxAndJacobian(i + 1, flux, cvstage, jac, false);
//		//Ui = -(GetJacobian(jac) - GetJacobian(jac).cwiseAbs()) / 2. * /*(GetEigenVector(cvstage, i + 1) - GetEigenVector(cvold_, i + 1))*/DdWijk[i] * (a[i + 1] + a[i]) / 2.;
//		rhsadd += Ui*1;
//
//		rhsadd -= D_SGS[i] * vol[i] * (GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))/*DdWijk[i]*/ * 1.;
//
//
//		/*rhsadd = 
//			(L_SGS[i] * (GetEigenVector(cvl_, i) - GetEigenVector(cvl, i)) + U_SGS[i] * (GetEigenVector(cvr_, i) - GetEigenVector(cvr, i))) * (a[i + 1] + a[i]) / 2.
//			- (L_SGS[i-1] * (GetEigenVector(cvl_, i-1) - GetEigenVector(cvl, i-1)) + U_SGS[i-1] * (GetEigenVector(cvr_, i-1) - GetEigenVector(cvr, i-1))) * (a[i - 1] + a[i]) / 2.
//			- D_SGS[i] * vol[i] * (GetEigenVector(cvstage, i) - GetEigenVector(cvold, i));*/
//
//		rhsstage[rks][0][i] = rhs[0][i] + rhsadd(0) * 1.;
//		rhsstage[rks][1][i] = rhs[1][i] + rhsadd(1) * 1.;
//		rhsstage[rks][2][i] = rhs[2][i] + rhsadd(2) * 1.;
//
//		adtv = ck[rks] * cfl * dt[i];
//
//		/*rhsstage[rks][0][i] = -(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))(0) * a[i] * vol[i] / adtv;
//		rhsstage[rks][1][i] = -(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))(1) * a[i] * vol[i] / adtv;
//		rhsstage[rks][2][i] = -(GetEigenVector(cvstage, i) - GetEigenVector(cvold_, i))(2) * a[i] * vol[i] / adtv;
//
//		for (int l = 0; l < rks; ++l)
//		{
//			rhsstage[rks][0][i] -= ak[rks][l] * rhsstage[l][0][i];
//			rhsstage[rks][1][i] -= ak[rks][l] * rhsstage[l][1][i];
//			rhsstage[rks][2][i] -= ak[rks][l] * rhsstage[l][2][i];
//		}
//		rhsstage[rks][0][i] /= ak[rks][rks];
//		rhsstage[rks][1][i] /= ak[rks][rks];
//		rhsstage[rks][2][i] /= ak[rks][rks];*/
//	}
//}

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

	Li = Lower(LeftFluxJac, cv_, Lind);
	Ui = Upper(RightFluxJac, cv_, Uind);
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
		R1 = ak[rks][rks] * GetEigenVector(rhs, i);
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

vector < vector < double > > Solver::MakeCV(vector < vector < double > >& fv_)
{
	vector < vector < double > > cv_(eq_num, vector < double > (imax));

	for (int i = 0; i < ib2 /*+ 1*/; ++i)
	{
		cv_[0][i] = fv_[RHO][i]/* * a[i]*/;

		cv_[1][i] = fv_[RHO][i] * fv_[U][i]/* * a[i]*/;

		cv_[2][i] = /*fv_[RHO][i] **/ (fv_[P][i] / (gamma - 1.) + 0.5 * fv_[RHO][i] * pow(fv_[U][i], 2)) /*/ fv_[RHO][i]*//* * a[i]*/;
	}

	return cv_;
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
	for (int i = 0; i < imax; ++i)
	{
		fv[RHO][i] = cv[0][i] / a[i];
		fv[U][i] = cv[1][i] / cv[0][i];
		fv[P][i] = p[i];
		//fv[P][i] = (gamma_ - 1.)* (cv[RHO_E_A][i] / a[i] - 0.5 * cv[RHO_U_A][i] / a[i] * cv[RHO_U_A][i] / a[i] * a[i] / cv[RHO_A][i]);
		fv[H][i] = gamma_ / (gamma_ - 1.) * fv[P][i] / fv[RHO][i] + 0.5 * pow(fv[U][i], 2);
		//fv[n][i] = cv[N_A][i] / a[i];
	}
}

void Solver::RhoUPH(vector < vector < double > >& cv_, vector < vector < double > >& fv_)
{
	double gamma_ = GetGamma();
	double rrho, rhou, rhoe;

	for (int i = 0; i < imax; ++i)
	{
		rrho = a[i] / cv_[RHO_A][i];
		rhou = cv_[RHO_U_A][i] / a[i];
		rhoe = cv_[RHO_E_A][i] / a[i];
		//p_[i] = (gamma_ - 1.) * (rhoe - 0.5 * rhou * rhou * rrho);

		fv_[RHO][i] = cv_[0][i] / a[i];
		fv_[U][i] = cv_[1][i] / cv_[0][i];
		fv_[P][i] = (gamma_ - 1.) * (rhoe - 0.5 * rhou * rhou * rrho);
		fv_[H][i] = gamma_ / (gamma_ - 1.) * fv_[P][i] / fv_[RHO][i] + 0.5 * pow(fv_[U][i], 2);
	}
}

void Solver::TimeSteps(bool local_time, double dt_)		// Blazek, Section 6.1.4, p. 173, "Determination of the maximum time step"
{
	double rho;
	double u;
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
		rho = cv[RHO_A][i] / a[i];
		u = cv[RHO_U_A][i] / cv[RHO_A][i];
		cs = sqrt(gamma * p[i] / rho);
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

	rho = cv_[RHO_A][i] / a[i];
	u = cv_[RHO_U_A][i] / cv_[RHO_A][i];
	rhoe = cv_[RHO_E_A][i] / a[i];
	cs = sqrt(gamma * (gamma - 1.) * (rhoe - 0.5 * rho * u * u) / rho);

	//return cs * sqrt(dx * dx + pow(a[i], 2)) + abs(u) * a[i];
	return (cs  + abs(u)) * a[i];
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

void Solver::SourceTerm(int i)
{
	double da;

	//for (int i = 1; i < ib2; ++i)
	//{

	//da = ((a[i + 1] - a[max(i, 0)]) * (x[i] - x[max(i - 1, 0)] + 1e-20)
	//	+ (a[i] - a[max(i - 1, 0)]) * (x[i + 1] - x[max(i, 0)] + 1e-20))
	//	/ ((x[i + 1] - x[max(i - 1, 0)]) + 1e-20);

		da = 0.5 * (a[i + 1] - a[i - 1]);
		rhs[RHO_U_A][i] = rhs[RHO_U_A][i] - p[i] * da;
	//}
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
		dr = cv[RHO_A][i] - cvold[RHO_A][i];
		//dn = cv[N_A][i] - cvold[N_A][i];

		avms = avms + cv[RHO_U_A][i];
		drho = drho + dr * dr;
		//dconc = dconc + dn * dn;
		if (abs(dr) >= drmax) {
			drmax = abs(dr);
			idrho = i;
		}
		rho = fv[RHO][i];		// cv[RHO_A][i] / a[i];
		u = fv[U][i];			// cv[RHO_U_A][i] / cv[RHO_A][i];
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
		rho = fv[RHO][i]; // cv[0][i] / a[i];
		u = fv[U][i]; // cv[1][i] / cv[0][i];
		temp = p[i] / (rgas * rho);
		c = sqrt(gamma * p[i] / rho);
		mach = /*abs*/(u) / c;
		fprintf(file, "%14.10lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
			x[i], a[i], rho, u, p[i], temp, mach, cv[RHO_U_A][i]);

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

//void Solver::FillJacobian(vector < vector < double > >& M_SGS, vector < double >& jac, double s)
void Solver::FillJacobian(MatrixXd& M_SGS, vector < double >& jac, double s)
{
	for (int i = 0; i < eq_num; ++i)
	{
		for (int j = 0; j < eq_num; ++j)
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
		else if (solver_s == "hlle") {

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
		/*if (config["steadiness"].as< string >() == "steady") {
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
		}*/

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
	if (solver_s == "hlle") {
		return new HLLE(sol_init);

	}
}

double Sign(double value)
{
	if (value > 0.) return 1.;
	if (value < 0.) return -1.;
	return 0.;
}

//adept::adouble aSign(adept::adouble value)
//{
//	if (value > 0.) return 1.;
//	if (value < 0.) return -1.;
//	return 0.;
//}

int count_(vector < int >& vec, int val)
{
	return count(vec.begin(), vec.end(), val);
}