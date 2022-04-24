// Central Differences Scheme

#include "HLLE.h"

HLLE::HLLE(sol_struct& sol_init_) : Solver(sol_init_)
{
	cout << "Solver is '" << solver_name << "'" << endl;

}

void HLLE::LRState(string var_)
{
	if (ls.size() == 0) {
		ls.resize(var_num);
		rs.resize(var_num);
		for (int i = 0; i < var_num; ++i) {
			ls[i].resize(imax, 0.);
			rs[i].resize(imax, 0.);
		}
	}

	double delt_l, delt_r;
	double* d1 = dummy.data() + 1;


	// first differences of rho, u, p

	for (int i = 0; i < ib2; ++i)
	{
		for (int eq = 0; eq < var_num; ++eq) {
			d1[eq * (imax + 1) + i] = fv[eq][i + 1] - fv[eq][i];
		}
	}
	for (int eq = 0; eq < var_num; ++eq) {
		d1[eq * (imax + 1) - 1] = d1[eq * (imax + 1) + 0];
		d1[eq * (imax + 1) + imax - 1] = d1[eq * (imax + 1) + ib2 - 1];
	}

	// left / right state

	for (int i = 0; i < ib2; ++i)
	{
		for (int eq = 0; eq < var_num; ++eq) {
			int offset = eq * (imax + 1);

			enum class flux_side {left, right};
			auto get_derivatives_dx = [&](int point, flux_side fs, double& dfdx1, double& dfdx2, double& dx) -> void
			{
				int i = fs == flux_side::left
					  ? point
					  : point + 1;

				double dx1 = (x[i] - x[			i - 1 > -1 ?
										i - 1
										  : 0]);
				double dx2 = (x[				i + 2 < ib2 ?
								i + 1 :
								ib2 - 1] - x[i]);

				double df1 = d1[i - 1 + offset];
				double df2 = d1[i + offset];

				dfdx1 = (dx1 > 0. ? df1 / dx1 : 0.);
				dfdx2 = (dx2 > 0. ? df2 / dx2 : 0.);
				dx = fs == flux_side::left
				   ? dx2
				   : dx1;
			};

			double dfdx1, dfdx2, dx;
			get_derivatives_dx(i, flux_side::left, dfdx1, dfdx2, dx);

			delt_l = MUSCL0(dfdx1, dfdx2, dx);
			delt_l *= dx;

			get_derivatives_dx(i, flux_side::right, dfdx1, dfdx2, dx);
			delt_r = MUSCL0(dfdx1, dfdx2, dx);
			delt_r *= dx;

			rs[eq][i] = fv[eq][i + 1] - 0.5 * delt_r;
			ls[eq][i] = fv[eq][i + 0] + 0.5 * delt_l;
		}
		rs[g_P][i] = max(rs[g_P][i], 1e-20);
		ls[g_P][i] = max(ls[g_P][i], 1e-20);
	}

}

double HLLE::MUSCL0(double a, double b, double dx)
{
	double eps = 0.09 * dx * 1e-10 + 1e-30;
	return (a * (b * b + eps) + b * (a * a + eps)) / (a * a + b * b + 2. * eps);
}

double HLLE::MUSCL0_new(double a, double b, double dx)
{
	double eps = 0.09 * dx * 1e-10 + 1e-30;
	return (a * (b * b + eps) + b * (a * a + eps)) / (a * a + b * b + 2. * eps);
}

double HLLE::MUSCL1_3(double a, double b, double dx)
{
	double eps = 0.09 * dx * 1e-10 + 1e-30;
	return (a * (b * b + 2. * eps) + b * (2. * a * a + eps)) / (2. * a * a + 2. * b * b - a * b + 3. * eps);
}

void HLLE::ComputeRHSandJacobian(bool NO_JAC)
{
	vector < double > fluxp(eq_num);
	vector < double > fluxn(eq_num);
	vector < double > source(eq_num);
	vector < double > jac(eq_num * eq_num);
	double si;
	//bool SGS(L_SGS.size() > 0);

	ResetDummy();
#ifdef USE_OMP
#pragma omp parallel for num_threads(MAX_THREAD_NUM) private(si, fluxp, fluxn, source)
#endif
	for (int i = 0; i < ib2; ++i) {
		si = 0.5 * (a[i] + a[i + 1]);
		GetPositiveFluxAndJacobian(i, fluxp, jac);
		if (/*SGS && */!time_expl && !NO_JAC)
			FillJacobian(L_SGS[i], jac, si);

		GetNegativeFluxAndJacobian(i, fluxn, jac);
		if (/*SGS && */!time_expl && !NO_JAC)
			FillJacobian(U_SGS[i], jac, si);

		SetFluxes(i, fluxp, fluxn);

		GetSourceAndJacobian(i, source, jac);
		if (/*SGS && */!time_expl && !NO_JAC)
			FillJacobian(D_SGS[i], jac, si);
	}
	SetRHS();
}

void HLLE::GetSourceAndJacobian(int i, vector < double >& y_val, vector < double >& jac)
{
	vector < double > x_val(eq_num);
	double gamma_ = GetGamma();

	vector<double*> fv_ref(var_num);
	for (int v = 0; v < var_num; ++v)
		fv_ref[v] = &fv[v][i];

	for (int eq = 0; eq < eq_num; ++eq)
		x_val[eq] = make_equation<double>(eq, Solver::equation::term_name::dt, fv_ref);

#ifdef USE_ADEPT
	vector<double_type> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<double_type> y(eq_num);
	ComputeSourceTerm(&x[0], &y[0], i);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
#else
	y_val.resize(eq_num);
	ComputeSourceTerm(&x_val[0], &y_val[0], i);
#endif
}

void HLLE::ComputeSourceTerm(const double_type* cv, double_type* source, int i)
{
	//using double_type;

	double da, dx;
	double gamma_ = GetGamma();

	vector<double_type> fv_a(var_num);
	vector<double_type*> fv_ref(var_num);
	for (int v = 0; v < var_num; ++v)
		fv_ref[v] = &fv_a[v];

	fill_fv_equations<double_type>(filling_type::common, fv_ref, i, true, cv);
	for (int eq = 0; eq < eq_num; ++eq)
		source[eq] = make_equation<double_type>(eq, Solver::equation::term_name::source, fv_ref);
}

void HLLE::SetRHS()
{
	// sum of fluxes = RHS
	vector<double*> fv_ref(var_num);
	vector<double> fv_a(var_num);
	for (int v = 0; v < var_num; ++v)
		fv_ref[v] = &fv_a[v];

	for (int i = 1; i < ib2; ++i)
	{
		fill_fv_equations<double>(filling_type::common, fv_ref, i);
		for (int eq = 0; eq < eq_num; ++eq) {
			rhs[eq][i] = dummy[eq * imax + i] - dummy[eq * imax + i - 1];
			rhs[eq][i] += make_equation<double>(eq, Solver::equation::term_name::source, fv_ref);		// Adding source term
		}
	}
}

void HLLE::GetPositiveFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac)
{
	vector < double > x_val(eq_num /** 3*/);	// 2 additional for gradient calculation
	double gamma_ = GetGamma();

	vector<vector<double*>> ls_vec(/*3*/ 1, vector<double*>(var_num));
	for (int var = 0; var < var_num /** 3*/; ++var)
	{
		int id = (var / var_num);
		if (id > 0)
			id = pow(-1, id);
		id = max(0, i + id);
		ls_vec[var / var_num][var % var_num] = &ls[var % var_num][id];
	}

	for (int eq = 0; eq < eq_num /** 3*/; ++eq)
		x_val[eq] = make_equation<double>(eq % eq_num, Solver::equation::term_name::dt, ls_vec[eq / eq_num]);

#ifdef USE_ADEPT
	vector<double_type> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<double_type> y(eq_num);
	ComputePositiveFlux(&x[0], &y[0], i);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
#else
	y_val.resize(eq_num);
	ComputePositiveFlux(&x_val[0], &y_val[0], i);
#endif
}

void HLLE::GetNegativeFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac)
{
	vector < double > x_val(eq_num /** 3*/);	// 2 additional for gradient calculation
	double gamma_ = GetGamma();

	vector<vector<double*>> rs_vec(/*3*/ 1, vector<double*>(var_num));
	for (int var = 0; var < var_num /** 3*/; ++var)
	{
		int id = (var / var_num);
		if (id > 0)
			id = pow(-1, id);
		id = max(0, i + id);
		rs_vec[var / var_num][var % var_num] = &rs[var % var_num][id];
	}

	for (int eq = 0; eq < eq_num /** 3*/; ++eq)
		x_val[eq] = make_equation<double>(eq % eq_num, Solver::equation::term_name::dt, rs_vec[eq / eq_num]);
	
#ifdef USE_ADEPT
	vector<double_type> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<double_type> y(eq_num);
	ComputeNegativeFlux(&x[0], &y[0], i);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
#else
	y_val.resize(eq_num);
	ComputeNegativeFlux(&x_val[0], &y_val[0], i);
#endif
}

void HLLE::ComputeFlux(const double_type* x_, double_type* fcav, int i, int direction)
{
	//using double_type;
	bool contact = solver_name == "hllc";
	
	double gamma_ = GetGamma();
	double sign_(Sign(direction));
	bool use_WAF = false && i != 0 && i != ib2 - 1;

	vector<double_type*> fv_ref(var_num);
	vector<double_type> fv_val(var_num);
	for (int var = 0; var < var_num; ++var)
		fv_ref[var] = &fv_val[var];

	fill_fv_equations<double_type>(filling_type::common, fv_ref, i, true, x_);

	vector<double_type*> fv_o(var_num);
	vector<double_type> state(var_num);
	for (int var = 0; var < var_num; ++var)
	{
		state[var] = direction > 0 ? rs[var][i] : ls[var][i];
		fv_o[var] = &state[var];
	}

	double_type c_l, c_r;
	double_type RT, uh, Hh, ch, SLm, SRp;

	vector<double_type*>& fv_r = direction > 0 ? fv_o : fv_ref;
	vector<double_type*>& fv_l = direction < 0 ? fv_o : fv_ref;

	c_l = sqrt(gamma_ * *fv_l[g_P] / *fv_l[g_RHO]);
	c_r = sqrt(gamma_ * *fv_r[g_P] / *fv_r[g_RHO]);

	RT = double_sqrt(*fv_r[g_RHO] / *fv_l[g_RHO]);
	uh = (*fv_l[g_U] + RT * *fv_r[g_U]) / (1. + RT);
	Hh = (*fv_l[g_H] + RT * *fv_r[g_H]) / (1. + RT);
	ch = sqrt((gamma_ - 1.) * (Hh - uh * uh / 2.) );

	SLm = double_min(*fv_l[g_U] - c_l, uh - ch);
	SRp = double_max(*fv_r[g_U] + c_r, uh + ch);

	double_type S_star = (*fv_r[g_P] - *fv_l[g_P] + *fv_l[g_RHO] * *fv_l[g_U] * (SLm - *fv_l[g_U]) - *fv_r[g_RHO] * *fv_r[g_U] * (SRp - *fv_r[g_U])) / (*fv_l[g_RHO] * (SLm - *fv_l[g_U]) - *fv_r[g_RHO] * (SRp - *fv_r[g_U]));
	vector<double_type> FL(eq_num), FL_s(eq_num), FHLLE(eq_num), FR_s(eq_num), FR(eq_num);

	//			FL case								FR case
	if ( (SLm >= 0. && direction > 0.) || (SRp <= 0. && direction < 0.) || use_WAF) {
		vector<double_type> fcav_copy{ construct_side_flux_array(fv_ref, i) };
		move(fcav_copy.begin(), fcav_copy.end(), fcav);
		if (use_WAF)
		{
			if (direction > 0.)
				FL.assign(fcav, fcav + eq_num);
			else
				FR.assign(fcav, fcav + eq_num);
		}
	}

	if (!contact)
	{
		// FHLLE
		if ((SLm <= 0. && SRp >= 0.) || use_WAF) {
			vector<double_type> fcav_copy{ construct_hlle_flux_array(fv_ref, SLm, SRp, direction, i) };
			move(fcav_copy.begin(), fcav_copy.end(), fcav);
			if (use_WAF)
			{
				FHLLE.assign(fcav, fcav + eq_num);
			}
		}
	}
	else
	{
		// FHLLC
		//							FL*	case										FR* case
		if ( ((SLm <= 0. && S_star >= 0. && direction > 0.) || ((S_star <= 0. && SRp >= 0. && direction < 0.) || use_WAF)) )
		{
			vector<double_type> fcav_copy{ construct_hllc_flux_array(fv_ref, SLm, SRp, S_star, direction, i) };
			move(fcav_copy.begin(), fcav_copy.end(), fcav);
			if (use_WAF)
			{
				if (direction > 0.)
					FL_s.assign(fcav, fcav + eq_num);
				else
					FR_s.assign(fcav, fcav + eq_num);
			}
		}
	}

	if (use_WAF)
	{
		double dt_dx = dt[i] / (x[i + 1] - x[i]);
		if (!contact)
		{
			vector < double_type > c_k(4);
			c_k[0] = -1.;
			c_k[3] = 1.;
			c_k[1] = min (max (dt_dx * SLm, c_k[0]), c_k[3]);
			c_k[2] = min (max (dt_dx * SRp, c_k[0]), c_k[3]);

			vector < double_type > beta_k(3);
			for (int k = 0; k < 3; ++k)
				beta_k[k] = 0.5 * (c_k[k + 1] - c_k[k]);

			for (int eq = 0; eq < eq_num; ++eq)
				fcav[eq] =   beta_k[0] * FL[eq]   + beta_k[1] * FHLLE[eq]   + beta_k[2] * FR[eq];
		}
		if (contact)
		{
			auto sign = [](double_type a) {if (a < 0.) return -1.; else if (a > 0.) return 1.; else return 0.; };
			vector < double_type > c_k(5);
			c_k[0] = -1.;
			c_k[4] = 1.;
			c_k[1] = min(max(dt_dx * SLm, c_k[0]), c_k[4]);
			c_k[2] = min(max(dt_dx * S_star, c_k[0]), c_k[4]);
			c_k[3] = min(max(dt_dx * SRp, c_k[0]), c_k[4]);

			bool debug = false;

			if (debug)
				cout << "r_k\ti = " << i << "\t";
			vector < double_type > r_k(5);
			for (int k = 0; k < r_k.size(); ++k)
			{
				double upw_rho = 0.;
				double loc_rho = 0.;

				if (c_k[k] > 0)
					upw_rho = fv[g_RHO][i] - fv[g_RHO][i - 1];
				else
					upw_rho = fv[g_RHO][i + 2] - fv[g_RHO][i + 1];

				loc_rho = fv[g_RHO][i + 1] - fv[g_RHO][i];
				if (loc_rho == 0.)
					r_k[k] = 0.;
				else
					r_k[k] = upw_rho / loc_rho;

				if (debug)
					cout << r_k[k] << "\t";
			}
			if (debug)
				cout << "\n";

			vector < double_type > fi_k(5);
			for (int k = 0; k < fi_k.size(); ++k)
			{
				if (r_k[k] <= 0.)
					fi_k[k] = 1.;
				else if (r_k[k] >= 0. && r_k[k] <= 0.5)
					fi_k[k] = 1. - 2. * (1. - abs(c_k[k])) * r_k[k];
				else if (r_k[k] >= 0.5 && r_k[k] <= 1.)
					fi_k[k] = abs(c_k[k]);
				else if (r_k[k] >= 1. && r_k[k] <= 2.)
					fi_k[k] = 1. - (1. - abs(c_k[k])) * r_k[k];
				else if (r_k[k] >= 2.)
					fi_k[k] = 2. * abs(c_k[k]) - 1.;
			}

			vector < double_type > beta_k(4);
			for (int k = 0; k < 4; ++k)
				beta_k[k] = 0.5 * (c_k[k + 1] - c_k[k]);

			if (false)
			{
				for (int eq = 0; eq < eq_num; ++eq)
					fcav[eq] = beta_k[0] * FL[eq] + beta_k[1] * FL_s[eq] + beta_k[2] * FR_s[eq] + beta_k[3] * FR[eq];
			}
			else
			{
				for (int eq = 0; eq < eq_num; ++eq)
					fcav[eq] = 0.5 * (FL[eq] + FR[eq]) - 0.5 * sign(c_k[1]) * fi_k[1] * (FL_s[eq] - FL[eq]) - 0.5 * sign(c_k[2]) * fi_k[2] * (FR_s[eq] - FL_s[eq]) - 0.5 * sign(c_k[3]) * fi_k[3] * (FR[eq] - FR_s[eq]);
			}
		}
	}
}

void HLLE::ComputePositiveFlux(const double_type* x, double_type* fcavp, int i)
{
	ComputeFlux(x, fcavp, i, 1);
}

void HLLE::ComputeNegativeFlux(const double_type* x, double_type* fcavn, int i)
{
	ComputeFlux(x, fcavn, i, -1);
}

void HLLE::SetFluxes(int i, vector < double >& fluxp, vector < double >& fluxn)
{
	double si = 0.5 * (a[i + 1] + a[i]);

	for (int eq = 0; eq < eq_num; ++eq) {
		//dummy[eq * imax + i] *= si;
		dummy[eq * imax + i] += 1. * (fluxp[eq] + fluxn[eq]) * si;
	}
}

const vector < int >& HLLE::GetDissFlag()
{
	return vector < int >();
}

const vector < double >& HLLE::GetDissBlend()
{
	return vector < double >();
}

void HLLE::Dissipation(double beta) {}

void HLLE::GetFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac, bool POS_NEG, bool simple) {}

vector<double_type> HLLE::construct_hlle_flux_array(const vector<double_type*>& vars, const double_type SLm, const double_type SRp, const double direction, const int i)
{
	vector<double_type> flux(eq_num);
	unsigned int eq = 0;

	for (const auto& equation : equations)
	{
		double_type term = (direction > 0. ? SRp : SLm) * direction;
		flux[eq] += make_equation<double_type>(eq, Solver::equation::term_name::dx, vars) * term;
		term = -SRp * SLm * direction;
		flux[eq] += make_equation<double_type>(eq, Solver::equation::term_name::dt, vars) * term;
		flux[eq] /= (SRp - SLm);
		++eq;
	}
	return flux;
}

vector<double_type> HLLE::construct_hllc_flux_array(const vector<double_type*>& vars, const double_type SLm, const double_type SRp, const double_type S_star, const double direction, const int i)
{
	vector<double_type> flux(eq_num);
	vector<double_type> D_star = { 0., 1., S_star };
	double_type S_K = (direction > 0. ? SLm : SRp);
	double_type p_star = S_K * (*vars[g_P] + *vars[g_RHO] * (S_K - *vars[g_U]) * (S_star - *vars[g_U]));

	unsigned int eq = 0;
	for (const auto& equation : equations)
	{
		double_type term = (direction > 0. ? SLm : SRp);
		flux[eq] += make_equation<double_type>(eq, Solver::equation::term_name::dt, vars) * term;
		term = -1.;
		flux[eq] += make_equation<double_type>(eq, Solver::equation::term_name::dx, vars) * term;
		flux[eq] *= S_star;
		flux[eq] += p_star * D_star[eq];
		
		flux[eq] /= (S_K - S_star);
		++eq;
	}
	return flux;
}