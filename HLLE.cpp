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

	double af, bf;
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
		rs[P][i] = max(rs[P][i], 1e-20);
		ls[P][i] = max(ls[P][i], 1e-20);
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

	vector<adept::adouble> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<adept::adouble> y(eq_num);
	ComputeSourceTerm(&x[0], &y[0], i);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
}

void HLLE::ComputeSourceTerm(const adept::adouble* cv, adept::adouble* source, int i)
{
	using adept::adouble;

	double da, dx;
	double gamma_ = GetGamma();

	vector<adouble> fv_a(var_num);
	vector<adept::adouble*> fv_ref(var_num);
	for (int v = 0; v < var_num; ++v)
		fv_ref[v] = &fv_a[v];

	fill_fv_equations<adept::adouble>(filling_type::common, fv_ref, i, true, cv);
	for (int eq = 0; eq < eq_num; ++eq)
		source[eq] = make_equation<adouble>(eq, Solver::equation::term_name::source, fv_ref);
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
			rhs[eq][i] += make_equation<double>(eq, Solver::equation::term_name::source, fv_ref);
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

	vector<adept::adouble> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<adept::adouble> y(eq_num);
	ComputePositiveFlux(&x[0], &y[0], i);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
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

	vector<adept::adouble> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<adept::adouble> y(eq_num);
	ComputeNegativeFlux(&x[0], &y[0], i);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
}

void HLLE::ComputeFlux(const adept::adouble* x_, adept::adouble* fcav, int i, int direction)
{
	using adept::adouble;
	bool contact = solver_name == "hllc";
	
	double gamma_ = GetGamma();
	double sign_(Sign(direction));
	bool use_WAF = false && i != 0 && i != ib2 - 1;

	vector<adouble*> fv_ref(var_num);
	vector<adouble> fv_val(var_num);
	for (int var = 0; var < FIELD_VAR_COUNT; ++var)
		fv_ref[var] = &fv_val[var];

	fill_fv_equations<adouble>(filling_type::common, fv_ref, i, true, x_);

	vector<adouble*> fv_o(var_num);
	vector<adouble> state(var_num);
	for (int var = 0; var < FIELD_VAR_COUNT; ++var)
	{
		state[var] = direction > 0 ? rs[var][i] : ls[var][i];
		fv_o[var] = &state[var];
	}

	adouble c_l, c_r;
	adouble RT, uh, Hh, ch, SLm, SRp;

	vector<adouble*>& fv_r = direction > 0 ? fv_o : fv_ref;
	vector<adouble*>& fv_l = direction < 0 ? fv_o : fv_ref;

	c_l = sqrt(gamma_ * *fv_l[P] / *fv_l[RHO]);
	c_r = sqrt(gamma_ * *fv_r[P] / *fv_r[RHO]);

	RT = adept::sqrt(*fv_r[RHO] / *fv_l[RHO]);
	uh = (*fv_l[U] + RT * *fv_r[U]) / (1. + RT);
	Hh = (*fv_l[H] + RT * *fv_r[H]) / (1. + RT);
	ch = sqrt((gamma_ - 1.) * (Hh - uh * uh / 2.) );

	SLm = adept::min(*fv_l[U] - c_l, uh - ch);
	SRp = adept::max(*fv_r[U] + c_r, uh + ch);

	adouble S_star = (*fv_r[P] - *fv_l[P] + *fv_l[RHO] * *fv_l[U] * (SLm - *fv_l[U]) - *fv_r[RHO] * *fv_r[U] * (SRp - *fv_r[U])) / (*fv_l[RHO] * (SLm - *fv_l[U]) - *fv_r[RHO] * (SRp - *fv_r[U]));
	vector<adept::adouble> FL(eq_num), FL_s(eq_num), FHLLE(eq_num), FR_s(eq_num), FR(eq_num);

	//			FL case								FR case
	if ( (SLm >= 0. && direction > 0.) || (SRp <= 0. && direction < 0.) || use_WAF) {
		vector<adept::adouble> fcav_copy = construct_side_flux_array(fv_ref, i);
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
			vector<adept::adouble> fcav_copy = construct_hlle_flux_array(fv_ref, SLm, SRp, direction, i);
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
			vector<adept::adouble> fcav_copy = construct_hllc_flux_array(fv_ref, SLm, SRp, S_star, direction, i);
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
			vector < adouble > c_k(4);
			c_k[0] = -1.;
			c_k[3] = 1.;
			c_k[1] = min (max (dt_dx * SLm, c_k[0]), c_k[3]);
			c_k[2] = min (max (dt_dx * SRp, c_k[0]), c_k[3]);

			vector < adouble > beta_k(3);
			for (int k = 0; k < 3; ++k)
				beta_k[k] = 0.5 * (c_k[k + 1] - c_k[k]);

			fcav[RHO_A] =   beta_k[0] * FL[RHO_A]   + beta_k[1] * FHLLE[RHO_A]   + beta_k[2] * FR[RHO_A];
			fcav[RHO_U_A] = beta_k[0] * FL[RHO_U_A] + beta_k[1] * FHLLE[RHO_U_A] + beta_k[2] * FR[RHO_U_A];
			fcav[RHO_E_A] = beta_k[0] * FL[RHO_E_A] + beta_k[1] * FHLLE[RHO_E_A] + beta_k[2] * FR[RHO_E_A];
		}
		if (contact)
		{
			auto sign = [](adouble a) {if (a < 0.) return -1.; else if (a > 0.) return 1.; else return 0.; };
			vector < adouble > c_k(5);
			c_k[0] = -1.;
			c_k[4] = 1.;
			c_k[1] = min(max(dt_dx * SLm, c_k[0]), c_k[4]);
			c_k[2] = min(max(dt_dx * S_star, c_k[0]), c_k[4]);
			c_k[3] = min(max(dt_dx * SRp, c_k[0]), c_k[4]);

			bool debug = false;

			if (debug)
				cout << "r_k\ti = " << i << "\t";
			vector < adouble > r_k(5);
			for (int k = 0; k < r_k.size(); ++k)
			{
				double upw_rho = 0.;
				double loc_rho = 0.;

				if (c_k[k] > 0)
					upw_rho = fv[RHO][i] - fv[RHO][i - 1];
				else
					upw_rho = fv[RHO][i + 2] - fv[RHO][i + 1];

				loc_rho = fv[RHO][i + 1] - fv[RHO][i];
				if (loc_rho == 0.)
					r_k[k] = 0.;
				else
					r_k[k] = upw_rho / loc_rho;

				if (debug)
					cout << r_k[k] << "\t";
			}
			if (debug)
				cout << "\n";

			vector < adouble > fi_k(5);
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

			vector < adouble > beta_k(4);
			for (int k = 0; k < 4; ++k)
				beta_k[k] = 0.5 * (c_k[k + 1] - c_k[k]);

			if (false)
			{
				fcav[RHO_A] = beta_k[0] * FL[RHO_A] + beta_k[1] * FL_s[RHO_A] + beta_k[2] * FR_s[RHO_A] + beta_k[3] * FR[RHO_A];
				fcav[RHO_U_A] = beta_k[0] * FL[RHO_U_A] + beta_k[1] * FL_s[RHO_U_A] + beta_k[2] * FR_s[RHO_U_A] + beta_k[3] * FR[RHO_U_A];
				fcav[RHO_E_A] = beta_k[0] * FL[RHO_E_A] + beta_k[1] * FL_s[RHO_E_A] + beta_k[2] * FR_s[RHO_E_A] + beta_k[3] * FR[RHO_E_A];
			}
			else
			{
				fcav[RHO_A] = 0.5 * (FL[RHO_A] + FR[RHO_A]) - 0.5 * sign(c_k[1]) * fi_k[1] * (FL_s[RHO_A] - FL[RHO_A]) - 0.5 * sign(c_k[2]) * fi_k[2] * (FR_s[RHO_A] - FL_s[RHO_A]) - 0.5 * sign(c_k[3]) * fi_k[3] * (FR[RHO_A] - FR_s[RHO_A]);
				fcav[RHO_U_A] = 0.5 * (FL[RHO_U_A] + FR[RHO_U_A]) - 0.5 * sign(c_k[1]) * fi_k[1] * (FL_s[RHO_U_A] - FL[RHO_U_A]) - 0.5 * sign(c_k[2]) * fi_k[2] * (FR_s[RHO_U_A] - FL_s[RHO_U_A]) - 0.5 * sign(c_k[3]) * fi_k[3] * (FR[RHO_U_A] - FR_s[RHO_U_A]);
				fcav[RHO_E_A] = 0.5 * (FL[RHO_E_A] + FR[RHO_E_A]) - 0.5 * sign(c_k[1]) * fi_k[1] * (FL_s[RHO_E_A] - FL[RHO_E_A]) - 0.5 * sign(c_k[2]) * fi_k[2] * (FR_s[RHO_E_A] - FL_s[RHO_E_A]) - 0.5 * sign(c_k[3]) * fi_k[3] * (FR[RHO_E_A] - FR_s[RHO_E_A]);
			}
		}
	}
}

void HLLE::ComputePositiveFlux(const adept::adouble* x, adept::adouble* fcavp, int i)
{
	ComputeFlux(x, fcavp, i, 1);
}

void HLLE::ComputeNegativeFlux(const adept::adouble* x, adept::adouble* fcavn, int i)
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

vector<adept::adouble> HLLE::construct_hlle_flux_array(const vector<adept::adouble*>& vars, const adept::adouble SLm, const adept::adouble SRp, const double direction, const int i)
{
	vector<adept::adouble> flux(eq_num);
	unsigned int eq = 0;

	for (const auto& equation : equations)
	{
		adept::adouble term = (direction > 0. ? SRp : SLm) * direction;
		flux[eq] += make_equation<adept::adouble>(eq, Solver::equation::term_name::dx, vars) * term;
		term = -SRp * SLm * direction;
		flux[eq] += make_equation<adept::adouble>(eq, Solver::equation::term_name::dt, vars) * term;
		flux[eq] /= (SRp - SLm);
		++eq;
	}
	return flux;
}

vector<adept::adouble> HLLE::construct_hllc_flux_array(const vector<adept::adouble*>& vars, const adept::adouble SLm, const adept::adouble SRp, const adept::adouble S_star, const double direction, const int i)
{
	vector<adept::adouble> flux(eq_num);
	vector<adept::adouble> D_star = { 0., 1., S_star };
	adept::adouble S_K = (direction > 0. ? SLm : SRp);
	adept::adouble p_star = S_K * (*vars[P] + *vars[RHO] * (S_K - *vars[U]) * (S_star - *vars[U]));

	unsigned int eq = 0;
	for (const auto& equation : equations)
	{
		adept::adouble term = (direction > 0. ? SLm : SRp);
		flux[eq] += make_equation<adept::adouble>(eq, Solver::equation::term_name::dt, vars) * term;
		term = -1.;
		flux[eq] += make_equation<adept::adouble>(eq, Solver::equation::term_name::dx, vars) * term;
		flux[eq] *= S_star;
		flux[eq] += p_star * D_star[eq];
		
		flux[eq] /= (S_K - S_star);
		++eq;
	}
	return flux;
}