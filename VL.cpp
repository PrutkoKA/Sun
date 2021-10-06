// Central Differences Scheme

#include "VL.h"

VL::VL(sol_struct& sol_init_) : Solver( sol_init_)
{
	solver_name = "vanleer";

	cout << "Solver is '" << solver_name << "'" << endl;

}

void VL::LRState()
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

	vector < double > delt_l(var_num, 0.);
	vector < double > delt_r(var_num, 0.);

	double* d1 = dummy.data() + 1;
	double* d2 = dummy.data() + 1 + (imax + 1);
	double* d3 = dummy.data() + 1 + (imax + 1) * 2;


	// first differences of rho, u, p

	for (int i = 0; i < ib2; ++i)
	{
		d1[i] = cv[0][i + 1] / a[i + 1] - cv[0][i] / a[i];
		d2[i] = cv[1][i + 1] / cv[0][i + 1] - cv[1][i] / cv[0][i];
		d3[i] = p[i + 1] - p[i];
	}
	d1[-1] = d1[0];
	d2[-1] = d2[0];
	d3[-1] = d3[0];
	d1[imax - 1] = d1[ib2 - 1];
	d2[imax - 1] = d2[ib2 - 1];
	d3[imax - 1] = d3[ib2 - 1];

	// left / right state

	for (int i = 0; i < ib2; ++i)
	{
		// delt[0] = 0.25 * (d1[i + 1] + d1[ i - 1]) *
		// 	  CUSPLimiter(d1[i + 1],  d1[i - 1]);

		delt_l[0] = MUSCL0(d1[i],     d1[i - 1], x[i + 1] - x[i]);
		delt_r[0] = MUSCL0(d1[i + 1], d1[i],     x[i + 1] - x[i]);

		delt_l[1] = MUSCL0(d2[i],     d2[i - 1], x[i + 1] - x[i]);
		delt_r[1] = MUSCL0(d2[i + 1], d2[i],     x[i + 1] - x[i]);

		delt_l[2] = MUSCL0(d3[i],     d3[i - 1], x[i + 1] - x[i]);
		delt_r[2] = MUSCL0(d3[i + 1], d3[i],     x[i + 1] - x[i]);

		// delt_r[0] = 0.25 * (d1[i + 1] + d1[ i - 1]) *
		// 	  CUSPLimiter(d1[i + 1],  d1[i - 1]);

		// delt[1] = 0.25 * (d2[i + 1] + d2[i - 1]) *
		// 	  CUSPLimiter(d2[i + 1],  d2[i - 1]);

		// delt[2] = 0.25 * (d3[i + 1] + d3[i - 1]) *
		// 	  CUSPLimiter(d3[i + 1],  d3[i - 1]);

		//rs[0][i] = cv[0][i + 1] / a[i + 1] - delt[0];

		rs[0][i] = cv[0][i + 1] / a[i + 1] - 0.5 * delt_r[0];
		ls[0][i] = cv[0][i]     / a[i]     + 0.5 * delt_l[0];

		rs[1][i] = cv[1][i + 1] / cv[0][i + 1] - 0.5 * delt_r[1];
		ls[1][i] = cv[1][i + 0] / cv[0][i + 0] + 0.5 * delt_l[1];

		rs[2][i] = p[i + 1] - 0.5 * delt_r[2];
		ls[2][i] = p[i + 0] + 0.5 * delt_l[2];

		// rs[1][i] = cv[1][i + 1] / cv[0][i + 1] - delt[1];
		// rs[2][i] = p[i + 1] - delt[2];

		// ls[0][i] = cv[0][i + 0] / a[i + 0] + delt[0];
		// ls[1][i] = cv[1][i + 0] / cv[0][i + 0] + delt[1];
		// ls[2][i] = p[i + 0] + delt[2];
	}

}

void VL::LRState(string var_)
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
			delt_l = MUSCL0(d1[eq * (imax + 1) + i],     d1[eq * (imax + 1) + i - 1], x[i + 1] - x[i]);
			delt_r = MUSCL0(d1[eq * (imax + 1) + i + 1], d1[eq * (imax + 1) + i],     x[i + 1] - x[i]);

			rs[eq][i] = fv[eq][i + 1] - 0.5 * delt_r;
			ls[eq][i] = fv[eq][i + 0] + 0.5 * delt_l;
		}
		rs[2][i] = max(rs[2][i], 1e-20);
		ls[2][i] = max(ls[2][i], 1e-20);
	}

}

void VL::LRState(vector < vector < double > >& cv_, vector < vector < double > >& ls_, vector < vector < double > >& rs_)
{
	vector < double > dummy_((imax + 1) * var_num, 0.);
	vector < vector < double > > fv_(var_num, vector < double >(imax, 0.));
	double af, bf;
	double delt_l, delt_r;
	double* d1 = dummy_.data() + 1;

	RhoUPH(cv_, fv_);


	// first differences of rho, u, p

	for (int i = 0; i < ib2; ++i)
	{
		for (int eq = 0; eq < var_num; ++eq) {
			d1[eq * (imax + 1) + i] = fv_[eq][i + 1] - fv_[eq][i];
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
			delt_l = MUSCL0(d1[eq * (imax + 1) + i], d1[eq * (imax + 1) + i - 1], x[i + 1] - x[i]);
			delt_r = MUSCL0(d1[eq * (imax + 1) + i + 1], d1[eq * (imax + 1) + i], x[i + 1] - x[i]);

			rs_[eq][i] = fv_[eq][i + 1] - 0.5 * delt_r;
			ls_[eq][i] = fv_[eq][i + 0] + 0.5 * delt_l;
		}
		//rs_[2][i] = max(rs_[2][i], 1e-20);
		//ls_[2][i] = max(ls_[2][i], 1e-20);
	}
	for (int eq = 0; eq < var_num; ++eq) {
		rs_[eq][ib2] = rs_[eq][ib2-1];
		ls_[eq][ib2] = ls_[eq][ib2 - 1];
	}
}

// CUSPLIM = original CUSP (SLIP) limiter (Eq. (4.121))

double VL::MUSCL0(double a, double b, double dx)
{
	double eps = 0.09 * dx * 1e-10 + 1e-30;
	return (a * (b * b + eps) + b * (a * a + eps)) / (a * a + b * b + 2. * eps);
}

double VL::MUSCL1_3(double a, double b, double dx)
{
	double eps = 0.09 * dx * 1e-10 + 1e-30;
	return (a * (b * b + 2. * eps) + b * (2. * a * a + eps)) / (2. * a * a + 2. * b * b - a * b + 3. * eps);
}

// dt("RhoUA")
// dx("RhoUU + p", "A")
// Q("p", "dA")

void VL::Fluxes()
{
	double ggm1;
	double si, rl, ul, pl, hl, qrl, rr, ur, pr, hr, qrr, rav, uav, pav, cav, machn, afac, bfac, h1;
	vector < double > fcav(eq_num, 0.), fdiss(eq_num, 0.);
	double gamma_ = GetGamma();

	ggm1 = gamma_ / (gamma_ - 1.);

//	dt_term("RhoUA");
//	dx_term("RhoUU+p");
//	RHS();

	for (int i = 0; i < ib2; ++i) {

		si = 0.5 * (a[i + 1] + a[i]);
		
		// average of left and right convective fluxes

		rl = ls[0][i];
		ul = ls[1][i];
		pl = ls[2][i];
		hl = ggm1 * pl / rl + 0.5 * ul * ul;
		qrl = ul * rl;

		rr = rs[0][i];
		ur = rs[1][i];
		pr = rs[2][i];
		hr = ggm1 * pr / rr + 0.5 * ur * ur;
		qrr = ur * rr;

		fcav[0] = qrl + qrr;
		fcav[1] = qrl * ul + qrr * ur + pl + pr;
		fcav[2] = qrl * hl + qrr * hr;

		// dissipative fluxes

		rav = 0.5 * (cv[0][i + 1] / a[i + 1] + cv[0][i] / a[i]);
		uav = 0.5 * (cv[1][i + 1] / cv[0][i + 1] + cv[1][i] / cv[0][i]);
		pav = 0.5 * (p[i+1] + p[i]);
		cav = sqrt(gamma_ * pav / rav);
		machn = (uav) / cav;

		if (machn >= 0. && machn < 1.) {
			h1 = 2. * machn - 1.;
			bfac = max(0., h1);
		} else if (machn > -1. && machn < 0.) {
			h1 = 2. * machn + 1.;
			bfac = min(0., h1);
		} else {
			bfac = Sign(h1);
		}

		afac = abs(uav) - bfac * uav;

		fdiss[0] = afac * (rr - rl) + bfac * (qrr - qrl);
		fdiss[1] = afac * (rr * ur - rl * ul) + bfac * (qrr * ur - qrl * ul + pr - pl);
		fdiss[2] = afac * (rr * hr - rl * hl) + bfac * (qrr * hr - qrl * hl);

		// total fluxes at i+1/2

		dummy[0 * imax + i] = 0.5 * (fcav[0] - fdiss[0]) * si;
		dummy[1 * imax + i] = 0.5 * (fcav[1] - fdiss[1]) * si;
		dummy[2 * imax + i] = 0.5 * (fcav[2] - fdiss[2]) * si;

	}

	// sum of fluxes = RHS

	for (int i = 1; i < ib2; ++i)
	{
		rhs[0][i] = dummy[0 * imax + i] - dummy[0 * imax + i - 1];
		rhs[1][i] = dummy[1 * imax + i] - dummy[1 * imax + i - 1];
		rhs[2][i] = dummy[2 * imax + i] - dummy[2 * imax + i - 1];

	}
}

void VL::RHS(int i)
{
	double ggm1;
	double si, rl, ul, pl, hl, qrl, rr, ur, pr, hr, qrr, rav, uav, pav, cav, machn, afac, bfac, h1;
	double cl, cr, machl, machr, MpL, MmR, M, f_mass_p, f_mass_m, u, c;
	vector < double > fcavm(eq_num, 0.), fcavp(eq_num, 0.), fdiss(eq_num, 0.);
	double gamma_ = GetGamma();
	double term;

	ggm1 = gamma_ / (gamma_ - 1.);

	for (int i = 0; i < ib2; ++i) {


		rl = ls[RHO][i];
		ul = ls[U][i];
		pl = ls[P][i];
		cl = sqrt(gamma_ * pl / rl);
		machl = (ul) / cl;

		rr = rs[RHO][i];
		ur = rs[U][i];
		pr = rs[P][i];
		cr = sqrt(gamma_ * pr / rr);
		machr = (ur) / cr;

		if (machl >= 1.) 
		{
			MpL = machl;
		} else if (machl <= -1.) {
			MpL = 0.;
		} else {
			MpL = 0.25 * pow(machl + 1., 2) *Sign(machl);
		}

		if (machr >= 1.) 
		{
			MmR = 0.;
		} else if (machl <= -1.) {
			MmR = machr;
		} else {
			MmR = 0.25 * pow(machr - 1., 2) *Sign(machr);
		}

		M = MpL + MmR;

		if (M >= 1.) {
			for (int eq = 0; eq < eq_num; ++eq) {
				fcavm[eq] = 0.;
				vector < vector < int > >& cur_dx = equations[eq].cur_dx;

				fcavp[eq] = 0.;		// Flux Convective AVerage
				for (int id = 0; id < cur_dx.size(); ++id)
				{
					term = 1.;
					for (auto var : cur_dx[id])
					{
						term *= ls[var][i];
					}
					fcavp[eq] += term;
				}
			}
		}
		else if (M <= -1.) {
			for (int eq = 0; eq < eq_num; ++eq) {
				fcavp[eq] = 0.;
				vector < vector < int > >& cur_dx = equations[eq].cur_dx;

				fcavm[eq] = 0.;		// Flux Convective AVerage
				for (int id = 0; id < cur_dx.size(); ++id)
				{
					term = 1.;
					for (auto var : cur_dx[id])
					{
						term *= rs[var][i];
					}
					fcavm[eq] += term;
				}
			}
		}
		else {
			f_mass_p = rl * cl * pow(machl + 1., 2) * 0.25;
			f_mass_m = -rr * cr * pow(machr - 1., 2) * 0.25;

			fcavp[0] = f_mass_p;
			fcavm[0] = f_mass_m;

			fcavp[1] = f_mass_p * ((-ul + 2. * cl) / gamma_ + ul);
			fcavm[1] = f_mass_m * ((-ur - 2. * cr) / gamma_ + ur);

			fcavp[2] = f_mass_p * (pow((gamma_ - 1.) * ul + 2. * cl, 2)) / (2. * (pow(gamma_, 2) - 1.));
			fcavm[2] = f_mass_m * (pow((gamma_ - 1.) * ur - 2. * cr, 2)) / (2. * (pow(gamma_, 2) - 1.));
		}

		 si = 0.5 * (a[i + 1] + a[i]);

		// total fluxes at i+1/2

		for (int eq = 0; eq < eq_num; ++eq) {
			dummy[eq * imax + i] = (fcavp[eq] + fcavm[eq]) *si;
		}
	}

	// sum of fluxes = RHS

	for (int i = 1; i < ib2; ++i)
	{
		for (int eq = 0; eq < eq_num; ++eq) {
			int dt_var = equations[eq].dt_var;

			rhs[dt_var][i] = dummy[eq * imax + i] - dummy[eq * imax + i - 1];
		}
	}
}

void VL::ComputeRHSandJacobian(bool NO_JAC)
{
	vector < double > fluxp(eq_num);
	vector < double > fluxn(eq_num);
	vector < double > source(eq_num);
	vector < double > jac(eq_num * eq_num);
	double si;
	bool simple = false;
	//bool SGS(L_SGS.size() > 0);

	for (int i = 0; i < ib2; ++i) {
		si = 0.5 * (a[i] + a[i + 1]);
		GetPositiveFluxAndJacobian(i, fluxp, jac, simple);
		if (/*SGS && */!time_expl && !NO_JAC)
			FillJacobian(L_SGS[i], jac, si);

		GetNegativeFluxAndJacobian(i, fluxn, jac, simple);
		if (/*SGS && */!time_expl && !NO_JAC)
			FillJacobian(U_SGS[i], jac, si);

		SetFluxes(i, fluxp, fluxn);

		GetSourceAndJacobian(i, source, jac);
		if (/*SGS && */!time_expl && !NO_JAC)
			FillJacobian(D_SGS[i], jac, si);
	}
	SetRHS();
}

void VL::GetSourceAndJacobian(int i, vector < double >& y_val, vector < double >& jac)
{
	vector < double > x_val(eq_num);
	double gamma_ = GetGamma();

	x_val[RHO_A] = fv[RHO][i];
	x_val[RHO_U_A] = fv[RHO][i] * fv[U][i];
	x_val[RHO_E_A] = (fv[P][i] / (gamma_ - 1.) + 0.5 * fv[RHO][i] * pow(fv[U][i], 2));

	vector<adept::adouble> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<adept::adouble> y(eq_num);
	ComputeSourceTerm(eq_num, &x[0], eq_num, &y[0], i);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
}

void VL::ComputeSourceTerm(int n, const adept::adouble* cv, int m, adept::adouble* source, int i)
{
	using adept::adouble;

	double da, dx;
	adouble r, u, p;
	double gamma_ = GetGamma();

	r = cv[RHO_A];
	u = cv[RHO_U_A] / cv[RHO_A];
	p = (gamma_ - 1.) * (cv[RHO_E_A] - 0.5 * pow(r * u, 2) / r);

	da = 0.5 * (a[i + 1] - a[max(i - 1, 0)]);
	dx = 0.5 * (x[i + 1] - x[max(i - 1, 0)]) + 1e-20;
	source[0] = 0.;
	source[1] = p * da / dx;
	source[2] = 0.;
	//rhs[RHO_U_A][i] = rhs[RHO_U_A][i] - p[i] * da;
}

void VL::GetPositiveFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac, bool simple)
{
	vector < double > x_val(eq_num);
	double gamma_ = GetGamma();

	x_val[RHO_A] = ls[RHO][i];
	x_val[RHO_U_A] = ls[RHO][i] * ls[U][i];
	x_val[RHO_E_A] = (ls[P][i] / (gamma_ - 1.) + 0.5 * ls[RHO][i] * pow(ls[U][i], 2));
	//adept::Stack stack;
	vector<adept::adouble> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<adept::adouble> y(eq_num);
	ComputePositiveFlux(eq_num, &x[0], eq_num, &y[0], i, simple);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
}

void VL::GetPositiveFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac)
{
	GetFluxAndJacobian(i, y_val, cv_, jac, true);
}

void VL::GetNegativeFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac)
{
	GetFluxAndJacobian(i, y_val, cv_, jac, false);
}

void VL::GetFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac, bool POS_NEG, bool simple)
{
	vector < double > x_val(eq_num);
	double gamma_ = GetGamma();

	x_val[RHO_A] = cv_[RHO_A][i];
	x_val[RHO_U_A] = cv_[RHO_U_A][i]; // *cv_[U][i];
	x_val[RHO_E_A] = cv_[RHO_E_A][i]; // (cv_[P][i] / (gamma_ - 1.) + 0.5 * cv_[RHO][i] * pow(cv_[U][i], 2));
	//adept::Stack stack;
	vector<adept::adouble> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<adept::adouble> y(eq_num);
	if (POS_NEG == true) {
		ComputePositiveFlux(eq_num, &x[0], eq_num, &y[0], i, simple);
	}
	else {
		ComputeNegativeFlux(eq_num, &x[0], eq_num, &y[0], i, simple);
	}

	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
}

void VL::GetNegativeFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac, bool simple)
{
	vector < double > x_val(eq_num);
	double gamma_ = GetGamma();

	x_val[RHO_A] = rs[RHO][i];
	x_val[RHO_U_A] = rs[RHO][i] * rs[U][i];
	x_val[RHO_E_A] = (rs[P][i] / (gamma_ - 1.) + 0.5 * rs[RHO][i] * pow(rs[U][i], 2));
	//adept::Stack stack;
	vector<adept::adouble> x(eq_num);
	adept::set_values(&x[0], eq_num, x_val.data());
	stack.new_recording();
	vector<adept::adouble> y(eq_num);
	ComputeNegativeFlux(eq_num, &x[0], eq_num, &y[0], i, simple);
	if (!time_expl) {
		stack.independent(&x[0], eq_num);
		stack.dependent(&y[0], eq_num);
		stack.jacobian(jac.data());
	}
	for (int iy = 0; iy < eq_num; ++iy)
		y_val[iy] = y[iy].value();
}

void VL::ComputeFlux(int n, const adept::adouble * x, int m, adept::adouble * fcav, int i, int direction, bool simple)
{
	using adept::adouble;

	adouble r, u, p_, h, c, mach, f_mass_p;
	double ro, uo, po, ho, co, macho, MpL, MmR, M;	//Other
	double gamma_ = GetGamma();
	double sign_(Sign(direction));

		r = x[RHO_A];
		u = x[RHO_U_A] / x[RHO_A];
		p_ = (gamma_ - 1.) * (x[RHO_E_A] - 0.5 * pow(r * u, 2) / r);
		h = gamma_ / (gamma_ - 1.) * p_ / r + 0.5 * pow(u, 2);
		c = sqrt(gamma_ * p_ / r);
		mach = (u) / c;

		ro = (direction > 0? rs[RHO][i] : ls[RHO][i]);
		uo = (direction > 0 ? rs[U][i] : ls[U][i]);
		po = (direction > 0 ? rs[P][i] : ls[P][i]);
		co = sqrt(gamma_ * po / ro);
		macho = (uo) / co;

		if (direction > 0) {
			if (mach >= 1.)
			{
				MpL = mach.value();
			}
			else if (mach <= -1.) {
				MpL = 0.;
			}
			else {
				MpL = 0.25 * pow(mach.value() + 1., 2) * Sign(mach.value());;// *aSign(machl);
			}

			if (macho >= 1.)
			{
				MmR = 0.;
			}
			else if (macho <= -1.) {
				MmR = macho;
			}
			else {
				MmR = 0.25 * pow(macho - 1., 2) * Sign(macho);// *aSign(machr);
			}
		}
		else {
			if (macho >= 1.)
			{
				MpL = macho;
			}
			else if (macho <= -1.) {
				MpL = 0.;
			}
			else {
				MpL = 0.25 * pow(macho + 1., 2) * Sign(macho);// *aSign(machl);
			}

			if (mach >= 1.)
			{
				MmR = 0.;
			}
			else if (mach <= -1.) {
				MmR = mach.value();
			}
			else {
				MmR = 0.25 * pow(mach.value() - 1., 2) * Sign(mach.value());// *aSign(machr);
			}
		}

		M = MpL + MmR;

		if (M * sign_ >= 1. || simple) {
			fcav[RHO_A] = r * u;
			fcav[RHO_U_A] = r * u * u + p_;
			fcav[RHO_E_A] = r * h * u;
		}
		else if (M * sign_ <= -1.) {
			for (int eq = 0; eq < eq_num; ++eq) {
				fcav[eq] = 0.;
			}
		}
		else {
			f_mass_p = sign_ * r * c * pow(mach + sign_ * 1., 2) * 0.25;
			fcav[0] = f_mass_p;
			fcav[1] = f_mass_p * ((-u + sign_ * 2. * c) / gamma_ + u);
			fcav[2] = f_mass_p * (pow((gamma_ - 1.) * u + sign_ * 2. * c, 2)) / (2. * (pow(gamma_, 2) - 1.));
		}
}

void VL::ComputePositiveFlux(int n, const adept::adouble* x, int m, adept::adouble* fcavp, int i, bool simple)
{
	ComputeFlux(n, x, m, fcavp, i, 1, simple);
}

void VL::ComputeNegativeFlux(int n, const adept::adouble* x, int m, adept::adouble* fcavm, int i, bool simple)
{
	ComputeFlux(n, x, m, fcavm, i, -1, simple);
}

void VL::SetFluxes(int i, vector < double >& fluxp, vector < double >& fluxn)
{
	double si = 0.5 * (a[i + 1] + a[i]);

	for (int eq = 0; eq < eq_num; ++eq) {
		dummy[eq * imax + i] = (fluxp[eq] + fluxn[eq]) * si;
	}
}

void VL::SetRHS()
{
	//double da;

	// sum of fluxes = RHS
	for (int i = 1; i < ib2; ++i)
	{
		//da = 0.5 * (a[i + 1] - a[i - 1]);
		for (int eq = 0; eq < eq_num; ++eq) {
			int dt_var = equations[eq].dt_var;

			rhs[dt_var][i] = dummy[eq * imax + i] - dummy[eq * imax + i - 1];
			/*if (dt_var == RHO_U_A)
				rhs[dt_var][i] = rhs[dt_var][i] - p[i] * da;*/
		}
		SourceTerm(i);
	}
}

const vector < int >& VL::GetDissFlag() 
{ 
	return vector < int >(); 
}

const vector < double >& VL::GetDissBlend() 
{ 
	return vector < double >(); 
}

void VL::Dissipation(double beta) {}
