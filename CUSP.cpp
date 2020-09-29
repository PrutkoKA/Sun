// Central Differences Scheme

#include "CUSP.h"

CUSP::CUSP(sol_struct& sol_init_) : Solver( sol_init_)
{
	solver_name = "cusp";

	cout << "Solver is '" << solver_name << "'" << endl;

}

void CUSP::LRState()
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
	vector < double > delt(var_num, 0.);
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
		delt[0] = 0.25 * (d1[i + 1] + d1[ i - 1]) *
			  CUSPLimiter(d1[i + 1],  d1[i - 1]);

		delt[1] = 0.25 * (d2[i + 1] + d2[i - 1]) *
			  CUSPLimiter(d2[i + 1],  d2[i - 1]);

		delt[2] = 0.25 * (d3[i + 1] + d3[i - 1]) *
			  CUSPLimiter(d3[i + 1],  d3[i - 1]);

		rs[0][i] = cv[0][i + 1] / a[i + 1] - delt[0];
		rs[1][i] = cv[1][i + 1] / cv[0][i + 1] - delt[1];
		rs[2][i] = p[i + 1] - delt[2];

		ls[0][i] = cv[0][i + 0] / a[i + 0] + delt[0];
		ls[1][i] = cv[1][i + 0] / cv[0][i + 0] + delt[1];
		ls[2][i] = p[i + 0] + delt[2];
	}

}

void CUSP::LRState(string var_)
{
	if (ls.size() == 0) {
		ls.resize(var_num);
		rs.resize(var_num);
		for (int i = 0; i < var_num; ++i) {
			ls[i].resize(imax, 0.);
			rs[i].resize(imax, 0.);
		}
	}

	double delt;
	double* d1 = dummy.data() + 1;

	// first differences of rho, u, p

	for (int i = 0; i < ib2; ++i)
	{
		for (int var = 0; var < var_num; ++var) {
			d1[var * (imax + 1) + i] = fv[var][i + 1] - fv[var][i];
		}
	}
	for (int var = 0; var < var_num; ++var) {
		d1[var * (imax + 1) - 1] = d1[var * (imax + 1) + 0];
		d1[var * (imax + 1) + imax - 1] = d1[var * (imax + 1) + ib2 - 1];
	}

	// left / right state

	for (int i = 0; i < ib2; ++i)
	{
		for (int var = 0; var < var_num; ++var) {
			delt = 0.25 * (d1[var * (imax + 1) + i + 1] + d1[var * (imax + 1) +  i - 1]) *
				  CUSPLimiter(d1[var * (imax + 1) + i + 1],  d1[var * (imax + 1) + i - 1]);

			rs[var][i] = fv[var][i + 1] - delt;
			ls[var][i] = fv[var][i + 0] + delt;
		}
	}

}

// CUSPLIM = original CUSP (SLIP) limiter (Eq. (4.121))

double CUSP::CUSPLimiter(double af, double bf)
{
	return 1. - pow( (af - bf) / (abs(af) + abs(bf) + 1e-20), 2.);
}

void CUSP::Fluxes()
{
	double ggm1;
	double si, rl, ul, pl, hl, qrl, rr, ur, pr, hr, qrr, rav, uav, pav, cav, machn, afac, bfac, h1;
	vector < double > fcav(eq_num, 0.), fdiss(eq_num, 0.);
	double gamma_ = GetGamma();

	ggm1 = gamma_ / (gamma_ - 1.);

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

void CUSP::RHS(int i)
{
	double ggm1;
	double si, rav, uav, pav, cav, machn, afac, bfac, h1;
	vector < double > fcav(eq_num, 0.), fdiss(eq_num, 0.);
	double gamma_ = GetGamma();
	double term_l;
	double term_r;

	ggm1 = gamma_ / (gamma_ - 1.);

	for (int i = 0; i < ib2; ++i) {

		si = 0.5 * (a[i + 1] + a[i]);
		
		// average of left and right convective fluxes
		for (int eq = 0; eq < eq_num; ++eq) {
			vector < vector < int > >& cur_dx = equations[eq].cur_dx;

			fcav[eq] = 0;		// Flux Convective AVerage
			for (int id = 0; id < cur_dx.size(); ++id) 
			{
				term_l = 1.;
				term_r = 1.;
				for (auto var : cur_dx[id]) 
				{
					term_l *= ls[var][i];
					term_r *= rs[var][i];
				}
				fcav[eq] += term_l + term_r;
			}
		}

		// dissipative fluxes

		rav = 0.5 * (fv[RHO][i + 1] + fv[RHO][i]);
		uav = 0.5 * (fv[U][i + 1] + fv[U][i]);
		pav = 0.5 * (fv[P][i+1] + fv[P][i]);
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

		for (int eq = 0; eq < eq_num; ++eq) {
			vector < vector < int > >& cur_dx = equations[eq].cur_dx;
			vector < vector < int > >& cur_dt = equations[eq].cur_dt;

			fdiss[eq] = 0.;
			for (int id = 0; id < cur_dt.size(); ++id) 
			{
				term_l = 1.;
				term_r = 1.;
				for (auto var : cur_dt[id]) 
				{
					term_l *= ls[var][i];
					term_r *= rs[var][i];
				}
				fdiss[eq] += afac * (term_r - term_l);
			}
			for (int id = 0; id < cur_dx.size(); ++id) 
			{
				term_l = 1.;
				term_r = 1.;
				for (auto var : cur_dx[id]) 
				{
					term_l *= ls[var][i];
					term_r *= rs[var][i];
				}
				fdiss[eq] += bfac * (term_r - term_l);
			}
		}

		// total fluxes at i+1/2

		for (int eq = 0; eq < eq_num; ++eq) {
			dummy[eq * imax + i] = 0.5 * (fcav[eq] - fdiss[eq]) * si;
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

const vector < int >& CUSP::GetDissFlag() 
{ 
	return vector < int >(); 
}

const vector < double >& CUSP::GetDissBlend() 
{ 
	return vector < double >(); 
}

void CUSP::Dissipation(double beta) {}
