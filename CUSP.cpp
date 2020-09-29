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
		ls.resize(3);
		rs.resize(3);
		for (int i = 0; i < 3; ++i) {
			ls[i].resize(imax, 0.);
			rs[i].resize(imax, 0.);
		}
	}

	double af, bf;
	double delt[3];
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

// CUSPLIM = original CUSP (SLIP) limiter (Eq. (4.121))

double CUSP::CUSPLimiter(double af, double bf)
{
	return 1. - pow( (af - bf) / (abs(af) + abs(bf) + 1e-20), 2.);
}

void CUSP::Fluxes()
{
	double ggm1;
	double si, rl, ul, pl, hl, qrl, rr, ur, pr, hr, qrr, rav, uav, pav, cav, machn, afac, bfac, h1;
	double fcav[3], fdiss[3];
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

const vector < int >& CUSP::GetDissFlag() 
{ 
	return vector < int >(); 
}

const vector < double >& CUSP::GetDissBlend() 
{ 
	return vector < double >(); 
}

void CUSP::Dissipation(double beta) {}
