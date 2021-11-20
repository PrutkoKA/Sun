// Central Differences Scheme

#include "CDS.h"

CDS::CDS(sol_struct& sol_init_, cds_struct& cds_init_) : Solver( sol_init_),
			vis2(cds_init_.vis2),
			vis4(1. / cds_init_.vis4),
			diss_blend(cds_init_.diss_blend),
			diss_flag(cds_init_.diss_flag) 
{
	cout << "Solver is '" << solver_name << "'" << endl;

	cout << "\tvis2: " << vis2 << endl;
	cout << "\tvis4: " << vis4 << endl;

	cout << "\tdiss_blend: ";
	for (auto blend_ : diss_blend)
		cout << blend_ << "   ";

	cout << endl;

	cout << "\tdiss_flag: ";
	for (auto flag_ : diss_flag)
		cout << flag_ << "   ";

	cout << endl;
}

const vector < int >& CDS::GetDissFlag()
{
	return diss_flag;
}

const vector < double >& CDS::GetDissBlend()
{
	return diss_blend;
}

void CDS::Dissipation(double beta)
{
	int im1, ip1, ip2, offset;
	double eval, pmax, eps2, eps4, beta1;

	double* dp = dummy.data();
	double* d = dummy.data() + imax;
	double fac, a_, b_;

	// double sign_ (GetInflowId() < GetOutflowId() ? 1. : -1.);

	// pressure sensor (divided second differences)
	for (int i = 1; i < ib2; ++i) 
	{
		dp[i] = abs((p[i + 1] - 2. * p[i] + p[i - 1]) /
					(p[i + 1] + 2. * p[i] + p[i - 1]));
	}
	dp[0] = dp[1];
	dp[imax - 1] = dp[ib2 - 1];

	// dissipation fluxes (at i+1/2)
	for (int i = 0; i < ib2 + 1; ++i)
	{
		im1 = max(i - 1, 0);
		ip1 = min(i + 1, imax - 1);
		ip2 = min(i + 2, imax - 1);

		eval = 0.5 * (vol[i] / dt[i] + vol[ip1] / dt[ip1]); // equation 4.56 - scaling factor
		pmax = max(dp[i], dp[ip1]);
		eps2 = eval * vis2 * pmax;
		eps4 = eval * vis4;
		eps4 = max(0., eps4 - eps2);
		for (int eq = 0; eq < eq_num; ++eq)
			d[eq * imax + i] = eps2 * (cv[eq][ip1] - cv[eq][i]) -
				eps4 * (cv[eq][ip2] - 3. * cv[eq][ip1] + 3. * cv[eq][i] - cv[eq][im1]);

	}
	
	beta1 = 1. - beta;
	offset = 1;
	for (int i = 1; i < ib2; ++i)
	{
		if (i == 1 || i == ib2 - 1 || true) {
			a_ = 1.;
			b_ = 1.;
		} else {
			a_ = 1. / (x[i] - x[i - 1]);
			b_ = 1. / (x[i + 1] - x[i]);
			fac = a_ + b_;
			a_ = a_ / fac * 2.;
			b_ = b_ / fac * 2.;
		}
		for (int eq = 0; eq < eq_num; ++eq)
			diss[eq][i] = offset * beta * (b_ * d[eq * imax + i] - a_ * d[eq * imax + i - offset]) + beta1 * diss[eq][i];

	}
}

void CDS::Fluxes()
{
	double si, rav, ruav, reav, pav, rhav, qs;
	double* f = dummy.data();
	int offset = 1;
	double fac, a_, b_;

	// flux term (average of variables)

	for (int i = 0; i < ib2 + 1*0; ++i)
	{
		si = 0.5 * (a[i + 1] + a[i]);
		rav = 0.5 * (cv[RHO_A][i + 1] / a[i + 1] + cv[RHO_A][i] / a[i]); //(fv[RHO][i + i] + fv[RHO][i]); // (cv[0][i + 1] / a[i + 1] + cv[0][i] / a[i]);
		ruav = 0.5 * (cv[RHO_U_A][i + 1] / a[i + 1] + cv[RHO_U_A][i] / a[i]);
		reav = 0.5 * (cv[RHO_E_A][i + 1] / a[i + 1] + cv[RHO_E_A][i] / a[i]);
		pav = 0.5 * (p[i + 1] + p[i]);
		rhav = reav + pav;
		qs = ruav * si / rav;

		f[0 * imax + i] = rav * qs;
		f[1 * imax + i] = ruav * qs + pav * si;
		f[2 * imax + i] = rhav * qs;
	}

	// flux + dissipation = RHS

	for (int i = 1; i < ib2; ++i)
	{
		if (i == 1 || i == ib2 - 1 || true) {
			a_ = 1.;
			b_ = 1.;
		} else {
			a_ = 1. / (x[i] - x[i - 1]);
			b_ = 1. / (x[i + 1] - x[i]);
			fac = a_ + b_;
			a_ = a_ / fac * 2.;
			b_ = b_ / fac * 2.;
		}
		for (int eq = 0; eq < eq_num; ++eq)
			rhs[eq][i] = offset * (b_ * f[eq * imax + i] - a_ * f[eq * imax + i - offset]) - diss[eq][i];

	}
}

void CDS::RHS(int i)
{
	double si, rav, ruav, reav, pav, rhav, qs;
	double* f = dummy.data();
	int offset = 1;
	double fac, a_, b_;
	double term_;
	bool U_IS_PRESENT;

	// flux term (average of variables)

	for (int i = 0; i < ib2 + 1*0; ++i)
	{
		si = 0.5 * (a[i + 1] + a[i]);

		for (int eq = 0; eq < eq_num; ++eq) {
			vector < vector < int > >& cur_dx = equations[eq].cur_dx;

			f[eq * imax + i] = 0;		// Flux Convective AVerage
			for (int id = 0; id < cur_dx.size(); ++id) 
			{
				term_ = 1.;
				for (auto var : cur_dx[id]) 
					term_ *= 0.5 * (fv[var][i + 1] + fv[var][i]);

				term_ *= si;
				f[eq * imax + i] += term_;
			}
		}
	}

	// flux + dissipation = RHS

	for (int i = 1; i < ib2; ++i)
	{
		if (i == 1 || i == ib2 - 1 || true) {
			a_ = 1.;
			b_ = 1.;
		}
		else {
			a_ = 1. / (x[i] - x[i - 1]);
			b_ = 1. / (x[i + 1] - x[i]);
			fac = a_ + b_;
			a_ = a_ / fac * 2.;
			b_ = b_ / fac * 2.;
		}
		for (int eq = 0; eq < eq_num; ++eq)
			rhs[eq][i] = offset * (b_ * f[eq * imax + i] - a_ * f[eq * imax + i - offset]) - diss[eq][i];

	}
}

void CDS::ComputeRHSandJacobian(bool NO_JAC)
{

}

void CDS::LRState()
{
	// Do nothing
}

void CDS::LRState(string var_)
{
	// Do nothing
}

void CDS::LRState(vector < vector < double > >& cv_, vector < vector < double > >& ls_, vector < vector < double > >& rs_)
{
	// Do nothing
}

void CDS::GetFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac, bool POS_NEG, bool simple)
{}