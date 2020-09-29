// Central Differences Scheme

#include "CDS.h"

CDS::CDS(sol_struct& sol_init_, cds_struct& cds_init_) : Solver( sol_init_),
			vis2(cds_init_.vis2),
			vis4(1. / cds_init_.vis4),
			diss_blend(cds_init_.diss_blend),
			diss_flag(cds_init_.diss_flag) 
{
	solver_name = "cds";

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
		// if (direction > 0 || true) {
			im1 = max(i - 1, 0);
			ip1 = min(i + 1, imax - 1);
			ip2 = min(i + 2, imax - 1);
		// } else {
		// 	im1 = min(i + 1, imax - 1);
		// 	ip1 = max(i - 1, 0);
		// 	ip2 = max(i - 2, 0);
		// }
		eval = 0.5 * (vol[i] / dt[i] + vol[ip1] / dt[ip1]); // equation 4.56 - scaling factor
		pmax = max(dp[i], dp[ip1]);
		eps2 = eval * vis2 * pmax;
		eps4 = eval * vis4;
		eps4 = max(0., eps4 - eps2);
		d[0 * imax + i] = eps2 * (cv[0][ip1] - cv[0][i]) -
				  eps4 * (cv[0][ip2] - 3. * cv[0][ip1] + 3. * cv[0][i] - cv[0][im1]);

		d[1 * imax + i] = eps2 * (cv[1][ip1] - cv[1][i]) -
				  eps4 * (cv[1][ip2] - 3. * cv[1][ip1] + 3. * cv[1][i] - cv[1][im1]);

		d[2 * imax + i] = eps2 * (cv[2][ip1] - cv[2][i]) -
				  eps4 * (cv[2][ip2] - 3. * cv[2][ip1] + 3. * cv[2][i] - cv[2][im1]);	
	}
	
	beta1 = 1. - beta;
	offset = 1;// * direction;
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
		diss[0][i] = offset*beta * (b_ * d[0 * imax + i] - a_ * d[0 * imax + i - offset]) + beta1 * diss[0][i];
		diss[1][i] = offset*beta * (b_ * d[1 * imax + i] - a_ * d[1 * imax + i - offset]) + beta1 * diss[1][i];
		diss[2][i] = offset*beta * (b_ * d[2 * imax + i] - a_ * d[2 * imax + i - offset]) + beta1 * diss[2][i];
	}
}

void CDS::Fluxes()
{
	double si, rav, ruav, reav, pav, rhav, qs;
	double* f = dummy.data();
	int offset = 1;// * direction;
	double fac, a_, b_;
	// double sign_ (GetInflowId() < GetOutflowId() ? 1. : -1.);

	// flux term (average of variables)

	for (int i = 0; i < ib2 + 1; ++i)
	{
		si = 0.5 * (a[i + 1] + a[i]);
		rav = 0.5 * (cv[0][i + 1] / a[i + 1] + cv[0][i] / a[i]);
		ruav = 0.5 * (cv[1][i + 1] / a[i + 1] + cv[1][i] / a[i]);
		reav = 0.5 * (cv[2][i + 1] / a[i + 1] + cv[2][i] / a[i]);
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
		rhs[0][i] = offset*(b_ * f[0 * imax + i] - a_ * f[0 * imax + i - offset]) - diss[0][i];
		rhs[1][i] = offset*(b_ * f[1 * imax + i] - a_ * f[1 * imax + i - offset]) - diss[1][i];
		rhs[2][i] = offset*(b_ * f[2 * imax + i] - a_ * f[2 * imax + i - offset]) - diss[2][i];
	}
}

void CDS::LRState()
{
	// Do nothing
}