#include "LavalTest.h"
#include "CDS.h"

//double gamma_;
//double cpgas;
//double p01;
//double t01;
//double p2;
//
//double gam1;
//double gap1;
//double rgas;
//double temp;
//double rho;
//double mach;
//double cs;
//double u;
//double mass;
//double e;
//
//vector < double > p;
//vector < vector < double > > cv;
//vector < vector < double > > cvold;
//vector < vector < double > > diss;
//vector < vector < double > > rhs;
//vector < vector < double > > rs;
//vector < double > dt;
//vector < double > ls;
//vector < double > dum;
//
//vector < double > l_sections;
//vector < double > l_coords;
//vector < double > l_volumes;
//
//double volref;
//double rhoref;
//double uref;
//double pref;
//double vis2(0.7), vis4(64.);	// artificial dissipation coefficient - k2   (central scheme)
//								// artificial dissipation coefficient - 1/k4 (central scheme)
//double cfl = 0.5;
//double epsirs = 0.8;      // coefficient of implicit residual smoothing
//
//double drho1;
//
//int iter;
//int imax;
//int ib2;
//int ncells;
//int maxiter = 50000;
//double convtol = 1e-5;
//int mxdum;
//
//vector < double > ark { 0.2500, 0.1667, 0.3750, 0.5000, 1.000 };		// stage coefficients
//vector < double > betrk { 1.00  , 0.00  , 0.56  , 0.00  , 0.44 };		// dissipation blending coeff.
//vector < int > ldiss { 1, 0, 1, 0, 1 };	// dissipation evaluation (1=yes)
//vector < int > lsmoo{ 1, 1, 1, 1, 1 };	// residual smoothing (1=yes)
//
//char kdissip = 'c';		// central scheme with artificial dissipation
//
//double drho;

void Laval2()	// Standard example as in Blazek (Central scheme)
{
	cout << endl << "  ---  Laval2 test  ---  " << endl;

	string file_name = "Input/solver.yml";

	Loop laval;
	Solver *central_s;

	central_s = CreateReadConfigFile(file_name);

	laval.ReadFile("Input/grid.txt");
	central_s->SetGrid(laval);
	central_s->ReadBoundaries("Input/boundary.yml");
	central_s->SetOutputFile("Output/output.txt");
	central_s->HideOutput();
	// central_s->InverseGrid();

	double convtol = central_s->GetTolerance();
	double drho = 1.;

	int maxiter = central_s->GetMaxIterNum();
	double gamma__ = central_s->GetGamma();
	double cpgas_ = central_s->GetCp();
	double t01_ = central_s->GetInflowT();
	double p01_ = central_s->GetInflowP();
	double p2_ = central_s->GetOutflowP();
	int id = central_s->GetInflowId();
	double section = central_s->GetSection()[(id > 0 ? id - 2 : id + 1)];

	double gam1_ = gamma__ - 1.;	// $\gamma_ - 1$
	double gap1_ = gamma__ + 1.;	// $\gamma_ + 1$
	double rgas_ = gam1_ * cpgas_ / gamma__;	// $r_{gas} = (\gamma_ - 1) \frac{c_p}{\gamma_}$
	double temp_ = t01_ * pow((p2_ / p01_), gam1_ / gamma__);		// $T = T_{01} \left( \frac{p_2} {p_{01}} \right)^\frac{\gamma_ - 1}{\gamma_}$
	double rho_ = p2_ / (rgas_ * temp_);		// $\rho = \frac{p_2} {r_{gas}T}$
	double mach_ = sqrt(2. * ((t01_ / temp_) - 1.) / gam1_);	// $ M = \sqrt{ 2 \frac{ T_{01} / T - 1}{\gamma_ - 1} } $
	double cs_ = sqrt(gamma__ * p2_ / rho_);		// $c_s = \sqrt{\frac{\gamma_ p_2}{\rho}}$
	double u_ = cs_ * mach_;		// $u = c_s M$
	double mass_ = rho_ * u_ * section;	// $ Q_M = \rho u S_{[1]}$
	double e_ = (cpgas_ - rgas_) * t01_;		// $E = (c_{p gas} - r_{gas}) T_{01}$

	// central_s->InverseGrid();
	central_s->InitFlow(rho_, mass_, e_, p2_);
	central_s->RefreshBoundaries();

	// convtol = 1e-5;
	for (int iter = 0; iter < maxiter && drho > convtol; ++iter)
	{
		drho = central_s->Solve();
	}

	central_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

void Laval3()	// Changing uniform grid (Central scheme)
{
	cout << endl << "  ---  Laval3 test  ---  " << endl;

	string file_name = "Input/solver.yml";

	Loop laval;
	Solver *central_s;
	vector < double > range_;
	vector < vector < double > > new_tab;

	central_s = CreateReadConfigFile(file_name);

	laval.ReadFile("Input/grid.txt");
	range_ = laval.MakeRange(0., 1., 0.1);
	range_.push_back(range_[range_.size() - 1]);
	range_.insert(range_.begin(), range_.begin(), range_.begin() + 1);
	new_tab = laval.NewTable("coordinate", range_);
	laval.SetData(new_tab);
	// laval.ShowData();

	central_s->SetGrid(laval);
	central_s->ReadBoundaries("Input/boundary.yml");
	central_s->SetOutputFile("Output/output2.txt");
	// central_s->HideOutput();
	// central_s->InverseGrid();

	double convtol = central_s->GetTolerance();
	double drho = 1.;

	int maxiter = central_s->GetMaxIterNum();
	double gamma__ = central_s->GetGamma();
	double cpgas_ = central_s->GetCp();
	double t01_ = central_s->GetInflowT();
	double p01_ = central_s->GetInflowP();
	double p2_ = central_s->GetOutflowP();
	int id = central_s->GetInflowId();
	double section = central_s->GetSection()[(id > 0 ? id - 2 : id + 1)];

	double gam1_ = gamma__ - 1.;	// $\gamma_ - 1$
	double gap1_ = gamma__ + 1.;	// $\gamma_ + 1$
	double rgas_ = gam1_ * cpgas_ / gamma__;	// $r_{gas} = (\gamma_ - 1) \frac{c_p}{\gamma_}$
	double temp_ = t01_ * pow((p2_ / p01_), gam1_ / gamma__);		// $T = T_{01} \left( \frac{p_2} {p_{01}} \right)^\frac{\gamma_ - 1}{\gamma_}$
	double rho_ = p2_ / (rgas_ * temp_);		// $\rho = \frac{p_2} {r_{gas}T}$
	double mach_ = sqrt(2. * ((t01_ / temp_) - 1.) / gam1_);	// $ M = \sqrt{ 2 \frac{ T_{01} / T - 1}{\gamma_ - 1} } $
	double cs_ = sqrt(gamma__ * p2_ / rho_);		// $c_s = \sqrt{\frac{\gamma_ p_2}{\rho}}$
	double u_ = cs_ * mach_;		// $u = c_s M$
	double mass_ = rho_ * u_ * section;	// $ Q_M = \rho u S_{[1]}$
	double e_ = (cpgas_ - rgas_) * t01_;		// $E = (c_{p gas} - r_{gas}) T_{01}$

	// central_s->InverseGrid();
	central_s->InitFlow(rho_, mass_, e_, p2_);
	central_s->RefreshBoundaries();

	// convtol = 1e-5;
	for (int iter = 0; iter < maxiter && drho > convtol; ++iter)
	{
		drho = central_s->Solve();
	}

	central_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

void Laval4()	// Examining non-uniform grid (Central scheme)
{
	cout << endl << "  ---  Laval4 test  ---  " << endl;

	string file_name = "Input/solver.yml";

	Loop laval;
	Solver *central_s;
	vector < double > range_(23, 0.);
	vector < vector < double > > new_tab;

	central_s = CreateReadConfigFile(file_name);

	laval.ReadFile("Input/grid.txt");
	
	double sum = 0.;
	double fac = 1.3;
	double h = 0.0115;
	range_[21] = 1.;
	range_[22] = 1.;
	for (int i = 2; i < 13; ++i) {
		range_[i] = range_[i-1] + h * pow(fac, double(12 - i));
		sum += range_[i] - range_[i-1];
		cout << i << "\t" << range_[i] << "\t" << sum << "\t" << pow(fac, double(12 - i)) << endl;
	}
	for (int i = 13; i < 21; ++i) {
		range_[i] = range_[i-1] + h * pow(fac, double(i - 13));
		sum += range_[i] - range_[i-1];
		cout << i << "\t" << range_[i] << "\t" << sum << "\t" << pow(fac, double(i - 13)) << endl;
	}

	new_tab = laval.NewTable("coordinate", range_);
	laval.SetData(new_tab);
//	laval.ShowData();

	central_s->SetGrid(laval);
	central_s->ReadBoundaries("Input/boundary.yml");
	central_s->SetOutputFile("Output/output3.txt");
	// central_s->HideOutput();
	// central_s->InverseGrid();

	double convtol = central_s->GetTolerance();
	double drho = 1.;

	int maxiter = central_s->GetMaxIterNum();
	double gamma__ = central_s->GetGamma();
	double cpgas_ = central_s->GetCp();
	double t01_ = central_s->GetInflowT();
	double p01_ = central_s->GetInflowP();
	double p2_ = central_s->GetOutflowP();
	int id = central_s->GetInflowId();
	double section = central_s->GetSection()[(id > 0 ? id - 2 : id + 1)];

	double gam1_ = gamma__ - 1.;	// $\gamma_ - 1$
	double gap1_ = gamma__ + 1.;	// $\gamma_ + 1$
	double rgas_ = gam1_ * cpgas_ / gamma__;	// $r_{gas} = (\gamma_ - 1) \frac{c_p}{\gamma_}$
	double temp_ = t01_ * pow((p2_ / p01_), gam1_ / gamma__);		// $T = T_{01} \left( \frac{p_2} {p_{01}} \right)^\frac{\gamma_ - 1}{\gamma_}$
	double rho_ = p2_ / (rgas_ * temp_);		// $\rho = \frac{p_2} {r_{gas}T}$
	double mach_ = sqrt(2. * ((t01_ / temp_) - 1.) / gam1_);	// $ M = \sqrt{ 2 \frac{ T_{01} / T - 1}{\gamma_ - 1} } $
	double cs_ = sqrt(gamma__ * p2_ / rho_);		// $c_s = \sqrt{\frac{\gamma_ p_2}{\rho}}$
	double u_ = cs_ * mach_;		// $u = c_s M$
	double mass_ = rho_ * u_ * section;	// $ Q_M = \rho u S_{[1]}$
	double e_ = (cpgas_ - rgas_) * t01_;		// $E = (c_{p gas} - r_{gas}) T_{01}$

	// central_s->InverseGrid();
	central_s->InitFlow(rho_, mass_, e_, p2_);
	central_s->RefreshBoundaries();

	// convtol = 1e-5;
	for (int iter = 0; iter < maxiter && drho > convtol; ++iter)
	{
		drho = central_s->Solve();
	}

	central_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

void Laval5()	// Standard example as in Blazek (CUSP)
{
	cout << endl << "  ---  Laval5 test  ---  " << endl;

	string file_name = "Input/cusp.yml";

	Loop laval;
	Solver *cusp_s;

	cusp_s = CreateReadConfigFile(file_name);

	laval.ReadFile("Input/grid.txt");
	cusp_s->SetGrid(laval);
	cusp_s->ReadBoundaries("Input/boundary.yml");
	cusp_s->SetOutputFile("Output/output4.txt");
	// cusp_s->HideOutput();
	// cusp_s->InverseGrid();

	double convtol = cusp_s->GetTolerance();
	double drho = 1.;

	int maxiter = cusp_s->GetMaxIterNum();
	double gamma__ = cusp_s->GetGamma();
	double cpgas_ = cusp_s->GetCp();
	double t01_ = cusp_s->GetInflowT();
	double p01_ = cusp_s->GetInflowP();
	double p2_ = cusp_s->GetOutflowP();
	int id = cusp_s->GetInflowId();
	double section = cusp_s->GetSection()[(id > 0 ? id - 2 : id + 1)];

	double gam1_ = gamma__ - 1.;	// $\gamma_ - 1$
	double gap1_ = gamma__ + 1.;	// $\gamma_ + 1$
	double rgas_ = gam1_ * cpgas_ / gamma__;	// $r_{gas} = (\gamma_ - 1) \frac{c_p}{\gamma_}$
	double temp_ = t01_ * pow((p2_ / p01_), gam1_ / gamma__);		// $T = T_{01} \left( \frac{p_2} {p_{01}} \right)^\frac{\gamma_ - 1}{\gamma_}$
	double rho_ = p2_ / (rgas_ * temp_);		// $\rho = \frac{p_2} {r_{gas}T}$
	double mach_ = sqrt(2. * ((t01_ / temp_) - 1.) / gam1_);	// $ M = \sqrt{ 2 \frac{ T_{01} / T - 1}{\gamma_ - 1} } $
	double cs_ = sqrt(gamma__ * p2_ / rho_);		// $c_s = \sqrt{\frac{\gamma_ p_2}{\rho}}$
	double u_ = cs_ * mach_;		// $u = c_s M$
	double mass_ = rho_ * u_ * section;	// $ Q_M = \rho u S_{[1]}$
	double e_ = (cpgas_ - rgas_) * t01_;		// $E = (c_{p gas} - r_{gas}) T_{01}$

	// cusp_s->InverseGrid();
	cusp_s->InitFlow(rho_, mass_, e_, p2_);
	cusp_s->RefreshBoundaries();

	// convtol = 1e-5;
	for (int iter = 0; iter < maxiter && drho > convtol; ++iter)
	{
		drho = cusp_s->Solve();
	}

	cusp_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

//void Laval()
//{
//	cout << endl << "  ---  Laval test  ---  " << endl;
//
//	Loop laval;
//
//	if (kdissip == 'c' || kdissip == 'C') vis4 = 1./vis4;
//
//	laval.ReadFile("Input/grid.txt");
//	//laval.ShowData(5);
//	l_sections = laval.GetSection();
//	l_coords = laval.GetCoordinates();
//	l_volumes.resize(l_coords.size(), 0.);
//
//	imax = laval.col_size;
//	ib2 = imax - 1;
//	ncells = imax - 3;
//	mxdum = imax * 4;
//
//	for (int i = 1; i < ib2; ++i)
//	{
//		//l_volumes[i] = l_sections[i - 1] + 2. * l_sections[i] + l_sections[i + 1];
//		l_volumes[i] = (l_coords[i + 1] - l_coords[i - 1]) / 8. *
//			(l_sections[i - 1] + 2. * l_sections[i] + l_sections[i + 1]);
//	}
//	l_volumes[0] = l_volumes[1];
//	l_volumes[imax - 1] = l_volumes[ib2 - 1];
//	laval.AddColumn("volume", l_volumes);
//	//laval.ShowData(5);
//
//	//  Problem parameters
//	gamma_ = 1.4;			// ratio of specific heats
//	cpgas = 1005.;		// specific heat coeff.at constant pressure [J / kgK]
//	p01 = 1e5;			// inlet total pressure [Pa]
//	t01 = 288.;			// inlet total temperature [K]
//	p2 = 0.37e5;			// outlet static pressure [Pa]
//
//	//  Initial flow parameters
//	gam1 = gamma_ - 1.;	// $\gamma_ - 1$
//	gap1 = gamma_ + 1.;	// $\gamma_ + 1$
//	rgas = gam1 * cpgas / gamma_;	// $r_{gas} = (\gamma_ - 1) \frac{c_p}{\gamma_}$
//	temp = t01 * pow((p2 / p01), gam1 / gamma_);		// $T = T_{01} \left( \frac{p_2} {p_{01}} \right)^\frac{\gamma_ - 1}{\gamma_}$
//	rho = p2 / (rgas * temp);		// $\rho = \frac{p_2} {r_{gas}T}$
//	mach = sqrt(2. * ((t01 / temp) - 1.) / gam1);	// $ M = \sqrt{ 2 \frac{ T_{01} / T - 1}{\gamma_ - 1} } $
//	cs = sqrt(gamma_ * p2 / rho);		// $c_s = \sqrt{\frac{\gamma_ p_2}{\rho}}$
//	u = cs * mach;		// $u = c_s M$
//	mass = rho * u * l_sections[1];	// $ Q_M = \rho u S_{[1]}$
//	e = (cpgas - rgas) * t01;		// $E = (c_{p gas} - r_{gas}) T_{01}$
//
//	dt.resize(imax, 0.);
//	ls.resize(imax, 0.);
//	dum.resize(mxdum, 0.);
//
//	cv.resize(3);
//	cvold.resize(3);
//	diss.resize(3);
//	rhs.resize(3);
//	rs.resize(3);
//	for (int i = 0; i < 3; ++i) {
//		cv[i].resize(imax, 0.);
//		cvold[i].resize(imax, 0.);
//		diss[i].resize(imax, 0.);
//		rhs[i].resize(imax, 0.);
//		rs[i].resize(imax, 0.);
//	}
//
//	p.resize(imax, 0.);
//
//	for (int i = 0; i < imax; ++i)
//	{
//		cv[0][i] = rho * l_sections[i];		//	$\rho S$
//		cv[1][i] = mass;					//	$\rho u S$
//		cv[2][i] = rho * e * l_sections[i];	//	$\rho E S$
//		p[i] = p2;
//	}
//
//	// limiter reference values
//	volref = 1.;			// characteristic length**2
//	rhoref = rho;
//	uref = u;
//	pref = p2;
//
//	drho1 = 0.;
//
//	bcond();
//
//	drho = 1.;
//
//	for (iter = 0; iter < maxiter && drho > convtol; ++iter)
//	{
//		solver(imax, ib2, mxdum, l_coords, l_sections, l_volumes, cv, p, cvold, diss, rhs, dt, ls, rs, dum);
//	}
//
//	Output(imax, ib2, l_coords, l_sections, cv, p);
//	Wsolut(imax, cv, p);
//
//
//	cout << endl << "   OK   " << endl;
//}
//
//void bcond()
//{
//	double cs2;
//	double c02;
//	double rinv;
//	double dis;
//	double cb;
//	double cc02;
//	double tb;
//	double pb;
//	double rhob;
//	double ub;
//
//	u = cv[1][1] / cv[0][1];		///< $u = \frac{(\rho u S)_{[1]}} {(\rho  S)_{[1]}}$
//	cs2 = gamma_ * p[1] * l_sections[1] / cv[0][1];	// Speed of sound $c_s^2 = \frac{\gamma_ p_{[1]} S_{[1]}} { \rho  S_{[1]}}$
//	c02 = cs2 + 0.5 * gam1 * u * u;			// Stagnation speed of sound$c_0^2 = c_s^2 + 0.5 (\gamma_ - 1) u^2$
//	rinv = u - 2. * sqrt(cs2) / gam1;		// Riemann invariant $R^-$, $r_{inv} = u - 2 \frac{\sqrt{c_s^2}}{\gamma_ - 1} $
//	dis = gap1 * c02 / (gam1 * rinv * rinv) - 0.5 * gam1;	// $dis = \frac{(\gamma_ + 1) c_0^2}{(\gamma_ - 1) r_{inv}^2}$
//	if (dis < 0)
//		dis = 1e-20;
//
//	cb = -rinv * (gam1 / gap1) * (1 + sqrt(dis));	// Boundary speed of sound $c_b = -r_{inv} \frac{\gamma_ - 1}{\gamma_ + 1} (1 + \sqrt{dis})$
//	cc02 = min(cb * cb / c02, 1.);		// It is like cb shoud be not greater than 1. // stagnation speed of sound, $cc_0^2 = min\left(\frac{cb^2}{c_0^2}, 2\right)$
//	tb = t01 * cc02;			// $T_b = T_{01} cc_0^2$
//	pb = p01 * pow((tb / t01), gamma_ / gam1);	// $p_b = p_{01} \left(\frac{T_b}{T_{01}}\right)^\frac{\gamma_}{\gamma_ - 1}$
//	rhob = pb / (rgas * tb);		//  $\rho_b = \frac{p_b}{r_{gas} T_b}$
////	cout << "cs2 = " << cs2 << endl;
////	cout << "c02 = " << c02 << endl;
////	cout << "l_sections[1] = " << l_sections[1] << endl;
////	cout << "p[1] = " << p[1] << endl;
////	cout << "cv[0][1] = " << cv[0][1] << endl;
////	cout << "u = " << u << endl;
////	cout << "rinv = " << rinv << endl;
////	cout << "dis = " << dis << endl;
////	cout << "cb = " << cb << endl;
////	cout << "tb = " << tb << endl;
////	cout << "pb = " << pb << endl;
////	cout << "rhob = " << rhob << endl;
//	ub = sqrt(2. * cpgas * (t01 - tb));		// $u_b = \sqrt{2 c_{p, gas} \left(T_{01} - T_b\right)}$
//
//	cv[0][0] = rhob * l_sections[1];		//	$\rho S = \rho_b S_{[1]}$
//	cv[1][0] = rhob * l_sections[1] * ub;					//	$\rho u S = \rho_b S_{[1]} u_b$
//	cv[2][0] = (pb / gam1 + 0.5 * rhob * ub * ub) * l_sections[1];	//	$\rho E S = \left(\frac{p_b}{\gamma_ - 1} + 0.5 \rho_b u_b^2\right) S_{[1]}$
//	p[0] = pb;
//
//	rho = cv[0][ib2 - 1] / l_sections[ib2 - 1];
//	u = cv[1][ib2 - 1] / cv[0][ib2 - 1];
//	cs = sqrt(gamma_ * p[ib2 - 1] / rho);
//
////	cout << cv[0][0] << "\t" << cv[1][0] << "\t" << cv[2][0] << "\t" << p[0] << endl;
//
//	// Now outflow
//	if (u >= cs) {		// supersonic flow
//		pb = p[ib2 - 1];
//		rhob = rho;
//		ub = u;
//	}
//	else {				//subsonic flow
//		pb = p2;
//		rhob = rho + (p2 - p[ib2 - 1]) / (cs * cs);		// $\rho_b = \rho + \frac{(p_2 - p_{[ib2 - 1]})}{c_s^2}$
//		ub = u - (p2 - p[ib2 - 1]) / (cs * rho);		// $u_b = u - \frac{(p_2 - p_{[ib2 - 1]})}{c_s \rho}$
//	}
//
//	cv[0][imax - 1] = rhob * l_sections[ib2 - 1];
//	cv[1][imax - 1] = rhob * ub * l_sections[ib2 - 1];
//	cv[2][imax - 1] = (pb / gam1 + 0.5 * rhob * ub * ub) * l_sections[ib2 - 1];
//	p[imax - 1] = pb;
//
////	cout << cv[0][imax - 1] << "\t" << cv[1][imax - 1] << "\t" << cv[2][imax - 1] << "\t" << p[imax - 1] << endl;
//}
//
//void solver(int imax, int ib2, int mxdum, vector < double >& l_coords, vector < double >& l_sections, vector < double >& l_volumes,
//	vector < vector < double > >& cv, vector < double >& p, vector < vector < double > >& cvold,
//	vector < vector < double > >& diss, vector < vector < double > >& rhs, vector < double >& dt, vector < double >& ls,
//	vector < vector < double > >& rs, vector < double >& dum)
//{
//	int i, irk, idim;
//	double fac, adtv, rrho, rhou, rhoe;
//
//	int nrk = 5; // number of stages
//
//	// Store previous solution; set dissipation = 0
//	for (int i = 0; i < imax; ++i)
//	{
//		cvold[0][i] = cv[0][i];
//		cvold[1][i] = cv[1][i];
//		cvold[2][i] = cv[2][i];
//		diss[0][i] = 0.;
//		diss[1][i] = 0.;
//		diss[2][i] = 0.;
//	}
//
//	// Calculate time step
//
//	Tstep(imax, ib2, l_coords, l_sections, l_volumes, cv, p, dt);
//
//	// loop over the R. - K.stages:
//
//	for (int irk = 0; irk < nrk; irk++)
//	{
//		if (kdissip == 'c') {
//			if (ldiss[irk] == 1)	// artificial dissipation
//				Dissip(imax, ib2, betrk[irk], l_volumes, cv, p, dt, dum.data(), (dum.data() + imax), diss);
//
//			// convective flux
//			Flux_cen(imax, ib2, l_sections, cv, p, diss, dum, rhs);
//		}
//
//		// source term
//
//		Srcterm(imax, ib2, l_sections, p, rhs);
//
//		// residual * time step
//
//		fac = ark[irk] * cfl;
//		for (int i = 1; i < ib2; ++i)
//		{
//			adtv = fac * dt[i] / l_volumes[i];
//			rhs[0][i] = adtv * rhs[0][i];
//			rhs[1][i] = adtv * rhs[1][i];
//			rhs[2][i] = adtv * rhs[2][i];
//		}
//
//		// implicit rediual smoothing
//
//		if (lsmoo[irk] > 0 && epsirs > 0.)
//			Irsmoo(imax, ib2, rhs, dum);
//
//		// update (conservative variables and pressure)
//
//		for (int i = 1; i < ib2; ++i)
//		{
//			cv[0][i] = cvold[0][i] - rhs[0][i];
//			cv[1][i] = cvold[1][i] - rhs[1][i];
//			cv[2][i] = cvold[2][i] - rhs[2][i];
//
//			rrho = l_sections[i] / cv[0][i];
//			rhou = cv[1][i] / l_sections[i];
//			rhoe = cv[2][i] / l_sections[i];
//			p[i] = gam1 * (rhoe - 0.5 * rhou * rhou *rrho);
//		}
//
//		// boundary conditions
//
//		bcond();
//
//	}
//
//	// print out convergence
//
//	Conver(imax, ib2, l_sections, cv, cvold, p);
//}
//
//// Local time stepping
//void Tstep(int imax, int ib2, vector < double >& l_coords, vector < double >& l_sections, vector < double >& l_volumes,
//			vector < vector < double > >&cv, vector < double >& p, vector < double >& dt)
//{
//	double rho;
//	double u;
//	double cs;
//	double dx;
//	double sprad;
//
//	for (int i = 1; i < ib2; ++i)
//	{
//		rho = cv[0][i] / l_sections[i];
//		u = cv[1][i] / cv[0][i];
//		cs = sqrt(gamma_ * p[i] / rho);
//		dx = 0.5 * (l_coords[i + 1] - l_coords[i - 1]);
//		sprad = cs * sqrt(dx * dx + pow(l_sections[i], 2)) + abs(u) * l_sections[i];
//		dt[i] = l_volumes[i] / sprad;
//	}
//	dt[0] = dt[1];
//	dt[imax - 1] = dt[ib2 - 1];
//}
//
///*!
//	Artificial dissipation flux
//*/
//void Dissip(int imax, int ib2, double beta, vector < double >& l_volumes, vector < vector < double > >& cv,
//			vector < double >& p, vector < double >& dt, double* dp, double* d, vector < vector < double > >& diss)
//{
//	int im1, ip1, ip2;
//	double eval, pmax, eps2, eps4, beta1;
//
//	// pressure sensor (divided second differences)
//	for (int i = 1; i < ib2; ++i)
//	{
//		dp[i] = abs((p[i + 1] - 2. * p[i] + p[i - 1]) /
//					(p[i + 1] + 2. * p[i] + p[i - 1]));
//	}
//	dp[0] = dp[1];
//	dp[imax - 1] = dp[ib2 - 1];
//
//	// dissipation fluxes (at i+1/2)
//	for (int i = 0; i < ib2; ++i)
//	{
//		im1 = max(i - 1, 0);
//		ip1 = 	  i + 1;
//		ip2 = min(i + 2, imax - 1);
//		eval = 0.5 * (l_volumes[i] / dt[i] + l_volumes[ip1] / dt[ip1]); // equation 4.56 - scaling factor
//		pmax = max(dp[i], dp[ip1]);
//		eps2 = eval * vis2 * pmax;
//		eps4 = eval * vis4;
//		eps4 = max(0., eps4 - eps2);
//		d[0 * imax + i] = eps2 * (cv[0][ip1] - cv[0][i]) -
//				  eps4 * (cv[0][ip2] - 3. * cv[0][ip1] + 3. * cv[0][i] - cv[0][im1]);
//
//		d[1 * imax + i] = eps2 * (cv[1][ip1] - cv[1][i]) -
//				  eps4 * (cv[1][ip2] - 3. * cv[1][ip1] + 3. * cv[1][i] - cv[1][im1]);
//
//		d[2 * imax + i] = eps2 * (cv[2][ip1] - cv[2][i]) -
//				  eps4 * (cv[2][ip2] - 3. * cv[2][ip1] + 3. * cv[2][i] - cv[2][im1]);
//	}
//
//	// dissipation term
//
////	cout << "l_volumes[1] = " << l_volumes[1] << endl;
////	cout << "dt[1] = " << dt[1] << endl;
////	cout << "eval = " << eval << endl;
////	cout << "eps2 = " << eps2 << endl;
////	cout << "eps4 = " << eps4 << endl;
////	cout << d[0 * imax + ib2 - 1] << "\t" << d[0 * imax + ib2 - 2] << endl;
////	cout << d[0 * imax + 1] << "\t" << d[0 * imax + 0] << endl;
////	cout << "---\t" << 	d[0 * imax + 0] << "\t" <<
////						cv[0][4] - 3. * cv[0][3] + 3. * cv[0][2] - cv[0][1] << "\t" <<
////						cv[0][1] - cv[0][0] << "\t" << cv[0][2] - cv[0][1] << endl;
//
//	beta1 = 1. - beta;
//	for (int i = 1; i < ib2; ++i)
//	{
//		diss[0][i] = beta * (d[0 * imax + i] - d[0 * imax + i - 1]) + beta1 * diss[0][i];
//		diss[1][i] = beta * (d[1 * imax + i] - d[1 * imax + i - 1]) + beta1 * diss[1][i];
//		diss[2][i] = beta * (d[2 * imax + i] - d[2 * imax + i - 1]) + beta1 * diss[2][i];
//	}
////	cout << "beta = " << beta << endl;
////	cout << "diss = " << diss[0][1] << endl;
//}
//
//void Flux_cen(int imax, int ib2, vector < double >& l_sections, vector < vector < double > >& cv,
//				vector < double >& p, vector < vector < double > >& diss, vector < double >& f,
//				vector < vector < double > >& rhs)
//{
//	double si, rav, ruav, reav, pav, rhav, qs;
//
//	// flux term (average of variables)
//
//	for (int i = 0; i < ib2; ++i)
//	{
//		si = 0.5 * (l_sections[i + 1] + l_sections[i]);
//		rav = 0.5 * (cv[0][i + 1] / l_sections[i + 1] + cv[0][i] / l_sections[i]);
//		ruav = 0.5 * (cv[1][i + 1] / l_sections[i + 1] + cv[1][i] / l_sections[i]);
//		reav = 0.5 * (cv[2][i + 1] / l_sections[i + 1] + cv[2][i] / l_sections[i]);
//		pav = 0.5 * (p[i + 1] + p[i]);
//		rhav = reav + pav;
//		qs = ruav * si / rav;
//
//		f[0 * imax + i] = rav * qs;
//		f[1 * imax + i] = ruav * qs + pav * si;
//		f[2 * imax + i] = rhav * qs;
//	}
//
//	// flux + dissipation = RHS
//
//	for (int i = 1; i < ib2; ++i)
//	{
//		rhs[0][i] = f[0 * imax + i] - f[0 * imax + i - 1] - diss[0][i];
//		rhs[1][i] = f[1 * imax + i] - f[1 * imax + i - 1] - diss[1][i];
//		rhs[2][i] = f[2 * imax + i] - f[2 * imax + i - 1] - diss[2][i];
//	}
//}
//
//// Didn't understand
//
//void Srcterm(int imax, int ib2, vector < double >& l_sections, vector < double >& p,
//				vector < vector < double > >& rhs)
//{
//	double da;
//
//	for (int i = 1; i < ib2; ++i)
//	{
//		da = 0.5 * (l_sections[i + 1] - l_sections[i - 1]);
//		rhs[1][i] = rhs[1][i] - p[i] * da;
//	}
//}
//
//void Irsmoo(int imax, int ib2, vector < vector < double > >& rhs, vector < double >& d)
//{
//	double eps2, t;
//	int i;
//
//	eps2 = 2. * epsirs + 1.;
//	d[0] = 0.;
//	rhs[0][0] = 0.;
//	rhs[1][0] = 0.;
//	rhs[2][0] = 0.;
//
//	// elimination step
//
//	for (int i = 1; i < ib2; ++i)
//	{
//		t = 1. / (eps2 - epsirs * d[i - 1]);
//		d[i] = t * epsirs;
//		rhs[0][i] = t * (rhs[0][i] + epsirs * rhs[0][i - 1]);
//		rhs[1][i] = t * (rhs[1][i] + epsirs * rhs[1][i - 1]);
//		rhs[2][i] = t * (rhs[2][i] + epsirs * rhs[2][i - 1]);
//	}
//
//	// backward substitution
//
//	i = ib2 - 1;
//	for (int ii = 2; ii < ib2; ++ii)
//	{
//		i = i - 1;
//		rhs[0][i] = rhs[0][i] + d[i] * rhs[0][i + 1];
//		rhs[1][i] = rhs[1][i] + d[i] * rhs[1][i + 1];
//		rhs[2][i] = rhs[2][i] + d[i] * rhs[2][i + 1];
//	}
//}
//
//void Conver(int max, int ib2, vector < double >& l_sections, vector < vector < double > >& cv,
//			vector < vector < double > >& cvold, vector < double >& p)
//{
//	int idrho, nsup;
//	double dr, drmax, rho, u, c, avms;
//
//	// get the change of density and the mass flow
//
//	drho = 0.;
//	drmax = -1e20;
//	avms = 0.;
//	nsup = 0;
//
//	for (int i = 1; i < ib2; ++i)
//	{
//		dr = cv[0][i] - cvold[0][i];
//		avms = avms + cv[1][i];
//		drho = drho + dr * dr;
//		if (abs(dr) >= drmax) {
//			drmax = abs(dr);
//			idrho = i;
//		}
//		rho = cv[0][i] / l_sections[i];
//		u = cv[1][i] / cv[0][i];
//		c = sqrt(gamma_* p[i] / rho);
//		if (u > c) nsup++;
//	}
//	avms = avms / double(ncells + 1);
//
//	// print convergence history
//
//	if (iter == 0) drho1 = sqrt(drho / double(ncells + 1));
//
//	drho = sqrt(drho / double(ncells + 1)) / drho1;
//
//	cout << "conver" << "\t";
//	cout << iter << "\t";
//	cout << log10(drho) << "\t";
//	cout << drmax << "\t";
//	cout << idrho << "\t";
//	cout << avms << "\t";
//	cout << cv[1][1] - cv[1][ib2 - 1] << "\t";
//	cout << nsup << "\n";
//
//}
//
//void Output(int imax, int ib2, vector < double >& l_coords, vector < double >& l_sections,
//			vector < vector < double > >& cv, vector < double >& p)
//{
//	double rho, u, temp, c, mach;
//
//	FILE* file;
//
//	file = fopen("Output/output.txt", "w");
//	fprintf(file, "x\tA\trho\tu\tp\tT\tM\tmass_flow\n");
//	for (int i = 1; i < ib2; ++i)
//	{
//		rho = cv[0][i] / l_sections[i];
//		u = cv[1][i] / cv[0][i];
//		temp = p[i] / (rgas * rho);
//		c = sqrt(gamma_ * p[i] / rho);
//		mach = u / c;
//		fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
//			l_coords[i], l_sections[i], rho, u, p[i], temp, mach, cv[1][i]);
//	}
//
//}
//
//void Wsolut(int imax, vector < vector < double > >& cv, vector < double >& p)
//{
//
//}
