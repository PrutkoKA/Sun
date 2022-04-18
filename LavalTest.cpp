#include "LavalTest.h"
#include <ctime>

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
	//central_s->HideOutput();
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
	central_s->RhoUPH();
	central_s->RefreshBoundaries();

	// convtol = 1e-5;
	for (int iter = 0; iter < maxiter && drho > convtol; ++iter)
	{
		drho = central_s->SolveExplImpl();
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
		drho = central_s->SolveExplImpl();
	}

	central_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

void Laval4()	// Examining non-uniform grid (Central scheme)
{
	cout << endl << "  ---  Laval4 test  ---  " << endl;

	//string file_name = "Input/solver.yml";
	string file_name = "Input/condensed.yml";

	Loop laval;
	Solver *central_s;
	const int range_size = 73;
	vector < double > range_(range_size, 0.);
	vector < vector < double > > new_tab;

	central_s = CreateReadConfigFile(file_name);

	//laval.ReadFile("Input/grid.txt");
	laval.ReadFile("Input/condensed_grid.txt");
	
	double sum = 0.;
	double fac = 1.04;
	double h = 0.; // 0.0115;
	
	const int middle = (range_size + 3) / 2;
	double sum_ = 0.;
	for (int i = 2; i < middle; ++i) {
		sum_ += pow(fac, double(middle - 1 - i));
	}
	h = 0.64 / sum_;
	range_[range_size - 2] = 1.;
	range_[range_size - 1] = 1.;
	
	for (int i = 2; i < middle; ++i) {
		range_[i] = range_[i-1] + h * pow(fac, double(middle - 1 - i));
		sum += range_[i] - range_[i-1];
		cout << i << "\t" << range_[i] << "\t" << sum << "\t" << pow(fac, double(middle - 1 - i)) << endl;
	}
	sum_ = 0.;
	for (int i = middle; i < range_size - 1; ++i) {
		sum_ += pow(fac, double(i - middle));
	}
	h = 0.36 / sum_;
	for (int i = middle; i < range_size - 2; ++i) {
		range_[i] = range_[i-1] + h * pow(fac, double(i - middle));
		sum += range_[i] - range_[i-1];
		cout << i << "\t" << range_[i] << "\t" << sum << "\t" << pow(fac, double(i - middle)) << endl;
	}

	new_tab = laval.NewTable("coordinate", range_);
	laval.SetData(new_tab);
//	laval.ShowData();

	central_s->SetGrid(laval);
	central_s->ReadBoundaries("Input/boundary.yml");
	central_s->SetOutputFile("Output/condensed_output.txt");
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
	central_s->RhoUPH();
	central_s->RefreshBoundaries();

	central_s->calculate_mass_matrix();
	central_s->fill_inverse_mass_matrix();
	//central_s->print_inversed_mass_matrix();

	// convtol = 1e-5;
	for (int iter = 0; iter < maxiter && drho > convtol; ++iter)
	{
		drho = central_s->SolveExplImpl();
	}

	central_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

void Laval6()	// Standard example as in Blazek (CUSP)
{
	cout << endl << "  ---  Laval6 test  ---  " << endl;

	string file_name = "Input/sod.yml";

	Loop laval;
	Solver *cusp_s;

	cusp_s = CreateReadConfigFile(file_name);

	laval.ReadFile("Input/sod.txt");


	/*vector < double > range_;
	vector < vector < double > > new_tab;
	range_ = laval.MakeRange(0., 1., 0.002);
	range_.push_back(range_[range_.size() - 1]);
	range_.insert(range_.begin(), range_.begin(), range_.begin() + 1);
	new_tab = laval.NewTable("coordinate", range_);
	laval.SetData(new_tab);*/

	cusp_s->SetGrid(laval);
	cusp_s->ReadBoundaries("Input/sod_boundary.yml");
	cusp_s->SetOutputFile("Output/output5.txt");
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
	//double gap1_ = gamma__ + 1.;	// $\gamma_ + 1$
	//double rgas_ = gam1_ * cpgas_ / gamma__;	// $r_{gas} = (\gamma_ - 1) \frac{c_p}{\gamma_}$
	//double temp_ = t01_ * pow((p2_ / p01_), gam1_ / gamma__);		// $T = T_{01} \left( \frac{p_2} {p_{01}} \right)^\frac{\gamma_ - 1}{\gamma_}$
	//double rho_ = p2_ / (rgas_ * temp_);		// $\rho = \frac{p_2} {r_{gas}T}$
	//double mach_ = sqrt(2. * ((t01_ / temp_) - 1.) / gam1_);	// $ M = \sqrt{ 2 \frac{ T_{01} / T - 1}{\gamma_ - 1} } $
	//double cs_ = sqrt(gamma__ * p2_ / rho_);		// $c_s = \sqrt{\frac{\gamma_ p_2}{\rho}}$
	//double u_ = cs_ * mach_;		// $u = c_s M$
	//double mass_ = rho_ * u_ * section;	// $ Q_M = \rho u S_{[1]}$
	//double e_ = (cpgas_ - rgas_) * t01_;		// $E = (c_{p gas} - r_{gas}) T_{01}$

	double rho_[2] = { 1, 0.125 };
	double u_[2] = { 0.75, 0 };
	double mass_[2];
	double p_[2] = { 1., 0.1 };
	double e_[2];

	for (int i = 0; i < 2; ++i) {
		e_[i] = p_[i] / gam1_ / rho_[i] + 0.5 * u_[i] * u_[i];
		mass_[i] = rho_[0] * u_[i];
	}

	// cusp_s->InverseGrid();
	cusp_s->InitFlow(rho_, mass_, e_, p_, 0.3);
	cusp_s->RhoUPH();
	cusp_s->RefreshBoundaries();

	double physDt_ = 0.2e-2;

	// convtol = 1e-5;
	//cusp_s->cvn = cusp_s->cv;
	cusp_s->cvnm1 = cusp_s->cv;
	//cusp_s->ForwardEuler(physDt_);
	cusp_s->iter = 0.;
	cusp_s->Global_Time = 0.;
	//for (int iter = 0; iter < maxiter && drho > convtol && false; ++iter)
	//{
	//	drho = cusp_s->Solve();
	//	//drho = cusp_s->SolveImplicit();
	//}
	cusp_s->cvn = cusp_s->cv;

	cusp_s->S.resize((cusp_s->ib2 - 1) * cusp_s->eq_num, (cusp_s->ib2 - 1) * cusp_s->eq_num);
	cusp_s->S.reserve(VectorXi::Constant((cusp_s->ib2 - 1) * cusp_s->eq_num, cusp_s->eq_num * cusp_s->eq_num));

	cusp_s->CalculateTimeSource(cusp_s->cvn, cusp_s->cvnm1, physDt_);
	cusp_s->iter = 0.;
	int iter = 0;
	double ttime = 0.;
	if (cusp_s->time_stepping == 0) {
		while (ttime < 0.2 / 1.) {
			//while (cusp_s->Global_Time < 0.2*50) {
			cusp_s->iter = 0;
			drho = 1.;
			for (int iter = 0; iter < maxiter && drho > convtol && true; ++iter)
			{
				drho = cusp_s->Solve(physDt_);
				//drho = cusp_s->SolveImplicit();
			}
			cusp_s->cvnm1 = cusp_s->cvn;
			cusp_s->cvn = cusp_s->cv;
			//if (cusp_s->iter == 1)
			//	physDt_ = cusp_s->Global_Time / cusp_s->GetCFL();

			//cout << cusp_s->Global_Time << endl;
			ttime += physDt_;
			iter += cusp_s->iter;
			cout << iter << endl;
		}
	}
	if (cusp_s->time_stepping == 1) {
		//cusp_s->Global_Time = 1e-20;
		physDt_ = 0.2e-0;
		//while (ttime < 0.2e-4) {
		while (cusp_s->Global_Time < 0.2*60) {
			//cusp_s->iter = 0;
			//drho = 1.;
			//for (int iter = 0; iter < maxiter && drho > convtol && true; ++iter)
			//{
				drho = cusp_s->Solve(physDt_);
				if (drho < -400.) {
					return;
				}
				//drho = cusp_s->SolveImplicit();
			//}
			cusp_s->cvnm1 = cusp_s->cvn;
			cusp_s->cvn = cusp_s->cv;
			if (cusp_s->iter == 1)
				physDt_ = cusp_s->Global_Time / cusp_s->GetCFL();

			cout << cusp_s->Global_Time << endl;
			//ttime += physDt_;
		}
	}
	
	cusp_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

void Laval5()	// Standard example as in Blazek (CUSP)
{
	cout << endl << "  ---  Laval5 test  ---  " << endl;

	string file_name = "Input/cusp.yml";

	Loop laval;
	Solver* cusp_s;

	cusp_s = CreateReadConfigFile(file_name);

	laval.ReadFile("Input/grid.txt");


	/*vector < double > range_;
	vector < vector < double > > new_tab;
	range_ = laval.MakeRange(0., 1., 0.002);
	range_.push_back(range_[range_.size() - 1]);
	range_.insert(range_.begin(), range_.begin(), range_.begin() + 1);
	new_tab = laval.NewTable("coordinate", range_);
	laval.SetData(new_tab);*/

	cusp_s->SetGrid(laval);
	cusp_s->ReadBoundaries("Input/boundary.yml");
	cusp_s->SetOutputFile("Output/output4.txt");
	// cusp_s->HideOutput();
	// cusp_s->InverseGrid();

	// Euler euqations
	cusp_s->SetEquation("mass", { "Rho", "*A" }, { "Rho", "*U", "*A" }, { "" }, cusp_s->vars, cusp_s->vars_o);	// RhoA, RhoUA
	cusp_s->SetEquation("impulse", { "Rho", "*U", "*A" }, { "Rho", "*U^2", "+p", "*A" }, { "-p", "*dA" }, cusp_s->vars, cusp_s->vars_o);		// RhoUA, (RhoUU+p)A
	cusp_s->SetEquation("energy", { "Rho", "*E", "*A" }, { "Rho", "*E", "+p", "*U", "*A" }, { "" }, cusp_s->vars, cusp_s->vars_o);		// RhoEA, (RhoEU+pU)A

	cusp_s->set_fv_equation(		// Rho = RhoA / A
		"Rho",
		{ "RhoA", "/A", "" }		// There is dummy for unambiguous conservation
	);
	cusp_s->set_fv_equation(		// E = RhoEA / RhoA
		"E",
		{ "RhoEA", "/RhoA", "" }		// There is dummy for unambiguous conservation
	);
	cusp_s->set_fv_equation(		// U = RhoUA / RhoA
		"U",
		{ "RhoUA", "/RhoA", "" }
	);
	cusp_s->set_fv_equation(		// p = (RhoEA / RhoA - 0.5U^2) * (gamma - 1) * RHO
		"p",
		{ "RhoEA", "/RhoA", "-0.5U^2", "*GAMMAM", "*Rho" }
	);
	cusp_s->set_fv_equation(		// H = gamma / (gamma - 1) * p / RHO + 0.5U^2
		"H",
		{ "GAMMA", "/GAMMAM", "*p", "/Rho", "+0.5U^2" }
	);
	cusp_s->set_fv_equation(		// A
		"A",
		{ "A" }
	);
	cusp_s->set_fv_equation(		// T = p * gamma / (gamma - 1) / Cp / rho
		"T",
		{ "p", "*GAMMA", "/GAMMAM", "/CP", "/Rho" }
	);

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
	cusp_s->RhoUPH();
	cusp_s->RefreshBoundaries();

	double physDt_ = 1e-1 * 0;

	// convtol = 1e-5;
	cusp_s->cvn = cusp_s->cv;
	cusp_s->cvnm1 = cusp_s->cv;

	cusp_s->S.resize((cusp_s->ib2 - 1) * cusp_s->eq_num, (cusp_s->ib2 - 1) * cusp_s->eq_num);
	cusp_s->S.reserve(VectorXi::Constant((cusp_s->ib2 - 1) * cusp_s->eq_num, cusp_s->eq_num * cusp_s->eq_num));

	cusp_s->CalculateTimeSource(cusp_s->cvn, cusp_s->cvnm1, physDt_);
	for (int iter = 0; iter < maxiter && drho > convtol; ++iter)
	{
		drho = cusp_s->Solve(physDt_);
		//drho = cusp_s->SolveImplicit();
	}
	cusp_s->cvnm1 = cusp_s->cvn;
	cusp_s->cvn = cusp_s->cv;

	cusp_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

void Laval7()	// Standard example as in Blazek (CUSP)
{
	cout << endl << "  ---  Laval6 test  ---  " << endl;

	string file_name = "Input/sod_explicit.yml";

	Loop laval;
	Solver* cusp_s;

	cusp_s = CreateReadConfigFile(file_name);

	laval.ReadFile("Input/sod.txt");


	/*vector < double > range_;
	vector < vector < double > > new_tab;
	range_ = laval.MakeRange(0., 1., 0.002);
	range_.push_back(range_[range_.size() - 1]);
	range_.insert(range_.begin(), range_.begin(), range_.begin() + 1);
	new_tab = laval.NewTable("coordinate", range_);
	laval.SetData(new_tab);*/

	cusp_s->SetGrid(laval);
	cusp_s->ReadBoundaries("Input/sod_boundary.yml");
	cusp_s->SetOutputFile("Output/sod_explicit_adaptive_grid.txt");
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
	//double gap1_ = gamma__ + 1.;	// $\gamma_ + 1$
	//double rgas_ = gam1_ * cpgas_ / gamma__;	// $r_{gas} = (\gamma_ - 1) \frac{c_p}{\gamma_}$
	//double temp_ = t01_ * pow((p2_ / p01_), gam1_ / gamma__);		// $T = T_{01} \left( \frac{p_2} {p_{01}} \right)^\frac{\gamma_ - 1}{\gamma_}$
	//double rho_ = p2_ / (rgas_ * temp_);		// $\rho = \frac{p_2} {r_{gas}T}$
	//double mach_ = sqrt(2. * ((t01_ / temp_) - 1.) / gam1_);	// $ M = \sqrt{ 2 \frac{ T_{01} / T - 1}{\gamma_ - 1} } $
	//double cs_ = sqrt(gamma__ * p2_ / rho_);		// $c_s = \sqrt{\frac{\gamma_ p_2}{\rho}}$
	//double u_ = cs_ * mach_;		// $u = c_s M$
	//double mass_ = rho_ * u_ * section;	// $ Q_M = \rho u S_{[1]}$
	//double e_ = (cpgas_ - rgas_) * t01_;		// $E = (c_{p gas} - r_{gas}) T_{01}$

	double rho_[2] = { 1, 0.125 };
	double u_[2] = { 0.75, 0. };
	//double rho_[2] = { 1, 1. };
	//double u_[2] = { -19.59745, -19.59745 };
	double mass_[2];
	double p_[2] = { 1., 0.1 };
	//double p_[2] = { 1000., 0.01 };
	double e_[2];

	for (int i = 0; i < 2; ++i) {
		e_[i] = p_[i] / gam1_ / rho_[i] + 0.5 * u_[i] * u_[i];
		mass_[i] = rho_[0] * u_[i];
	}

	// cusp_s->InverseGrid();
	double shock_pos = 0.3;
	cusp_s->InitFlowAG(rho_, mass_, e_, p_, shock_pos);

	FILE* file;
	file = fopen("Output/RemeshedCoords.txt", "w");
	fclose(file);
	cusp_s->grid.PrintWholeRow("Output/RemeshedCoords.txt", cusp_s->grid.GetCoordinates());

	// Parameters to adjust mesh were used
	for (int i = 0; i < 100 && cusp_s->remesh; ++i) {
		double relax_coef = i < 300 ? 0.5
			: i < 500 ? 0.2
			: i < 1700 ? 0.05
			: i < 2000 ? 0.05 
			: i < 5000 ? 0.02 : 0.005;
		cusp_s->AdjustMesh(rho_, mass_, e_, p_, shock_pos, 1.);
		cusp_s->grid.PrintWholeRow("Output/RemeshedCoords.txt", cusp_s->grid.GetCoordinates());
	}
	cusp_s->InitFlowAG(rho_, mass_, e_, p_, shock_pos);
	cusp_s->CalculateVolumes();
	cusp_s->grid.SetRow("volume", cusp_s->vol);
	//cusp_s->grid.RefreshR("R");
	//cusp_s->grid.RefreshN("n");
	cusp_s->grid.PrintTable("Output/RemeshedGrid.txt");
	//return;

	cusp_s->RhoUPH();
	cusp_s->RefreshBoundaries();

	double fac = 1. / 1e1;
	double physDt_ = 0.2e-2 * fac / 1.;
	//double physDt_ = 1e-10;

	// convtol = 1e-5;
	//cusp_s->cvn = cusp_s->cv;
	cusp_s->cvnm1 = cusp_s->cv;
	//cusp_s->ForwardEuler(physDt_);
	cusp_s->iter = 0.;
	cusp_s->Global_Time = 0.;
	//for (int iter = 0; iter < maxiter && drho > convtol && false; ++iter)
	//{
	//	drho = cusp_s->Solve();
	//	//drho = cusp_s->SolveImplicit();
	//}
	cusp_s->cvn = cusp_s->cv;

	// Must make mass matrix
	// Must remesh cvn, cvnm1, cvold for appropriate dual time stepping
	// Try relaxation (with neighbours) of conc method (several cycles)
	// Mooving mesh accounting?

	cusp_s->S.resize((cusp_s->ib2 - 1) * cusp_s->eq_num, (cusp_s->ib2 - 1) * cusp_s->eq_num);
	cusp_s->S.reserve(VectorXi::Constant((cusp_s->ib2 - 1) * cusp_s->eq_num, cusp_s->eq_num * cusp_s->eq_num));

	//cusp_s->CalculateTimeSource(cusp_s->cvn, cusp_s->cvnm1, physDt_);
	cusp_s->grid.AddColumn("old_coords", cusp_s->grid.GetValues("coordinate"));
	cusp_s->iter = 0.;
	int iter = 0;
	double ttime = 0.;
	if (cusp_s->time_stepping == 0) {
		//while (ttime < 0.2 / 1. * 1.) {
		while (ttime < 0.2e-2 * fac * 1000. || ttime < 0.075) {
			//while (cusp_s->Global_Time < 0.2*50) {
			cusp_s->iter = 0;
			drho = 1.;
			cusp_s->calculate_mass_matrix();
			cusp_s->fill_inverse_mass_matrix();

			vector <double> very_old_coords = cusp_s->grid.GetCoordinates();
			vector <vector <vector <double> > > cvss;
			cvss.push_back(vector <vector <double> >());
			for (int var = 0; var < cusp_s->CONS_VAR_COUNT; ++var)
				cvss[0].push_back(cusp_s->cvn[var]);
			cvss.push_back(vector <vector <double> >());
			for (int var = 0; var < cusp_s->CONS_VAR_COUNT; ++var)
				cvss[1].push_back(cusp_s->cvnm1[var]);

			cusp_s->grid.SetRow("old_coords", cusp_s->grid.GetValues("coordinate"));
			cusp_s->cvn_old = cusp_s->cvn;
			cusp_s->cvnm1_old = cusp_s->cvnm1;


			for (int iter = 0; iter < maxiter && drho > convtol && true; ++iter)
			{
				drho = cusp_s->Solve(physDt_);
			}

			vector <double> old_coords = cusp_s->grid.GetCoordinates();

			auto remesh = [&](int count, bool use_func, double tau) -> void {
				if (cusp_s->remesh && true) {

					// New adaptive grid
					for (int i = 0; i < count; ++i) {
						for (int var = 0; var < cusp_s->CONS_VAR_COUNT; ++var)
							cusp_s->grid.SetRow(cusp_s->c_var_name[var], cusp_s->cv[var]);

						vector< vector <double> > functions;;
						vector< double > Fs;

						if (use_func)
						{
							functions.push_back(cusp_s->cvn[cusp_s->RHO_A]);
							Fs.push_back(1.);
						}

						cusp_s->grid.CalculateResolution(1., 1., cusp_s->c_var_name[cusp_s->RHO_A], "coordinate", functions, Fs);
						cusp_s->grid.CalculateConcentration(1., "coordinate");

						for (int var = 0; var < cusp_s->CONS_VAR_COUNT; ++var)
							cusp_s->grid.SetRow(cusp_s->c_var_name[var], cusp_s->cv[var]);
						cusp_s->grid.SetRow("coordinate", old_coords);

						cusp_s->x = cusp_s->grid.RefineMesh(physDt_, tau);

						for (int var = 0; var < cusp_s->CONS_VAR_COUNT; ++var)
							cusp_s->cv[var] = cusp_s->grid.GetValues(cusp_s->c_var_name[var]);

						cusp_s->RefreshBoundaries();						// Refresh boundary conditions
						cusp_s->RhoUPH();

						cusp_s->CalculateVolumes();
						cusp_s->grid.SetRow("volume", cusp_s->vol);
					}
					// Now we will reevaluate cvn, cvnm1, cvold
					vector < vector < vector < double > >* > cvs{ &cusp_s->cvn, &cusp_s->cvnm1/*, &cusp_s->cvold*/ };

					if (count > 0)
					{
						int i = 0;
						for (auto it : cvs)
						{
							cusp_s->grid.SetRow("coordinate", very_old_coords);
							for (int var = 0; var < cusp_s->CONS_VAR_COUNT; ++var)
								cusp_s->grid.SetRow(cusp_s->c_var_name[var], cvss[i][0]);
							++i;

							vector < string > ignore;
							vector < vector < double > > new_tab;
							ignore.push_back(cusp_s->grid.TYPE_COL);
							new_tab = cusp_s->grid.NewTable("coordinate", cusp_s->x, ignore, false);
							cusp_s->grid.SetData(new_tab);

							for (int var = 0; var < cusp_s->CONS_VAR_COUNT; ++var)
								(*it)[var] = cusp_s->grid.GetValues(cusp_s->c_var_name[var]);
						}

						cusp_s->grid.SetRow("coordinate", cusp_s->x);
						for (int var = 0; var < cusp_s->CONS_VAR_COUNT; ++var)
							cusp_s->grid.SetRow(cusp_s->c_var_name[var], cusp_s->cv[var]);

						// Checking cvn transformation
						/*cusp_s->cv[cusp_s->RHO_A] = cusp_s->cvn[cusp_s->RHO_A];
						cusp_s->cv[cusp_s->RHO_U_A] = cusp_s->cvn[cusp_s->RHO_U_A];
						cusp_s->cv[cusp_s->RHO_E_A] = cusp_s->cvn[cusp_s->RHO_E_A];
						cusp_s->RhoUPH();*/

						cusp_s->calculate_mass_matrix();
						cusp_s->fill_inverse_mass_matrix();
					}
				}
			};

			//remesh(0, false, 1e-1);
			remesh(0, false, 1e15);

			cusp_s->cvnm1 = cusp_s->cvn;
			cusp_s->cvn = cusp_s->cv;
			//if (cusp_s->iter == 1)
			//	physDt_ = cusp_s->Global_Time / cusp_s->GetCFL();

			//cout << cusp_s->Global_Time << endl;
			ttime += physDt_;
			iter += cusp_s->iter;
			cout << iter << "\t" << ttime << endl;

			//if (!cusp_s->steadiness)
				//cusp_s->PrintResult();
		}
	}
	if (cusp_s->time_stepping == 1) {
		//cusp_s->Global_Time = 1e-20;
		physDt_ = 0.2e-0;
		//while (ttime < 0.2e-4) {
		while (cusp_s->Global_Time < 0.2 * 60) {
			//cusp_s->iter = 0;
			//drho = 1.;
			//for (int iter = 0; iter < maxiter && drho > convtol && true; ++iter)
			//{
			drho = cusp_s->Solve(physDt_);
			if (drho < -400.) {
				return;
			}
			//drho = cusp_s->SolveImplicit();
		//}
			cusp_s->cvnm1 = cusp_s->cvn;
			cusp_s->cvn = cusp_s->cv;
			if (cusp_s->iter == 1)
				physDt_ = cusp_s->Global_Time / cusp_s->GetCFL();

			cout << cusp_s->Global_Time << endl;
			//ttime += physDt_;
		}
	}

	cusp_s->PrintResult();

	cout << endl << "   OK   " << endl;
}

void unsteady_sod_test(const string &output_file, const string& yml_file, const double end_time_)	// Standard example as in Blazek (CUSP)
{
	cout << endl << "  ---  unsteady sod test  ---  " << endl;

	string file_name = yml_file;

	Loop sod;
	Solver* hllc_s;

	hllc_s = CreateReadConfigFile(file_name);

	hllc_s->HideOutput();

	sod.ReadFile("Input/sod_mesh_test.txt");

	hllc_s->SetGrid(sod);
	hllc_s->ReadBoundaries("Input/sod_boundary_test.yml");
	hllc_s->SetOutputFile(output_file);

	// Euler euqations
	hllc_s->SetEquation("mass", { "Rho", "*A" }, { "Rho", "*U", "*A" }, { "" }, hllc_s->vars, hllc_s->vars_o);	// RhoA, RhoUA
	hllc_s->SetEquation("impulse", { "Rho", "*U", "*A" }, { "Rho", "*U^2", "+p", "*A" }, { "-p", "*dA" }, hllc_s->vars, hllc_s->vars_o);		// RhoUA, (RhoUU+p)A
	hllc_s->SetEquation("energy", { "Rho", "*E", "*A" }, { "Rho", "*E", "+p", "*U", "*A" }, { "" }, hllc_s->vars, hllc_s->vars_o);		// RhoEA, (RhoEU+pU)A

	hllc_s->set_fv_equation(		// Rho = RhoA / A
		"Rho",
		{ "RhoA", "/A", "" }		// There is dummy for unambiguous conservation
	);
	hllc_s->set_fv_equation(		// E = RhoEA / RhoA
		"E",
		{ "RhoEA", "/RhoA", "" }		// There is dummy for unambiguous conservation
	);
	hllc_s->set_fv_equation(		// U = RhoUA / RhoA
		"U",
		{ "RhoUA", "/RhoA", "" }
	);
	hllc_s->set_fv_equation(		// p = (RhoEA / RhoA - 0.5U^2) * (gamma - 1) * RHO
		"p",
		{ "RhoEA", "/RhoA", "-0.5U^2", "*GAMMAM", "*Rho" }
	);
	hllc_s->set_fv_equation(		// H = gamma / (gamma - 1) * p / RHO + 0.5U^2
		"H",
		{ "GAMMA", "/GAMMAM", "*p", "/Rho", "+0.5U^2" }
	);
	hllc_s->set_fv_equation(		// A
		"A",
		{ "A" }
	);
	hllc_s->set_fv_equation(		// T = p * gamma / (gamma - 1) / Cp / rho
		"T",
		{ "p", "*GAMMA", "/GAMMAM", "/CP", "/Rho" }
	);

	double convtol = hllc_s->GetTolerance();
	double drho = 1.;

	int maxiter = hllc_s->GetMaxIterNum();
	double gamma__ = hllc_s->GetGamma();
	double cpgas_ = hllc_s->GetCp();
	double t01_ = hllc_s->GetInflowT();
	double p01_ = hllc_s->GetInflowP();
	double p2_ = hllc_s->GetOutflowP();
	int id = hllc_s->GetInflowId();
	double section = hllc_s->GetSection()[(id > 0 ? id - 2 : id + 1)];

	double gam1_ = gamma__ - 1.;	// $\gamma_ - 1$

	double rho_[2] = { 1, 0.125 };
	double u_[2] = { 0.75, 0. };

	double mass_[2];
	double p_[2] = { 1., 0.1 };
	double e_[2];

	for (int i = 0; i < 2; ++i) {
		e_[i] = p_[i] / gam1_ / rho_[i] + 0.5 * u_[i] * u_[i];
		mass_[i] = rho_[0] * u_[i];
	}

	double shock_pos = 0.3;
	hllc_s->InitFlowAG(rho_, mass_, e_, p_, shock_pos);

	// Parameters to adjust mesh were used
	for (int i = 0; i < 100 && hllc_s->remesh; ++i) {
		hllc_s->AdjustMesh(rho_, mass_, e_, p_, shock_pos, 1.);
	}
	hllc_s->InitFlowAG(rho_, mass_, e_, p_, shock_pos);
	hllc_s->CalculateVolumes();
	hllc_s->grid.SetRow("volume", hllc_s->vol);

	hllc_s->RhoUPH();
	hllc_s->RefreshBoundaries();

	double fac = 1. / 1e1;
	double physDt_ = 0.2e-2 * fac / 1.;

	hllc_s->cvnm1 = hllc_s->cv;
	hllc_s->iter = 0.;
	hllc_s->Global_Time = 0.;
	hllc_s->cvn = hllc_s->cv;

	hllc_s->S.resize((hllc_s->ib2 - 1) * hllc_s->eq_num, (hllc_s->ib2 - 1) * hllc_s->eq_num);
	hllc_s->S.reserve(VectorXi::Constant((hllc_s->ib2 - 1) * hllc_s->eq_num, hllc_s->eq_num * hllc_s->eq_num));

	hllc_s->grid.AddColumn("old_coords", hllc_s->grid.GetValues("coordinate"));
	hllc_s->iter = 0.;

	int iter = 0;
	double ttime = 0.;

	double start_time = clock();
	if (hllc_s->time_stepping == 0) {
		while (ttime < end_time_) {
			hllc_s->iter = 0;
			drho = 1.;
			hllc_s->calculate_mass_matrix();
			hllc_s->fill_inverse_mass_matrix();

			// For good remeshing. Interpolation of values from old grig to new one.
			vector <double> very_old_coords = hllc_s->grid.GetCoordinates();
			vector <vector <vector <double> > > cvss;
			cvss.push_back(vector <vector <double> >());
			for (int var = 0; var < hllc_s->CONS_VAR_COUNT; ++var)
				cvss[0].push_back(hllc_s->cvn[var]);
			cvss.push_back(vector <vector <double> >());
			for (int var = 0; var < hllc_s->CONS_VAR_COUNT; ++var)
				cvss[1].push_back(hllc_s->cvnm1[var]);

			hllc_s->grid.SetRow("old_coords", hllc_s->grid.GetValues("coordinate"));
			hllc_s->cvn_old = hllc_s->cvn;
			hllc_s->cvnm1_old = hllc_s->cvnm1;

			for (int iter = 0; iter < maxiter && drho > convtol && true; ++iter)
			{
				drho = hllc_s->Solve(physDt_);
			}

			vector <double> old_coords = hllc_s->grid.GetCoordinates();

			auto remesh = [&](int count, bool use_func, double tau) -> void {
				if (hllc_s->remesh) {

					// New adaptive grid
					for (int i = 0; i < count; ++i) {
						for (int var = 0; var < hllc_s->CONS_VAR_COUNT; ++var)
							hllc_s->grid.SetRow(hllc_s->c_var_name[var], hllc_s->cv[var]);

						vector< vector <double> > functions;;
						vector< double > Fs;

						if (use_func)
						{
							functions.push_back(hllc_s->cvn[hllc_s->RHO_A]);
							Fs.push_back(1.);
						}

						hllc_s->grid.CalculateResolution(1., 1., hllc_s->c_var_name[hllc_s->RHO_A], "coordinate", functions, Fs);
						hllc_s->grid.CalculateConcentration(1., "coordinate");

						for (int var = 0; var < hllc_s->CONS_VAR_COUNT; ++var)
							hllc_s->grid.SetRow(hllc_s->c_var_name[var], hllc_s->cv[var]);
						hllc_s->grid.SetRow("coordinate", old_coords);

						hllc_s->x = hllc_s->grid.RefineMesh(physDt_, tau);

						for (int var = 0; var < hllc_s->CONS_VAR_COUNT; ++var)
							hllc_s->cv[var] = hllc_s->grid.GetValues(hllc_s->c_var_name[var]);

						hllc_s->RefreshBoundaries();						// Refresh boundary conditions
						hllc_s->RhoUPH();

						hllc_s->CalculateVolumes();
						hllc_s->grid.SetRow("volume", hllc_s->vol);
					}
					// Now we will reevaluate cvn, cvnm1, cvold
					vector < vector < vector < double > >* > cvs{ &hllc_s->cvn, &hllc_s->cvnm1/*, &hllc_s->cvold*/ };

					if (count > 0)
					{
						int i = 0;
						for (auto it : cvs)
						{
							hllc_s->grid.SetRow("coordinate", very_old_coords);
							for (int var = 0; var < hllc_s->CONS_VAR_COUNT; ++var)
								hllc_s->grid.SetRow(hllc_s->c_var_name[var], cvss[i][0]);
							++i;

							vector < string > ignore;
							vector < vector < double > > new_tab;
							ignore.push_back(hllc_s->grid.TYPE_COL);
							new_tab = hllc_s->grid.NewTable("coordinate", hllc_s->x, ignore, false);
							hllc_s->grid.SetData(new_tab);

							for (int var = 0; var < hllc_s->CONS_VAR_COUNT; ++var)
								(*it)[var] = hllc_s->grid.GetValues(hllc_s->c_var_name[var]);
						}

						hllc_s->grid.SetRow("coordinate", hllc_s->x);
						for (int var = 0; var < hllc_s->CONS_VAR_COUNT; ++var)
							hllc_s->grid.SetRow(hllc_s->c_var_name[var], hllc_s->cv[var]);

						hllc_s->calculate_mass_matrix();
						hllc_s->fill_inverse_mass_matrix();
					}
				}
			};

			remesh(0, false, 1e-2);

			hllc_s->cvnm1 = hllc_s->cvn;
			hllc_s->cvn = hllc_s->cv;

			ttime += physDt_;
			iter += hllc_s->iter;
			cout << iter << "\t" << ttime << endl;
		}
	}
	cout << "Execution time: " << (clock() - start_time) / CLOCKS_PER_SEC << "sec." << endl;
	hllc_s->PrintResult();
	hllc_s->deactivate_adept_stack();

	cout << endl << "   OK   " << endl;
}

void Rad_function(const std::vector<std::vector<double>>& params, const std::map < std::string, int >& var_name_ids, const std::vector<std::string>& names, int i, double& result)
{
	const int T_id = var_name_ids.at(names[0]);
	const double& T = params[T_id][i];
	constexpr double C = 3e-22;
	if (T < 2e4)
	{
		result = C * pow((0.5e4 * T), 3);
	}
	else if (T <= 2e5)
	{
		result = C;
	}
	else
	{
		result = C / sqrt(0.5e5 * T) + 2e-23 * sqrt(1e-8 * T);
	}
}

void AdjustMesh(Solver* solver, Loop& loop)
{
	
	vector<vector<double>> funcs(1, solver->grid.GetValues("T"));
	constexpr double X_max = 12500000.;
	constexpr double Rho_max = 1.0793596511952E-10;
	solver->grid.CalculateResolution(X_max, Rho_max, solver->c_var_name[solver->RHO_A], "coordinate", funcs, { 1056100. });
	solver->grid.CalculateConcentration(X_max, "coordinate");

	solver->x = solver->grid.RefineMesh(1., 10., 1.);
	vector < double > ran1(solver->x);

	solver->SetGrid(loop);

	vector < string > ignore(1, loop.TYPE_COL);
	vector < vector < double > > new_tab = solver->grid.NewTable("coordinate", ran1, ignore, false);
	solver->grid.SetData(new_tab);

	init_sun_atmosphere(solver);

	for (unsigned int var = 0; var < solver->CONS_VAR_COUNT; ++var)
		solver->cv[var] = solver->grid.GetValues(solver->c_var_name[var]);

	solver->RefreshBoundaries();						// Refresh boundary conditions
	solver->RhoUPH();
}

void init_sun_atmosphere(Solver* solver)
{
	cout << "Initializing flow" << endl;

	int eq_num = solver->eq_num;
	int imax = solver->imax;
	int var_num = solver->var_num;

	solver->cv.resize(eq_num);
	solver->cvold.resize(eq_num);
	solver->diss.resize(eq_num);
	solver->rhs.resize(eq_num);
	solver->Q_star.resize(eq_num);
	for (int i = 0; i < eq_num; ++i) {
		solver->cv[i].resize(imax, 0.);
		solver->cvold[i].resize(imax, 0.);
		solver->diss[i].resize(imax, 0.);
		solver->rhs[i].resize(imax, 0.);
		solver->Q_star[i].resize(imax, 0.);
	}
	solver->p.resize(imax, 0.);
	solver->dt.resize(imax, 0.);
	solver->dummy.resize((imax + 1) * (var_num + 1), 0.);

	const vector<double>& rho = solver->grid.GetValues("Rho");
	const vector<double>& Temp = solver->grid.GetValues("T");
	constexpr double k_Planc = 6.626e-34;
	constexpr double m_proton_rep = 5.9786e+26;
	for (int i = 0; i < imax; ++i)
	{
		solver->cv[solver->RHO_A][i] = rho[i] * solver->a[i];		//	$\rho S$
		solver->cv[solver->RHO_U_A][i] = 0.;					//	Initial velocity iz zero, I think $\rho u S$
		solver->p[i] = 2. * m_proton_rep * rho[i] * k_Planc * Temp[i];
		solver->cv[solver->RHO_E_A][i] = 1.5 * solver->p[i] * solver->a[i];	//	$\rho E S$
	}

	for (unsigned int var = 0; var < solver->CONS_VAR_COUNT; ++var)
		solver->grid.AddColumn(solver->c_var_name[var], solver->cv[var]);

	solver->grid.CalculateResolution(1., 1., solver->c_var_name[solver->RHO_A], "coordinate");
	solver->grid.CalculateConcentration(1., "coordinate");

	if (!solver->time_expl) {
		solver->L_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
		solver->U_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
		solver->D_SGS.resize(imax + 1, MatrixXd(eq_num, eq_num));
	}
}

void loop_foot_point(const string& output_file, const string& yml_file, const double end_time_)	// Sun loop
{
	cout << endl << "  ---  Loop foot point heating test  ---  " << endl;

	string file_name = yml_file;

	Loop loop;
	Solver* hll;

	hll = CreateReadConfigFile(file_name);

	hll->HideOutput();

	loop.ReadFile("Input/sun_loop_mesh.txt");

	// Make even grid
	//vector < vector < double > > new_tab;
	//vector < double > ran1 = loop.MakeRange(0., 12500000., 50000.);
	//vector < string > ignore;
	//ignore.push_back(loop.TYPE_COL);
	//new_tab = loop.NewTable("coordinate", ran1, ignore, true);
	//loop.PrintTable("Output/sum_mesh.txt", new_tab);

	hll->SetGrid(loop);
	hll->ReadBoundaries("Input/sun_loop_boundary.yml");
	hll->SetOutputFile(output_file);

	// Euler euqations
	hll->SetEquation("mass", { "Rho", "*A" }, { "Rho", "*U", "*A" }, { "" }, hll->vars, hll->vars_o);	// RhoA, RhoUA
	hll->SetEquation("impulse", { "Rho", "*U", "*A" }, { "Rho", "*U^2", "+p", "*A" }, { "Rho", "*g" }, hll->vars, hll->vars_o);
	hll->SetEquation("energy", { "Rho", "*E", "*A" }, { "Rho", "*E", "+p", "*U", "-FT", "*A" }, { "n^2", "*RadFunc" }, hll->vars, hll->vars_o);

	hll->set_fv_equation(		// Rho = RhoA / A
		"Rho",
		{ "RhoA", "/A", "" }		// There is dummy for unambiguous conservation
	);
	hll->set_fv_equation(		// n = 5.9786e+26 * Rho
		"n",
		{ "5.9786e+26Rho" }
	);
	hll->set_fv_equation(		// E = RhoEA / RhoA
		"E",
		{ "RhoEA", "/RhoA", "" }		// There is dummy for unambiguous conservation
	);
	hll->set_fv_equation(		// U = RhoUA / RhoA
		"U",
		{ "RhoUA", "/RhoA", "" }
	);
	hll->set_fv_equation(		// p = (RhoEA / RhoA - 0.5U^2) * (gamma - 1) * RHO
		"p",
		{ "0.66666667Rho", "*E", "" }
	);
	hll->set_fv_equation(		// H = gamma / (gamma - 1) * p / RHO + 0.5U^2
		"H",
		{ "GAMMA", "/GAMMAM", "*p", "/Rho", "+0.5U^2" }
	);
	hll->set_fv_equation(		// A
		"A",
		{ "A" }
	);
	hll->set_fv_equation(		// T = p * gamma / (gamma - 1) / Cp / rho
		"T",
		{ "p", "/2n", "/k" }
	);
	hll->set_fv_equation(		// FT = 10e-6 * T^2.5 * dT
		"FT",
		{ "1e-6T^2.5", "*dT", "" }
	);
	hll->AddDelayedFvEquation({ "FT" });

	hll->set_function("RadFunc", Rad_function, hll->fv, hll->vars_o, { "T" });
	//hll->set_fv_equation(		// Lam_a = 3.75e-11 * T^3
	//	"Rad",
	//	{ "RadFunc" }
	//);

	double convtol = hll->GetTolerance();
	double drho = 1.;
	int maxiter = hll->GetMaxIterNum();

	init_sun_atmosphere(hll);

	// Parameters to adjust mesh were used
	if (hll->remesh)
		for (int i = 0; i < 200; ++i) {
			AdjustMesh(hll, loop);
		}

	hll->CalculateVolumes();
	hll->grid.SetRow("volume", hll->vol);

	hll->grid.PrintTable("Output/sun_table.txt");

	return;

	hll->RhoUPH();
	hll->RefreshBoundaries();

	double fac = 1. / 1e1;
	double physDt_ = 0.2e-2 * fac / 1.;

	hll->cvnm1 = hll->cv;
	hll->iter = 0.;
	hll->Global_Time = 0.;
	hll->cvn = hll->cv;

	hll->S.resize((hll->ib2 - 1) * hll->eq_num, (hll->ib2 - 1) * hll->eq_num);
	hll->S.reserve(VectorXi::Constant((hll->ib2 - 1) * hll->eq_num, hll->eq_num * hll->eq_num));

	hll->iter = 0.;

	int iter = 0;
	double ttime = 0.;

	double start_time = clock();
	if (hll->time_stepping == 0) {
		while (ttime < end_time_) {
			hll->iter = 0;
			drho = 1.;
			hll->calculate_mass_matrix();
			hll->fill_inverse_mass_matrix();

			hll->cvn_old = hll->cvn;
			hll->cvnm1_old = hll->cvnm1;

			for (int iter = 0; iter < maxiter && drho > convtol && true; ++iter)
			{
				drho = hll->Solve(physDt_);
			}

			hll->cvnm1 = hll->cvn;
			hll->cvn = hll->cv;

			ttime += physDt_;
			iter += hll->iter;
			cout << iter << "\t" << ttime << endl;
		}
	}
	cout << "Execution time: " << (clock() - start_time) / CLOCKS_PER_SEC << "sec." << endl;
	hll->PrintResult();
	hll->deactivate_adept_stack();

	cout << endl << "   OK   " << endl;
}

void compare_files(const string &first_file, const string &second_file, const double eps)
{
	FILE* file1;
	FILE* file2;
	char str_c_1[20] = "\0", str_c_2[20] = "\0";
	string str_s_1, str_s_2;
	double value_1 = 0., value_2 = 0.;

	cout << "Opening first file " << "\"" << first_file.c_str() << "\"..." << endl;
	if ((file1 = fopen(first_file.c_str(), "r")) == NULL) {
		cout << "File \"" << first_file.c_str() << "\" wasn't opened." << endl;
		return;
	}

	cout << "Opening second file " << "\"" << second_file.c_str() << "\"..." << endl;
	if ((file2 = fopen(second_file.c_str(), "r")) == NULL) {
		cout << "File \"" << second_file.c_str() << "\" wasn't opened." << endl;
		return;
	}

	vector <string> header;
	int col_num = 0;
	bool is_data = false;
	int data_count = 0;
	int diff_count = 0;

	while (!feof(file1)) {
		if (feof(file1))
		{
			cout << "file1 fscanf feof" << endl;
			break;
		}
		else
			fscanf(file1, "%s", str_c_1, sizeof(str_c_1));

		if (feof(file2))
		{
			cout << "file2 fscanf feof" << endl;
			break;
		}
		else
			fscanf(file2, "%s", str_c_2, sizeof(str_c_2));

		str_s_1 = str_c_1;
		str_s_2 = str_c_2;
		if (isalpha(str_s_1[0]))
			header.push_back(str_s_1);
		else if (!is_data)
		{
			col_num = header.size();
			is_data = true;
		}

		if (str_s_1 != str_s_2)
		{
			sscanf(str_s_1.c_str(), "%le", &value_1);
			sscanf(str_s_2.c_str(), "%le", &value_2);
			double diff = (value_1 - value_2) * 2. / (value_1 + value_2);
			if (fabs(diff) > eps)
			{
				++diff_count;
				int line = data_count / col_num;
				int data_type = data_count % col_num;
				cout << "Difference: line " << line << "; " << header[data_type] << " : " << value_1 << "\t" << value_2 << " (" << diff << ")" << endl;
				//return;
			}
		}
		if (is_data)
			++data_count;
	}
	if (diff_count == 0)
		cout << "Files coincide" << endl;
	cout << "Comparison is finished" << endl;
}