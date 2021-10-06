#include "LavalTest.h"

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

	string file_name = "Input/sod.yml";

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
	cusp_s->InitFlowAG(rho_, mass_, e_, p_, 0.3);

	for (int i = 0; i < 1; ++i) {
		cusp_s->AdjustMesh(rho_, mass_, e_, p_, 0.3);
	}
	cusp_s->CalculateVolumes();
	cusp_s->grid.SetRow("volume", cusp_s->vol);

	cusp_s->RhoUPH();
	cusp_s->RefreshBoundaries();

	//double physDt_ = 0.2e-2;
	double physDt_ = 1e-10;

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
	cusp_s->iter = 0.;
	int iter = 0;
	double ttime = 0.;
	if (cusp_s->time_stepping == 0) {
		//while (ttime < 0.2 / 1. * 1.) {
		while (ttime < physDt_) {
			//while (cusp_s->Global_Time < 0.2*50) {
			cusp_s->iter = 0;
			drho = 1.;
			for (int iter = 0; iter < maxiter && drho > convtol && true; ++iter)
			{
				drho = cusp_s->Solve(physDt_);
				//drho = cusp_s->SolveImplicit();
				//break;
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