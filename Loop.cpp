#include "Loop.h"
#include <numeric>

vector< double > Loop::GetCoordinates()
{
	return GetValues(COORD_COL);
}

vector< double > Loop::GetSection()
{
	return GetValues(SECT_COL);
}

void Loop::SetTime(double time_)
{
	time = time_;
}

void Loop::IncrementTime(double dt)
{
	time += dt;
}

double Loop::GetTime()
{
	return time;
}

void Loop::Refresh(string col_name, function<double(vector < double >, int)> foo, vector < double > parameters)
{
	for (int i = 0; i < col_size; i++)
		Data[columns[col_name]][i] = foo(parameters, i);
}

void Loop::Refresh(string col_name, DataBase db)
{
	double val = db.GetFromPoint(TIME_COL, time)[db.columns[col_name]];

	SetRow(col_name, val);
}

void Loop::RefreshX(string col_name)
{
	SetRow(col_name, x);
}

void Loop::RefreshR(string col_name)
{
	SetRow(col_name, R);
}

void Loop::RefreshN(string col_name)
{
	SetRow(col_name, n);
}

void Loop::CalculateResolution(double X, double F, string function_name, string x_name, const vector < vector <double> >& functions, const  vector <double>& Fs)
{
	if (R.size() == 0)
		R.resize(col_size, 0.);

	vector < double > f = GetValues(function_name);
	vector < double > x = GetValues(x_name);

	for (int i = 1; i < col_size - 1; ++i) 
	{
		/*if (x[i + 1] - x[i] > 0.) {
			R[i] = sqrt(1. + pow(X / F * (f[i + 1] - f[i]) / (x[i + 1] - x[i]), 2));
		}
		else {
			R[i] = R[i - 1];
		}*/
		if (x[i + 1] - x[i - 1] > 0.) {
			//R[i] = sqrt(1. + pow(X / F * (f[i + 1] - f[i - 1]) / (x[i + 1] - x[i - 1]), 2));
			double factor;
			R[i] = 0.;
			factor = 0.;
			if (x[i + 1] - x[i] > 0.) {
				double sum = pow(X / F * (f[i + 1] - f[i]) / (x[i + 1] - x[i]), 2);
				for (int f_i = 0; f_i < functions.size(); ++f_i)
					sum += pow(X / Fs[f_i] * (functions[f_i][i + 1] - functions[f_i][i]) / (x[i + 1] - x[i]), 2);

				R[i] += sqrt(1. + sum);
				factor++;
			}
			if (x[i] - x[i - 1] > 0.) {
				double sum = pow(X / F * (f[i] - f[i - 1]) / (x[i] - x[i - 1]), 2);
				for (int f_i = 0; f_i < functions.size(); ++f_i)
					sum += pow(X / Fs[f_i] * (functions[f_i][i] - functions[f_i][i - 1]) / (x[i] - x[i - 1]), 2);

				R[i] += sqrt(1. + sum);
				factor++;
			}
			R[i] /= factor;
		}
		else {
			R[i] = R[i - 1];
		}
	}
	R[0] = R[1];
	R[col_size - 1] = R[col_size - 2];
}

void Loop::CalculateConcentration(double X, string x_name)
{
	if (n.size() == 0)
		n.resize(col_size, 0.);

	vector < double > x = GetValues(x_name);

	for (int i = 2; i < col_size - 2; ++i)
	{
		/*if (x[i + 1] - x[i] > 0.) {
			n[i] = X / (x[i + 1] - x[i]);
		}
		else {
			n[i] = n[i - 1];
		}*/
		if (x[i + 1] - x[i - 1] > 0.) {
			double factor(i == 1 || i == col_size - 2 ? 1. : 2.);
			//n[i] = factor * X / (x[i + 1] - x[i - 1]);
			n[i] = 0.;
			factor = 0.;
			if (x[i + 1] - x[i] > 0.) {
				n[i] += X / (x[i + 1] - x[i]);
				factor++;
			}
			if (x[i] - x[i - 1] > 0.) {
				n[i] += X / (x[i] - x[i - 1]);
				factor++;
			}
			n[i] /= factor;
		}
		else {
			n[i] = n[i - 1];
		}
	}
	n[1] = n[2];
	n[0] = n[1];
	n[col_size - 2] = n[col_size - 3];
	n[col_size - 1] = n[col_size - 2];

	if (TotalSum == 0.) {
		for (int i = 1; i < col_size; ++i) {
			//TotalSum += 0.5 * (x[i] - x[i - 1]) * (n[i] + n[i - 1]);
			if (i == 1 || i == col_size - 1) {
				TotalSum += FindArea(n.data() + i - 1, x.data() + i - 1, 0);
			}
			else {
				TotalSum += FindArea(n.data() + i - 1, x.data() + i - 1, 1);
			}
		}
	}
}

void Loop::CalculateConcentrationWave()
{
	if (n_w.size() == 0)
		n_w.resize(col_size, 0.);

	for (int i = 1; i < col_size - 1; ++i)
	{
		n_w[i] = n[i] - alpha * (alpha + 1.) * (n[i + 1] - 2. * n[i] + n[i - 1]);
	}
	n_w[0] = n_w[1];
	n_w[col_size - 1] = n_w[col_size - 2];
}

SparseMatrix< double > Loop::ILU_0(SparseMatrix< double >& SM)
{
	SparseMatrix< double > LU;
	int k_start, k_finish, j_start, j_finish;
	int size_ = SM.rows();

	//LU.reserve(VectorXi::Constant((ib2 - 1) * eq_num, eq_num * eq_num));
	LU = SM;

	// There is more to optimize
	for (int i = 1; i < size_; ++i) {
		k_start = max(0, i - 2);
		k_finish = i;
		for (int k = k_start; k < k_finish; ++k) {
			LU.coeffRef(i, k) = LU.coeff(i, k) / LU.coeff(k, k);
			j_start = k + 1;
			j_finish = min(size_, i + 1 + 1);
			for (int j = j_start; j < j_finish; ++j) {
				LU.coeffRef(i, j) = LU.coeff(i, j) - LU.coeff(i, k) * LU.coeff(k, j);		// <- Here for instance. We don't need compute if LU(k, j) is zero (No, it's OK)
			}
		}
	}
	LU.isCompressed();

	return LU;
}

Vector<Scalar, Dynamic> Loop::SolveFromLU(SparseMatrix<double>& LU, Vector<Scalar, Dynamic>& b)
{
	int size_ = LU.rows();
	Vector<Scalar, Dynamic> y(size_), x_(size_);
	int j_start, j_finish;

	y(0) = b(0);
	for (int i = 1; i < size_; ++i)
	{
		j_start = max(0, i - 2);
		y(i) = b(i);
		for (int j = j_start; j < i; ++j) {
			y(i) -= y(j) * LU.coeff(i, j);
		}
	}
	x_(size_ - 1) = y(size_ - 1) / LU.coeff(size_ - 1, size_ - 1);
	for (int i = size_ - 1 - 1; i >= 0; --i)
	{
		j_finish = min(size_, i + 1 + 1);

		x_(i) = y(i);
		for (int j = i + 1; j < j_finish; ++j) {
			x_(i) -= x_(j) * LU.coeff(i, j);
		}
		x_(i) /= LU.coeff(i, i);
	}
	return x_;
}

vector < double > Loop::Redistribute(string col_name, vector < double >& conc)
{
	vector< double > x_ = GetValues(col_name);
	vector< double > c_ = GetConcentration();
	vector< double > new_x(x_.size());
	double X = x_[x_.size() - 1] - x_[0];
	double relax = 100;

	//x_.erase(x_.begin());
	//x_.pop_back();
	conc.insert(conc.begin(), conc[0]);
	conc.push_back(conc[conc.size() - 1]);

	new_x[0] = x_[0];
	new_x[1] = x_[1];
	conc[1] = conc[2];
	conc[0] = conc[1];
	conc[col_size - 2] = conc[col_size - 3];
	conc[col_size - 1] = conc[col_size - 2];
	//for (int i = 2; i < col_size - 2; ++i)
	//{
	//	new_x[i] = x_[i] + abs(X / conc[i] - X / c_[i]) * (-(X / conc[i + 1] - X / conc[i]) + (X / conc[i] - X / conc[i - 1])) * relax;
	//	//new_x[i] = x_[i] + abs(X / conc[i] - X / c_[i]) * (-(X / conc[i + 1] - X / c_[i + 1]) + (X / conc[i - 1] - X / c_[i - 1])) /*/ pow(double(x_.size() - 3), 2)*/ * relax;
	//	//new_x[i] = x_[i] + abs(X / conc[i] - X / c_[i]) * ( ((conc[i + 1] - c_[i + 1]) - (conc[i - 1] - c_[i - 1]) >= 0. ? 1. : -1. )) /*/ pow(double(x_.size() - 3), 2)*/ * relax;
	//}

	double S_star = TotalSum / (col_size - 3.);
	double S_left = S_star;
	double sec_area;
	double dx, dx1, dx2;
	int id = 2;
	double a, b, c, D;
	double conc_cur;
	int searchMethod = 2;
	int areaMethod = 0;
	bool PointsAreFilled = false;

	if (searchMethod == 2) {
		areaMethod = 1;
	}
	else {
		areaMethod = 0;
	}

	for (int i = 1; i < col_size - 1 && !PointsAreFilled; ++i)
	{
		//sec_area = 0.5 * (conc[i + 1] + conc[i]) * (x_[i + 1] - x_[i]);
		//FindArea(conc.data() + i, x_.data() + i, 1);
		sec_area = FindArea(conc.data() + i, x_.data() + i, areaMethod);
		if (sec_area < S_left) {
			S_left -= sec_area;
			continue;
		}
		else {
			double slope = (conc[i + 1] - conc[i]) / (x_[i + 1] - x_[i]);

			if (searchMethod == 0) {
				slope = 0.;
			}
			conc_cur = conc[i];

			SetNewPoint(S_left, sec_area, conc_cur, searchMethod, new_x, x_[i], id, slope, conc.data() + i, x_.data() + i);
			if (id >= col_size - 2) {
				PointsAreFilled = true;
				break;
			}
			S_left = S_star;
			while (sec_area >= S_left) {
				SetNewPoint(S_left, sec_area, conc_cur, searchMethod, new_x, new_x[id - 1], id, slope, conc.data() + i, x_.data() + i);
				if (id >= col_size - 2) {
					PointsAreFilled = true;
					break;
				}
			}
			S_left -= sec_area;
		}
	}

	//new_x[0] = x_[0];
	//new_x[1] = x_[1];
	new_x[col_size - 2] = x_[col_size - 2];
	new_x[col_size - 1] = x_[col_size - 1];

	double coef = 0.2;
	double coef2 = 0.2;
	vector < double > x2 = new_x;

	// Relaxation 1
	//for (int i = 2; i < col_size - 2; ++i)
	//{
	//	new_x[i] = coef * new_x[i] + (1. - coef) * x_[i];
	//	//new_x[i] = coef * x2[i] + (1. - coef) * 0.25 * (x2[i - 1] + x2[i + 1] + x_[i - 1] + x_[i + 1]);
	//}

	// Relaxation 2
	/*for (int i = 3; i < col_size - 3; ++i) {
		if ((new_x[i + 1] - new_x[i]) / (new_x[i] - new_x[i - 1]) > (alpha + 1.) / alpha) {
			new_x[i + 1] += (max(new_x[i] + (alpha + 1.) / alpha * (new_x[i] - new_x[i - 1]), new_x[i]) - new_x[i + 1]) * coef2;
		}
		if ((new_x[i + 1] - new_x[i]) / (new_x[i] - new_x[i - 1]) < alpha / (alpha + 1.)) {
			new_x[i + 1] += (min(new_x[i] + alpha / (alpha + 1.) * (new_x[i] - new_x[i - 1]), new_x[i + 2]) - new_x[i + 1]) * coef2;
		}
	}
	for (int i = col_size - 4; i <= 3; ++i) {
		if ((new_x[i] - new_x[i - 1]) / (new_x[i + 1] - new_x[i]) > (alpha + 1.) / alpha) {
			new_x[i - 1] += (min(new_x[i] - (alpha + 1.) / alpha * (new_x[i + 1] - new_x[i]), new_x[i]) - new_x[i - 1]) * coef2;
		}
		if ((new_x[i] - new_x[i - 1]) / (new_x[i + 1] - new_x[i]) < alpha / (alpha + 1.)) {
			new_x[i - 1] += (max(new_x[i] - alpha / (alpha + 1.) * (new_x[i + 1] - new_x[i]), new_x[i - 2]) - new_x[i - 1]) * coef2;
		}
	}*/

	// Relaxation 3
	//double coef3 = 0.5;
	//for (int times = 0; time < 5; ++time) {
	//	for (int i = 2; i < col_size - 2; ++i)
	//	{
	//		//new_x[i] = coef * new_x[i] + (1. - coef) * x_[i];
	//		new_x[i] = coef3 * x2[i] + (1. - coef3) * 0.5 * (x2[i - 1] + x2[i + 1]);
	//	}
	//	x2 = new_x;
	//}

	/*for (int i = 2; i < col_size - 2; ++i) {
		c_[i] = (conc[i] - alpha * (alpha + 1.) * (conc[i + 1] - 2. * conc[i] + conc[i - 1])) / R[i];
	}
	c_ = R;*/


	return new_x;
}

double Loop::FindArea(double* conc, double* x_, int method)
{
	if (method == 0 || x_[-1] == x_[0] || x_[1] == x_[2]) {
		return 0.5 * (conc[1] + conc[0]) * (x_[1] - x_[0]);
	}
	if (method == 1) {
		double x1, x2, x3;
		double n1, n2, n3;
		double a, b, c;
		double det;
		double S = 0.;
		
		for (int i = 0; i < 2; ++i) {
			x1 = x_[-1 + i];
			x2 = x_[0 + i];
			x3 = x_[1 + i];
			n1 = conc[-1 + i];
			n2 = conc[0 + i];
			n3 = conc[1 + i];

			det = (x2 - x3) * (x1 * (x1 - x2 - x3) + x2 * x3);
			a = ((x2 - x3) * n1 + (x3 - x1) * n2 + (x1 - x2) * n3) / det;
			b = ((x3 * x3 - x2 * x2) * n1 + (x1 * x1 - x3 * x3) * n2 + (x2 * x2 - x1 * x1) * n3) / det;
			c = ((x2 * x2 * x3 - x3 * x3 * x2) * n1 + (x3 * x3 * x1 - x1 * x1 * x3) * n2 + (x1 * x1 * x2 - x2 * x2 * x1) * n3) / det;

			S += FindCubicFunctionValue(a / 3., b / 2., c, 0., x_[1]) - FindCubicFunctionValue(a / 3., b / 2., c, 0., x_[0]);
		}
		S /= 2.;

		return S;
	}
}

double Loop::FindCubicFunctionValue(double a, double b, double c, double d, double x)
{
	return a * pow(x, 3) + b * pow(x, 2) + c * x + d;
}

void Loop::SetNewPoint(double S_left, double& sec_area, double& conc_cur, int searchMethod, vector < double >& new_x, double prev_x, int& id, double slope, double* conc, double* x_)
{
	double dx;

	//dxSearch(S_left, conc_cur, 2, slope, prev_x, conc, x_);
	dx = dxSearch(S_left, conc_cur, searchMethod, slope, prev_x, conc, x_);
	if (id >= 3 && dx > (alpha + 1.) / alpha * (new_x[id - 1] - new_x[id - 2])) {
		dx = (alpha + 1.) / alpha * (new_x[id - 1] - new_x[id - 2]);
	}
	/*if (id >= 3 && dx < alpha / (alpha + 1.) * (new_x[id - 1] - new_x[id - 2])) {
		dx = alpha / (alpha + 1.) * (new_x[id - 1] - new_x[id - 2]);
	}*/
	conc_cur = conc_cur + dx * slope;
	new_x[id] = prev_x + dx;
	id++;
	sec_area -= S_left;
}

double Loop::dxSearch(double S_left, double conc_, int method, double slope, double prev_x, double* conc, double* x_)
{
	if (method == 0 || slope == 0) {
		return S_left / conc_;
	}
	if (method == 1 || x_[-1] == x_[0] || x_[1] == x_[2]) {
		double a, b, c, D, dx1, dx2;
		a = slope / 2.;
		b = conc_;		// is always positive
		c = -S_left;
		D = b * b - 4. * a * c;		// must be positive
		dx1 = (-b + sqrt(D)) / (2. * a);	// if a is positive - this is the only option
		dx2 = (-b - sqrt(D)) / (2. * a);	// if a is negative both variants are option
		if (a > 0.) {
			return dx1;
		}
		else {
			if (dx1 < 0.) {
				return dx2;
			}
			else {
				return dx1;
			}
			//cout << "Think!" << endl;
		}
	}
	if (method == 2) {
		double x1, x2, x3;
		double n1, n2, n3;
		double a = 0., b = 0., c = 0., d = 0.;
		double det;
		double S = 0.;
		double p, q, Q;
		double alpha, beta;
		double x_ans;

		for (int i = 0; i < 2; ++i) {
			x1 = x_[-1 + i];
			x2 = x_[0 + i];
			x3 = x_[1 + i];
			n1 = conc[-1 + i];
			n2 = conc[0 + i];
			n3 = conc[1 + i];

			det = (x2 - x3) * (x1 * (x1 - x2 - x3) + x2 * x3);
			a += ((x2 - x3) * n1 + (x3 - x1) * n2 + (x1 - x2) * n3) / det / 3. * 0.5;
			b += ((x3 * x3 - x2 * x2) * n1 + (x1 * x1 - x3 * x3) * n2 + (x2 * x2 - x1 * x1) * n3) / det / 2. * 0.5;
			c += ((x2 * x2 * x3 - x3 * x3 * x2) * n1 + (x3 * x3 * x1 - x1 * x1 * x3) * n2 + (x1 * x1 * x2 - x2 * x2 * x1) * n3) / det * 0.5;
		}

		d = -FindCubicFunctionValue(a, b, c, 0, prev_x) - S_left;

		p = (3. * a * c - b * b) / (3. * a * a);
		q = (2. * pow(b, 3) - 9. * a * b * c + 27 * a * a * d) / (27. * pow(a, 3));
		Q = pow(p / 3., 3) + pow(q / 2., 2);
		//Q = max(Q, 0.);
		if (Q >= 0.) {
			alpha = -q / 2. + sqrt(Q);
			beta = -q / 2. - sqrt(Q);
			alpha = cbrt(alpha);
			beta = cbrt(beta);

			x_ans = alpha + beta - b / (3. * a);
		}
		else {
			double Re = -q / 2.;
			double Im = -sqrt(abs(Q));
			double r = sqrt(Re * Re + Im * Im);
			double phi = atan(Im / (Re + 1e-40));

			alpha = cbrt(r) * (cos(phi / 3.) + sin(phi / 3.));
			beta  = cbrt(r) * (cos(phi / 3.) - sin(phi / 3.));

			double x_3ans[3];
			x_3ans[0] = alpha + beta - b / (3. * a);
			x_3ans[1] = -(alpha + beta) / 2. + (alpha - beta) / 2. * sqrt(3) - b / (3. * a);
			x_3ans[2] = -(alpha + beta) / 2. - (alpha - beta) / 2. * sqrt(3) - b / (3. * a);
			double dx_min = 1e10;
			int id_dx_min = 0;

			for (int j = 0; j < 3; ++j) {
				if (abs(x_3ans[j] - prev_x) < dx_min && x_3ans[j] - prev_x > 0.) {
					dx_min = x_3ans[j] - prev_x;
					id_dx_min = j;
				}
			}

			x_ans = x_3ans[id_dx_min];
		}

		return x_ans - prev_x;

			//S += FindCubicFunctionValue(a, b, c, d, x_[1]) - FindCubicFunctionValue(a, b, c, d, x_[0]);
	}
}

vector< double > Loop::RefineMesh(double dt_, double tau_, double relax_coef, vector <string> ignore)
{
	double tau = tau_;
	double d_t = dt_;
	double t_frac = 1. + tau / d_t;
	double coef;
	int off = 2.;
	int j;

	CalculateConcentrationWave();
	//n_w = n_w;

	Vector<Scalar, Dynamic> b(col_size - off);
	Vector<Scalar, Dynamic> NewConc, check;
	S.resize(col_size - off, col_size - off);
	S.reserve(VectorXi::Constant((col_size - off), 4));
	//S.setZero();
	for (int i = 0; i < col_size - off; ++i) {
		j = i + off / 2.;
		if (i > 1 && i < col_size - off - 2) {
			coef = alpha * (alpha + 1.);
			S.insert(i, i - 2) = - coef / R[j - 1] * t_frac;
			S.insert(i, i - 1) = ((2. * coef + 1.) / R[j - 1] + coef / R[j]) * t_frac;
			S.insert(i, i + 0) = (-coef / R[j - 1] - (2. * coef + 1.) / R[j]) * t_frac;
			S.insert(i, i + 1) = coef / R[j] * t_frac;
			b(i) = tau / d_t * (n_w[j - 1] / R[j - 1] - n_w[j] / R[j]);
		} 
		else if (i == 0 || i == col_size - off - 1) {
			S.insert(i, i) = 1.;
			b(i) = GetConcentration()[j];
		} 
		else if (i == 1 || i == col_size - off - 1 - 1) {
			if (i == 1) {
				S.insert(i, i - 1) = -1.;
			}
			else {
				S.insert(i, i + 1) = -1.;
			}
			S.insert(i, i) = 1.;
			b(i) = 0.;
		}
	}
	S.makeCompressed();

	//MatrixXd A (S);

	SparseMatrix<double> LU = ILU_0(S);
	NewConc = SolveFromLU(LU, b);		// check 9-10-11 iterations..

	//BiCGSTAB<SparseMatrix<double> > solver;

	//SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
	//solver.analyzePattern(S);
	//solver.factorize(S);
	//
	//solver.compute(S);

	//NewConc = solver.solve(b);
	/*NewConc = A.colPivHouseholderQr().solve(b);*/


	//check = S * NewConc;
	//cout << (check - b).sum() << "\t" << b.sum();
	vector < double > x_ = GetValues("coordinate");;
	vector < double > vec;
	vec.insert(vec.end(), &NewConc.data()[0], &NewConc.data()[col_size - off]);

	// Appropriate coordinates redistribution.
	vector <double> result_dx(vec.size() - 1);
	//std::transform(vec.begin(), vec.end() - 1, result_dx.begin(), [](double n) -> double {return 1. / n;  });
	result_dx[0] = 1. / vec[0];
	std::transform(vec.begin() + 1, vec.end() - 1, result_dx.begin(), result_dx.begin() + 1, [](double n, double prev_dx) -> double {return 1. / (2. * n - 1. / prev_dx); });
	double x_sum = std::accumulate(result_dx.begin(), result_dx.end(), 0.);
	std::transform(vec.begin(), vec.end(), vec.begin(), [x_sum](double n) -> double {return n * x_sum; });

	result_dx[0] = 1. / vec[0];
	std::transform(vec.begin() + 1, vec.end() - 1, result_dx.begin(), result_dx.begin() + 1, [](double n, double prev_dx) -> double {return 1. / (2. * n - 1. / prev_dx); });
	//x_sum = std::accumulate(result_dx.begin(), result_dx.end(), 0.);

	vec.insert(vec.begin(), 0.);
	//std::transform(vec.begin() + 1, vec.end() - 1, vec.begin(), vec.begin() + 1, [](double n, double prev_x) -> double {return prev_x + 1. / n; });
	std::transform(result_dx.begin(), result_dx.end() - 1, vec.begin(), vec.begin() + 1, [](double dx, double prev_x) -> double {return prev_x + dx; });

	std::transform(vec.begin(), vec.end(), vec.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, x_.back() - x_[0]));
	if (fabs(x_[0]) > 1e-10)
		std::transform(vec.begin(), vec.end(), vec.begin(), std::bind(std::plus<double>(), std::placeholders::_1, x_[0]));

	vec.insert(vec.begin(), x_[0]);
	vec[vec.size() - 2] = x_.back();
	vec.back() = x_.back();

	//double TS = 0.;
	//for (int i = 1; i < col_size - off; ++i) {
	//	//TS += 0.5 * (x_[i + off/2.] - x_[i + off/2. - 1]) * (vec[i] + vec[i - 1]);
	//	if (i <= 1 || i >= col_size - 1 - off) {
	//		TS += FindArea(vec.data() + i - 1, x_.data() + i + int(off) / 2 - 1, 0);
	//	}
	//	else {
	//		TS += FindArea(vec.data() + i - 1, x_.data() + i + int(off) / 2 - 1, 1);
	//	}
	//}
	//for (int i = 0; i < col_size - off; ++i) {
	//	vec[i] /= TS / TotalSum;
	//}

	//vec = Redistribute("coordinate", vec);

	//coef = relax_coef;
	//for (int i = 0; i < vec.size(); ++i)
	//{
	//	vec[i] = coef * vec[i] + (1. - coef) * x_[i];
	//}

	//vector < string > ignore;
	vector < vector < double > > new_tab;
	ignore.push_back(TYPE_COL);
	new_tab = NewTable("coordinate", vec, ignore, false);
	SetData(new_tab);
	/*CalculateResolution(1., 1., "rho", "coordinate");
	CalculateConcentration(1., "coordinate");*/

	return vec;
}

void Loop::FindMesh(double dt, string x_name)
{
	n_w_old = n_w;

	//n = n_w;

	double sum = 0.;
	double R_integral = 0.;
	vector < vector < double > > matrix;

	if (matrix.size() == 0) {
		matrix.resize(col_size);
		for (int i = 0; i < col_size; ++i)
		{
			matrix[i].resize(col_size + 1, 0.);
		}
	}
	else {
		for (int i = 0; i < col_size; ++i)
		{
			matrix[i].assign(col_size + 1, 0.);
		}
	}
	

	x = GetValues(x_name);

	for (int i = 0; i < col_size - 1; ++i)
	{
		//sum += (x[i + 1] - x[i]) * (n[i + 1] + n[i]) / 2.;
		sum += (x[i + 1] - x[i]) * n[i];
		R_integral += (x[i + 1] - x[i]) * R[i];
	}
	cout << sum << endl;

	double sum2 = 0.;

	double delta_t = 0.01;
	double tau = 0.001;

	double aap1 = alpha * (alpha + 1.);
	for (int i = 2; i < col_size - 2; ++i)
	{
		/*matrix[i - 2][i - 1] = -R[i - 1] / R[i] * alpha*(alpha + 1.) - 2.*alpha*(alpha + 1.);
		matrix[i - 2][i] = R[i - 1] / R[i] * (1. + 2.*alpha*(alpha + 1.)) + alpha * (alpha + 1.);
		matrix[i - 2][i + 1] = -R[i - 1] / R[i] * alpha*(alpha + 1.);*/
		matrix[i][i - 2] = (-aap1 / R[i - 1]) * (1. + tau / delta_t);
		matrix[i][i - 1] = (((2. * aap1 + 1.) / R[i - 1] + aap1 / R[i])) * (1. + tau / delta_t);
		matrix[i][i + 0] = (-(aap1 / R[i - 1] + (2. * aap1 + 1.) / R[i])) * (1. + tau / delta_t);
		matrix[i][i + 1] = (aap1 / R[i]) * (1. + tau / delta_t);
	}
	matrix[0][0] = 1.;
	matrix[1][1] = 1.;
	matrix[col_size - 2][col_size - 2] = 1.;
	matrix[col_size - 1][col_size - 1] = 1.;
	matrix[0][col_size] = n[0];
	matrix[1][col_size] = n[1];
	matrix[col_size - 2][col_size] = n[col_size - 2];
	matrix[col_size - 1][col_size] = n[col_size - 1];
	for (int i = 2; i < col_size - 2; ++i)
	{
		matrix[i][col_size] = (n[i - 1] / R[i - 1] - n[i] / R[i]) * tau / delta_t;
	}
	gauss(matrix, n);

	n[1] = n[2];
	n[0] = n[1];
	n[col_size - 2] = n[col_size - 3];
	n[col_size - 1] = n[col_size - 2];

	//SmoothN(0.95);

	for (int i = 1; i < col_size; ++i)
	{
		//n[i] = n[i - 1] * R[i] / R[i - 1];
		sum2 += (x[i] - x[i - 1]) * n[i - 1];
	}
	cout << " - " << sum2 << endl;
	for (int i = 0; i < col_size; ++i)
	{
		n[i] = n[i] / sum2 * sum;
	}

	vector < double > new_x(col_size, 0.);
	double sum3 = 0.;
	double area, slope;
	double a, b, c, D, x1, x2, b_prev;
	int id = 0;
	double smooth = 0.00;

	new_x[id++] = x[0];
	double area2 = 0.;
	// Integrating

	//for (int i = 0; i < col_size - 1; ++i)
	//{
	//	new_x[i + 1] = new_x[i] + 1 / n[i];
	//}

	for (int i = 0; i < col_size - 1 && true; ++i)
	{
		if (id == col_size) {
			break;
		}
		slope = (n[i + 1] - n[i]) / (x[i + 1] - x[i]);
		//area = (x[i + 1] - x[i]) * (n[i + 1] + n[i]) / 2.;
		area = (x[i + 1] - x[i]) * n[i];
		area2 += (x[i + 1] - x[i]) * R[i];
		//cout << i << "\t" << area2 / R_integral << endl;
		if (sum3 + area > 1.) {
			a = slope / 2.;
			a = 0.;
			/*if (abs(a) < 1e-10)
				a = 1e-10;*/

			if (a == 0.) {
				x1 = (1. - sum3) / n[i];
				x1 = x[i] + x1;
				new_x[id++] = x1;
				b = n[i] + smooth * slope * (new_x[id - 1] - x[i]);

				area -= (1. - sum3);
				area *= b / n[i];
				if (area < 1.) {
					sum3 = area;
				} else {
					//area -= (1. - sum3);
					sum3 = 0.;
					while (sum3 + area >= 1.) {
						//x1 = 1. / n[i];
						x1 = 1. / b;
						x1 = new_x[id - 1] + x1;
						if (id == col_size) break;
						new_x[id++] = x1;
						if (id == col_size) break;

						area -= 1.;

						b_prev = b;
						b = n[i] + smooth * slope * (new_x[id - 1] - x[i]);
						area *= b / b_prev;

						//sum3 = 0.;
					}
					sum3 = area;
				}
			}
			else {
				b = n[i];
				c = sum3 - 1.;
				D = b * b - 4 * a * c;
				x1 = (-b + sqrt(D)) / 2. / a;
				x2 = (-b - sqrt(D)) / 2. / a;
				if (x1 < 0.) x1 = x[col_size - 1];
				if (x2 < 0.) x2 = x[col_size - 1];

				x1 = min(x1, x2);
				x1 = x[i] + x1;
				new_x[id++] = x1;
				if (sum3 + area < 2.) {
					sum3 = sum3 + area - 1.;
				} else {
					area -= (1. - sum3);
					sum3 = 0.;
					while (sum3 + area >= 1.) {
						b = n[i] + slope * (new_x[id - 1] - x[i]);
						c = -1.;
						D = b * b - 4 * a * c;
						x1 = (-b + sqrt(D)) / 2. / a;
						x2 = (-b - sqrt(D)) / 2. / a;
						if (x1 < 0.) x1 = x[col_size - 1];
						if (x2 < 0.) x2 = x[col_size - 1];

						x1 = min(x1, x2);
						/*if (x1 > x[i + 1]) {

						}*/
						x1 = new_x[id - 1] + x1;
						if (id == col_size) break;
						new_x[id++] = x1;
						if (id == col_size) break;

						area -= 1.;
						//sum3 = 0.;
					}
					sum3 = area;
				}
			}
		} else {
			sum3 += area;
		}
	}
	new_x[col_size - 1] = x[col_size - 1];
	x = new_x;
	Data[columns[x_name]] = x;
}

void Loop::SmoothN(double coef)
{
	vector < double > smoothed_n(col_size);

	coef = min(max(0., coef), 1.);
	double side_coef = (1. - coef) / 2.;

	for (int i = 1 ; i < col_size - 1; ++i)
	{
		smoothed_n[i] = (pow(n[i - 1], 1) + pow(n[i + 1], 1)) /*/ (n[i - 1] + n[i + 1])*/ * side_coef + n[i] * coef;
	}
	smoothed_n[0] = n[0];
	smoothed_n[col_size - 1] = n[col_size - 1];

	n = smoothed_n;
}

// New smooting algorithm
const vector < double > Loop::smooth_least_square(const vector < double > x, const vector < double > f, unsigned int poly_degree, unsigned int half_count_points)
{
	if (poly_degree == 0)
		return f;

	int vec_size = x.size();
	vector <double> smoothed_f(vec_size);

	MatrixXf A(2 * half_count_points + 1, poly_degree + 1);
	VectorXf f_vec(2 * half_count_points + 1);
	VectorXf poly_coefs(poly_degree + 1);
	for (int row = 0; row < 2 * half_count_points + 1; ++row)
		A.coeffRef(row, poly_degree) = 1.;

	for (int i = 0; i < vec_size; ++i)
	{
		int A_row = 0;
		for (int row = i - int(half_count_points); row <= i + int(half_count_points); ++row)
		{
			int item = min(max(0, row), vec_size - 1);
			for (int A_col = int(poly_degree) - 1; A_col >= 0; --A_col)
				A.coeffRef(A_row, A_col) = x[item] * A(A_row, A_col + 1);

			f_vec.coeffRef(A_row) = f[item];
			++A_row;
		}
		poly_coefs = A.bdcSvd(ComputeThinU | ComputeThinV).solve(f_vec);
		for (unsigned int j = 0; j < poly_degree - 1; ++j)
			smoothed_f[i] += pow(x[i], poly_degree - j) * poly_coefs(j);
		smoothed_f[i] += x[i] * poly_coefs(poly_degree - 1) + poly_coefs(poly_degree);
	}

	return smoothed_f;
}

void Loop::AddCollumnR()
{
	AddColumn("R", R);
}

void Loop::AddCollumnN()
{
	AddColumn("n", n);
}

vector < double > Loop::GetConcentration()
{
	return n;
}
vector < double > Loop::GetResolution()
{
	return R;
}

double Loop::GetAlpha()
{
	return alpha;
}

int gauss(vector < vector < double > >& a, vector < double >& ans)
{
	int n = a.size();
	int m = a[0].size() - 1;
	double EPS = 1e-5;
	int INF = -1;

	vector < int > where(m, -1);
	for (int col = 0, row = 0; col < m && row < n; ++col)
	{
		int sel = row;
		for (int i = row; i < n; ++i)
			if (abs(a[i][col]) > abs(a[sel][col]))
				sel = i;
		if (abs(a[sel][col]) < EPS)
			continue;

		for (int i = col; i <= m; ++i)
			swap(a[sel][i], a[row][i]);

		where[col] = row;

		for (int i = 0; i < n; ++i)
			if (i != row) {
				double c = a[i][col] / a[row][col];
				for (int j = col; j <= m; ++j)
					a[i][j] -= a[row][j] * c;
			}
		++row;
	}

	ans.assign(m, 0.);
	for (int i = 0; i < m; ++i)
		if (where[i] != -1)
			ans[i] = a[where[i]][m] / a[where[i]][i];

	for (int i = 0; i < n; ++i)
	{
		double sum = 0.;

		for (int j = 0; j < m; ++j)
			sum += ans[j] * a[i][j];

		if (abs(sum - a[i][m]) > EPS)
			return 0.;

	}

	for (int i = 0; i < m; ++i)
		if (where[i] == -1)
			return INF;

	return 1;
}