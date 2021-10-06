// Sun.cpp: определяет точку входа для консольного приложения.
//

//#include "stdafx.h"
#include <omp.h>
#include <functional>
#include "BaseTest.h"
#include "TimeTest.h"
#include "LavalTest.h"
#include "GridTest.h"
//#include "adept/adept_source.h"
//#include "adept/algorithm.h"
#include "adept.h"

/*!
	\brief Main procedure

	Doing something

*/
void algo(int n, const adept::adouble* x, int m, adept::adouble* y);
//double algorithm_ad(const double x_val[2], // Input values
//	double* Y_ad,          // Input-output adjoint
//	double x_ad[2]);

int main()
{
	//int in;
	////MatrixReplacement A;
	////SparseMatrix<double> A(4, 4);
	//MatrixXd A(4, 4);
	//Eigen::VectorXd b(4), x;
	//b << 1, 1, 1, 1;
	//A << 35, 0.5, 0.7, 4,
	//	2, 25, 0.04, 1,
	//	1, -2, 25, 8,
	//	8, -5, 6, 46;
	//// solve Ax = b using CG with matrix-free version:
	//Eigen::ConjugateGradient < /*MatrixReplacement*/ MatrixXd , Eigen::Lower | Eigen::Upper, MyPreconditioner<double> > cg;
	//Eigen::VectorXd invdiag(4);
	//invdiag << 1. / 3., 1. / 4., 1. / 4., 1. / 3.;
	//cg.setMaxIterations(1000.);
	//cg.setTolerance(1e-10);
	//cg.preconditioner().setInvDiag(invdiag);
	//cg.compute(A);
	//x = cg.solve(b);
	//std::cout << "#iterations: " << cg.iterations() << std::endl;
	//std::cout << "estimated error: " << cg.error() << std::endl;
	//cout << x << endl;
	//cin >> in;

	// Base();

	// Time();

	 //Laval2();

	// Laval3();

	 Laval4();

	//Laval5();

	//Laval6();		// Unsteady

	//Laval7();		// Unsteady + Frid

	//Grid();

	//SodGrid();

	/*double x[2] = { 2.0, 3.0 };
	double y_ad = 1.0;
	double x_ad[2];
	double y = algorithm_ad(x, &y_ad, x_ad);
	std::cout << "x[0] = " << x[0] << "\n"
		<< "x[1] = " << x[1] << "\n"
		<< "y    = " << y << "\n"
		<< "y_ad = " << y_ad << "\n"
		<< "x_ad[0]=" << x_ad[0] << "\n"
		<< "x_ad[1]=" << x_ad[1] << "\n";*/

	//double x_val[3];
	//double y_val[3];
	//double jac[9];

	//x_val[0] = 1.;
	//x_val[1] = 2.;
	//x_val[2] = 3.;
	//adept::Stack stack;

	//cout << "Start" << endl;
	//for (int i = 0; i < 1000; ++i)
	//{
	//	x_val[2] = i;
	//	//adept::Stack stack;
	//	vector<adept::adouble> x(3);
	//	adept::set_values(&x[0], 3, x_val);
	//	/*x.resize(3);
	//	x[0] = 1.;
	//	x[1] = 1.;
	//	x[2] = 1.;*/

	//	stack.new_recording();
	//	vector<adept::adouble> y(3);
	//	algo(3, &x[0], 3, &y[0]);
	//	stack.independent(&x[0], 3);
	//	stack.dependent(&y[0], 3);
	//	stack.jacobian(jac);
	//	//stack.clear_dependents();
	//	//stack.clear_independents();
	//	for (int iy = 0; iy < 3; ++iy)
	//		y_val[iy] = y[iy].value();

	//	cout << jac[8] << "\t";
	//}
	//cout << "\nEnd" << endl;

    return 0;
}

void algo(int n, const adept::adouble* x, int m, adept::adouble* y)
{
	//adept::aVector y(3);
	y[0] = x[0];
	y[1] = x[1] * x[1];
	//if (x[2] < 500) {
	//	y[2] = x[2] * x[2] * x[2];
	//}
	//else {
		y[2] = x[2] * x[2] * x[2] * x[2];
	//}
}

//double algorithm_ad(const double x_val[2], // Input values
//	double* Y_ad,          // Input-output adjoint
//	double x_ad[2]) {      // Output adjoint
//	using namespace adept;                   // Import Stack and adouble from adept
//	Stack stack;                             // Where differential information is stored
//	adouble x[2] = { x_val[0], x_val[1] };     // Initialize adouble inputs
//	stack.new_recording();                   // Start recording derivatives
//	adouble Y = algorithm(x);                // Version overloaded for adouble args
//	Y.set_gradient(*Y_ad);                   // Load the input-output adjoint
//	stack.reverse();                         // Run the adjoint algorithm
//	x_ad[0] = x[0].get_gradient();           // Extract the output adjoint for x[0]
//	x_ad[1] = x[1].get_gradient();           //   ...and x[1]
//	*Y_ad = Y.get_gradient();              // Input-output adjoint has changed too
//	return Y.value();                        // Return result of simple computation
//}