// Central Differences Scheme
#pragma once

//using DiagonalFunc = std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int,
//	std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >, std::function < MatrixXd(std::vector < MatrixXd >&, std::vector < std::vector < double > >&, int) >) >;

#include "Solver.h"
//#include "adept/adept_source.h"

class VL : public Solver {
public:
	VL(sol_struct& sol_init_);
	virtual void LRState();
	virtual void LRState(string var_);
	virtual void LRState(vector < vector < double > >& cv_, vector < vector < double > >& ls_, vector < vector < double > >& rs_);
	virtual void Fluxes();
	double MUSCL0(double a, double b, double dx);
	double MUSCL1_3(double a, double b, double dx);
	//double CUSPLimiter(double af, double bf);
	// virtual double Solve();
	virtual const vector < int >& GetDissFlag();
	virtual const vector < double >& GetDissBlend();
	virtual void Dissipation(double beta);

	virtual void RHS(int i);
	virtual void ComputeRHSandJacobian(bool NO_JAC = false);
	//void ComputeRightLeftState(const adept::aVector& x, adept::aVector& LS, adept::aVector& RS);
	void GetSourceAndJacobian(int i, vector < double >& y_val, vector < double >& jac);
	void ComputeSourceTerm(int n, const adept::adouble* cv, int m, adept::adouble* source, int i);
	void GetPositiveFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac, bool simple = false);
	void GetPositiveFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac);
	void GetNegativeFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac, bool simple = false);
	void GetNegativeFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac);
	virtual void GetFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac, bool POS_NEG, bool simple = false);
	void ComputeFlux(int n, const adept::adouble* x, int m, adept::adouble* fcavp, int i, int direction = 0, bool simple = false);
	void ComputePositiveFlux(int n, const adept::adouble* x, int m, adept::adouble* fcavp, int i, bool simple = false);
	void ComputeNegativeFlux(int n, const adept::adouble* x, int m, adept::adouble* fcavp, int i, bool simple = false);
	void SetFluxes(int i, vector < double >& fluxp, vector < double >& fluxn);
	void SetRHS();
	// virtual void Fluxes();

private:
	vector < vector < double > > ls;
	vector < vector < double > > rs;
	adept::Stack stack;
	// vector < double > diss_blend;		///< dissipation blending coeffs for different stages
	// vector < int > diss_flag;			///< dissipation evaluation (1=yes)

	// double vis2;						///< artificial dissipation coefficient - k2   (central scheme)
	// double vis4;						///< artificial dissipation coefficient - 1/k4 (central scheme)
};
