#pragma once

#include "Solver.h"

class HLLE : public Solver {
public:
	HLLE(sol_struct& sol_init_);
	virtual void LRState() {};
	virtual void LRState(string var_);
	virtual void LRState(vector < vector < double > >& cv_, vector < vector < double > >& ls_, vector < vector < double > >& rs_) {};
	virtual void Fluxes() {};
	double MUSCL0(double a, double b, double dx);
	double MUSCL0_new(double a, double b, double dx);
	double MUSCL1_3(double a, double b, double dx);
	//virtual double Solve() { return 0.; };
	virtual const vector < int >& GetDissFlag();
	virtual const vector < double >& GetDissBlend();
	virtual void Dissipation(double beta);

	virtual void RHS(int i) {};
	virtual void ComputeRHSandJacobian(bool NO_JAC = false);
	void GetPositiveFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac);
	void GetNegativeFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac);
	void ComputePositiveFlux(int n, const adept::adouble* x, int m, adept::adouble* fcavp, int i);
	void ComputeNegativeFlux(int n, const adept::adouble* x, int m, adept::adouble* fcavp, int i);
	void ComputeFlux(int n, const adept::adouble* x, int m, adept::adouble* fcav, int i, int direction = 0);
	void SetFluxes(int i, vector < double >& fluxp, vector < double >& fluxn);
	void GetSourceAndJacobian(int i, vector < double >& y_val, vector < double >& jac);
	void ComputeSourceTerm(int n, const adept::adouble* cv, int m, adept::adouble* source, int i);
	void SetRHS();

	virtual void GetFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac, bool POS_NEG, bool simple = false);

private:
	vector < vector < double > > ls;
	vector < vector < double > > rs;
	adept::Stack stack;
};