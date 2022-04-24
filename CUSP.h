#pragma once

#include "Solver.h"
//#include "adept/adept_source.h"

// struct cds_struct
// {
// 	vector < double > diss_blend;
// 	vector < int > diss_flag;

// 	double vis2;
// 	double vis4;
// };

class CUSP : public Solver {
public:
	CUSP(sol_struct& sol_init_);
	virtual void LRState();
	virtual void LRState(string var_);
	virtual void LRState(vector < vector < double > >& cv_, vector < vector < double > >& ls_, vector < vector < double > >& rs_);
	virtual void Fluxes();
	double CUSPLimiter(double af, double bf);
	// virtual double Solve();
	virtual const vector < int >& GetDissFlag();
	virtual const vector < double >& GetDissBlend();
	virtual void Dissipation(double beta);

	virtual void RHS(int i);
	virtual void ComputeRHSandJacobian(bool NO_JAC = false);
	void GetPositiveFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac);
	void GetNegativeFluxAndJacobian(int i, vector < double >& y_val, vector < double >& jac);
	void ComputePositiveFlux(int n, const double_type* x, int m, double_type* fcavp, int i);
	void ComputeNegativeFlux(int n, const double_type* x, int m, double_type* fcavp, int i);
	void ComputeFlux(int n, const double_type* x, int m, double_type* fcav, int i, int direction = 0);
	void SetFluxes(int i, vector < double >& fluxp, vector < double >& fluxn);
	void GetSourceAndJacobian(int i, vector < double >& y_val, vector < double >& jac);
	void ComputeSourceTerm(int n, const double_type* cv, int m, double_type* source, int i);
	void SetRHS();

	virtual void GetFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac, bool POS_NEG, bool simple = false);
	// virtual void Fluxes();
	virtual void deactivate_adept_stack()
	{
		stack.deactivate();
	};

private:
	vector < vector < double > > ls;
	vector < vector < double > > rs;
	adept::Stack stack;
	// vector < double > diss_blend;		///< dissipation blending coeffs for different stages
	// vector < int > diss_flag;			///< dissipation evaluation (1=yes)

	// double vis2;						///< artificial dissipation coefficient - k2   (central scheme)
	// double vis4;						///< artificial dissipation coefficient - 1/k4 (central scheme)
};
