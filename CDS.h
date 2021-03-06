// Central Differences Scheme
#pragma once

#include "Solver.h"
//#include "adept/adept_source.h"

struct cds_struct
{
	vector < double > diss_blend;
	vector < int > diss_flag;

	double vis2;
	double vis4;
};

class CDS : public Solver {
public:
	CDS(sol_struct& sol_init_, cds_struct& cds_init_);
	// virtual double Solve();
	virtual const vector < int >& GetDissFlag();
	virtual const vector < double >& GetDissBlend();
	virtual void Dissipation(double beta);
	virtual void LRState();
	virtual void LRState(string var_);
	virtual void LRState(vector < vector < double > >& cv_, vector < vector < double > >& ls_, vector < vector < double > >& rs_);
	virtual void Fluxes();
	virtual void RHS(int i);
	virtual void ComputeRHSandJacobian(bool NO_JAC = false);
	virtual void GetFluxAndJacobian(int i, vector < double >& y_val, vector < vector < double > >& cv_, vector < double >& jac, bool POS_NEG, bool simple = false);
	virtual void deactivate_adept_stack()
	{
		// do nothing
	};

private:
	vector < double > diss_blend;		///< dissipation blending coeffs for different stages
	vector < int > diss_flag;			///< dissipation evaluation (1=yes)

	double vis2;						///< artificial dissipation coefficient - k2   (central scheme)
	double vis4;						///< artificial dissipation coefficient - 1/k4 (central scheme)
};
