// Central Differences Scheme
#pragma once

#include "Solver.h"

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
	virtual void Fluxes();

private:
	vector < double > diss_blend;		///< dissipation blending coeffs for different stages
	vector < int > diss_flag;			///< dissipation evaluation (1=yes)

	double vis2;						///< artificial dissipation coefficient - k2   (central scheme)
	double vis4;						///< artificial dissipation coefficient - 1/k4 (central scheme)
};
