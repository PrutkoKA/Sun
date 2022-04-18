#pragma once

//#include "Loop.h"
#include "Solver.h"

//double gamma;
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
//
//double volref;
//double rhoref;
//double uref;
//double pref;

void Laval2();

void Laval3();

void Laval4();

void Laval5();

void Laval6();

void Laval7();

void unsteady_sod_test(const string &output_file, const string& yml_file, const double end_time_ = 0.02);

void init_sun_atmosphere(Solver* solver);
void loop_foot_point(const string& output_file, const string& yml_file, const double end_time_);	// Sun loop

void compare_files(const string& first_file, const string& second_file, const double eps = 1e-6);

void Laval();

void bcond();

void solver(int imax, int ib2, int mxdum, vector < double >& l_coords, vector < double >& l_sections, vector < double >& l_volumes, 
			vector < vector < double > >& cv, vector < double >& p, vector < vector < double > >& cvold,
			vector < vector < double > >& diss, vector < vector < double > >& rhs, vector < double >& dt, vector < double >& ls, 
			vector < vector < double > >& rs, vector < double >& dum);

void Tstep(int imax, int ib2, vector < double >& l_coords, vector < double >& l_sections, vector < double >& l_volumes,
	vector < vector < double > >&cv, vector < double >& p, vector < double >& dt);

void Dissip(int imax, int ib2, double beta, vector < double >& l_volumes, vector < vector < double > >& cv,
	vector < double >& p, vector < double >& dt, double* dp, double* d, vector < vector < double > >& diss);

void Flux_cen(int imax, int ib2, vector < double >& l_sections, vector < vector < double > >& cv, 
				vector < double >& p, vector < vector < double > >& diss, vector < double >& f, 
				vector < vector < double > >& rhs);

void Srcterm(int imax, int ib2, vector < double >& l_sections, vector < double >& p, 
				vector < vector < double > >& rhs);

void Irsmoo(int imax, int ib2, vector < vector < double > >& rhs, vector < double >& d);

void Conver(int max, int ib2, vector < double >& l_sections, vector < vector < double > >& cv, 
			vector < vector < double > >& cvold, vector < double >& p);

void Output(int imax, int ib2, vector < double >& l_coords, vector < double >& l_sections, 
			vector < vector < double > >& cv, vector < double >& p);

void Wsolut(int imax, vector < vector < double > >& cv, vector < double >& p);
