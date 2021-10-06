#include "GridTest.h"

void Grid()
{
	cout << endl << "  ---  Grid test  ---  " << endl;

	Loop grid;

	grid.ReadFile("Input/mesh.txt");
	grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("x"));

	vector < double > parameters, coords;

	coords = grid.GetValues("x");
	parameters.push_back(0.);
	parameters.insert(parameters.end(), coords.begin(), coords.end());

	grid.Refresh("value", sharp_foo, parameters);

	//grid.PrintColumn("Output/value.txt", grid.GetValues("value"));

	grid.CalculateResolution(1., 1., "value", "x");
	grid.CalculateConcentration(1., "x");
	grid.AddCollumnR();
	grid.AddCollumnN();

	grid.PrintTable("Output/new_mesh_0.txt");

	grid.CalculateConcentrationWave();
	grid.FindMesh(0.1, "x");
	//grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("x"));
	//grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("x"));

	for (int i = 0; i < 100; ++i) {
		cout << i + 1 << "  ";

		grid.RefreshX("x");

		coords = grid.GetValues("x");
		parameters.clear();
		parameters.push_back(0.);
		parameters.insert(parameters.end(), coords.begin(), coords.end());
		grid.Refresh("value", sharp_foo, parameters);
		//grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("value"));
		grid.CalculateResolution(1., 1., "value", "x");
		grid.CalculateConcentration(1., "x");

		grid.RefreshR("R");
		grid.RefreshN("n");

		char num[4];
		_itoa(i + 1, num, 10);
		grid.PrintTable(string("Output/new_mesh_") + num + ".txt");

		grid.CalculateConcentrationWave();
		grid.FindMesh(0.1, "x");
		//grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("x"));
	}

}

double sharp_foo(vector < double > parameters, int i)
{
	return 0.5 * (1. + tanh(1e3 * (parameters[i + 1] - 0.4))) * exp(-pow((parameters[i + 1] - 0.4) / 0.2, 2.));
}

double sharp_rho(vector < double > parameters, int i)
{
	if (parameters[i] <= 0.3) {
		return 1.;
	}
	else {
		return 0.125;
	}
}

void SodGrid()
{
	cout << endl << "  ---  SodGrid test  ---  " << endl;

	double dt = 0.1;

	Loop grid;

	grid.ReadFile("Input/sod_start.txt");
	grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("x"));

	vector < double > parameters, coords;

	coords = grid.GetValues("x");
	parameters.push_back(0.);
	parameters.insert(parameters.end(), coords.begin(), coords.end());

	grid.Refresh("rho", sharp_rho, parameters);

	//grid.PrintColumn("Output/value.txt", grid.GetValues("value"));

	grid.CalculateResolution(1., 1., "rho", "x");
	grid.CalculateConcentration(1., "x");
	grid.AddCollumnR();
	grid.AddCollumnN();

	grid.PrintTable("Output/new_mesh_0.txt");

	grid.CalculateConcentrationWave();
	grid.FindMesh(dt, "x");
	//grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("x"));
	//grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("x"));

	for (int i = 0; i < 100; ++i) {
		cout << i + 1 << "  ";

		grid.RefreshX("x");

		coords = grid.GetValues("x");
		parameters.clear();
		parameters.push_back(0.);
		parameters.insert(parameters.end(), coords.begin(), coords.end());
		grid.Refresh("rho", sharp_rho, parameters);
		//grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("value"));
		grid.CalculateResolution(1., 1., "rho", "x");
		grid.CalculateConcentration(1., "x");

		grid.RefreshR("R");
		grid.RefreshN("n");

		char num[4];
		_itoa(i + 1, num, 10);
		//grid.PrintTable(string("Output/new_mesh_") + num + ".txt");

		grid.CalculateConcentrationWave();
		grid.FindMesh(dt, "x");
		//grid.PrintColumn("Output/GridMoving.txt", grid.GetValues("x"));
	}
	grid.PrintTable(string("Output/new_mesh_sod") + ".txt");
}