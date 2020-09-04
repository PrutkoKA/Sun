// Sun.cpp: определяет точку входа для консольного приложения.
//

//#include "stdafx.h"
#include "Atmosphere.h"
#include "Loop.h"
#include <omp.h>
#include <functional>

using namespace std::placeholders;

double section_foo(vector < double > parameters, int i);

int main()
{
	Atmosphere new_atm;
	Loop new_loop;
	DataBase section;

	valarray< double > new_tab;
	string file_name = "Output\\Inter.txt";

	new_atm.ReadFile("Input\\atm.txt");
	new_atm.ShowData(5);
	//new_atm.ShowHeader();
	//new_atm.PrintHeader(file_name);

	vector < double > ran1 = new_atm.MakeRange(1., 4.95, 0.05);
	vector < double > ran2 = new_atm.MakeRange(5., 13.46, .25);
	vector < string > ignore;

	ignore.push_back(new_atm.TYPE_COL);
	ran1.insert(ran1.end(), ran2.begin(), ran2.end());
	new_tab = new_atm.NewTable("height", ran1, ignore, true);
	//new_tab = new_atm.NewTable("height", 0.2, true);

	new_atm.PrintTable(file_name, new_tab);
	new_atm.ShowData(new_tab, 5);
	new_atm.SetData(new_tab);
	new_atm.ShowData(5);

	section.ReadFile("Input\\section.txt");

	auto sec_foo = bind(&section_foo, _1, _2);
	new_loop.ReadFile("Input\\loop.txt");
	new_loop.ShowData();
	new_loop.SetTime(0.);

	vector < double > parameters, coords;

	coords = new_loop.GetCoordinates();
	parameters.push_back(new_loop.GetTime());
	parameters.insert(parameters.end(), coords.begin(), coords.end());

	//new_loop.Refresh("section", sec_foo, parameters);
	new_loop.Refresh("section", section);
	new_loop.PrintColumn("Output\\Section.txt", new_loop.GetSection());
	for (int i = 0; i < 5; i++)
	{
		new_loop.IncrementTime(0.1);
		parameters[0] = new_loop.GetTime();

		//new_loop.Refresh("section", sec_foo, parameters);
		new_loop.Refresh("section", section);
		new_loop.PrintColumn("Output\\Section.txt", new_loop.GetSection());
	}
	

	/*cout << "Time " << new_loop.GetTime() << endl;
	cout << "section" << endl;
	new_loop.ShowColumn("section");*/


	

	/*cout << "Time " << new_loop.GetTime() << endl;
	cout << "section" << endl;
	new_loop.ShowColumn("section");*/

	/*cout << "coordinates" << endl;
	for (auto coord : new_loop.GetCoordinates())
		cout << coord << endl;*/
	//for (double h = new_atm.GetHeights()[0]; h < new_atm.GetHeights()[new_atm.GetHeights().size()-1]; h += 0.1) {
	//	//new_atm.ShowRow(new_atm.GetFromPoint("height", 3.5));
	//	new_atm.PrintRow(file_name, new_atm.GetFromPoint("height", h));
	//}
	/*cout << "heights" << endl;
	for (auto height : new_atm.GetHeights())
		cout << height << endl;*/

	/*cout << "temperatures" << endl;
	for (auto temperature : new_atm.GetTemperatures())
		cout << temperature << endl;

	cout << "some" << endl;
	for (auto some : new_atm.GetValues("some"))
		cout << some << endl;

	cout << "density" << endl;
	for (auto density : new_atm.GetValues("density"))
		cout << density << endl;*/

    return 0;
}

double section_foo(vector < double > parameters, int i)
{
	return 2. + sin(parameters[0] * 10. - parameters[i + 1] * 2.); // +cos(parameters[i + 1] * 2.);
}