#include "TimeTest.h"

void Time()
{
	cout << endl << "  ---  Time test  ---  " << endl;

	Loop new_loop;
	DataBase section;

	section.ReadFile("Input/section.txt");

	auto sec_foo = bind(&section_foo, _1, _2);
	new_loop.ReadFile("Input/loop.txt");
	new_loop.ShowData();
	new_loop.SetTime(0.);

	vector < double > parameters, coords;

	coords = new_loop.GetCoordinates();
	parameters.push_back(new_loop.GetTime());
	parameters.insert(parameters.end(), coords.begin(), coords.end());

	new_loop.Refresh("section", sec_foo, parameters);
	//new_loop.Refresh("section", section);
	new_loop.PrintColumn("Output/Section.txt", new_loop.GetSection());
	for (int i = 0; i < 5; i++)
	{
		new_loop.IncrementTime(0.1);
		parameters[0] = new_loop.GetTime();

		new_loop.Refresh("section", sec_foo, parameters);
		//new_loop.Refresh("section", section);
		new_loop.PrintColumn("Output/Section.txt", new_loop.GetSection());
	}

	cout << endl << "   OK" << endl;
}

double section_foo(vector < double > parameters, int i)
{
	return 2. + sin(parameters[0] * 10. - parameters[i + 1] * 2.); // +cos(parameters[i + 1] * 2.);
}