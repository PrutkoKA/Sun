#include "BaseTest.h"

void Base()
{
	cout << endl << "  ---  Base functions test  ---  " << endl;

	Atmosphere new_atm;

	//valarray< double > new_tab;
	vector < vector < double > > new_tab;
	string file_name = "Output\\Inter.txt";

	new_atm.ReadFile("Input/atm.txt");
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

	cout << endl << "   OK   " << endl;
}