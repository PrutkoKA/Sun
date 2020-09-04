#pragma once
#include <functional>
#include "DataBase.h"

class Loop : public DataBase {
public:

	vector< double > GetCoordinates();
	vector< double > GetSection();
	void SetTime(double time_);
	void IncrementTime(double dt);
	double GetTime();
	void Refresh(string col_name, function<double(vector < double >, int)> foo, vector < double > parameter);
	void Refresh(string col_name, DataBase db);
private:
	const string COORD_COL = "coordinate";
	const string SECT_COL = "section";
	const string TIME_COL = "time";
	double time;
};
