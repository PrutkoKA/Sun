#include "Loop.h"

vector< double > Loop::GetCoordinates()
{
	return GetValues(COORD_COL);
}

vector< double > Loop::GetSection()
{
	return GetValues(SECT_COL);
}

void Loop::SetTime(double time_)
{
	time = time_;
}

void Loop::IncrementTime(double dt)
{
	time += dt;
}

double Loop::GetTime()
{
	return time;
}

void Loop::Refresh(string col_name, function<double(vector < double >, int)> foo, vector < double > parameters)
{
	for (int i = 0; i < col_size; i++)
		//Data_[columns[col_name] + i * col_num] = foo(parameters, i);
		Data[columns[col_name]][i] = foo(parameters, i);
}

void Loop::Refresh(string col_name, DataBase db)
{
	double val = db.GetFromPoint(TIME_COL, time)[db.columns[col_name]];

	SetRow(col_name, val);
}