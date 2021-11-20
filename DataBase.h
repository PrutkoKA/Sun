#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <valarray>
#include <string>
//#include <conio.h>
#include <algorithm>


using namespace std;

//struct AtmData {
//	double height;
//	double temperature;
//	double density;
//	// Other parameters;
//};

/*!
	\brief Main database class

	Description of entities

*/
class DataBase {
public:
	vector < vector < double > > Data;
	//valarray < double > Data_;
	map< string, int > columns;
	int col_num;	///< Number of columns
	int col_size;	///< Number of elements in a column
	const string TYPE_COL = "type";

	int ReadFile(string file_name);
	void SetData(vector < vector < double > > table);
	void SetData(vector < vector < double > > table, map< string, int > columns_);
	void SetRow(string col_name, double value);
	void SetRow(string col_name, vector < double > values);
	void ShowData(int lines = -1);

	void ShowData(vector < vector < double > > table, int lines = -1);
	int PrintHeader(string file_name);
	int PrintRow(string file_name, valarray< double > row);
	int PrintRow(string file_name, vector< double > row);
	int PrintWholeRow(string file_name, valarray< double > row);
	int PrintWholeRow(string file_name, vector< double > row);
	int PrintColumn(string file_name, vector< double > column);
	int PrintTable(string file_name);
	int PrintTable(string file_name, vector < vector < double > > table);
	void ShowHeader();
	void ShowRow(int row_num);
	void ShowRow(int row_num, vector < vector < double > > table);
	void ShowColumn(string col_name);
	void ShowColumn(vector< double > column);
	vector< double > GetValues(string col_name);
	bool ColumnExists(string col_name) { return !(columns.find(col_name) == columns.end()); };

	vector < vector < double > > GetTable();
	map< string, int > GetColumnNames();

	vector < vector < double > > NewTable(string col_name, double step,
		vector < string > ignore_cols = vector < string >(), bool INCLUDE_POINTS = false);

	vector < vector < double > > NewTable(string col_name, vector < double > range,
		vector < string > ignore_cols = vector < string >(), bool INCLUDE_POINTS = false);

	valarray < double > GetRow(int row_n);
	valarray < double > GetRow(int row_n, vector < vector < double > > table);
	int AddColumn(string col_name, vector < double > new_column);
	valarray< double > GetFromPoint(string col_name, double value);
	int FindIndexInCol(string col_name, double value);
	vector < double > MakeRange(double start, double end, double step);

	void InverseData(string col_name);
private:
	/*const string TYPE_COL = "type";*/

	//vector < double > Data;
	/*valarray < double > Data_;*/
	/*map< string, int > columns;*/
	/*int col_num, col_size;*/
};
