#include "DataBase.h"

int DataBase::ReadFile(string file_name)
{
	FILE* file;
	char str_c[12];
	string str_s;
	double value;

	col_num = 0;
	col_size = 0;

	cout << "Opening file " << "\"" << file_name.c_str() << "\"..." << endl; 
	if ((file = fopen(file_name.c_str(), "r")) == NULL) {
		cout << "File \"" << file_name.c_str() << "\" wasn't opened." << endl;
		return -1;
	}
	cout << "File was opened successfuly." << endl;

	str_s = " ";
	while (str_s != "//") {
		fscanf(file, "%s", str_c, sizeof(str_c));
		str_s = str_c;
		if (str_s != "//")
			columns.insert(pair< string, int >(str_s, col_num++) );
	}
	Data.resize(columns.size());
	while (!feof(file)) {
		col_size++;
		for (int i = 0; i < col_num; i++) {
			fscanf(file, "%lf", &value, sizeof(value));
			Data[i].push_back(value);
		}
	}

	return 0;
}

void DataBase::SetData(vector < vector < double > > table)
{
	col_size = table[0].size();
	Data = table;
}

void DataBase::SetData(vector < vector < double > > table, map< string, int > columns_)
{
	columns = columns_;
	col_num = columns.size();
	col_size = table[0].size();
	Data = table;
}

void DataBase::SetRow(string col_name, double value)
{
	vector < double > values(col_size, value);

	SetRow(col_name, values);
}

void DataBase::SetRow(string col_name, vector < double > values)		// Column may be???
{
	Data[columns[col_name]] = values;
}

void DataBase::ShowData(int lines)
{
	if (lines < -1) {
		cout << "Incorrect number of lines to show" << endl;
		return;
	}
	ShowHeader();
	for (int i = 0; i < (lines == -1 ? col_size : min(col_size, lines)); i++) {
		ShowRow(i);
	}
}

void DataBase::ShowData(vector < vector < double > > table, int lines)
{
	if (lines < -1) {
		cout << "Incorrect number of lines to show" << endl;
		return;
	}
	ShowHeader();
	for (int i = 0; i < (lines == -1 ? table[0].size() / col_num : min(int(table[0].size()) / col_num, lines)); i++) {
		ShowRow(i, table);
	}
}

int DataBase::PrintHeader(string file_name)
{
	FILE* file;

	file = fopen(file_name.c_str(), "w");
	if (!file) {
		cout << "File \"" << file_name.c_str() << "\" wasn't opened." << endl;
		return -1;
	}

	for (auto col = columns.begin(); col != columns.end(); ++col) {
		fprintf(file, "%s\t", (col->first).c_str());
	}
	fprintf(file, "\n");
	fclose(file);

	return 0;
}

int DataBase::PrintRow(string file_name, valarray< double > row)
{
	FILE* file;

	file = fopen(file_name.c_str(), "a");
	if (!file) {
		cout << "File \"" << file_name.c_str() << "\" wasn't opened." << endl;
		return -1;
	}

	for (auto col = columns.begin(); col != columns.end(); ++col) {
		fprintf(file, "%lf\t", row[col->second]);
	}
	fprintf(file, "\n");
	fclose(file);

	return 0;
}

int DataBase::PrintRow(string file_name, vector< double > row)
{
	return PrintRow(file_name, valarray< double >(row.data(), row.size()));
}

int DataBase::PrintWholeRow(string file_name, vector< double > row)
{
	return PrintWholeRow(file_name, valarray< double >(row.data(), row.size()));
}

int DataBase::PrintWholeRow(string file_name, valarray< double > row)
{
	FILE* file;

	file = fopen(file_name.c_str(), "a");
	if (!file) {
		cout << "File \"" << file_name.c_str() << "\" wasn't opened." << endl;
		return -1;
	}

	for (unsigned int i = 0; i < row.size(); ++i) {
		fprintf(file, "%lf\t", row[i]);
	}
	fprintf(file, "\n");
	fclose(file);

	return 0;
}

int DataBase::PrintColumn(string file_name, vector< double > column)
{
	FILE* file;

	file = fopen(file_name.c_str(), "a");
	if (!file) {
		cout << "File \"" << file_name.c_str() << "\" wasn't opened." << endl;
		return -1;
	}
	for (auto col : column) {
		fprintf(file, "%lf\n", col);
	}
	fclose(file);
}

int DataBase::PrintTable(string file_name)
{
	if (PrintTable(file_name, Data) != 0)
		return -1;

	return 0;
}

int DataBase::PrintTable(string file_name, vector < vector < double > > table)
{
	if (PrintHeader(file_name) != 0)
		return -1;

	for (int i = 0; i < table[0].size(); i++)
		if (PrintRow(file_name, GetRow(i, table)) != 0)
			return -1;

	return 0;
}

void DataBase::ShowHeader()
{
	for (auto col = columns.begin(); col != columns.end(); ++col) {
		cout << col->first << "\t";
	}
	cout << endl;
}

void DataBase::ShowRow(int row_num)
{
	for (auto col = columns.begin(); col != columns.end(); ++col) {
		cout << Data[col->second][row_num] << "\t";
	}
	cout << endl;
}

void DataBase::ShowRow(int row_num, vector < vector < double > > table)
{
	for (auto col = columns.begin(); col != columns.end(); ++col) {
		cout << table[col->second][row_num] << "\t";
	}
	cout << endl;
}

void DataBase::ShowColumn(string col_name)
{
	ShowColumn(GetValues(col_name));
}

void DataBase::ShowColumn(vector< double > column)
{
	for (auto col : column)
		cout << col << endl;
}

vector< double > DataBase::GetValues(string col_name)
{
	if (columns.find(col_name) == columns.end())
	{
		cout << "There is no column \"" << col_name << "\"." << endl;
		return vector< double >();
	}

	return Data[columns[col_name]];
}

vector < vector < double > > DataBase::GetTable()
{
	return Data;
}

map< string, int > DataBase::GetColumnNames()
{
	return columns;
}

vector < vector < double > > DataBase::NewTable(string col_name, double step, vector < string > ignore_cols, bool INCLUDE_POINTS)
{
	double init_val, fin_val;
	vector < double > range_;

	init_val = GetValues(col_name)[0];
	fin_val = GetValues(col_name)[GetValues(col_name).size() - 1];
	range_ = MakeRange(init_val, fin_val, step);
	return NewTable(col_name, range_, ignore_cols, INCLUDE_POINTS);
}

vector < vector < double > > DataBase::NewTable(string col_name, vector < double > range, vector < string > ignore_cols, bool INCLUDE_POINTS)
{
	vector < vector < double > > temp_v(col_num);
	
	valarray < double > temp_va(col_num);
	int previous_index = 0;
	int index;
	double point_val = GetRow(0)[columns[col_name]];
	
	double val;

	for (int i = 0; i < range.size(); i ++) {
		val = range[i];
		if (INCLUDE_POINTS) {
			index = FindIndexInCol(col_name, val);

			if (index != previous_index) {// && val != point_val) {
				point_val = GetRow(index)[columns[col_name]];
				if (val != point_val) {
					val = point_val;
					i--;
				}
			}
			previous_index = index;
		}
		temp_va = GetFromPoint(col_name, val);
		if (ignore_cols.size() > 0) {
			//index = FindIndexInCol(col_name, val);
			index = i;
			for (auto col : ignore_cols) {
				temp_va[columns[col]] = GetRow(index)[columns[col]];
			}
		}
		for (auto col = columns.begin(); col != columns.end(); ++col)
		{
			temp_v[col->second].push_back(temp_va[col->second]);
		}
	}

	return temp_v;
}

valarray < double > DataBase::GetRow(int row_n)
{
	vector < double > row;

	row.resize(col_num);
	for (auto col = columns.begin(); col != columns.end(); ++col) {
		row[col->second] = Data[col->second][row_n];
	}

	return valarray< double >(row.data(), row.size());
}

valarray < double > DataBase::GetRow(int row_n, vector < vector < double > > table)
{
	vector < double > row;

	row.resize(col_num);
	for (auto col = columns.begin(); col != columns.end(); ++col) {
		row[col->second] = table[col->second][row_n];
	}

	return valarray< double >(row.data(), row.size());
}

int DataBase::AddColumn(string col_name, vector < double > new_column)
{
	if (new_column.size() != col_size)
	{
		cout << "Sizes of columns don't match" << endl;
		return -1;
	}
	columns[col_name] = col_num;
	Data.push_back(vector < double >());
	Data[col_num++] = new_column;

	return 0;
}

valarray< double > DataBase::GetFromPoint(string col_name, double value)
{
	int index = FindIndexInCol(col_name, value);
	valarray < double > im_row, i_row, ip_row, ipp_row, result, result_l, result_r;

	valarray < double > a0(col_num), a1(col_num), a2(col_num);
	double blend;

	i_row = GetRow(index);				//i row
	if (i_row[columns[TYPE_COL]] == -1) {
		i_row[columns[col_name]] = value;
		return i_row;
	}

	if (i_row[columns[col_name]] == value) {
		return i_row;
	}

	if (index + 1 < col_size) {
		ip_row = GetRow(index + 1);		//i+1 row
	}
	else {
		//Something
		return i_row;
	}
	if (i_row[columns[TYPE_COL]] == 1 && ip_row[columns[TYPE_COL]] == 1)
	{
		blend = (value - i_row[columns[col_name]]) / (ip_row[columns[col_name]] - i_row[columns[col_name]]);
		if (index > 0) {
			im_row = GetRow(index - 1);	//i-1 row
			if (i_row[columns[col_name]] != im_row[columns[col_name]]) {
					a2 = (ip_row - im_row) /
					((ip_row[columns[col_name]] - im_row[columns[col_name]]) * (ip_row[columns[col_name]] - i_row[columns[col_name]])) -
					(i_row - im_row) /
					((i_row[columns[col_name]] - im_row[columns[col_name]]) * (ip_row[columns[col_name]] - i_row[columns[col_name]]));
				a1 = (i_row - im_row) / (i_row[columns[col_name]] - im_row[columns[col_name]]) -
					a2 * (i_row[columns[col_name]] + im_row[columns[col_name]]);
				a0 = im_row - a1 * im_row[columns[col_name]] - a2 * pow(im_row[columns[col_name]], 2);
				result_l = a0 + a1 * value + a2 * pow(value, 2);
			} else {
				result_l = valarray< double >(0., col_num);
				blend = 1;
			}
		} else {
			result_l = valarray< double >(0., col_num);
			blend = 1;
		}
		if (index < col_size - 2) {
			ipp_row = GetRow(index + 2);	//i+2 row
			if (ipp_row[columns[col_name]] != ip_row[columns[col_name]]) {
				a2 = (ipp_row - i_row) /
					((ipp_row[columns[col_name]] - i_row[columns[col_name]]) * (ipp_row[columns[col_name]] - ip_row[columns[col_name]])) -
					(ip_row - i_row) /
					((ip_row[columns[col_name]] - i_row[columns[col_name]]) * (ipp_row[columns[col_name]] - ip_row[columns[col_name]]));
				a1 = (ip_row - i_row) / (ip_row[columns[col_name]] - i_row[columns[col_name]]) -
					a2 * (ip_row[columns[col_name]] + i_row[columns[col_name]]);
				a0 = i_row - a1 * i_row[columns[col_name]] - a2 * pow(i_row[columns[col_name]], 2);
				result_r = a0 + a1 * value + a2 * pow(value, 2);
			} else {
				result_r = valarray< double > (0., col_num);
				blend = 0;
			}
		} else {
			result_r = valarray< double > (0., col_num);
			blend = 0;
		}
		result = (1. - blend) * result_l + blend * result_r;
	}
	else {
		a2 = 0;
		a1 = (ip_row - i_row) / (ip_row[columns[col_name]] - i_row[columns[col_name]]);
		a0 = i_row - a1 * i_row[columns[col_name]];
		result = a0 + a1 * value + a2 * pow(value, 2);
	}
	
	return result;
}

int DataBase::FindIndexInCol(string col_name, double value)
{
	if (value <= Data[columns[col_name]][0])
		return 0;

	if (value >= Data[columns[col_name]][col_size - 1])
		return col_size - 1;

	int index = 0;
	for (auto col_val : Data[columns[col_name]])
	{
		if (col_val > value)
			return index - 1;

		index++;
	}

	return -1;
}

vector < double > DataBase::MakeRange(double start, double end, double step)
{
	vector < double > temp;

	temp.resize((end - start) / step + 1);
	for (int i = 0; i < temp.size(); i++)
		temp[i] = start + i * step;

	return temp;
}

void DataBase::InverseData(string col_name)
{
	double swap;
	for (auto col = columns.begin(); col != columns.end(); ++col)
	{
		if (col->first == col_name)
			continue;

		for (int i = 0; i < col_size / 2; ++i)
		{
			swap = Data[col->second][i];
			Data[col->second][i] = Data[col->second][col_size - i - 1];
			Data[col->second][col_size - i - 1] = swap;
		}
	}
}
