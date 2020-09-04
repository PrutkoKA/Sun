#pragma once
#include "DataBase.h"

class Atmosphere : public DataBase {
public:

	vector< double > GetHeights();
	vector< double > GetTemperatures();
private:
	const string HEIGHT_COL = "height";
	const string TEMP_COL = "temperature";
};