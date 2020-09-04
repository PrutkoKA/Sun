#include "Atmosphere.h"

vector< double > Atmosphere::GetHeights()
{
	return GetValues(HEIGHT_COL);
}

vector< double > Atmosphere::GetTemperatures()
{
	return GetValues(TEMP_COL);
}