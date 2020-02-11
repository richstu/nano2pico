//Check if your run is inJSON by calling bool inJSON(VRunLumi,run,lumiblock) in the event loop..
//ported from babymaker

#ifndef H_IN_JSON
#define H_IN_JSON

#include <string>
#include <vector>

/**
 * MakeVRunLumi - reads golden JSON and creates a matrix of integers needed for inJSON to determine 
 * if a given event is deemed good
 * input is either "golden201x" or the path to the golden cert file
 */
std::vector< std::vector<int> > MakeVRunLumi(std::string input);

/**
 * inJSON - returns true if an event in a given run and lumiblock is deemed good by the golden JSON 
 * and false otherwise
 */
bool inJSON(std::vector< std::vector<int> > VVRunLumi, int Run, int LumiBlock);

#endif
