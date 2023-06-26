#include<iostream>

//#include"Driver.h"
#include <vector>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <unordered_set>


#ifndef PAIRDRIVER_PREPROCESS_H
#define PAIRDRIVER_PREPROCESS_H

#endif //PAIRDRIVER_PREPROCESS_H

void StringSplit(string long_str, string sep, vector<string>& vec_str);

unordered_set<string> ReadMutGenes(string mutationPath);
unordered_set<string> ReadExpGenes(string expNormalPath, string expTumorPath);
unordered_set<string> ReadNetGenes(string networkPath, map<int, string> mapNCBIToOfficial);





