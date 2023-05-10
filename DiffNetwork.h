#include<iostream>

#include"Driver.h"
#include <vector>
#include <ctime>
#include <algorithm>
#include <numeric>
#ifndef PAIRDRIVER_DIFFNETWORK_H
#define PAIRDRIVER_DIFFNETWORK_H

#endif //PAIRDRIVER_DIFFNETWORK_H

struct SolfThreshold {
    map<int, vector<double>> corrMatrix;
    vector<double> connectivity;
};
struct SolfThreshold1 {
    vector<vector<double>> corrMatrix;
    vector<double> connectivity;
};
struct Expression {
    map<string, map<string, double>> m;
    vector<string> p;
};
struct Diff {
    map<int, vector<double>> m;
    vector<string> g;
    vector<string> bp;
};
struct Net {
    map<string, double> geneMaxWeight;
    map<string, double> geneWeight;
    vector<string> genesInNet;
};
struct PPI {
    int Id;
    double Id_weight;
    PPI(int _Id, double _Id_weight) {
        this->Id = _Id;
        this->Id_weight = _Id_weight;
    }
};

Expression ReadExpressionMatrix(string matrix_path);
Expression ReadExpressionMatrixWithPreGenes(string matrix_path, unordered_set<string> preGenes);
Expression ReadExpressionMatrixForRandom(string matrix_path, float keepRate);
void ReadPreSamples(string expTumorPath, string expNormalPath);
//void ReadExpressionMatrix1(string expTumorPath, string expNormalPath, unordered_set<string> preGenes);
map<int, vector<string>> ReadDiffMatrix(string matrix_path);
map<int, vector<double>> GetDoubleDiffMatrix(map<int ,vector<string>> diffMatrixRaw);
vector<string> GetMajorGenes(map<int ,vector<string>> diffMatrixRaw);
map<int, vector<double>> GetCorrMatrix (map<int, vector<double>> diffMatrix, int k);
vector<double> GetConnectivity (map<int, vector<double>> corrMatrix);
double InnerConnectValue(vector<double> vec1, int pos1, vector<double> vec2, int pos2);
double OrdinaryLeastSquare(const vector<double>& x, const vector<double>& y);
SolfThreshold GetOptimalCorrMatrix (map<int, vector<double>> diffMatrix, float threshold=0.9, int maxPower=6);
SolfThreshold1 GetOptimalCorrMatrix1 (map<int, vector<double>> diffMatrix, float threshold=0.90, int maxPower=6);
//map<string, map<string, double>> GetFinalMatrix (map<int, vector<double>> corrMatrix, vector<double> connectivity, vector<string> majorGenes, map<string, int> mapOfficialToNCBI);
void WriteDiffNetwork(string writeDNPath, map<string, map<string, double>> finalMatrix, vector<string> majorGenes);
//Diff GetDiffMatrix(string mutationPath, string netPath1, string geneIDPath, string expTumorPath, string expNormalPath, unordered_set<string> preGenes);
void WriteDiffMatrix(string diffMatrixPath, map<int, vector<double>> diffMatrix, vector<string> patients, vector<string> genes);