#include <iostream>
#include <vector>
#include "Driver.h"
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <time.h>
#include "Preprocess.h"
#include <mutex>
#include <thread>

std::mutex mtx;  // 定义互斥锁

int Rand(int i) { return rand()%i;}

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
//struct Diff {
//    map<int, vector<double>> m;
//    vector<string> g;
//    vector<string> bp;
//};
struct Net {
    map<string, double> geneMaxWeight;
    map<string, double> geneWeight;
    vector<string> genesInNet;

};



map<string, string> GetMapIdName(string geneIDPath) {
    map<string, string> mapIdName;
    ifstream ifs;
    ifs.open(geneIDPath.c_str());
    string str_temp;
    vector<string> vec_string;
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        mapIdName[vec_string[1]] = vec_string[0];
    }
    return mapIdName;
}

Expression ReadExpressionMatrix(string matrix_path) {
    Expression exp;
    ifstream ifs;
    ifs.open(matrix_path.c_str());
    string str_temp;
    vector<string> vec_string, patients;
    map<string, map<string, double>> expMatrix;

    getline(ifs, str_temp);
    StringSplit(str_temp, STR_SEPTABLE, vec_string);
    for (int i = 0; i < vec_string.size() - 1; i++) patients.push_back(vec_string[i].substr(0, 12));
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        vector<double> vec_double;
        map<string, double> mapRow;
        for (int i = 1; i < vec_string.size() - 1; i++)  {
            double value = log2(stof(vec_string[i]) + 1);
            mapRow[patients[i - 1]] = value;
            vec_double.push_back(value);
        }
        expMatrix[vec_string[0]] = mapRow;

    }
    exp.m = expMatrix;
    exp.p = patients;
    return exp;
}

Expression ReadExpressionMatrixWithPreGenes(string matrix_path, unordered_set<string> preGenes) {
    Expression exp;
    ifstream ifs;
    ifs.open(matrix_path.c_str());
    string str_temp;
    vector<string> vec_string, patients;
    map<string, map<string, double>> expMatrix;

    getline(ifs, str_temp);
    StringSplit(str_temp, STR_SEPTABLE, vec_string);
    for (int i = 0; i < vec_string.size() - 1; i++) patients.push_back(vec_string[i].substr(0, 12));
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        vector<double> vec_double;
        map<string, double> mapRow;
//        auto iter = preGenes.find(vec_string[0]);
        if (preGenes.count(vec_string[0])) {
            for (int i = 1; i < vec_string.size() - 1; i++)  {
                double value = log2(stof(vec_string[i]) + 1);
                mapRow[patients[i - 1]] = value;
                vec_double.push_back(value);
            }
            expMatrix[vec_string[0]] = mapRow;
        }
    }
    exp.m = expMatrix;
    exp.p = patients;
    return exp;
}

Expression ReadExpressionMatrixForRandom(string matrix_path, float keepRate) {
    Expression exp;
    ifstream ifs;
    ifs.open(matrix_path.c_str());
    string str_temp;
    vector<string> vec_string, patients, rawPatients;
    map<string, map<string, double>> expMatrix;

    getline(ifs, str_temp);
    StringSplit(str_temp, STR_SEPTABLE, vec_string);
    for (int i = 0; i < vec_string.size() - 1; i++) {
        rawPatients.push_back(vec_string[i].substr(0, 12));
        patients.push_back(vec_string[i].substr(0, 12));
    }

    srand(unsigned(time(0)));
    random_shuffle(patients.begin(), patients.end(), Rand);
    vector<string> randomPatients;
    for (int i = 0; i < (int)(patients.size() * keepRate); i++) {
        randomPatients.push_back(patients[i]);
    }
    vector<int> posRandomPatients;
    for (auto & patient: randomPatients) {
        auto iter = find(rawPatients.begin(), rawPatients.end(), patient);
        posRandomPatients.push_back(distance(rawPatients.begin(), iter) + 1);
    }
    std::sort(posRandomPatients.begin(), posRandomPatients.end());
    for (auto &pos:posRandomPatients) cout << pos - 1 << ' ';
    cout << endl;
    vector<string> cisPatients;
    for (auto & pos:posRandomPatients) {
        cisPatients.push_back(rawPatients[pos - 1]);
    }

    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        vector<double> vec_double;
        map<string, double> mapRow;

        for (int i = 1; i < vec_string.size() - 1; i++)  {
            auto iter = find(posRandomPatients.begin(), posRandomPatients.end(), i);
            if (iter != posRandomPatients.end()) {
                double value = log2(stof(vec_string[i]) + 1);
                mapRow[patients[i - 1]] = value;
                vec_double.push_back(value);
            }
        }
        expMatrix[vec_string[0]] = mapRow;

    }
    exp.m = expMatrix;
    exp.p = cisPatients;
    return exp;
}

//void ReadPreSamples(string expTumorPath, string expNormalPath, vector<int>& posTumor, vector<int>& posNormal) {
//    ifstream ifs, ifs1;
//
//    ifs.open(expTumorPath.c_str());
//    string str_temp;
//    vector<string> vec_string, tumorPatients, normalPatients;
//    // 读取第一行样本名称；
//    getline(ifs, str_temp);
//    StringSplit(str_temp, STR_SEPTABLE, vec_string);
//    for (int i = 0; i < vec_string.size() - 1; i++) tumorPatients.push_back(vec_string[i].substr(0, 12));
//
//    ifs1.open(expNormalPath.c_str());
//    // 读取第一行样本名称；
//    getline(ifs1, str_temp);
//    StringSplit(str_temp, STR_SEPTABLE, vec_string);
//    for (int i = 0; i < vec_string.size() - 1; i++) normalPatients.push_back(vec_string[i].substr(0, 12));
//
//    for (int i = 0; i < normalPatients.size(); i++) {
//        auto iter = find(tumorPatients.begin(), tumorPatients.end(), normalPatients[i]);
//        if (iter != tumorPatients.end()) {
//            posTumor.push_back(distance(tumorPatients.begin(), iter));
//            posNormal.push_back(i);
//        }
//    }
//}

//void ReadExpressionMatrix1(string expTumorPath, string expNormalPath, unordered_set<string> preGenes) {
//    // 初始化矩阵；
//    vector<int> posTumor, posNormal;
//    ReadPreSamples(expTumorPath, expNormalPath, posTumor, posNormal);
//    Matrix tumorM(preGenes.size(), posNormal.size());
//    // 读取Tumor；
//    ifstream ifs;
//    ifs.open(expTumorPath.c_str());
//    string str_temp;
//    vector<string> vec_string, patients;
//    // 跳过第一行；
//    getline(ifs, str_temp);
//    // 开始获取矩阵数据;
//    int row = 0;
//    while (getline(ifs, str_temp)) {
//        StringSplit(str_temp, STR_SEPTABLE, vec_string);
//        if (preGenes.count(vec_string[0])) {
//            for (int i = 1; i < vec_string.size() - 1; i++)  {
//                tumorM(row, i - 1) = log2(stof(vec_string[i]) + 1);
//                for (int j = 0; j < posTumor.size(); j++) {
//                    tumorM(row, j) = log2(stof(vec_string[posTumor[j] + 1]) + 1);
//                }
//            }
//            row += 1;
//            if (row % 1000 == 0) cout << row << endl;
//        }
//    }
//}

map<int, vector<string>> ReadDiffMatrix(string matrix_path) {
    ifstream ifs;
    ifs.open(matrix_path.c_str());
    string str_temp;
    vector<string> vec_string;
    map<int, vector<string>> diffMatrix;
    int bin_index = 0;
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        diffMatrix[bin_index] = vec_string;
        bin_index ++;
    }
    return diffMatrix;
}

vector<string> GetMutationGenes(string mutationPath) {
    ifstream ifs;
    ifs.open(mutationPath.c_str());
    string str_temp;
    vector<string> vec_string, mutationGenes;
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        if (vec_string[0].substr(0, 4) == "TCGA") continue;
        vector<double> vec_float;
        for (int i = 1; i < vec_string.size() - 1; i++) vec_float.push_back(stoi(vec_string[i]));
        if (count(vec_float.begin(), vec_float.end(), 1)) mutationGenes.push_back(vec_string[0]);
    }
    return mutationGenes;
}

vector<string> vectors_intersection(vector<string> v1, vector<string> v2) {
    vector<string> v;
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));//求交集
    return v;
}


Net GetGeneMaxWeight(string netPath1, string geneIDPath) {
    Net net;
    map<string, string> mapIdName = GetMapIdName(geneIDPath);
    ifstream ifs;
    ifs.open(netPath1.c_str());
    string str_temp;
    vector<string> vec_string, netGenes1, netGenes2;
    map<string, double> mapGeneMaxWeight, mapGeneWeight;
    ofstream ofs;
    ofs.open("C:/Users/63121/CLionProjects/PairDriver/GenesMaxWeight.txt");
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        // MaxWeight;
        if (mapGeneMaxWeight.find(mapIdName[vec_string[0]]) != mapGeneMaxWeight.end()) {
            if (atof(vec_string[2].c_str()) > mapGeneMaxWeight[mapIdName[vec_string[0]]]) mapGeneMaxWeight[mapIdName[vec_string[0]]] = atof(vec_string[2].c_str());
        } else {
            mapGeneMaxWeight[mapIdName[vec_string[0]]] = atof(vec_string[2].c_str());
            netGenes1.push_back(mapIdName[vec_string[0]]);
        }
        if (mapGeneMaxWeight.find(mapIdName[vec_string[1]]) != mapGeneMaxWeight.end()) {
            if (atof(vec_string[2].c_str()) > mapGeneMaxWeight[mapIdName[vec_string[0]]]) mapGeneMaxWeight[mapIdName[vec_string[1]]] = atof(vec_string[2].c_str());
        } else {
            mapGeneMaxWeight[mapIdName[vec_string[1]]] = atof(vec_string[2].c_str());
            netGenes1.push_back(mapIdName[vec_string[1]]);
        }
        // Weight;
        if (mapGeneWeight.find(mapIdName[vec_string[0]]) != mapGeneWeight.end()) {
            mapGeneWeight[mapIdName[vec_string[0]]] += atof(vec_string[2].c_str());
        } else {
            mapGeneWeight[mapIdName[vec_string[0]]] = atof(vec_string[2].c_str());
        }
        if (mapGeneWeight.find(mapIdName[vec_string[1]]) != mapGeneWeight.end()) {
            mapGeneWeight[mapIdName[vec_string[1]]] += atof(vec_string[2].c_str());
        } else {
            mapGeneWeight[mapIdName[vec_string[1]]] = atof(vec_string[2].c_str());
        }
    }
    ifs.close();
    net.geneMaxWeight = mapGeneMaxWeight; net.genesInNet = netGenes1; net.geneWeight = mapGeneWeight;
    return net;
}

vector<string> vectors_set_union(vector<string> v1, vector<string> v2) {
    vector<string> v;
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));//求交集
    return v;
}

vector<double> CalScaledVector(vector<double> vec) {
    double mean = 0, sd = 0;
    for (auto & i : vec) mean += i;
    mean = mean / vec.size();
    for (auto & i : vec) sd += pow((i - mean), 2);
    sd = sqrt(sd / (vec.size() - 1));
    for (auto & i : vec) i = (i - mean) / sd;
    return vec;
}

vector<double> absVector(vector<double> vec) {
    for (auto & i : vec) i = abs(i);
    return vec;
}

vector<double> log2Vector(vector<double> vec) {
    for (auto & i : vec) i = log2(i);
    return vec;
}

bool cmp(pair<string, double> a, pair<string, double> b) {
    return a.second > b.second;
}

vector<string> GetMajorGenes(map<string, double> mapGeneMaxWeight, map<string, double> mapGeneWeight,map<string, double> geneDiffScore, vector<string> majorGeneCandidates, float t1 = 0.3, float t2 = 0, float t3 = 0.04) {
    vector<pair<string, double>> vecNetMaxWeight, vecNetWeight, vecDiff;
    int num = majorGeneCandidates.size();
    vector<string> vec1, vec2, vec3, vecUnion;
    // Max Weight;
    for (auto & i : mapGeneMaxWeight) {
        if (std::find(majorGeneCandidates.begin(), majorGeneCandidates.end(), i.first) != majorGeneCandidates.end()) {
            vecNetMaxWeight.push_back(pair<string, double> (i.first, i.second));
        }
    }
    // Sum of weight;
    for (auto & i : mapGeneWeight) {
        if (std::find(majorGeneCandidates.begin(), majorGeneCandidates.end(), i.first) != majorGeneCandidates.end()) {
            vecNetWeight.push_back(pair<string, double> (i.first, i.second));
        }
    }
    // Diff score;
    for (auto & i : geneDiffScore) {
        if (std::find(majorGeneCandidates.begin(), majorGeneCandidates.end(), i.first) != majorGeneCandidates.end()) {
            vecDiff.push_back(pair<string, double> (i.first, i.second));
        }
    }
    sort(vecNetMaxWeight.begin(), vecNetMaxWeight.end(), cmp);
    sort(vecNetWeight.begin(), vecNetWeight.end(), cmp);
    sort(vecDiff.begin(), vecDiff.end(), cmp);
    vector<pair<string, double>>::iterator iter;
    for (int i = 0; i < ((int)(num * t1) - 1); i++) vec1.push_back(vecNetMaxWeight[i].first);
    for (int i = 0; i < ((int)(num * t2) - 1); i++) vec2.push_back(vecNetWeight[i].first);
    for (int i = 0; i < ((int)(num * t3) - 1); i++) vec3.push_back(vecDiff[i].first);
    vecUnion = vectors_set_union(vec1, vec2);
    return vectors_set_union(vecUnion, vec3);
}

double Pearson(vector<double> &inst1, vector<double> &inst2) {
    if(inst1.size() != inst2.size()) {
        std::cout<<"the size of the vectors is not the same\n";
        return 0;
    }
    size_t n = inst1.size();
    double pearson = n * inner_product(inst1.begin(), inst1.end(), inst2.begin(), 0.0) -
                     accumulate(inst1.begin(), inst1.end(), 0.0)*accumulate(inst2.begin(), inst2.end(), 0.0);
    double temp1= n * inner_product(inst1.begin(), inst1.end(), inst1.begin(), 0.0) -
                  pow(accumulate(inst1.begin(), inst1.end(), 0.0), 2.0);
    double temp2 = n * inner_product(inst2.begin(), inst2.end(), inst2.begin(), 0.0) -
                   pow(accumulate(inst2.begin(), inst2.end(), 0.0), 2.0);
    temp1=sqrt(temp1);
    temp2=sqrt(temp2);
    pearson=pearson/(temp1*temp2);
    return pearson;
}

vector<string> GetMajorGeneCandidates(map<string, map<string, double>> normalMatrix, map<string, map<string, double>> tumorMatrix, vector<string> mutationGenes, vector<string> netGenes) {
    vector<string> majorGenesCandidates;
    for (auto & i : normalMatrix) {
        string geneTemp = i.first;
        // normal vector;
        int flag = 0;
        for (auto & j : i.second) {
            if (j.second == 0) {
                flag = 1;
                break;
            }
        }
        for (auto & j : tumorMatrix[geneTemp]) {
            if (j.second == 0) {
                flag = 1;
                break;
            }
        }
        if (flag == 1) {
            continue;
        } else {
            if (((count(mutationGenes.begin(), mutationGenes.end(), geneTemp)) != 0) && ((count(netGenes.begin(), netGenes.end(), geneTemp)) != 0)) {
                majorGenesCandidates.push_back(geneTemp);
            }
        }
    }
//    cout << "Candidates: " << majorGenesCandidates.size() << endl;
    return majorGenesCandidates;
}



Driver::Diff Driver::GetDiffMatrix(string mutationPath, string netPath1, string geneIDPath, string expTumorPath, string expNormalPath, unordered_set<string> preGenes) {
    cout << "====================================== STEP 1: Filter Major Genes =====================================" << endl;
    cout << "====                                                                                               ====" << endl;
    Diff diff;
    vector<string> genes;
    Net net = GetGeneMaxWeight(netPath1, geneIDPath);
    vector<string> mutationGenes = GetMutationGenes(mutationPath);
    // 读取Tumor/Normal表达矩阵；
    Expression expNormal = ReadExpressionMatrixWithPreGenes(expNormalPath, preGenes);
    Expression expTumor = ReadExpressionMatrixWithPreGenes(expTumorPath, preGenes);
//    Expression expNormal = ReadExpressionMatrixForRandom(expNormalPath, 0.8);
    cout << "====                                  Expression Tumor Size: " << expTumor.p.size() << "                                  ====" << endl;
    cout << "====                                  Expression Normal Size: " << expNormal.p.size() << "                                  ====" << endl;
    vector<string> majorGenesCandidates = GetMajorGeneCandidates(expNormal.m, expTumor.m, mutationGenes, net.genesInNet);

    // Get tumor/normal expression matrix of major gene candidates;
    map<string, map<string, double>> expNormalC, expTumorC;
    for (auto & pair : expNormal.m) {
        auto iter = std::find(majorGenesCandidates.begin(), majorGenesCandidates.end(), pair.first);
        if (iter != majorGenesCandidates.end()) {
            expNormalC[pair.first] = pair.second;
            expTumorC[pair.first] = expTumor.m[pair.first];
        }
    }
    // Filter major genes and get the differential expression matrix;
    vector<string> bothPatients = vectors_intersection(expTumor.p, expNormal.p);
    map<string, vector<double>> diffExpMatrix;
    map<string, double> geneDiffScoreEU;
//    map<string, double> geneDiffScorePearson;
    for (auto & pair : expTumorC) {
        vector<double> tempVec, tempVecT, tempVecN;
        for (auto & i : bothPatients) {
            auto iterString = pair.second.find(i);
            if (iterString != pair.second.end()) {
                tempVec.push_back(pair.second[i] - expNormalC[pair.first][i]);
            }
        }
        double diffScore = 0;
        for (auto & i : tempVec) diffScore += pow(i, 2);
        geneDiffScoreEU[pair.first] = sqrt(diffScore);

//        for (auto & i : tempVec) diffScore += abs(i);
//        vector<double> absVecTemp = absVector(tempVec);
//        diffScore = *(max_element(absVecTemp.begin(), absVecTemp.end()));
//        geneDiffScoreEU[pair.first] = diffScore;
        diffExpMatrix[pair.first] = tempVec;
    }
    vector<string> majorGenes = GetMajorGenes(net.geneMaxWeight, net.geneWeight, geneDiffScoreEU, majorGenesCandidates);
    map<int, vector<double>> finalDiffExpMatrix;
    int row = 0;
    for (auto & i : diffExpMatrix) {
        vector<double> vecNormed = absVector(CalScaledVector(i.second)); // ****
        auto iter = std::find(majorGenes.begin(), majorGenes.end(), i.first);
        if (iter != majorGenes.end()) {
            genes.push_back(i.first);
            finalDiffExpMatrix[row] = vecNormed;
            double expScore = 0;
            for (auto & value:vecNormed) expScore += value;
            expScore = expScore / vecNormed.size();
            mapIdExp[mapOfficialToNCBI[i.first]] = expScore;
            row++;
        }
    }
    diff.m = finalDiffExpMatrix; diff.g = genes; diff.bp = bothPatients;
    cout << "====                                   Major Genes' Number: " << majorGenes.size() << "                                   ====" << endl;
    cout << "====                                                                                               ====" << endl;
    cout << "=======================================================================================================" << endl;
    cout << endl;
//    cout << "Major genes is filtered, number is " << majorGenes.size() << "." << endl;
    return diff;
}

int findPosVector(vector <double> input , double number) {
    vector<double>::iterator iter=std::find(input.begin(),input.end(),number);
    if(iter == input.end())
    {
        return -1;
    } else{
        return distance(input.begin(),iter);
    }
}

double Spearman(vector<double>& vec1, vector<double>& vec2) {
    int n = vec1.size();
    vector<double> rankVec1, rankVec2, copyVec1, copyVec2;
    copyVec1.assign(vec1.begin(), vec1.end());
    copyVec2.assign(vec2.begin(), vec2.end());
    sort(copyVec1.rbegin(), copyVec1.rend());
    sort(copyVec2.rbegin(), copyVec2.rend());
    for (auto & i : vec1) rankVec1.push_back(findPosVector(copyVec1, i) + 1);
    for (auto & i : vec2) rankVec2.push_back(findPosVector(copyVec2, i) + 1);
    double value = 0;
    for (int i = 0; i < n; i++) {
        value += pow((rankVec1[i] - rankVec2[i]), 2);
    }
    value = 1 - (6 * value) / (n * (pow(n, 2) - 1));
    return value;
}

//map<int, vector<double>> GetDoubleDiffMatrix(map<int ,vector<string>> diffMatrixRaw) {
//    map<int, vector<double>> diffMatrix;
//    int row = 0;
//    for (auto & i : diffMatrixRaw) {
//        if (i.second[0].substr(0,4) == "TCGA") {
//            continue;
//        }
//        vector<double> tempVecDouble;
//        for (int j = 1; j < i.second.size(); j++) {
//            tempVecDouble.push_back(atof(i.second[j].c_str()));
//        }
//        diffMatrix[row] = tempVecDouble;
//        row++;
//    }
//    return diffMatrix;
//}
//
//vector<string> GetMajorGenes(map<int ,vector<string>> diffMatrixRaw) {
//    vector<string> majorGenes;
//    for (auto & i : diffMatrixRaw) {
//        majorGenes.push_back(i.second[0]);
//    }
//    majorGenes.erase(majorGenes.begin());
//    return majorGenes;
//}

//map<int, vector<double>> GetCorrMatrix (map<int, vector<double>> diffMatrix, int k) {
//    map<int, vector<double>> corrMatrix;
//    for (int i = 0; i < diffMatrix.size(); i++) {
//        vector<double> tempVecDouble;
//        for (int j = 0; j < diffMatrix.size(); j++) {
//            double correlation = (Spearman(diffMatrix[i], diffMatrix[j]) + Pearson(diffMatrix[i], diffMatrix[j])) / 2;
//            if (correlation > 0) {
//                tempVecDouble.push_back(pow(correlation, k));
//            } else {
//                tempVecDouble.push_back(0);
//            }
//        }
//        corrMatrix[i] = tempVecDouble;
//    }
//    return corrMatrix;
//}

map<int, vector<double>> GetCorrMatrix (map<int, vector<double>> diffMatrix) {
    map<int, vector<double>> corrMatrix;
    for (int i = 0; i < diffMatrix.size(); i++) {
        vector<double> tempVecDouble;
        for (int j = 0; j < diffMatrix.size(); j++) {
            double correlation = Pearson(diffMatrix[i], diffMatrix[j]);
            if (correlation > 0) {
                tempVecDouble.push_back(correlation);
            } else {
                tempVecDouble.push_back(0);
            }
        }
        corrMatrix[i] = tempVecDouble;
    }
    return corrMatrix;
}

vector<vector<double>> GetCorrMatrix1 (map<int, vector<double>> diffMatrix) {
    vector<vector<double>> corrMatrix;
    for (int i = 0; i < diffMatrix.size(); i++) {
        vector<double> tempVecDouble;
        for (int j = 0; j < diffMatrix.size(); j++) {
            double correlation = Pearson(diffMatrix[i], diffMatrix[j]);
            if (correlation > 0) {
                tempVecDouble.push_back(correlation);
            } else {
                tempVecDouble.push_back(0);
            }
        }
        corrMatrix.push_back(tempVecDouble);
    }
    return corrMatrix;
}

//map<int, map<int, double>> GetCorrMatrix (map<int, vector<double>> diffMatrix) {
//    map<int, map<int, double>> corrMatrix;
//    for (int i = 0; i < diffMatrix.size(); i++) {
//        for (int j = i; j < diffMatrix.size(); j++) {
//            double correlation = (Spearman(diffMatrix[i], diffMatrix[j]) + Pearson(diffMatrix[i], diffMatrix[j])) / 2;
//            if (i != j) {
//                corrMatrix[i][j] = correlation;
//                corrMatrix[j][i] = correlation;
//            } else corrMatrix[i][j] = correlation;
//        }
//    }
//    return corrMatrix;
//}

map<int, vector<double>> GetPowerCorrMatrix (map<int, vector<double>> corrMatrix, int power) {
    for (auto & i : corrMatrix) {
        for (auto & j : i.second) {
            j = pow(j, power);
        }
    }
    return corrMatrix;
}

vector<vector<double>> GetPowerCorrMatrix1 (vector<vector<double>> corrMatrix, int power) {
    for (auto & i : corrMatrix) {
        for (auto & j : i) {
            j = pow(j, power);
        }
    }
    return corrMatrix;
}

vector<double> GetConnectivity (map<int, vector<double>> corrMatrix) {
    vector<double> connectivity;
    for (auto & i : corrMatrix) {
        double sum=0;
        for (auto & j : i.second) {
            sum += j;
        }
        connectivity.push_back(sum - 1);
    }
    return connectivity;
}

vector<double> GetConnectivity1 (vector<vector<double>> corrMatrix) {
    vector<double> connectivity;
    for (auto & i : corrMatrix) {
        double sum=0;
        for (auto & j : i) {
            sum += j;
        }
        connectivity.push_back(sum - 1);
    }
    return connectivity;
}

double InnerConnectValue(vector<double> vec1, int pos1, vector<double> vec2, int pos2) {
    vec1[pos1] = 0; vec2[pos2] = 0;
    double value = 0;
    for (int i = 0; i < vec1.size(); i++) {
        value += vec1[i] * vec2[i];
    }
    return value;
}

double OrdinaryLeastSquare(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    double value1 = 0 ,sumX = 0, sumY = 0, value2 = 0;
    for (int i = 0; i < n; i++) {
        value1 += x[i] * y[i];
    }
    for (auto & i : x) sumX += i;
    for (auto & i : y) sumY += i;
    for (auto & i : x) value2 += pow(i, 2);
    double slope = (n * value1 - sumX * sumY) / (n * value2 - pow(sumX, 2));
    double intercept = (sumY / n) - slope * (sumX / n);
//    cout << slope << " " << intercept << endl;
    double SST = 0, SSE = 0;
    for (int i = 0; i < n; i++) {
        SST += pow((y[i] - (sumY / n)), 2);
    }
    for (int i = 0; i < n; i++) {
        SSE += pow((y[i] - (slope * x[i] + intercept)), 2);
    }
    double r2 = 1 - (SSE / SST);
    return r2;
}

void GetCorrMatrixForEachPower (map<int, vector<double>> initCorrMatrix, vector<double>& r2List, int power) {
    mtx.lock();  // 获取互斥锁
    map<int, vector<double>> corrMatrix = GetPowerCorrMatrix(initCorrMatrix, power);
    vector<double> connectivity = GetConnectivity(corrMatrix);
    double min = *min_element(connectivity.begin(), connectivity.end());
    double max = *max_element(connectivity.begin(), connectivity.end());
    int intervalNum = 10;
    double step = (max - min) / intervalNum;
    vector<double> countVec(intervalNum, 0);
    vector<double> x;
    for (auto & k : connectivity) {
        for (int i = 0; i < intervalNum; i++) {
            if (k < min + (i+1) * step & k >= min + i * step) countVec[i]++;
        }
    }
    countVec[intervalNum-1]++;
    for (int i = 0; i < intervalNum; i++) {
        x.push_back(min + i * step + step / 2);
    }
    for (auto & i : countVec) {
        i = i / connectivity.size();
    }
    for (auto & i : countVec) i = log10(i);
    for (auto & i : x) i = log10(i);
    double r2 = OrdinaryLeastSquare(x, countVec);
//    if (std::isnan(static_cast<double>(r2))) r2List[power - 1] = 0;
//    else r2List[power - 1] = r2;

    mtx.unlock();  // 释放互斥锁
}

void GetCorrMatrixForEachPower1 (vector<vector<double>> initCorrMatrix, vector<double>& r2List, int power) {
    mtx.lock();  // 获取互斥锁
    vector<vector<double>> corrMatrix = GetPowerCorrMatrix1(initCorrMatrix, power);
    vector<double> connectivity = GetConnectivity1(corrMatrix);
    double min = *min_element(connectivity.begin(), connectivity.end());
    double max = *max_element(connectivity.begin(), connectivity.end());
    int intervalNum = 10;
    double step = (max - min) / intervalNum;
    vector<double> countVec(intervalNum, 0);
    vector<double> x;
    for (auto & k : connectivity) {
        for (int i = 0; i < intervalNum; i++) {
            if (k < min + (i+1) * step & k >= min + i * step) countVec[i]++;
        }
    }
    countVec[intervalNum-1]++;
    for (int i = 0; i < intervalNum; i++) {
        x.push_back(min + i * step + step / 2);
    }
    for (auto & i : countVec) {
        i = i / connectivity.size();
    }
    for (auto & i : countVec) i = log10(i);
    for (auto & i : x) i = log10(i);
    double r2 = OrdinaryLeastSquare(x, countVec);
//    if (isnan(r2)) r2List[power - 1] = 0;
//    else r2List[power - 1] = r2;

    mtx.unlock();  // 释放互斥锁
}

SolfThreshold GetOptimalCorrMatrix (map<int, vector<double>> diffMatrix, float threshold=0.90, int maxPower=6) {
    cout << "================================ STEP 2: Generate Differential Network ================================" << endl;
    cout << "====                                                                                               ====" << endl;
//    cout << "Starting get the optimal correlation matrix." << endl;
    map<int, vector<double>> initCorrMatrix = GetCorrMatrix(diffMatrix);
//    cout << "Init correlation matrix is generated." << endl;
    SolfThreshold finalOutput;
    vector<double> r2List = {0, 0, 0, 0, 0, 0};
    std::vector<std::thread> threads;
    for (int power = 1; power < maxPower + 1; power++) threads.push_back(std::thread(GetCorrMatrixForEachPower, initCorrMatrix, std::ref(r2List), power));
    for (auto& t : threads) t.join();

    int flag = 0;
    for (int i = 0; i < r2List.size(); i++) {
        if (r2List[i] > threshold) {
            finalOutput.corrMatrix = GetPowerCorrMatrix(GetCorrMatrix(diffMatrix), i + 1);
            finalOutput.connectivity = GetConnectivity(finalOutput.corrMatrix);
            flag = 1;
        }
    }
    if (flag == 0) {
        vector<double>::iterator maxValue = max_element(r2List.begin(), r2List.end());
        int maxPos = distance(r2List.begin(), maxValue) + 1;
        finalOutput.corrMatrix = GetPowerCorrMatrix(GetCorrMatrix(diffMatrix), maxPos);
        finalOutput.connectivity = GetConnectivity(finalOutput.corrMatrix);
    }
    cout << "====                                              Done                                             ====" << endl;
    cout << "====                                                                                               ====" << endl;
    cout << "=======================================================================================================" << endl;
    cout << endl;
//    cout << "The optimal power does not exist across the threshold, we use the power " << maxPos << " that has max r2 " << *maxValue << "." << endl;
    return finalOutput;
}


SolfThreshold1 GetOptimalCorrMatrix1 (map<int, vector<double>> diffMatrix, float threshold=0.90, int maxPower=6) {
    cout << "===================== STEP 2: Compute Correlation Matrix | Generate Diff Network ======================" << endl;
//    cout << "Starting get the optimal correlation matrix." << endl;
    vector<vector<double>> initCorrMatrix1 = GetCorrMatrix1(diffMatrix);
//    cout << "Init correlation matrix is generated." << endl;
    SolfThreshold1 finalOutput1;
    vector<double> r2List = {0, 0, 0, 0, 0, 0};
    std::vector<std::thread> threads;
    for (int power = 1; power < maxPower + 1; power++) threads.push_back(std::thread(GetCorrMatrixForEachPower1, initCorrMatrix1, std::ref(r2List), power));
    for (auto& t : threads) t.join();

    int flag = 0;
    for (int i = 0; i < r2List.size(); i++) {
        if (r2List[i] > threshold) {
            finalOutput1.corrMatrix = GetPowerCorrMatrix1(GetCorrMatrix1(diffMatrix), i + 1);
            finalOutput1.connectivity = GetConnectivity1(finalOutput1.corrMatrix);
            flag = 1;
        }
    }
    if (flag == 0) {
        vector<double>::iterator maxValue = max_element(r2List.begin(), r2List.end());
        int maxPos = distance(r2List.begin(), maxValue) + 1;
        finalOutput1.corrMatrix = GetPowerCorrMatrix1(GetCorrMatrix1(diffMatrix), maxPos);
        finalOutput1.connectivity = GetConnectivity1(finalOutput1.corrMatrix);
    }
//    cout << "The optimal power does not exist across the threshold, we use the power " << maxPos << " that has max r2 " << *maxValue << "." << endl;
    return finalOutput1;
}



void Driver::GetFinalMatrix (map<int, vector<double>> corrMatrix, vector<double> connectivity, vector<string> majorGenes) {
    double maxValue = 0;
    for (int i = 0; i < corrMatrix.size(); i++) {
        auto iter_i = mapOfficialToNCBI.find(majorGenes[i]);
        if (iter_i == mapOfficialToNCBI.end()) continue;
        map<string, double> tempMap;
        list<PPI> pid_i;
        for (int j = 0; j < corrMatrix.size(); j++) {
            auto iter_j = mapOfficialToNCBI.find(majorGenes[j]);
            if (iter_j == mapOfficialToNCBI.end()) continue;
            double value;
            value = (InnerConnectValue(corrMatrix[i], i, corrMatrix[j], j) + corrMatrix[i][j]) /
                    (min(connectivity[i], connectivity[j]) + 1 - corrMatrix[i][j]);
            if (value > maxValue) maxValue = value;

            if (j >= i+1) {
                tempMap[majorGenes[j]] = value;
                pid_i.push_back(PPI(iter_j->second, value));
            }
            if (j < i){
                pid_i.push_back(PPI(iter_j->second, value));
            }
        }
//        cout << pid_i.size() << endl;
        ExpNetWork[iter_i->second] = pid_i;
    }
//    cout << "Optimal correlation matrix is generated." << endl;
}

void Driver::GetFinalMatrix1 (vector<vector<double>> corrMatrix, vector<double> connectivity, vector<string> majorGenes) {
    double maxValue = 0;
    for (int i = 0; i < corrMatrix.size(); i++) {
        auto iter_i = mapOfficialToNCBI.find(majorGenes[i]);
        if (iter_i == mapOfficialToNCBI.end()) continue;
        map<string, double> tempMap;
        list<PPI> pid_i;
        for (int j = 0; j < corrMatrix.size(); j++) {
            auto iter_j = mapOfficialToNCBI.find(majorGenes[j]);
            if (iter_j == mapOfficialToNCBI.end()) continue;
            double value;
            value = (InnerConnectValue(corrMatrix[i], i, corrMatrix[j], j) + corrMatrix[i][j]) /
                    (min(connectivity[i], connectivity[j]) + 1 - corrMatrix[i][j]);
            if (value > maxValue) maxValue = value;

            if (j >= i+1) {
                tempMap[majorGenes[j]] = value;
                pid_i.push_back(PPI(iter_j->second, value));
            }
            if (j < i){
                pid_i.push_back(PPI(iter_j->second, value));
            }
        }
//        cout << pid_i.size() << endl;
        ExpNetWork[iter_i->second] = pid_i;
    }
    cout << "Optimal correlation matrix is generated." << endl;
}

void WriteDiffMatrix(string diffMatrixPath, map<int, vector<double>> diffMatrix, vector<string> patients, vector<string> genes) {
    ifstream m_ifs;
    ofstream m_ofs;
    m_ofs.open(diffMatrixPath);
    for (auto & i : patients) m_ofs << i << "\t";
    m_ofs << "\n";
    for (int i = 0; i < diffMatrix.size(); i++) {
        m_ofs << genes[i] << "\t";
        for (auto & j : diffMatrix[i]) m_ofs << setprecision(16) << j << "\t";
        m_ofs << "\n";
    }
}

void WriteDiffNetwork(string writeDNPath, map<string, map<string, double>> finalMatrix, vector<string> majorGenes) {
    ifstream m_ifs;
    ofstream m_ofs;
    m_ofs.open(writeDNPath);
    for (auto & i : finalMatrix) {
        for (auto & j : i.second) {
            m_ofs << i.first << "\t" << j.first << '\t' << j.second << "\n";
        }
    }
    m_ofs.close();
}


