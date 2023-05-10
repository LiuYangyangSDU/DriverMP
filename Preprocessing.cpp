#include<iostream>
#include"Driver.h"
#include <vector>
#include <unordered_set>
//#include <cstdlib>

void StringSplit(string long_str, string sep, vector<string>& vec_str) {
    vec_str.clear();
    int n_end = 0;
    int n_start = 0;
    while (n_end != string::npos) {
        n_end = long_str.find(sep, n_start);
        if(n_end == string::npos) {
            vec_str.push_back(long_str.substr(n_start, long_str.length()-n_start));
        } else {
            vec_str.push_back(long_str.substr(n_start, n_end-n_start));
        }
        n_start = n_end + sep.length();
    }
}

void Driver::mapGeneID(string geneIDPath) {
    ifstream ifs;
    ifs.open(geneIDPath.c_str());
    string str_temp;
    vector<string> vec_string;
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        mapNCBIToOfficial[atoi(vec_string[1].c_str())] = vec_string[0];
        mapOfficialToNCBI[vec_string[0]] = atoi(vec_string[1].c_str());
    }
}

unordered_set<string> ReadMutGenes(string mutationPath) {
    ifstream ifs;
    ifs.open(mutationPath.c_str());
    ifs.imbue(std::locale("en_US.UTF-8"));
    try {
        if (!ifs.is_open()) throw "Failed to read mutation data, please check the input path.";
    } catch (const char* message) {
        std::cerr << "Error: " << message << std::endl;
        exit(1);
    }
    string str_temp;
    vector<string> vec_string;
    unordered_set<string> mutationGenes;
    getline(ifs, str_temp); // 跳过第一行；
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        mutationGenes.insert(vec_string[0]);
    }
    return mutationGenes;
}

unordered_set<string> ReadExpGenes(string expNormalPath, string expTumorPath) {
    ifstream ifs, ifs1;
    ifs.open(expNormalPath.c_str());
    ifs1.open(expTumorPath.c_str());
    try {
        if (!ifs.is_open() or !ifs1.is_open()) throw "Failed to read gene expression data, please check the input path.";
    } catch (const char* message) {
        std::cerr << "Error: " << message << std::endl;
        exit(1);
    }
    string str_temp;
    vector<string> vec_string;
    unordered_set<string> expressionGenes;
    getline(ifs, str_temp); // 跳过第一行；
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        expressionGenes.insert(vec_string[0]);
    }
    return expressionGenes;
}

unordered_set<string> Driver::ReadNetGenes(string networkPath) {
    MinWeight = 1.0;
    ifstream ifs;
    ifs.open(networkPath.c_str());
    try {
        if (!ifs.is_open()) throw "Failed to read PPI network data, please check the input path.";
    } catch (const char* message) {
        std::cerr << "Error: " << message << std::endl;
        exit(1);
    }
    string str_temp;
    vector<string> vec_string;
    unordered_set<string> networkGeneSet;
    while (getline(ifs, str_temp)) {
        StringSplit(str_temp, STR_SEPTABLE, vec_string);
        int gene1 = atoi(vec_string[0].c_str());
        int gene2 = atoi(vec_string[1].c_str());
        double weight = atof(vec_string[2].c_str());
        if (mapNCBIToOfficial.find(gene1) != mapNCBIToOfficial.end() and mapNCBIToOfficial.find(gene2) != mapNCBIToOfficial.end()) {
            networkGeneSet.insert(mapNCBIToOfficial[gene1]);
            networkGeneSet.insert(mapNCBIToOfficial[gene2]);
        }

        auto iter1 = mapNCBIToOfficial.find(gene1);
        auto iter2 = mapNCBIToOfficial.find(gene2);
        if (iter1 != mapNCBIToOfficial.end() and iter2 != mapNCBIToOfficial.end()) {
            vector<int> pairGenes = {gene1, gene2};
            FunctionalNet.emplace_back(pairGenes);
            if (weight < MinWeight) MinWeight = weight;
            // 对第一个基因；
            auto iter_1 = NetWork.find(gene1);
            if (iter_1 != NetWork.end()) iter_1->second.emplace_back(PPI(gene2, weight));
            else {
                list<PPI> listTemp;
                listTemp.emplace_back(PPI(gene2, weight));
                NetWork[gene1] = listTemp;
            }
            // 对第二个基因；
            auto iter_2 = NetWork.find(gene2);
            if (iter_2 != NetWork.end()) iter_2->second.emplace_back(PPI(gene1, weight));
            else {
                list<PPI> listTemp;
                listTemp.emplace_back(PPI(gene1, weight));
                NetWork[gene2] = listTemp;
            }

        }
    }
//    cout << "Edges number in PPI networks: " << FunctionalNet.size() << ", Nodes number: " << NetWork.size() << endl;
    return networkGeneSet;
}

unordered_set<string> Driver::GetPreGenes(string mutationPath, string expNormalPath, string expTumorPath, string networkPath) {
    unordered_set<string> mutGenes, expGenes, netGenes;
    mutGenes = ReadMutGenes(mutationPath);
    expGenes = ReadExpGenes(expNormalPath, expTumorPath);
    netGenes = ReadNetGenes(networkPath);

    std::unordered_set<string> result;
    for (const auto& elem : netGenes) {
        if (mutGenes.find(elem) != mutGenes.end() && expGenes.find(elem) != expGenes.end()) {
            result.insert(elem);
        }
    }
    try {
        if (result.size() == 0) throw "Please input formatted data.";
    } catch (const char* message) {
        std::cerr << "Error: " << message << std::endl;
        exit(1);
    }
    return result;
}

