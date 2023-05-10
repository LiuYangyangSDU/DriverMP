#pragma once
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<list>
#include<algorithm>
#include<cmath>
#include<math.h>
#include<fstream>
#include<ostream>
#include<string.h>
#include <unordered_set>
using namespace std;
#define STR_SEPTABLE "\t"
#ifndef _DRIVER_
#define _DRIVER_

class Driver {
    typedef  struct PPI {
        int Id;
        double Id_weight;
        PPI(int _Id, double _Id_weight) {
            this->Id = _Id;
            this->Id_weight = _Id_weight;
        }
    }PPI;

public:
    struct Diff {
        map<int, vector<double>> m;
        vector<string> g;
        vector<string> bp;
    };

    //MutationScore.cpp
    void ReadGeneID(string referencePath);//读取基因名称对应ID
    //void ReadGeneExpMIF();//读取由MaxMIF算出的差异基因的表达MIF
    bool ReadGenomicMatrix(string mutationPath);//读取突变矩阵
    void MutScore();//计算突变分数

    //PairMIF.cpp
    void ReadCancerGene(string referencePath);//读取参考数据集
//    bool ReadNetWork(string netPath);//读取PPI网络
//    void ReadGeneExp(string expPath);//读取基因表达差异的平均值
//    bool ReadGeneExpNetWork(string netPath);//读取由WCGNA方法算出来的，基因差异表达网络
    void PairGeneMutation();//计算成对基因的突变分数乘积
    map<vector<int>, double> ComputeDCIS(vector<string> majorGenes);//计算成对基因的pairMIF值
    void ComputePairMIF1();
    //SingleMIF.cpp
    void InitialMaxMIF();//计算用MaxMIF算法计算每个基因的MaxMIF值
    void OutDegree();//计算PPI网络中每个基因的度
    void EdgeWeight();//计算PPI网络中每个基因连边的权重
    void ComputeSingleMIF_degree(map<vector<int>, double> pairIdMIF);//用度对每对基因进行拆分
    void OutputRank(string savePath);//输出最终排名
    void ComputeDCISCircle(int gene1, int gene2, int index);
    vector<vector<int>> GetRequiredPairGenes(vector<string> majorGenes);
    void GetFinalMatrix (map<int, vector<double>> corrMatrix, vector<double> connectivity, vector<string> majorGenes);
    void GetFinalMatrix1 (vector<vector<double>> corrMatrix, vector<double> connectivity, vector<string> majorGenes);
    void mapGeneID(string geneIDPath);
    unordered_set<string> GetPreGenes(string mutationPath, string expNormalPath, string expTumorPath, string networkPath);
    unordered_set<string> ReadNetGenes(string networkPath);
    Diff GetDiffMatrix(string mutationPath, string netPath1, string geneIDPath, string expTumorPath, string expNormalPath, unordered_set<string> preGenes);


private:
    void StringSplit(string longstr, string sep, vector<string>& vectstr)//自定义分割函数
    {//分割函数
        vectstr.clear();
        int nend = 0;
        int nbegin = 0;
        while (nend != string::npos) {
            nend = longstr.find(sep, nbegin);
            if (nend == string::npos)
                vectstr.push_back(longstr.substr(nbegin, longstr.length() - nbegin));
            else vectstr.push_back(longstr.substr(nbegin, nend - nbegin));
            nbegin = nend + sep.length();
        }
    }

    static bool cmp(const pair<int, double>& a, const pair<int, double>& b) {
        return a.second > b.second;
    }

private:
    ifstream m_ifs;
    ofstream m_ofs;
    map<string, double> mapGeneMut;
    vector<int> m_vecTempInt;
    int intSampleNumber;
    string strGeneName;
    vector<string> vecGeneName;
    vector<string>::iterator iterVecString;
    vector<int>::iterator iterVecInt;
    map<int, vector<string> > mapSampleMutationGene;
    map<int, vector<string> >::iterator iterIntVecString;
    map<string, vector<int> > mapGeneMutationSample;
    map<string, vector<int> >::iterator iterStringVecInt;
    map<string, int>::iterator iterStringInt;
    vector<string> vecUnMutationGene;
    double minScore;
    map<int, int> mapSampleTotal;
    map<string, double> mapGeneMutation;
    map<int, double> mapIdMutation;
    vector<pair<string, double> >  vecMutationScore;
    map<string, double>::iterator iterStringDouble;
    vector<pair<string, double> >::iterator  vecPairMutationScore;
    set<int> CancerGene[1];
    int intGeneId, intGeneId1;
    double MinWeight;
    vector<vector<int> > FunctionalNet;
    map<int, list<PPI> > NetWork;
    map<int, list<PPI> >::iterator iterPPI;
    map<int, double>  mapIdExp;
    map<int, list<PPI> > ExpNetWork;
    int nodeInt, nodeInt1;
    map<int, double>::iterator iterIntDouble, iterIntDouble1;
    map<int, double> mapGeneDegree, mapGeneMaxWeight ,mapGeneSumWeight;
    list<PPI>::iterator iterListPPI, iterListPPI1;
    map<vector<int>, double> PairGeneMut;
    map<vector<int>, double>::iterator iterMapVecIntDouble, iterMapVecIntDouble1;
    map<double, int> mapSingleMIFID;
    map<double, string> mapSingleMIFGene;
    vector<pair<int, double> > vecIdSingleMIF;
    vector<double> vecMaxExpPairMIF;
    vector<double> pairMIF;
    map<int, string> mapNCBIToOfficial;
    map<string, int> mapOfficialToNCBI;
};



#endif