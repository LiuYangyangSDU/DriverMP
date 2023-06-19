#include<iostream>
#include"Driver.h"
#include"DiffNetwork.h"
#include <vector>
#include "Preprocess.h"
#include <unordered_set>
#include <chrono>
#include <cstring>
#define STR_SEPTABLE "\t"

using namespace std;

int threadNum = 1;

int main(int argc, char* argv[]) {
    string geneIDPath = "DriverMP_example//geneID.txt";
    string referencePath = "DriverMP_example//CGC.txt";
    string mutationPath;
    string expTumorPath;
    string expNormalPath;
    string netPath;
    string savePath = "Output//output.txt";

    if (argc == 5) {
        mutationPath = argv[1];
        expTumorPath = argv[2];
        expNormalPath = argv[3];
        netPath = argv[4];
    }
    else if (argc == 6){
        mutationPath = argv[1];
        expTumorPath = argv[2];
        expNormalPath = argv[3];
        netPath = argv[4];
        threadNum = stoi(argv[5]);
    }
    else if (argc == 2) {
        string arg = argv[1];
        if (arg == "-help" || arg == "-h") {
            std::cout << "Usage：" << std::endl;
            std::cout << "DriverMP [Mutation data] [Tumor expression data] [Normal expression data] [PPI network] [Thread Number]" << std::endl;
            std::cout << "[Mutation data]: Path to format-compliant non-silent somatic mutation data;" << std::endl;
            std::cout << "[Tumor expression data]: Path to RNA-Seq data (FPKM normalized) for specific cancer samples;" << std::endl;
            std::cout << "[Normal expression data]: Path to RNA-Seq data (FPKM normalized) of normal samples corresponding to the specific cancer;" << std::endl;
            std::cout << "[PPI network data]: Path to Protein-Protein Interaction Networks;" << std::endl;
            std::cout << "[Thread Number]: Number of cores for multithreading, the default is 1;" << std::endl;
            std::cout << "-v, -version Show program version number" << std::endl;
            // 添加更多的帮助信息
            return 0;
        } else if (arg == "-v" || arg == "-version") {
            std::cout << "Version of DriverMP：1.0" << std::endl;
            return 0;
        } else {
            cout << "Wrong number of parameters, please try again." << endl;
            return 1;
        }
    }
    else {
        cout << "Wrong number of parameters, please try again." << endl;
        return 1;
    }

    // Preprocessing the original genes;
    Driver D;
    D.mapGeneID(geneIDPath);
    unordered_set<string> preGenes = D.GetPreGenes(mutationPath, expNormalPath, expTumorPath, netPath);

    // Get the differential expression matrix;
    Driver::Diff diff = D.GetDiffMatrix(mutationPath, netPath, geneIDPath, expTumorPath, expNormalPath, preGenes);
    vector<string> majorGenes = diff.g;

    // Generate the differential expression network;
    SolfThreshold ast = GetOptimalCorrMatrix(diff.m);
    D.GetFinalMatrix(ast.corrMatrix, ast.connectivity, majorGenes);

    // Main body of DriverMP;
    D.ReadGenomicMatrix(mutationPath);
    D.MutScore();
    D.ReadCancerGene(referencePath);
    D.PairGeneMutation();
    map<vector<int>, double> pairIdMIF = D.ComputeDCIS(majorGenes, threadNum);
    D.OutDegree();
    D.ComputeSingleMIF_degree(pairIdMIF);
    D.OutputRank(savePath);

    return 0;
}
