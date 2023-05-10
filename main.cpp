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

int main(int argc, char* argv[]) {
    string geneIDPath = "Data//geneID.txt";
    string referencePath = "Data//CGC.txt";
    string mutationPath;
    string expTumorPath;
    string expNormalPath;
    string netPath;
    string savePath = "Output//output.txt";

    // default: BRCA
    if (argc == 1) {
        mutationPath = "Data//BRCA";
        expTumorPath = "Data//BRCA_TCGA_FPKM_tumor.txt";
        expNormalPath = "Data//BRCA_TCGA_FPKM_normal.txt";
        netPath = "Data//HumanNet";
    }
    if (argc == 5) {
        if (strcmp(argv[4], "S") == 0) {
            mutationPath = argv[1];
            expTumorPath = argv[2];
            expNormalPath = argv[3];
            netPath = "Data//STRINGv10";
        }
        else if (strcmp(argv[4], "H") == 0) {
            mutationPath = argv[1];
            expTumorPath = argv[2];
            expNormalPath = argv[3];
            netPath = "Data//HumanNet";
        }
        else {
            mutationPath = argv[1];
            expTumorPath = argv[2];
            expNormalPath = argv[3];
            netPath = argv[4];
        }
    }
    if (argc > 1 and argc < 5){
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
    map<vector<int>, double> pairIdMIF = D.ComputeDCIS(majorGenes);
    D.OutDegree();
    D.ComputeSingleMIF_degree(pairIdMIF);
    D.OutputRank(savePath);

    return 0;
}
