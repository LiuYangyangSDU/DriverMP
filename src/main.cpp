#include<iostream>
#include"Driver.h"
#include"DiffNetwork.h"
#include <vector>
#include "Preprocess.h"
#include <unordered_set>
#include <chrono>
#include <sstream>
#include <cstring>
#include <getopt.h>
#define STR_SEPTABLE "\t"
using namespace std;

const char * mut = "";
const char * tumor_exp = "";
const char * normal_exp = "";
const char * ppi = "";
bool Help = false;
bool Version=false;
const char * threads="1";
//const char * version="The current version is: v.1.0";
const char * short_opt="m:u:n:p:ht:v";
int threadNum = 1;

struct option long_opt[]={
        {"mut",1,0,'m'},
        {"tumor_exp",1,0,'u'},
        {"normal_exp",1,0,'n'},
        {"ppi",1,0,'p'},
        {"help",0,0,'h'},
        {"threads",1,0,'t'},
        {"version",0,0,'v'},
        {0,0,0,0}
};

string  usage(){
    stringstream use_info;
    use_info<<endl<<"================================================================"<<endl;
    use_info<<"Usage:"<<"\n"<<endl;
    use_info<<"  DriverMP [options] -m <mutation_file> -t <tumor_expression_file> -n <normal_expression_file> -p <PPI_network_file>"<<endl;
    use_info<<""<<endl;
    use_info<<"** Required **"<<"\n"<<endl
            <<"  --mut/-m <string>            "<<": Path to format-compliant non-silent somatic mutation data;"<<"\n"<<endl
            <<"  --tumor_exp/-u <string>      "<<": Path to RNA-Seq data (FPKM normalized) for specific cancer samples;"<<"\n"<<endl
            <<"  --normal_exp/-n <string>     "<<": Path to RNA-Seq data (FPKM normalized) of normal samples corresponding to the specific cancer;"<<"\n"<<endl
            <<"  --ppi/-p <string>            "<<": Path to Protein-Protein Interaction Networks."<<"\n"<<endl
            <<"** Optional **"<<"\n"<<endl
            <<"  --threads/-t <int>           "<<": Number of threads to launch, default: 1;"<<"\n"<<endl
            <<"  --version/-v                 "<<": Show current version of DriverMP;"<<"\n"<<endl
            <<"  --help/-h                    "<<": help infomation;"<<"\n"<<endl
            <<"** Typical commands **"<<"\n"<<endl
            <<"A example of DriverMP  might be:"<<"\n"<<endl
            <<"  ./DriverMP -m DriverMP_example//Mutation_data.txt -u DriverMP_example//Gene_expresstion_tumor.txt -n DriverMP_example//Gene_expression_normal.txt -p DriverMP_example//HumanNet"<<"\n"<<endl;
    use_info<<"================================================================"<<endl;
    return use_info.str();
}

string print_version() {
    stringstream use_info;
    use_info<<endl<<"================================================================"<<endl;
    use_info<<"The current version is: v.1.0"<<endl;
    use_info<<"================================================================"<<endl;
    return use_info.str();
}

int parse_options(int argc,char*argv[]){
    int option_index=0;
    int next_option;
    while((next_option=getopt_long(argc,argv,short_opt,long_opt,&option_index)) != -1){
        switch(next_option){
            case 'm':
                    mut=optarg;
                    break;
            case 'u':
                    tumor_exp=optarg;
                    break;
            case 'n':
                    normal_exp=optarg;
                    break;
            case 'p':
                    ppi=optarg;
                    break;
                case 'h':
                    Help=true;
                    break;
            case 'v':
            Version=true;
            break;
                case 't':
                    threads=optarg;
                    break;
                case '?':
                    cout<<usage();
                    exit(1);
        }
    }
    if(Help)
    {
        cout<<usage()<<endl;
        exit(1);
    }
    if(Version)
    {
        cout<<print_version()<<endl;
        exit(1);
    }
    if(mut == "")
    {
        cout<<"Error : -m option needs an argument!! "<<endl;
        cout<<usage();
        exit(1);
    }
    if(tumor_exp == "")
    {
        cout<<"Error : -u option needs an argument!! "<<endl;
        cout<<usage();
        exit(1);
    }
    if(normal_exp == "")
    {
        cout<<"Error : -n option needs an argument!! "<<endl;
        cout<<usage();
        exit(1);
    }
    if(ppi == "")
    {
        cout<<"Error : -p option needs an argument!! "<<endl;
        cout<<usage();
        exit(1);
    }

    return 0;
}


int main(int argc, char* argv[]) {
    parse_options(argc,argv);
    string geneIDPath = "DriverMP_example//geneID.txt";
    string referencePath = "DriverMP_example//CGC.txt";
    string mutationPath;
    string expTumorPath;
    string expNormalPath;
    string netPath;
    string savePath = "Output//output.txt";

    // Preprocessing the original genes;
    Driver D;
    D.mapGeneID(geneIDPath);
    unordered_set<string> preGenes = D.GetPreGenes(mut, normal_exp, tumor_exp, ppi);

    // Get the differential expression matrix;
    Driver::Diff diff = D.GetDiffMatrix(mut, ppi, geneIDPath, tumor_exp, normal_exp, preGenes);
    vector<string> majorGenes = diff.g;

    // Generate the differential expression network;
    SolfThreshold ast = GetOptimalCorrMatrix(diff.m);
    D.GetFinalMatrix(ast.corrMatrix, ast.connectivity, majorGenes);

    // Main body of DriverMP;
    D.ReadGenomicMatrix(mut);
    D.MutScore();
    D.ReadCancerGene(referencePath);
    D.PairGeneMutation();
    map<vector<int>, double> pairIdMIF = D.ComputeDCIS(majorGenes, std::stoi(threads));
    D.OutDegree();
    D.ComputeSingleMIF_degree(pairIdMIF);
    D.OutputRank(savePath);

    return 0;
}
