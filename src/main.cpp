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
#include <unistd.h>
#include <sys/stat.h>
#define STR_SEPTABLE "\t"
using namespace std;

const char * mut = "";
const char * tumor_exp = "";
const char * normal_exp = "";
const char * ppi = "";
bool Help = false;
bool Version=false;
const char * output = "DriverMP_results";
const char * threads = "1";
//const char * version="The current version is: v.1.0";
const char * short_opt = "m:u:n:p:ho:t:v";
int threadNum = 1;

struct option long_opt[]={
        {"mut",1,0,'m'},
        {"tumor_exp",1,0,'u'},
        {"normal_exp",1,0,'n'},
        {"ppi",1,0,'p'},
        {"help",0,0,'h'},
        {"output",1,0,'o'},
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
            <<"  --tumor_exp/-u <string>      "<<": Path to gene expression data (FPKM normalized) for specific cancer samples;"<<"\n"<<endl
            <<"  --normal_exp/-n <string>     "<<": Path to gene expression (FPKM normalized) of normal samples corresponding to the specific cancer;"<<"\n"<<endl
            <<"  --ppi/-p <string>            "<<": Path to Protein-Protein Interaction Networks."<<"\n"<<endl
            <<"** Optional **"<<"\n"<<endl
            <<"  --output/-o <string>         "<<": Output path, default: ./DriverMP_results/output.txt;"<<"\n"<<endl
            <<"  --threads/-t <int>           "<<": Number of threads to launch, default: 1;"<<"\n"<<endl
            <<"  --version/-v                 "<<": Show current version of DriverMP;"<<"\n"<<endl
            <<"  --help/-h                    "<<": help infomation;"<<"\n"<<endl
            <<"** Typical commands **"<<"\n"<<endl
            <<"A example of DriverMP  might be:"<<"\n"<<endl
            <<"  ./DriverMP DriverMP_example/Mutation_data.txt DriverMP_example/Gene_expresstion_tumor.txt DriverMP_example/Gene_expression_normal.txt DriverMP_example/HumanNet"<<"\n"<<endl;
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
                case 'o':
                    output=optarg;
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
    // 检查是否存在DriverMP_results文件夹，如果存在则删除该目录及其内容。
    char RM[100000];
    if( (access( output, 0 )) != -1 )
    {
        cout<<"[WARNING] : "<<output<<" exists already. It will be overwritten."<<endl;
        FILE * stream;
        strcpy(RM,"rm -r ");
        strcat(RM,output);
        stream = popen(RM,"r");
        pclose(stream);
        if (mkdir(output, 0777) != 0) std::cerr << "Failed to create folder." << std::endl;
    } else {
        if (mkdir(output, 0777) != 0) std::cerr << "Failed to create folder." << std::endl;
    }
    string geneIDPath = "DriverMP_example//geneID.txt";
    string referencePath = "DriverMP_example//Combination.txt";
    string mutationPath;
    string expTumorPath;
    string expNormalPath;
    string netPath;
    string savePath = string(output) + "//output.txt";

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
