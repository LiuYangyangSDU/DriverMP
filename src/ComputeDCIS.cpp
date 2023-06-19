#include"Driver.h"
#include <thread>
#include <mutex>
#include <omp.h>



void Driver::ReadCancerGene(string referencePath) {
    string path[5] = {referencePath };
    ifstream ifs;
    string strTemp;
    int geneID;
    for (int i = 0; i < 1; ++i) {
        ifs.open(path[i].c_str());
        while (getline(ifs, strTemp)) {
            geneID = atoi(strTemp.c_str());
            CancerGene[i].insert(geneID);
        }
        ifs.close();
    }
}

//bool Driver::ReadNetWork(string netPath) {
//    MinWeight = 1.0;
//    m_ifs.open(netPath.c_str());
//    if (!m_ifs.is_open()) {
//        cout << "\n" << "Fail to open NetWork data: " << netPath << "！\n"
//             << "please check if your PPI network file path is correct and try again later!" << "\n";
//        return false;
//    }
//
//    while (getline(m_ifs, m_strTemp)) {
//        StringSplit(m_strTemp, STR_SEPTABLE, m_vecTempString);
//        for (int i = 0; i < 2; i++) {
//            m_vecTempInt.push_back(atoi(m_vecTempString[i].c_str()));
//        }
//        FunctionalNet.push_back(m_vecTempInt);//Node-to-Node
//        m_vecTempInt.clear();
//
//
//        intGeneId = atoi(m_vecTempString[0].c_str());
//        intGeneId1 = atoi(m_vecTempString[1].c_str());
//        GeneGeneWeight = atof(m_vecTempString[2].c_str());
//        if (GeneGeneWeight < MinWeight) MinWeight = GeneGeneWeight;
//
//        iterPPI = NetWork.find(intGeneId);//NetWork是存了每个基因节点对应相连的其他基因数目
//        if (iterPPI != NetWork.end()) {
//            list<PPI>& temp1 = NetWork[intGeneId];
//            temp1.push_back(PPI(intGeneId1, GeneGeneWeight));
//        }
//        else {
//            list<PPI> pid_w11;
//            pid_w11.push_back(PPI(intGeneId1, GeneGeneWeight));
//            NetWork[intGeneId] = pid_w11;
//        }
//
//
//        iterPPI1 = NetWork.find(intGeneId1);
//        if (iterPPI1 != NetWork.end()) {
//            list<PPI>& temp2 = NetWork[intGeneId1];
//            temp2.push_back(PPI(intGeneId, GeneGeneWeight));
//        }
//        else {
//            list<PPI> pid_w12;
//            pid_w12.push_back(PPI(intGeneId, GeneGeneWeight));
//            NetWork[intGeneId1] = pid_w12;
//        }
//    }
//    m_ifs.close();
//    cout << "Edges number in PPI networks: " << FunctionalNet.size() << ", Nodes number: " << NetWork.size() << endl;
//    return false;
//}


//void Driver::ReadGeneExp(string expPath) {
//    MinGeneExp = 1000000.0;
//    m_ifs.open(expPath.c_str());
//    if (!m_ifs.is_open()) {
//        cout << "Fail to open gene expression file, please check the file path." << endl;
//    }
//    int count = 0;
//    while (getline(m_ifs, m_strTemp)) {
//        count = count + 1;
//        if (count == 1) { continue;}
//        StringSplit(m_strTemp, STR_SEPTABLE, m_vecTempString);
//        // Calculate the expression score for major genes.
//        doubleGeneExp = 0;
//        for (auto iter : m_vecTempString) {
//            doubleGeneExp += atof(iter.c_str());
//        }
//        doubleGeneExp = doubleGeneExp / (m_vecTempString.size() - 1);
//        if (doubleGeneExp < MinGeneExp) { MinGeneExp = doubleGeneExp; }
//        strGeneName = m_vecTempString[0].c_str();
////        mapGeneExp[strGeneName] = doubleGeneExp;
//
//        map<string, int>::iterator it = mapOfficialToNCBI.find(strGeneName);
//        if (it != mapOfficialToNCBI.end()) {
//            mapIdExp.insert(pair<int, double>(it->second, doubleGeneExp));
//        }
//    }
//    m_ifs.close();
//}


//bool Driver::ReadGeneExpNetWork(string netPath) {
//    MinDist = 100;
//    m_ifs.open(netPath.c_str());
//    if (!m_ifs.is_open()) {
//        cout << "\n" << "Fail to open GeneExpNetWork data: " << netPath << "！\n"
//             << "please check if your input is correct and try again later!" << "\n";
//        return false;
//    }
//
//    while (getline(m_ifs, m_strTemp)) {
//        StringSplit(m_strTemp, STR_SEPTABLE, m_vecTempString);
//        strGeneName = m_vecTempString[0].c_str();
//        strGeneName1 = m_vecTempString[1].c_str();
//        iterStringInt = mapGeneId.find(strGeneName);
//        iterStringInt1 = mapGeneId.find(strGeneName1);
//        if ((iterStringInt == mapGeneId.end()) | (iterStringInt1 == mapGeneId.end())) {
//            continue;
//        }
//        intGeneId = iterStringInt->second;
//        intGeneId1 = iterStringInt1->second;
//        GeneGeneDist = atof(m_vecTempString[2].c_str());
//
//        // 对第一个基因；
//        if (ExpNetWorkNew.find(intGeneId) == ExpNetWorkNew.end()) {
//            map<int, double> tempMap;
//            tempMap[intGeneId1] = GeneGeneWeight;
//            ExpNetWorkNew[intGeneId] = tempMap;
//        } else {
//            ExpNetWorkNew[intGeneId][intGeneId1] = GeneGeneWeight;
//        }
//        // 对第二个基因；
//        if (ExpNetWorkNew.find(intGeneId1) == ExpNetWorkNew.end()) {
//            map<int, double> tempMap;
//            tempMap[intGeneId] = GeneGeneWeight;
//            ExpNetWorkNew[intGeneId1] = tempMap;
//        } else {
//            ExpNetWorkNew[intGeneId1][intGeneId] = GeneGeneWeight;
//        }
//
//        // list<PPI>形式：第一个基因；
//        iterPPI = ExpNetWork.find(intGeneId);//NetWork是存了每个基因节点对应相连的其他基因数目
//        if (iterPPI != ExpNetWork.end()) {
//            list<PPI>& temp1 = ExpNetWork[intGeneId];
//            temp1.push_back(PPI(intGeneId1, GeneGeneDist));
//        }
//        else {
//            list<PPI> pid_w11;
//            pid_w11.push_back(PPI(intGeneId1, GeneGeneDist));
//            ExpNetWork[intGeneId] = pid_w11;
//        }
//        // 对第二个基因；
//        iterPPI = ExpNetWork.find(intGeneId1);
//        if (iterPPI != ExpNetWork.end()) {
//            list<PPI>& temp2 = ExpNetWork[intGeneId1];
//            temp2.push_back(PPI(intGeneId, GeneGeneDist));
//        }
//        else {
//            list<PPI> pid_w12;
//            pid_w12.push_back(PPI(intGeneId, GeneGeneDist));
//            ExpNetWork[intGeneId1] = pid_w12;
//        }
//
//    }
//    cout << ExpNetWork.size() << endl;
////    cout << ExpNetWorkNew.size() << endl;
//    m_ifs.close();
//    return false;
//}

void Driver::ComputeDCISCircle(int gene1, int gene2, int index) {
//    lock_guard<mutex> lock(m);
    // 两个基因的突变分数和表达分数；
    double nodeIntMut = mapIdMutation[gene1];
    double nodeInt1Mut = mapIdMutation[gene2];
    double nodeIntExp = mapIdExp[gene1];
    double nodeInt1Exp = mapIdExp[gene2];

    // 获得成对基因在G_PPI的局部网络结构；
    auto  iterPPI_1 = NetWork.find(gene1);
    auto  iterPPI_2 = NetWork.find(gene2);
    list<PPI> ListPPI, ListPPI1, mapTransList;
    if (iterPPI_1 != NetWork.end() && iterPPI_2 != NetWork.end()) {
        ListPPI = iterPPI_1->second;
        ListPPI1 = iterPPI_2->second;
    }
//    double maxWeight = 0;
    list<PPI>::iterator iterListPPI_1, iterListPPI_2;
    map<int, double> listTransMap;
    for (auto & link:ListPPI) listTransMap.insert(pair<int, double>(link.Id, link.Id_weight));
    for (auto & link:ListPPI1) {
        auto iter = listTransMap.find(link.Id);
        if (iter == listTransMap.end()) listTransMap.insert(pair<int, double>(link.Id, link.Id_weight));
        else listTransMap[link.Id] = max(link.Id_weight, iter->second);
    }
    for (auto & iter : listTransMap) mapTransList.emplace_back(iter.first, iter.second);

    //融合双节点后的基因表达网络
    auto  iterExpPPI = ExpNetWork.find(gene1);
    auto  iterExpPPI1 = ExpNetWork.find(gene2);
    list<PPI> ListExpPPI, ListExpPPI1, mapExpTransList;
    if (iterExpPPI != ExpNetWork.end() and iterExpPPI1 != ExpNetWork.end()) {
        ListExpPPI = iterExpPPI->second;
        ListExpPPI1 = iterExpPPI1->second;
    }
    list<PPI>::iterator iterListExpPPI, iterListExpPPI1;
//    double minDist = 0;
    map<int, double> listExpTransMap;
    for (auto & link:ListExpPPI) listExpTransMap.insert(pair<int, double>(link.Id, link.Id_weight));
    for (auto & link:ListExpPPI1) {
        auto iter = listExpTransMap.find(link.Id);
        if (iter == listExpTransMap.end()) listExpTransMap.insert(pair<int, double>(link.Id, link.Id_weight));
        else listExpTransMap[iter->first] = max(iter->second, link.Id_weight);
    }
    for (auto & iter : listExpTransMap) {
        mapExpTransList.emplace_back(iter.first, iter.second);
    }

    //计算PPI网络中成对基因节点的pairMIF值
    double internalWeight = 0;
    auto iterIntListPPI = NetWork.find(gene1);
    if (iterIntListPPI != NetWork.end()) {
        list<PPI> listPPITemp = iterIntListPPI->second;
        list<PPI>::iterator  iterListPPITemp;
        for (iterListPPITemp = listPPITemp.begin(); iterListPPITemp != listPPITemp.end(); ++iterListPPITemp) {
            if (iterListPPITemp->Id == gene2) {
                internalWeight = iterListPPITemp->Id_weight;
            }
        }
    }
    // 计算G_PPI中AS值；
    double MaxMutPairMIF = 0;
    // 如果该成对节点非孤立；
    if (!mapTransList.empty()) {
        for (auto & link:mapTransList) {
            double linkMutation = mapIdMutation[link.Id];
            double linkWeight = link.Id_weight;
            double TempScore = (nodeIntMut * nodeInt1Mut * pow(internalWeight, 2)) * linkMutation * pow(linkWeight, 2);
            if (TempScore > MaxMutPairMIF) MaxMutPairMIF = TempScore;
        }

    }
    else MaxMutPairMIF = (nodeIntMut * nodeInt1Mut * pow(internalWeight, 2)) * minScore * pow(MinWeight, 2);
    listExpTransMap.clear();
    mapTransList.clear();

    double MaxExpPairMIF = 0;
    if (!mapExpTransList.empty()) {//融合后不是独立节点
        vector<double> vecMutExp;
        for (auto & link:mapExpTransList) {
            double linkExp = mapIdExp[link.Id]; double linkDist = link.Id_weight;
            double TempExpScore = (max(nodeIntExp, nodeInt1Exp) * linkExp) / (pow((1 - linkDist), 2));
            if (TempExpScore > MaxExpPairMIF) MaxExpPairMIF = TempExpScore;
        }

        listExpTransMap.clear();
        mapExpTransList.clear();

        double FinalScore = MaxMutPairMIF * MaxExpPairMIF;
//        std::mutex mtx;
//        lock_guard<mutex> lock(mtx);
        pairMIF[index] = FinalScore;

    }
}

vector<vector<int>> Driver::GetRequiredPairGenes(vector<string> majorGenes) {
    vector<vector<int>> pairGenes;
    for (auto & triple : FunctionalNet) {
        int gene1 = triple[0];
        int gene2 = triple[1];

        // 判断两个基因是否为重要基因；
        auto findIter = find(majorGenes.begin(), majorGenes.end(), mapNCBIToOfficial[gene1]);
        auto findIter1 = find(majorGenes.begin(), majorGenes.end(), mapNCBIToOfficial[gene2]);
        if (findIter == majorGenes.end() or findIter1 == majorGenes.end()) continue;
        // 判断两个基因是否在表达网络中连接；
        int flag = 0;
//        if (!ExpNetWorkNew[gene1].count(gene2)) continue;
        for (auto & i : ExpNetWork[gene1]) {
            if (i.Id == gene2) flag = 1;
        }
        if (flag == 0) continue;

        vector<int> vecGeneInt;
        vecGeneInt.push_back(gene1);
        vecGeneInt.push_back(gene2);
        pairGenes.push_back(vecGeneInt);
    }
//    cout << pairGenes.size() << endl;
    return pairGenes;
}


map<vector<int>, double> Driver::ComputeDCIS(vector<string> majorGenes, int threadNum) {
    omp_set_num_threads(threadNum);
    cout << "===================================== STEP 4: Compute paired DCISs ====================================" << endl;
    cout << "====                                                                                               ====" << endl;
    // calculating the mutation score and network of pair genes
//    cout << "Computing DCIS for paired genes, it may take a few minutes." << "\n";
    map<vector<int>, double> pairIdMIF;
    vector<vector<int>> pairGenes = GetRequiredPairGenes(std::move(majorGenes));

    for (int i=0; i < pairGenes.size(); i++) pairMIF.push_back(0);

    // 可以运行；
    #pragma omp parallel for
    for (int index = 0; index < pairGenes.size(); index++) {
        ComputeDCISCircle(pairGenes[index][0], pairGenes[index][1], index);
//        int thread_id = omp_get_thread_num();
//        std::cout << "Current thread ID: " << thread_id << std::endl;
    }

    for (int i = 0; i < pairGenes.size(); i++) pairIdMIF[pairGenes[i]] = pairMIF[i];
//    cout << "The DCISs for all paired genes are calculated." << endl;
    cout << "====                                              Done                                             ====" << endl;
    cout << "====                                                                                               ====" << endl;
    cout << "=======================================================================================================" << endl;
    cout << endl;
    return pairIdMIF;
}


void Driver::PairGeneMutation() {

    for (unsigned int i = 0; i < FunctionalNet.size(); i++) {

        nodeInt = FunctionalNet[i][0];
        nodeInt1 = FunctionalNet[i][1];


        // Computing the mutation score for the paired genes.
        iterIntDouble = mapIdMutation.find(nodeInt);
        iterIntDouble1 = mapIdMutation.find(nodeInt1);

        vector<int> vecTempInt;
        vecTempInt.push_back(nodeInt);
        vecTempInt.push_back(nodeInt1);
        if (iterIntDouble == mapIdMutation.end() && iterIntDouble1 == mapIdMutation.end()) {

            PairGeneMut.insert(pair<vector<int>, double>(vecTempInt, minScore * minScore));
        }
        else if (iterIntDouble != mapIdMutation.end() && iterIntDouble1 == mapIdMutation.end()) {

            PairGeneMut.insert(pair<vector<int>, double>(vecTempInt, iterIntDouble->second * minScore));
        }
        else if (iterIntDouble == mapIdMutation.end() && iterIntDouble1 != mapIdMutation.end()) {

            PairGeneMut.insert(pair<vector<int>, double>(vecTempInt, minScore * iterIntDouble1->second));
        }
        else {
            PairGeneMut.insert(pair<vector<int>, double>(vecTempInt, iterIntDouble->second * iterIntDouble1->second));
        }

    }
}
