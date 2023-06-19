#include"Driver.h"
void Driver::OutDegree() {
    for (iterPPI = NetWork.begin(); iterPPI != NetWork.end(); ++iterPPI) {
        nodeInt = iterPPI->first;
        list<PPI> listPPITemp = iterPPI->second;
        double maxWeight = 0, Degree = 0, sumWeight = 0;
        for (iterListPPI = listPPITemp.begin(); iterListPPI != listPPITemp.end(); ++iterListPPI) {
            if ((iterListPPI->Id_weight) > maxWeight) maxWeight = iterListPPI->Id_weight;
            sumWeight += iterListPPI->Id_weight;
        }
        Degree = iterPPI->second.size();
        mapGeneSumWeight.insert(pair<int, double>(nodeInt, sumWeight));
        if (Degree > 0) mapGeneDegree.insert(pair<int, double>(nodeInt, Degree));
        if (maxWeight > 0) mapGeneMaxWeight.insert(pair<int, double>(nodeInt, maxWeight));
    }
}

void Driver::ComputeSingleMIF_degree(map<vector<int>, double> pairIdMIF) {
    // 获取diff网络其他的成对节点；
    cout << "===================================== STEP 5: Separate paired DCISs ===================================" << endl;
    cout << "====                                                                                               ====" << endl;
    map<vector<int>, double> mapOtherPairIntScore;
    auto minPosExp = min_element(vecMaxExpPairMIF.begin(), vecMaxExpPairMIF.end());
    for (iterPPI = NetWork.begin(); iterPPI != NetWork.end(); ++iterPPI)
    {
        double internalWeight = 0, pairMIF = 0, SingleMIF = 0;
        intGeneId = iterPPI->first;
        auto ListPPITemp = iterPPI->second;
        for (iterListPPI = ListPPITemp.begin(); iterListPPI != ListPPITemp.end(); ++iterListPPI) {
            internalWeight = iterListPPI->Id_weight;
            intGeneId1 = iterListPPI->Id;

            vector<int> vecGeneID, vecGeneID1;
            vecGeneID.push_back(intGeneId);
            vecGeneID.push_back(intGeneId1);
            iterMapVecIntDouble = pairIdMIF.find(vecGeneID);
            vecGeneID1.push_back(intGeneId1);
            vecGeneID1.push_back(intGeneId);
            iterMapVecIntDouble1 = pairIdMIF.find(vecGeneID1);
            if (iterMapVecIntDouble != pairIdMIF.end() | iterMapVecIntDouble1 != pairIdMIF.end()) {
                continue;
            }

        }
    }
    for (iterPPI = NetWork.begin(); iterPPI != NetWork.end(); ++iterPPI)
    {
        double internalWeight = 0, pairMIF = 0, SingleMIF = 0;
        intGeneId = iterPPI->first;
        auto ListPPITemp = iterPPI->second;
        for (iterListPPI = ListPPITemp.begin(); iterListPPI != ListPPITemp.end(); ++iterListPPI) {
            internalWeight = iterListPPI->Id_weight;
            intGeneId1 = iterListPPI->Id;
            map<vector<int>, double>::iterator iterMapVecIntDouble3, iterMapVecIntDouble4;

            vector<int> vecGeneID, vecGeneID1;
            vecGeneID.push_back(intGeneId);
            vecGeneID.push_back(intGeneId1);
            iterMapVecIntDouble = pairIdMIF.find(vecGeneID);
            vecGeneID1.push_back(intGeneId1);
            vecGeneID1.push_back(intGeneId);
            iterMapVecIntDouble1 = pairIdMIF.find(vecGeneID1);
            if (iterMapVecIntDouble != pairIdMIF.end() | iterMapVecIntDouble1 != pairIdMIF.end()) {
                if (iterMapVecIntDouble != pairIdMIF.end()) { pairMIF = iterMapVecIntDouble->second; }
                if (iterMapVecIntDouble1 != pairIdMIF.end()) { pairMIF = iterMapVecIntDouble1->second; }
            }
            double degree = 1.0 * mapGeneSumWeight[intGeneId] / 1.0 * (mapGeneSumWeight[intGeneId] + mapGeneDegree[intGeneId1]);
            double tempScore = degree * pairMIF;
            if (tempScore > SingleMIF) SingleMIF = tempScore;
        }
        m_ofs << intGeneId << "\t" << SingleMIF << "\n";
        mapSingleMIFID.insert(pair<double, int>(SingleMIF, intGeneId));
        mapSingleMIFGene.insert(pair<double, string>(SingleMIF, strGeneName));

    }
    cout << "====                                              Done                                             ====" << endl;
    cout << "====                                                                                               ====" << endl;
    cout << "=======================================================================================================" << endl;
    cout << endl;
//    cout << "All MME score are seperated." << endl;
    m_ofs.close();
}

void Driver::OutputRank(string savePath) {
    map<double, int >::iterator iterDoubleInt;
    for (iterDoubleInt = mapSingleMIFID.begin(); iterDoubleInt != mapSingleMIFID.end(); ++iterDoubleInt) {
        vecIdSingleMIF.push_back(pair<int, double>(iterDoubleInt->second, iterDoubleInt->first));
    }
    sort(vecIdSingleMIF.begin(), vecIdSingleMIF.end(), cmp);
    m_ofs.open(savePath.c_str());
    m_ofs << "Rank" << "\t" << "GeneID" << "\t" << "GeneName" << "\t" <<
          "CGC" << "\n";
    int rank = 0, count100 = 0, count200 = 0, count300 = 0, count400 = 0, count500 = 0;
    vector<pair<int, double> >::iterator ite;
    for (ite = vecIdSingleMIF.begin(); ite != vecIdSingleMIF.end(); ++ite) {
        rank++;
        m_ofs << rank << "\t" << ite->first << "\t";
        if (mapNCBIToOfficial.end() != mapNCBIToOfficial.find(ite->first)) {
            m_ofs << mapNCBIToOfficial[ite->first] << "\t";
        } else {
            m_ofs << "NULL" << "\t";
        }
        for (int j = 0; j < 1; j++) {
            int flag1 = CancerGene[j].count(ite->first);
            if (flag1 == 1) {
                m_ofs << "Y" << "\t";
                if (rank < 100 && j == 0) count100++;
                else if (rank < 200 && j == 0) count200++;
                else if (rank < 300 && j == 0) count300++;
                else if (rank < 400 && j == 0) count400++;
                else if (rank < 501 && j == 0) count500++;
                else continue;
            }
            else m_ofs << "N" << "\t";
        }
        m_ofs << "\n";
    }
    cout << "=======================================================================================================" << endl;
    cout << "====                                                                                               ====" << endl;
    cout << "====                                            All Done                                           ====" << endl;
    cout << "====                                                                                               ====" << endl;
    cout << "=======================================================================================================" << endl;
//    cout << count100 << endl;
//    cout << count200 << endl;
//    cout << count300 << endl;
//    cout << count400 << endl;
//    cout << count500 << endl;
//    cout << count100 + count200 + count300 + count400 + count500 << endl;
//    cout << '\n';
    m_ofs.close();
}
