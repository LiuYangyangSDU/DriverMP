#include"Driver.h"
//void Driver::ReadGeneID(string geneIDPath) {//读取GeneId数据
//    ifstream ifs;
//    string str_temp;
//    vector<string> vecTempStr;
//    int geneID;
//    ifs.open(geneIDPath.c_str());
//    if (!ifs.is_open()) {
//        cout << "Fail to open geneID" << endl;
//    }
//    while (getline(ifs, str_temp)) {
//        StringSplit(str_temp, STR_SEPTABLE, vecTempStr);
//        geneID = atoi(vecTempStr[1].c_str());
//        mapGeneId[vecTempStr[0]] = geneID;
//        mapIdGene[geneID] = vecTempStr[0];
//    }
//    ifs.close();
//}


bool Driver::ReadGenomicMatrix(string mutationPath) {//读取突变矩阵
    cout << "==================================== STEP 3: Compute Mutation score ===================================" << endl;
    cout << "====                                                                                               ====" << endl;
    ifstream ifs;
    string strTemp;
    vector<string> vecStringTemp;
    ifs.open(mutationPath.c_str());
    if (!ifs.is_open()) {
        cout << "\n\n" << "Fail to open the MutaionPath" << mutationPath << "!\n"
             << " please check if your input is correct and try again later! " << "\n";
        return false;
    }

    getline(ifs, strTemp);
    StringSplit(strTemp, STR_SEPTABLE, vecStringTemp);
    intSampleNumber = vecStringTemp.size() - 1;//获取样本个数
    if (intSampleNumber == 1) {
        cout << "\n" << "Type:Mutation frequancy data..." << "\n";
        while (getline(ifs, strTemp)) {
            StringSplit(strTemp, STR_SEPTABLE, vecStringTemp);
            strGeneName = vecStringTemp[0].c_str();
            float GeneMut = atof(vecStringTemp[1].c_str());
            if (GeneMut != 0.0f) {//GeneMut不为0.0f表示不为0，即发生突变的基因
                mapGeneMut[strGeneName] = GeneMut;
            }
        }
        ifs.close();
        return true;
    }

//    cout << "Somatic mutation matrix's sample number: " << intSampleNumber << "\n";
    cout << "====                                 Mutation Patients' Number: " << intSampleNumber << "                                ====" << endl;
    vector<int> vecIntTemp;
    string s1 = "1", s2 = "1.0";
    while (getline(ifs, strTemp)) {
        StringSplit(strTemp, STR_SEPTABLE, vecStringTemp);
        strGeneName = vecStringTemp[0].c_str();
        for (int i = 0; i < intSampleNumber; ++i) {
            if (vecStringTemp[i + 1].c_str() == s1 | vecStringTemp[i + 1] == s2) {
                vecIntTemp.push_back(i);//存放基因突变的样本
                iterIntVecString = mapSampleMutationGene.find(i);//mapSampleGene存放样本i中有哪些基因发生突变
                if (iterIntVecString == mapSampleMutationGene.end()) {//如果样本i不在mapSampleGene中就插入样本i
                    vecGeneName.push_back(strGeneName);
                    mapSampleMutationGene[i] = vecGeneName;
                    vecGeneName.clear();
                }
                else {
                    mapSampleMutationGene[i].push_back(strGeneName);//如果样本i已经存在，就插入突变的基因
                }
            }
        }
        if (!vecIntTemp.empty()) {
            mapGeneMutationSample[strGeneName] = vecIntTemp;//记录基因j在哪些样本中发生了突变
            vecIntTemp.clear();
        }
        else {
            iterStringInt = mapOfficialToNCBI.find(strGeneName);
            if (iterStringInt != mapOfficialToNCBI.end()) vecUnMutationGene.push_back(strGeneName);
        }
    }
    ifs.close();
    cout << "====                                                                                               ====" << endl;
    cout << "=======================================================================================================" << endl;
    cout << endl;
    return false;
}


void Driver::MutScore() {

    for (iterIntVecString = mapSampleMutationGene.begin(); iterIntVecString != mapSampleMutationGene.end(); ++iterIntVecString) {
        int total = 0;
        int sampleID = iterIntVecString->first;
        vecGeneName = iterIntVecString->second;
        for (iterVecString = vecGeneName.begin(); iterVecString != vecGeneName.end(); ++iterVecString) {
            total = total + mapGeneMutationSample[(*iterVecString)].size();
        }
        mapSampleTotal[sampleID] = total;//一个样本中所有基因突变总数之和
    }

    for (iterStringVecInt = mapGeneMutationSample.begin(); iterStringVecInt != mapGeneMutationSample.end(); ++iterStringVecInt) {
        strGeneName = iterStringVecInt->first;
        m_vecTempInt = iterStringVecInt->second;

        double tempScore = 0;
        for (iterVecInt = m_vecTempInt.begin(); iterVecInt != m_vecTempInt.end(); ++iterVecInt) {
            int sampleID = *iterVecInt;
            tempScore = tempScore + (double)m_vecTempInt.size() / mapSampleTotal[sampleID];
        }
        mapGeneMutation[strGeneName] = tempScore;//每个突变基因的mut分数

    }

    for (iterStringDouble = mapGeneMutation.begin(); iterStringDouble != mapGeneMutation.end(); ++iterStringDouble) {
        vecMutationScore.push_back(pair<string, double>(iterStringDouble->first, iterStringDouble->second));
    }
    sort(vecMutationScore.begin(), vecMutationScore.end());

    minScore = vecMutationScore.front().second;


    for (iterVecString = vecUnMutationGene.begin(); iterVecString != vecUnMutationGene.end(); ++iterVecString) {
        strGeneName = (*iterVecString);
        vecMutationScore.push_back(pair<string, double>(strGeneName, minScore));
        mapGeneMutation[strGeneName]=minScore;
    }

    for (vecPairMutationScore = vecMutationScore.begin(); vecPairMutationScore != vecMutationScore.end(); ++vecPairMutationScore) {
        strGeneName = vecPairMutationScore->first;
        iterStringInt = mapOfficialToNCBI.find(strGeneName);
        if (iterStringInt != mapOfficialToNCBI.end())//end终止判断
            mapIdMutation.insert(pair<int, double>(iterStringInt->second, vecPairMutationScore->second));
        continue;
    }
}