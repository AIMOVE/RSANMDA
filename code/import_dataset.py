import numpy as np
import random

def import_dis_mir_vec(data_path):
    print('\nDisease comprehensive similarity ...\n')
    DSS1 = np.loadtxt(data_path + 'DSS1.txt')
    DSS2 = np.loadtxt(data_path + 'DSS2.txt')
    DSS = (DSS1 + DSS2) / 2
    DGS = np.loadtxt(data_path + 'DGS.txt')
    ID = np.zeros(shape = (DSS.shape[0], DSS.shape[1]))
    for i in range(DSS.shape[0]):
        for j in range(DSS.shape[1]):
            if DSS[i][j] == 0:
                ID[i][j] = DGS[i][j]
            else:
                ID[i][j] = DSS[i][j]
    
    print('\nmiRNA comprehensive similarity ...\n')
    MFS = np.loadtxt(data_path + 'MFS.txt')
    MSS = np.loadtxt(data_path + 'MSS.txt')
    MFSS = (MFS + MSS) / 2
    MGS = np.loadtxt(data_path + 'MGS.txt')
    IM = np.zeros(shape = (MFSS.shape[0], MFSS.shape[1]))
    for i in range(MFSS.shape[0]):
        for j in range(MFSS.shape[1]):
            if MFSS[i][j] == 0:
                IM[i][j] = MGS[i][j]
            else:
                IM[i][j] = MFSS[i][j]
                
    print('\nconstruct positive pairs and unlabelled pairs ...\n')
    A = np.zeros(shape = (DSS.shape[0], MFSS.shape[1]))
    asso_file =  xlrd.open_workbook(data_path + 'Known disease-miRNA association number.xlsx')
    asso_pairs = asso_file.sheets()[0]
    for i in range(5430):
        asso = asso_pairs.row_values(i)
        m = int(asso[0])
        n = int(asso[1])
        A[n-1, m-1] = 1
        
    known=[]
    unknown=[]                         
    for x in range(383):
        for y in range(495):
            if A[x,y]==0:                 
                unknown.append((x,y))
            else:
                known.append((x,y))
    
    posi_list = []
    for i in range(5430):
        posi = ID[known[i][0],:].tolist() + IM[known[i][1],:].tolist() + [1, 0]
        posi_list.append(posi)
    unlabelled_list = []
    for j in range(184155):
        unlabelled=ID[unknown[j][0],:].tolist()+IM[unknown[j][1],:].tolist() + [0, 1]
        unlabelled_list.append(unlabelled)
    
    random.shuffle(posi_list)
    random.shuffle(unlabelled_list)
    return posi_list, unlabelled_list


# Divide the same number of unlabeled samples as the positive samples
def spliting_unlabelled_data(unlabelled_data, test_num):
    unlabelled_train_data = unlabelled_data[test_num:]
    unlabelled_cv_data = unlabelled_data[:test_num]
    return unlabelled_train_data, unlabelled_cv_data
