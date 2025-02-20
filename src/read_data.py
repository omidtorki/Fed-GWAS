
import random
import numpy as np
import math 


def read_feature(num_ind, num_snp, num_cov):
    geno_data = {}
    pheno_data = {}
    cov_data = {}
    count = 0
    num_ind_total = 569
    precentage = 1 # 100%
    index = math.ceil(num_ind_total*precentage)
    chose = np.zeros(num_ind_total)
    chose[0:index] = 1
    np.random.shuffle(chose)
    with open("../data/geno.txt", "r") as f:
        for _ in range(num_ind_total):
            line = f.readline()
            if chose[_] == 0:
              continue
            count += 1
            key = "user_" + str(count)
            line = line.split()[:num_snp]
            line = np.array(line, dtype=np.float64)
            geno_data[key] = line
    
    print(f"number of individual: {len(geno_data)}")
    print(f"number of snp: {geno_data['user_1'].shape}")        
    
    
    count = 0
    with open("../data/pheno.txt", "r") as f:
        for _ in range(num_ind_total):
            line = f.readline()
            if chose[_] == 0:
              continue
            count += 1
            key = "user_" + str(count)
            #line = line.split()[1:]
            line = np.array(line, dtype=np.float64)
            pheno_data[key] = line
    print("len pheno: ", len(pheno_data))        
    
    intercept = np.ones(1)
    count = 0
    with open("../data/cov.txt", "r") as f:
        for _ in range(num_ind_total):
            line = f.readline()
            if chose[_] == 0:
              continue
            count += 1
            key = "user_" + str(count)
            line = line.split()[:num_cov]
            line = np.array(line, dtype=np.float64)
            line = np.insert(line, 0, intercept)
            cov_data[key] = line
            

    print(f"number of covariate: {cov_data['user_1'].shape}")
    return geno_data, pheno_data, cov_data





        
