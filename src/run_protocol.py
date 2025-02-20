from run_gwas import *
from read_data import *
import time
import sys
import math
import random
import numpy as np


def restore_PRG_keys(num_ind, num_Miner):
    
    #Restore PRG keys Analyzer
    keys = []
    with open("../key/PRG_keys_Analyzer.bin", 'rb') as file:
        while True:
            key = file.read(16)
            if not key:
                break
            integer_key = int.from_bytes(key, byteorder='big')
            keys.append(integer_key)
    PRG_keys_Analyzer = np.reshape( keys, num_ind)
    
    
    
    
    #Restore PRG keys Miner
    keys = []
    with open("../key/PRG_keys_Miner.bin", 'rb') as file:
        while True:
            key = file.read(16)
            if not key:
                break
            integer_key = int.from_bytes(key, byteorder='big')
            keys.append(integer_key)
    PRG_keys_Miner_temp = np.reshape( keys, (num_Miner,math.ceil(num_ind/2)))
    PRG_keys_Miner = np.zeros((num_Miner, num_ind))
    PRG_keys_Miner_sign = np.zeros((num_Miner, num_ind))
    PRG_keys_residue = np.zeros(num_Miner)
    lst_ind = list(range(num_ind))
    for i in range(num_Miner):
        for j in range(math.floor(num_ind/2)):
            index1 = lst_ind.pop(random.randrange(len(lst_ind)))
            index2 = lst_ind.pop(random.randrange(len(lst_ind)))
            PRG_keys_Miner[i, index1] = PRG_keys_Miner_temp[i, j]
            PRG_keys_Miner_sign[i, index1] = 0
            PRG_keys_Miner[i, index2] = PRG_keys_Miner_temp[i, j]
            PRG_keys_Miner_sign[i, index2] = 1
        if num_ind % 2 == 1:   # Number of individual (Genome-owner) is odd
            index = lst_ind.pop(random.randrange(len(lst_ind)))
            PRG_keys_Miner[i, index] = PRG_keys_Miner_temp[i, j+1]
            PRG_keys_Miner_sign[i, index] = 0
            PRG_keys_residue[i] = PRG_keys_Miner_temp[i, j+1]
        lst_ind = list(range(num_ind))

        
    if num_ind % 2 == 1:
      PRG_keys_residue = list(PRG_keys_residue)
    else:
      PRG_keys_residue = []
             
    
    return PRG_keys_Analyzer, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_residue 
 



def run_protocol(log_file):
    param_file = open("../param/param.txt", "r")
    num_ind = int(param_file.readline().split()[1])
    num_snp = int(param_file.readline().split()[1])
    num_cov = int(param_file.readline().split()[1])
    num_Miner = int(param_file.readline().split()[1])
    param_file.close()
    geno_data, pheno_data, cov_data = read_feature(num_ind, num_snp, num_cov)
    PRG_keys_Analyzer, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_residue  = restore_PRG_keys(num_ind, num_Miner)
    print("PRG_keys_Analyzer.shape", PRG_keys_Analyzer.shape, "  ", "PRG_keys_Miner.shape", PRG_keys_Miner.shape, "  ", "PRG_keys_residue.shape", len(PRG_keys_residue))
    z_score, pval = FL_parallel_gwas(geno_data, pheno_data, cov_data, PRG_keys_Analyzer, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_residue, log_file)
    f = open("../out/pval.txt", "w")
    for i in range(pval.shape[0]):
        line = str(z_score[i]) + " " + str(pval[i]) + "\n"
        f.write(line)
    
    f.close()

