from run_protocol import *
from generate_key import *
import math

def generate_param_file(num_ind, num_snp, num_cov, num_Miner, num_runtime_iteration):
    with open ('../param/param.txt', 'w') as file:
        line = "Num_Ind " + str(num_ind) + "\n"
        file.write(line)
        line = "Num_SNP " + str(num_snp) + "\n"
        file.write(line)
        line = "Num_Cov " + str(num_cov) + "\n"
        file.write(line)
        line = "Num_Miner " + str(num_Miner) + "\n"
        file.write(line)
        line = "Runtime_iteration " + str(num_runtime_iteration) + "\n"
        file.write(line)

num_cov = 9
num_runtime_iteration = 10
num_Miner = 3
precentage = 1

'''
   precentage specifies the percentage of the dataset that you want to use in the Geos process.
   Use the same values ​​for precentage in this file and the read_data file.
'''

num_ind = math.ceil(precentage*569)
lst_num_ind = [num_ind]
lst_num_snp = [21]

for i in lst_num_ind:
    for j in lst_num_snp:
        generate_param_file(i,j,num_cov,num_Miner,num_runtime_iteration)
        logfile_name = "../log/log_snp_" + str(j) +"_miner_" + str(num_Miner)
        log_file = open(logfile_name, 'w')
        key_generation_runtime(log_file)
        generate_PRG_keys(i, num_Miner)
        run_protocol(log_file)
        log_file.close() 

print("GWAS process successfully completed")

