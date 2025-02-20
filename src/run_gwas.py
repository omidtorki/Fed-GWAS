import numpy as np
import math 
import scipy.stats as stats
from Pyfhel import Pyfhel
from PRG import *
import time



def fisherscoring(cov_data, pheno_data, bfv_encoder, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_residue, itter, log_file):
    fisherscoring_runtime = -1
    
    k = cov_data['user_1'].shape[0]
    global_beta = np.zeros((k,1))
    X_T_W_X_global = np.zeros((k,k))
    X_T_W_v_global = np.zeros((k,1))
    

    
    for t in range(itter):
        print(f"fisherscoring iteration {t+1}")

        for key in cov_data:
            start = time.time()
            local_parameter[key][2] = np.log(local_parameter[key][0]/(1-local_parameter[key][0]))+ (np.reshape(pheno_data[key], (1,1))-local_parameter[key][0])/local_parameter[key][1]
            X_T_W_local = np.matmul(np.reshape(cov_data[key], (k,1)) , local_parameter[key][1])
            X_T_W_X_local = np.matmul(np.reshape(X_T_W_local, (k,1)), np.reshape(cov_data[key], (1,k)))
            X_T_W_v_local = np.matmul(X_T_W_local , local_parameter[key][2])

            

            X_T_W_X_local_enc = secure_send(X_T_W_X_local, PRG_keys_Miner, PRG_keys_Miner_sign, [], bfv_encoder, int(key.split('_')[1]), t, False)
            X_T_W_v_local_enc = secure_send(X_T_W_v_local, PRG_keys_Miner, PRG_keys_Miner_sign, [], bfv_encoder, int(key.split('_')[1]), t, False)
            end = time.time()
            new_runtime = end - start
            
            '''
            in practice the below code is executed. However, python libraries does not impelement "+" operation for two encoded messeages.
            if len(X_T_W_X_global_enc):
                for i in range(len(X_T_W_X_global_enc)):
                    X_T_W_X_global_enc[i] += X_T_W_X_local_enc[i]
                for i in range(len(X_T_W_v_global_enc)):
                    X_T_W_v_global_enc[i] += X_T_W_v_local_enc[i]
            else:
                X_T_W_X_global_enc = X_T_W_X_local_enc.copy()
                X_T_W_v_global_enc = X_T_W_v_local_enc.copy()
            '''
                
            if key == 'user_1':
                X_T_W_X_global += secure_aggregate(X_T_W_X_local_enc,  bfv_encoder, X_T_W_X_global.shape[0], X_T_W_X_global.shape[1], 'Miner', [], PRG_keys_residue, t, False)
                X_T_W_v_global += secure_aggregate(X_T_W_v_local_enc,  bfv_encoder, X_T_W_v_global.shape[0], X_T_W_v_global.shape[1], 'Miner', [], PRG_keys_residue, t, False)
            else:
                X_T_W_X_global += secure_aggregate(X_T_W_X_local_enc,  bfv_encoder, X_T_W_X_global.shape[0], X_T_W_X_global.shape[1], 'Miner', [], [], t, False)
                X_T_W_v_global += secure_aggregate(X_T_W_v_local_enc,  bfv_encoder, X_T_W_v_global.shape[0], X_T_W_v_global.shape[1], 'Miner', [], [], t, False)
        X_T_W_X_global_temp = convert_2D(X_T_W_X_global)
        X_T_W_v_global_temp = convert_2D(X_T_W_v_global)
        X_T_W_X_global_temp = convert_to_float(X_T_W_X_global_temp)
        X_T_W_v_global_temp = convert_to_float(X_T_W_v_global_temp)
        X_T_W_X_global = convert_1D(X_T_W_X_global_temp, X_T_W_X_global.shape[0], X_T_W_X_global.shape[1])
        X_T_W_v_global = convert_1D(X_T_W_v_global_temp, X_T_W_v_global.shape[0], X_T_W_v_global.shape[1])


        global_beta = np.matmul(np.linalg.inv(X_T_W_X_global), X_T_W_v_global)

        X_T_W_X_global = np.zeros((k,k))
        X_T_W_v_global = np.zeros((k,1))
        
        for key in cov_data:
            party_id = int(key.split('_')[1])
            start = time.time()
            local_parameter[key][0] = np.matmul(np.reshape(cov_data[key], (1,k)), global_beta)
            local_parameter[key][0] = 1/(1 + np.exp(-local_parameter[key][0]))
            local_parameter[key][1] = local_parameter[key][0] * (1 - local_parameter[key][0])
            end = time.time()
            new_runtime += end - start
            fisherscoring_runtime = max(fisherscoring_runtime, new_runtime)

        line = "fisherscoring iteration " + str(t) + ":" + str(fisherscoring_runtime) + "\n"
        log_file.write(line)
        fisherscoring_runtime = -1
            
    return global_beta
            
                





def parallel_logistic_regression(geno_data, cov_data, pheno_data, bfv_encoder, PRG_keys_Analyzer, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_residue, itter, log_file):
    parallel_regression_runtime = -1
    global_beta = fisherscoring(cov_data, pheno_data, bfv_encoder, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_residue, itter, log_file)
    print('complete fisherscoring')
    k = cov_data['user_1'].shape[0]
    m = geno_data['user_1'].shape[0]
    X_T_W_X_global = np.zeros((k,k))
    X_T_W_v_global = np.zeros((k,1))
    X_T_W_S_global = np.zeros((k,m))
    c_global = np.zeros((m,1))
    d_global = np.zeros((m,1))
    X_T_W_X_global_enc = []
    X_T_W_v_global_enc = []
    X_T_W_S_global_enc = []
    c_global_enc = []
    d_global_enc = []
    print("parallel regression phase 1")
    for key in cov_data:
        party_id = int(key.split('_')[1])
        start = time.time()
        local_parameter[key][2] = np.log(local_parameter[key][0]/(1-local_parameter[key][0]))+ (np.reshape(pheno_data[key], (1,1))-local_parameter[key][0])/local_parameter[key][1]
        X_T_W_local = np.matmul(np.reshape(cov_data[key], (k,1)), local_parameter[key][1])
        X_T_W_X_local = np.matmul(X_T_W_local, np.reshape(cov_data[key], (1,k)))
        X_T_W_v_local = np.matmul(X_T_W_local, local_parameter[key][2])
        X_T_W_S_local = np.matmul(X_T_W_local, np.reshape(geno_data[key], (1,m)))

        X_T_W_X_local_enc = secure_send(X_T_W_X_local, PRG_keys_Miner, PRG_keys_Miner_sign, [], bfv_encoder, int(key.split('_')[1]), itter+1, False)
        X_T_W_v_local_enc = secure_send(X_T_W_v_local, PRG_keys_Miner, PRG_keys_Miner_sign, [], bfv_encoder, int(key.split('_')[1]), itter+1, False)
        X_T_W_S_local_enc = secure_send(X_T_W_S_local, PRG_keys_Miner, PRG_keys_Miner_sign, [], bfv_encoder, int(key.split('_')[1]), itter+1, False)
        end = time.time()
        parallel_regression_runtime = max(parallel_regression_runtime, (end - start))

        if key == 'user_1':
            X_T_W_X_global += secure_aggregate(X_T_W_X_local_enc, bfv_encoder, X_T_W_X_global.shape[0], X_T_W_X_global.shape[1], 'Miner', [], PRG_keys_residue, itter+1, False)
            X_T_W_v_global += secure_aggregate(X_T_W_v_local_enc, bfv_encoder, X_T_W_v_global.shape[0], X_T_W_v_global.shape[1], 'Miner', [], PRG_keys_residue, itter+1, False)
            X_T_W_S_global += secure_aggregate(X_T_W_S_local_enc, bfv_encoder, X_T_W_S_global.shape[0], X_T_W_S_global.shape[1], 'Miner', [], PRG_keys_residue, itter+1, False)
        else:
            X_T_W_X_global += secure_aggregate(X_T_W_X_local_enc, bfv_encoder, X_T_W_X_global.shape[0], X_T_W_X_global.shape[1], 'Miner', [], [], itter+1, False)
            X_T_W_v_global += secure_aggregate(X_T_W_v_local_enc, bfv_encoder, X_T_W_v_global.shape[0], X_T_W_v_global.shape[1], 'Miner', [], [], itter+1, False)
            X_T_W_S_global += secure_aggregate(X_T_W_S_local_enc, bfv_encoder, X_T_W_S_global.shape[0], X_T_W_S_global.shape[1], 'Miner', [], [], itter+1, False)
        
        if party_id %50 == 0:
                print(f"{party_id} individual complete...")

    X_T_W_X_global_temp = convert_2D(X_T_W_X_global)
    X_T_W_v_global_temp = convert_2D(X_T_W_v_global)
    X_T_W_S_global_temp = convert_2D(X_T_W_S_global)
    X_T_W_X_global_temp = convert_to_float(X_T_W_X_global_temp)
    X_T_W_v_global_temp = convert_to_float(X_T_W_v_global_temp)
    X_T_W_S_global_temp = convert_to_float(X_T_W_S_global_temp)
    X_T_W_X_global = convert_1D(X_T_W_X_global_temp, X_T_W_X_global.shape[0], X_T_W_X_global.shape[1])
    X_T_W_v_global = convert_1D(X_T_W_v_global_temp, X_T_W_v_global.shape[0], X_T_W_v_global.shape[1])
    X_T_W_S_global = convert_1D(X_T_W_S_global_temp, X_T_W_S_global.shape[0], X_T_W_S_global.shape[1])



    line = "parallel regression phase 1 :" + str(parallel_regression_runtime) + "\n"
    log_file.write(line)
    parallel_regression_runtime = -1
    
    temp1 = np.matmul(np.linalg.inv(X_T_W_X_global), X_T_W_S_global)
    temp2 = np.matmul(np.linalg.inv(X_T_W_X_global), X_T_W_v_global)
    print("parallel regression phase 2")
    for key in cov_data:
        party_id = int(key.split('_')[1])
        start = time.time()
        S_star = np.reshape(geno_data[key], (1,m)) - np.matmul(np.reshape(cov_data[key], (1,k)), temp1)
        v_star = local_parameter[key][2] - np.matmul(np.reshape(cov_data[key], (1,k)), temp2)
        
        c_local = np.matmul(S_star.T , local_parameter[key][1])
        c_local = np.matmul(c_local , v_star)
        d_local = np.matmul(S_star.T , local_parameter[key][1])
        #d_local = np.matmul(d_local, S_star)
        for i in range(m):
            d_local[i,0] = d_local[i,0] * S_star[0,i]

        c_local_enc = secure_send(c_local, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_Analyzer, bfv_encoder, int(key.split('_')[1]), itter+2, True)
        d_local_enc = secure_send(d_local, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_Analyzer, bfv_encoder, int(key.split('_')[1]), itter+2, True)
        end = time.time()
        parallel_regression_runtime = max(parallel_regression_runtime, (end - start))
        if key == 'user_1':

            c_global += secure_aggregate(c_local_enc, bfv_encoder, c_global.shape[0], c_global.shape[1], 'Miner', [], PRG_keys_residue, itter+2, False)
            d_global += secure_aggregate(d_local_enc, bfv_encoder, d_global.shape[0], d_global.shape[1], 'Miner', [], PRG_keys_residue, itter+2, False)
        else:
            c_global += secure_aggregate(c_local_enc, bfv_encoder, c_global.shape[0], c_global.shape[1], 'Miner', [], [], itter+2, False)
            d_global += secure_aggregate(d_local_enc, bfv_encoder, d_global.shape[0], d_global.shape[1], 'Miner', [], [], itter+2, False)
        if party_id %50 == 0:
                print(f"{party_id} individual complete...")


    line = "parallel regression phase 2 :" + str(parallel_regression_runtime) + "\n"
    log_file.write(line)
    c_global_enc = encode_model_param(c_global, bfv_encoder)
    d_global_enc = encode_model_param(d_global, bfv_encoder)

    analyzer_runtime = 0

    c_global, runtime = secure_aggregate(c_global_enc, bfv_encoder, c_global.shape[0], c_global.shape[1], 'Analyzer', PRG_keys_Analyzer, [], itter+2, True)
    analyzer_runtime += runtime
    d_global, runtime = secure_aggregate(d_global_enc, bfv_encoder, d_global.shape[0], d_global.shape[1], 'Analyzer', PRG_keys_Analyzer, [], itter+2, True)
    analyzer_runtime += runtime
  
    start = time.time()
    z_score = np.zeros(m)
    p_val = np.zeros(m)
    for i in range(m):
        z_score[i] = c_global[i,0]/math.sqrt(abs(d_global[i,0]))
        z_score[i] = - abs(z_score[i])

        p_val[i] = 2 * stats.norm.cdf(z_score[i])
    '''
    # encrypt p-value vector for send to smart contract
    lst_enc = []
    slot = 16384
    batch_num = math.ceil(p_val.shape[0]/slot)

    for i in range(batch_num-1):
        plain = p_val[slot*i : slot*i + slot]
        cipher = ckks_public_key.encrypt(plain)
        lst_enc.append(cipher)
    plain = p_val[slot*(batch_num-1) :]
    cipher = ckks_public_key.encrypt(plain)
    lst_enc.append(cipher)
    
    end = time.time()
    analyzer_runtime += end - start
    line = "derandomize, compute pval by Analyzer and then encrypt :" + str(analyzer_runtime) + "\n"
    log_file.write(line)
    '''

        
    return z_score, p_val
        
        
        


def FL_parallel_gwas(geno_data, pheno_data, cov_data, PRG_keys_Analyzer, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_residue, log_file):
    
    p = np.full((1,1), 0.5)
    w = 0.25*np.identity(1)
    v = np.zeros((1,1))
    for key in geno_data:
        local_parameter[key] = [np.copy(p), np.copy(w), np.copy(v)]
    itter = 4
    z_score, pval = parallel_logistic_regression(geno_data, cov_data, pheno_data, bfv_encoder, PRG_keys_Analyzer, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_residue, itter, log_file)
    return z_score, pval
    
 
    
    
    
local_parameter = {}


bfv_encoder = Pyfhel()
bfv_encoder.load_context("../key/bfv_encoder")
ckks_public_key = Pyfhel()
ckks_private_key = Pyfhel()
ckks_public_key.load_context("../key/ckks_context")
ckks_private_key.load_context("../key/ckks_context")
ckks_public_key.load_public_key("../key/ckks_pub.key")
ckks_private_key.load_public_key("../key/ckks_pub.key")
ckks_private_key.load_secret_key("../key/ckks_sec.key")


