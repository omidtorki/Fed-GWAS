import numpy as np
import random
import math
import time
from Pyfhel import Pyfhel



def secure_send(model_param, PRG_keys_Miner, PRG_keys_Miner_sign, PRG_keys_Analyzer, bfv_encoder, party_id,  _round, final_flag):
    q = 144115188075855873
    if model_param.ndim == 2:
        model_param = convert_2D(model_param)
    
    model_param = convert_to_int(model_param)
    perturbe_vector = pseudo_random_generator(PRG_keys_Miner, PRG_keys_Miner_sign, party_id, _round, len(model_param))
    model_param += perturbe_vector
    if final_flag == True:
        key = PRG_keys_Analyzer[party_id-1]
        random.seed(key)
        perturbe_vector = [random.randint(0, q-1) for _ in range(model_param.shape[0])]
        perturbe_vector = np.array(perturbe_vector)
        model_param += perturbe_vector
    
    for i in range(model_param.shape[0]):
        model_param[i] = model_param[i] % q 
    lst_enc = []
    slot = 32768
    batch_num = math.ceil(model_param.shape[0]/slot)

    for i in range(batch_num-1):
        plain = model_param[slot*i : slot*i + slot]
        cipher = bfv_encoder.encode(plain)
        lst_enc.append(cipher)
    plain = model_param[slot*(batch_num-1) :]
    cipher = bfv_encoder.encode(plain)
    lst_enc.append(cipher)

    return lst_enc




def secure_aggregate(enc_lst, bfv_encoder, shape1, shape2, invoker, PRG_keys_Analyzer, PRG_keys_residue, _round, to_float):
    q = 144115188075855873
    plain = np.zeros(0)
    if invoker == 'Miner':
        for i in range(len(enc_lst)):
            plain_tmp = bfv_encoder.decode(enc_lst[i])
            plain = np.concatenate([plain, plain_tmp])
        
        #print(plain.shape)
        if len(PRG_keys_residue):
            for key in PRG_keys_residue:
                key = key + _round
                random.seed(key)
                random_vector = [random.randint(0, q-1) for _ in range(plain.shape[0])]   # 65537 == 0
                random_vector = np.array(random_vector)
                for r in range(random_vector.shape[0]):
                    random_vector[r] = (q - random_vector[r]) % q
                plain += random_vector
            for i in range(plain.shape[0]):
                plain[i] %= q
        if to_float == True:
            plain = convert_to_float(plain)
            aggregate_res = convert_1D(plain, shape1, shape2)
            return aggregate_res
        aggregate_res = convert_1D(plain, shape1, shape2)
        return aggregate_res
            
    elif invoker == 'Analyzer':
        Analyzer_aggregation_runtime = 0
        start = time.time()
        for i in range(len(enc_lst)):
            plain_tmp = bfv_encoder.decode(enc_lst[i])
            plain = np.concatenate([plain, plain_tmp])

        end = time.time()
        Analyzer_aggregation_runtime += end - start
        # generate of random (below task is a preprocessing task)
        vector_size = plain.shape[0]
        perturbe_vector = np.zeros(vector_size)
        for key in PRG_keys_Analyzer:
            random.seed(key)
            random_vector = [random.randint(0, q-1) for _ in range(vector_size)]   # 65537 == 0
            random_vector = np.array(random_vector)
            perturbe_vector += random_vector
        for i in range(perturbe_vector.shape[0]):
            perturbe_vector[i] = (q - perturbe_vector[i]) % q
        start = time.time()
        plain += perturbe_vector
        for i in range(plain.shape[0]):
            plain[i] %= q
        plain = convert_to_float(plain)
        aggregate_res = convert_1D(plain, shape1, shape2)
        end = time.time()
        Analyzer_aggregation_runtime += end - start 
        return aggregate_res, Analyzer_aggregation_runtime
        


def pseudo_random_generator(PRG_keys_Miner, PRG_keys_Miner_sign, party_id, _round, vector_size):
    perturbe_vector = np.zeros(vector_size, dtype = int)
    q = 144115188075855873
    for i in range(PRG_keys_Miner.shape[0]):
        key = PRG_keys_Miner[i,party_id-1]
        key = key + _round
        random.seed(key)
        random_vector = [random.randint(0, q-1) for _ in range(vector_size)]   # 65537 == 0
        random_vector = np.array(random_vector)
        if PRG_keys_Miner_sign[i,party_id-1] == 1:
            for r in range(random_vector.shape[0]):
                random_vector[r] = (q - random_vector[r]) % q
        perturbe_vector += random_vector

        
    return perturbe_vector
    


def encode_model_param(model_param, bfv_encoder):
    q = 144115188075855873
    for i in range(model_param.shape[0]):
        model_param[i] %= q   
    if model_param.ndim == 2:
        model_param = convert_2D(model_param)
    lst_enc = []
    slot = 32768
    batch_num = math.ceil(len(model_param)/slot)

    for i in range(batch_num-1):
        plain = model_param[slot*i : slot*i + slot]
        cipher = bfv_encoder.encode(plain)
        lst_enc.append(cipher)
    plain = model_param[slot*(batch_num-1) :]
    cipher = bfv_encoder.encode(plain)
    lst_enc.append(cipher)

    return lst_enc




def convert_2D(arr):
    lst_1D = np.zeros(arr.shape[0]*arr.shape[1])
    index = 0
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            lst_1D[index] = (arr[i,j])
            index += 1
            
    return lst_1D
    
def convert_1D(lst, shape1, shape2):
    arr = np.zeros((shape1, shape2))
    index = 0
    for i in range(shape1):
        for j in range(shape2):
            arr[i,j] = lst[index]
            index += 1
            
    return arr
    

def convert_to_int(float_numbers):
    f = 37
    _2ef = 2**f  # 2^f : for faster computation
    q = 144115188075855873
    int_numbers = np.zeros(float_numbers.shape,dtype=int)
    for i in range(int_numbers.shape[0]):
        int_numbers[i] = math.floor(float_numbers[i]*(_2ef))% q
    return int_numbers

def convert_to_float(int_numbers):
    k = 57
    f = 37
    q = 144115188075855873
    for i in range(int_numbers.shape[0]):
      int_numbers[i] = int_numbers[i] % q
    _2e_f = 2**(-f) # 2^(-f)
    _2ek_1 = 2**(k-1) # 2^(k-1)
    float_numbers = np.zeros(int_numbers.shape)
    for i in range(int_numbers.shape[0]):
        if int_numbers[i] >= 0 and int_numbers[i] < _2ek_1:
            float_numbers[i] = int_numbers[i]*(_2e_f)
        elif int_numbers[i] >= q - _2ek_1 and int_numbers[i] <= q-1:
            float_numbers[i] = -(q-int_numbers[i])*(_2e_f)
    return float_numbers



