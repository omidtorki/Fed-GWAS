import os
import sys
import random
import math
import numpy as np
from tinyec import registry
import secrets
from Pyfhel import Pyfhel
import time
curve = registry.get_curve('brainpoolP256r1')
# ECC key runtime ---> generate PRG key, encrypt it for send to blockchain and decrypt


def compress_point(point):
    return hex(point.x) + hex(point.y % 2)[2:]

def ecc_calc_encryption_keys(pubKey, random_number):
    ciphertextPubKey = pubKey + (curve.g * random_number)
    shared_PRG_Key = pubKey + (curve.g * random_number)
    return (shared_PRG_Key, ciphertextPubKey)

def ecc_calc_decryption_key(privKey, ciphertextPubKey):
    shared_PRG_Key = ciphertextPubKey + (privKey * curve.g)
    return shared_PRG_Key



def key_generation_runtime(log_file):
    param_file = open("../param/param.txt", "r")
    num_ind = int(param_file.readline().split()[1])
    param_file.readline()
    param_file.readline()
    num_Miner = int(param_file.readline().split()[1])
    runtime_iteration = int(param_file.readline().split()[1])
    param_file.close()
    # Generate ECC key pair
    start = time.time()
    for _ in range(runtime_iteration):
        privKey = secrets.randbelow(curve.field.n)
        pubKey = privKey * curve.g
    end = time.time()
    ecc_key_generation_runtime = (end - start)/runtime_iteration
    line = "ECC public/private key generation :" + str(ecc_key_generation_runtime) + "\n"
    log_file.write(line)
    
    #prg_encryption_runtime = 0
    prg_decryption_runtime = 0
    for _ in range(runtime_iteration):
       #start = time.time()
        random_number = int.from_bytes(secrets.token_bytes(16), 'big')
        (encryptKey, ciphertextPubKey) = ecc_calc_encryption_keys(pubKey, random_number)
        #end = time.time()
        #prg_encryption_runtime += end - start
        start = time.time()
        decryptKey = ecc_calc_decryption_key(privKey, ciphertextPubKey)
        end = time.time()
        prg_decryption_runtime += end - start
        
    #prg_encryption_runtime /= runtime_iteration
    #prg_encryption_runtime *= (num_ind-1)
    prg_decryption_runtime /= runtime_iteration
    prg_decryption_runtime *= (num_Miner)
    #line = "PRG key encryption :" + str(prg_encryption_runtime) + "\n"
    #log_file.write(line)
    line = "PRG key decryption :" + str(prg_decryption_runtime) + "\n"
    log_file.write(line)
    
    






 
def generate_PRG_keys(num_ind, num_Miner):
    
    keys = []
    # generate num_Miner * ceil(num_ind/2) PRG-keys for n individuals.
    for _ in range(math.ceil(num_ind/2)):
        for _ in range(num_Miner):
            random_bytes = os.urandom(16)  # 16 bytes = 128 bits
            keys.append(random_bytes)
    with open("../key/PRG_keys_Miner.bin", 'wb') as file:
        for key in keys:
            file.write(key)
    keys = []
    for _ in range(num_ind):
        random_bytes = os.urandom(16)
        keys.append(random_bytes)
    with open("../key/PRG_keys_Analyzer.bin", 'wb') as file:
        for key in keys:
            file.write(key)
        
 

 
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
 


def generate_ckks_key():
    HE = Pyfhel()           
    ckks_params = {
    'scheme': 'CKKS', 'n': 2**15, 'scale': 2**40, 'qi_sizes': [60, 60]  
}
    HE.contextGen(**ckks_params)
    HE.keyGen()

    HE.save_context("../key/ckks_context")
    HE.save_public_key("../key/ckks_pub.key")
    HE.save_secret_key("../key/ckks_sec.key")


def generate_bfv_encoder():
    HE = Pyfhel() 
    bfv_params = {'scheme': 'BFV', 'n': 2**15, 't_bits': 60}
    HE.contextGen(**bfv_params)
    HE.save_context("../key/bfv_encoder")



if __name__ == "__main__":
    
    param_file = open("../param/param.txt", "r")
    num_ind = int(param_file.readline().split()[1])
    param_file.close()
    
    # Parse command line arguments
    function_name = sys.argv[1]

 
    # Call the appropriate function based on the provided name
    if function_name == 'generate_PRG_keys':
        result = generate_PRG_keys(num_ind)
    elif function_name == 'restore_PRG_keys':
        result = restore_PRG_keys(num_ind)
    elif function_name == 'runtime':
        result = key_generation_runtime(num_ind)
    elif function_name == 'generate_ckks_key':
        generate_ckks_key()
    elif function_name == 'generate_bfv_encoder':
        generate_bfv_encoder()
    else:
        print(f"Error: Function '{function_name}' not found.")
        
 


