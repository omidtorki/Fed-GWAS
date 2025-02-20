import numpy as np
import math 
import scipy.stats as stats



num_snp = 21
num_party = 569
num_cov = 9




def fisherscoring(x, y, itter):
    n = x.shape[0]
    k = x.shape[1]
    beta = np.zeros(k)
    p = np.full(n, 0.5)
    w = 0.25*np.identity(n)
    v = np.zeros(n)
    for t in range(itter):
        v = np.log(p/(1-p))+ (y-p)/np.diag(w)
        x_trans = np.transpose(x)
        beta = np.matmul(x_trans, w)
        temp = beta
        beta = np.matmul(beta, x)
        beta = np.linalg.inv(beta)
        beta = np.matmul(beta, temp)
        beta = np.matmul(beta, v)
        p = np.matmul(x,beta)
        p = 1/(1 + np.exp(-p))
        temp2 = np.multiply(p, 1-p)
        w= np.diag(temp2)
        
        
        
    return beta, p, w




def parallel_logistic_regression(s, x, y, itter):
    n= x.shape[0]
    k = x.shape[1]
    m = s.shape[1]
    beta , p , w = fisherscoring(x, y, itter)
    v = np.log(p/(1-p))+ (y-p)/np.diag(w)
    x_trans = np.transpose(x)
    temp1 = np.matmul(x_trans, w)
    temp2 = np.matmul(temp1, x)
    temp2 = np.linalg.inv(temp2)
    temp1 = np.matmul(temp2, temp1)
    temp1 = np.matmul(x, temp1)
    s_star = np.matmul(temp1, s)
    s_star = s - s_star
    v_star = np.matmul(temp1, v)
    v_star = v - v_star
    temp = np.matmul(np.transpose(s_star), w)
    c = np.matmul(temp, v_star)
    #**************************
    #d = np.diag(np.matmul(temp, s_star))
    del beta,p,w,temp1,temp2,v_star
    d = np.zeros(num_snp)
    for i in range(num_snp):
        d[i] = np.dot(temp[i,:], s_star[:,i])
    z_score = np.zeros(m)
    p_val = np.zeros(m)
    for i in range(m):
        z_score[i] = c[i]/math.sqrt(d[i])
        z_score[i] = - abs(z_score[i])
        #p_val[i] = 1- stats.norm.cdf(abs(z_score[i])) + stats.norm.cdf(-abs(z_score[i]))
        p_val[i] = 2 * stats.norm.cdf(z_score[i])

    return z_score, p_val




pheno = np.zeros(num_party, dtype=np.float64)
cov = np.zeros((num_party, num_cov), dtype=np.float64)
geno = np.zeros((num_party, num_snp), dtype=np.float64)

f = open("../data/pheno.txt", "r")
index = 0
for line in f:
    lst = line
    lst2 = np.array(lst, dtype=float)
    pheno[index] = lst2
    index += 1
    
f.close()


f = open("../data/cov.txt", "r")
index = 0
for line in f:
    lst = line.split()[:num_cov]
    lst2 = np.array(lst, dtype=float)
    cov[index,:] = lst2
    index += 1
    
f.close()

f = open("../data/geno.txt", "r")
index = 0
for line in f:
    lst = line.split()[:num_snp]
    lst2 = np.array(lst, dtype=float)
    geno[index,:] = lst2
    index += 1
    
f.close()
print("start...")

intercept = np.ones(num_party)

cov = np.insert(cov, 0, intercept, axis=1)

z_score , pval = parallel_logistic_regression(geno, cov, pheno, 4)

f = open("../out/centralized_pval.txt", "w")
for i in range(num_snp):
    line = str(z_score[i]) + " " + str(pval[i]) + "\n"
    f.write(line)
    
f.close()
print("complete...")



