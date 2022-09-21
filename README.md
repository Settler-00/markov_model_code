# markov_model_code
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 09:34:58 2022

@author: zWX1091054
"""

from scipy.optimize import root
import numpy as np
import math

def SUM1(a, b, c):
    s = 0
    while a <= b:
        s += c ** (a - 1)
        a += 1
    return s

def SUM2(a, b, c):
    s = 0
    if a > b or a < 0 or b < 0:
        return s
    else:
        while a <= b:
            s += c ** (b - a)
            a += 1
        return s

def f1(x, CW_min, CW_max, m, n):
    m_max = 60    
    
    f0 = (1 - x[0]) ** n
    for i in range(1, m + 1):
        f0 = f0 * (1 - x[i])
    
    tau_temp = np.zeros(shape = (1, m + 1))
    
    for k in range(0, m + 1):
        FAI = 0
        g = f0 / (1 - x[k])
        p = 1 - f0 / (1 - x[k])
        a = 0 
        c = 0
        while a <= m_max:
            if a <= CW_max[k] - CW_min[k]:
                c = a + CW_min[k]
            else:
                c = CW_max[k]    
            b = 1
            if b > 2 ** c - 1:
                FAI += 0
            else:
                while b <= 2 ** c - 1:
                    FAI += (2 ** c - b) * p ** a / 2 ** c
                    b += 1
            a += 1
        b_steady = ((1 - p ** (m_max + 1)) / (1 - p) + FAI / g) ** -1
        tau_temp[0, k] = (1 - p ** (m_max + 1)) / (1 - p) * b_steady
    
    eqs = []
    for j in range(0, m + 1):
        eqs.append(x[j] - tau_temp[0, j])
    return eqs

def ste(P, max_iter = 1000, tol = 1e-8):
    n_components = len(P)
    pi0 = np.array([1 / n_components] * n_components)
    for j in range(max_iter):
        pi1 = np.dot(P, pi0)
        if np.sum(np.abs(pi0 - pi1)) < tol:
            break
        pi0 = pi1
    return pi0
    
AIFSN = [3, 3]
slot = 9
SIFS = 16
DIFS = 32
ACK = 57 * 8 / 6
RTS = 64
CTS = 44
EIFS = ACK + DIFS + SIFS
AIFS = 3 * slot + SIFS

AGGR_TIME = 1000
PREAMBLE_TIME = 5
RATE_ont = 2144
AGGR_ont =  int(RATE_ont * (AGGR_TIME - PREAMBLE_TIME) / (1500 * 8))
DATA_ont = 1500 * 8 * AGGR_ont
Ts_ont = DATA_ont / RATE_ont + AIFS + ACK

RATE_sta = 2144
DATA_sta = 40 * 8
Ts_sta = DATA_sta / RATE_sta + AIFS

Ts_ap = 0.5 * (Ts_ont + Ts_sta)

Ts_in = Ts_ont
Ts1 = Ts_in
Tc = AGGR_TIME + EIFS + AIFS

CW_in = 4
NUM_in = 1
W = 3
list2 = []
list_beg = []
list_end = []
for i in range(0, W + 1):
    beg = int(i * (W + 2 - i) + i * (i - 1) / 2)
    end = int((i + 1) * (W + 1 - i) + i * (i + 1) / 2 - 1)
    list_beg.append(beg)
    list_end.append(end)

for CW_ont in range(5, 6):
    for CW_ap in range(1, 11):
        for CW_sta in range(5, 6):
            list0 = [] 
            list1 = []
            list_tau_in = [0, 0, 0, 0, 0, 0, 0]
            list_f = [0, 0, 0, 0, 0, 0, 0]
            list_P1 = [0, 0, 0, 0, 0, 0, 0]
            list_Ts2 = [0, 0, 0, 0, 0, 0, 0]
            list_DATA_down = [0, 0, 0, 0, 0, 0, 0]
            flag = [0, 0, 0, 0, 0, 0, 0]
            for j in range(0, W + 1):
                for i in range(0, W - j + 1):
                    if i ==0 and j == 0:
                        z = 0
                        if flag[z] == 0:
                            CW_min = [CW_in, CW_ont]
                            CW_max = [CW_in, CW_ont]
                            result = root(f1, [0.5, 0.5], tol = 1e-10, args = (CW_min, CW_max, 1, NUM_in))
                            tau_in = result.x[0]
                            tau_ont = result.x[1]
                            f = (1 - tau_ont) * (1 - tau_in) ** NUM_in
                            P1 = tau_ont * f / (1 - tau_ont)
                            Ts2 = Ts_ont
                            DATA_down = 0
                            list_tau_in[z] = tau_in
                            list_f[z] = f
                            list_P1[z] = P1
                            list_Ts2[z] = Ts2
                            list_DATA_down[z] = DATA_down
                            flag[z] = 1
                        else:
                            tau_in = list_tau_in[z] 
                            f = list_f[z]
                            P1 = list_P1[z]
                            Ts2 = list_Ts2[z]
                            DATA_down = list_DATA_down[z]
                    elif i == W and j == 0:
                        z = 1
                        if flag[z] == 0:
                            CW_min = [CW_in, CW_ap]
                            CW_max = [CW_in, CW_ap]
                            result = root(f1, [0.5, 0.5], tol = 1e-10, args = (CW_min, CW_max, 1, NUM_in))
                            tau_in = result.x[0]
                            tau_ap = result.x[1]
                            f = (1 - tau_ap) * (1 - tau_in) ** NUM_in
                            P1 = tau_ap * f / (1 - tau_ap)
                            Ts2 = Ts_ap
                            DATA_down = 0.5 * DATA_ont
                            list_tau_in[z] = tau_in
                            list_f[z] = f
                            list_P1[z] = P1
                            list_Ts2[z] = Ts2
                            list_DATA_down[z] = DATA_down
                            flag[z] = 1
                        else:
                            tau_in = list_tau_in[z] 
                            f = list_f[z]
                            P1 = list_P1[z]
                            Ts2 = list_Ts2[z]
                            DATA_down = list_DATA_down[z]
                    elif i == 0 and j == W:
                        z = 2
                        if flag[z] == 0:
                            CW_min = [CW_in, CW_sta]
                            CW_max = [CW_in, CW_sta]
                            result = root(f1, [0.5, 0.5], tol = 1e-10, args = (CW_min, CW_max, 1, NUM_in))
                            tau_in = result.x[0]
                            tau_sta = result.x[1]
                            f = (1 - tau_sta) * (1 - tau_in) ** NUM_in
                            P1 = tau_sta * f / (1 - tau_sta)
                            Ts2 = Ts_sta
                            DATA_down = 0
                            list_tau_in[z] = tau_in
                            list_f[z] = f
                            list_P1[z] = P1
                            list_Ts2[z] = Ts2
                            list_DATA_down[z] = DATA_down
                            flag[z] = 1
                        else:
                            tau_in = list_tau_in[z] 
                            f = list_f[z]
                            P1 = list_P1[z]
                            Ts2 = list_Ts2[z]
                            DATA_down = list_DATA_down[z]
                    elif i > 0 and i < W and j == 0:
                        z = 3
                        if flag[z] == 0:
                            CW_min = [CW_in, CW_ont, CW_ap]
                            CW_max = [CW_in, CW_ont, CW_ap]
                            result = root(f1, [0.5, 0.5, 0.5], tol = 1e-10, args = (CW_min, CW_max, 2, NUM_in))
                            tau_in = result.x[0]
                            tau_ont = result.x[1]
                            tau_ap = result.x[2]
                            f = (1 - tau_ont) * (1 - tau_ap) * (1 - tau_in) ** NUM_in
                            P1 = tau_ont * f / (1 - tau_ont) + tau_ap * f / (1 - tau_ap)
                            PO1 = tau_ont * (1 - tau_ap) / (P1 / (1 - tau_in) ** NUM_in)
                            PA1 = 0.5 * tau_ap * (1 - tau_ont) / (P1 / (1 - tau_in) ** NUM_in)
                            Ts2 = PO1 * Ts_ont + 2 * PA1 * Ts_ap
                            DATA_down = PA1 * DATA_ont
                            list_tau_in[z] = tau_in
                            list_f[z] = f
                            list_P1[z] = P1
                            list_Ts2[z] = Ts2
                            list_DATA_down[z] = DATA_down
                            flag[z] = 1
                        else:
                            tau_in = list_tau_in[z] 
                            f = list_f[z]
                            P1 = list_P1[z]
                            Ts2 = list_Ts2[z]
                            DATA_down = list_DATA_down[z]
                    elif i == 0 and j > 0 and j < W:
                        z = 4
                        if flag[z] == 0:
                            CW_min = [CW_in, CW_ont, CW_sta]
                            CW_max = [CW_in, CW_ont, CW_sta]
                            result = root(f1, [0.5, 0.5, 0.5], tol = 1e-10, args = (CW_min, CW_max, 2, NUM_in))
                            tau_in = result.x[0]
                            tau_ont = result.x[1]
                            tau_sta = result.x[2]
                            f = (1 - tau_ont) * (1 - tau_sta) * (1 - tau_in) ** NUM_in
                            P1 = tau_ont * f / (1 - tau_ont) + tau_sta * f / (1 - tau_sta)
                            PO2 = tau_ont * (1 - tau_sta) / (P1 / (1 - tau_in) ** NUM_in)
                            PS2 = tau_sta * (1 - tau_ont) / (P1 / (1 - tau_in) ** NUM_in)
                            Ts2 = PO2 * Ts_ont + PS2 * Ts_sta
                            DATA_down = 0
                            list_tau_in[z] = tau_in
                            list_f[z] = f
                            list_P1[z] = P1
                            list_Ts2[z] = Ts2
                            list_DATA_down[z] = DATA_down
                            flag[z] = 1
                        else:
                            tau_in = list_tau_in[z] 
                            f = list_f[z]
                            P1 = list_P1[z]
                            Ts2 = list_Ts2[z]
                            DATA_down = list_DATA_down[z]
                    elif i != 0 and i != W and i + j == W:
                        z = 5
                        if flag[z] == 0:
                            CW_min = [CW_in, CW_ap, CW_sta]
                            CW_max = [CW_in, CW_ap, CW_sta]
                            result = root(f1, [0.5, 0.5, 0.5], tol = 1e-10, args = (CW_min, CW_max, 2, NUM_in))
                            tau_in = result.x[0]
                            tau_ap = result.x[1]
                            tau_sta = result.x[2]
                            f = (1 - tau_ap) * (1 - tau_sta) * (1 - tau_in) ** NUM_in
                            P1 = tau_ap * f / (1 - tau_ap) + tau_sta * f / (1 - tau_sta)
                            PA2 = 0.5 * tau_ap * (1 - tau_sta) / (P1 / (1 - tau_in) ** NUM_in)
                            PS1 = tau_sta * (1 - tau_ap) / (P1 / (1 - tau_in) ** NUM_in)
                            Ts2 = 2 * PA2 * Ts_ap + PS1 * Ts_sta
                            DATA_down = PA2 * DATA_ont
                            list_tau_in[z] = tau_in
                            list_f[z] = f
                            list_P1[z] = P1
                            list_Ts2[z] = Ts2
                            list_DATA_down[z] = DATA_down
                            flag[z] = 1
                        else:
                            tau_in = list_tau_in[z] 
                            f = list_f[z]
                            P1 = list_P1[z]
                            Ts2 = list_Ts2[z]
                            DATA_down = list_DATA_down[z]
                    else:
                        z = 6
                        if flag[z] == 0:
                            CW_min = [CW_in, CW_ont, CW_ap, CW_sta]
                            CW_max = [CW_in, CW_ont, CW_ap, CW_sta]
                            result = root(f1, [0.5, 0.5, 0.5, 0.5], tol = 1e-10, args = (CW_min, CW_max, 3, NUM_in))
                            tau_in = result.x[0]
                            tau_ont = result.x[1]
                            tau_ap = result.x[2]
                            tau_sta = result.x[3]
                            f = (1 - tau_ont) * (1 - tau_ap) * (1 - tau_sta) * (1 - tau_in) ** NUM_in
                            P1 = tau_ont * f / (1 - tau_ont) + tau_ap * f / (1 - tau_ap) + tau_sta * f / (1 - tau_sta)
                            PO3 = tau_ont * (1 - tau_ap) * (1 - tau_sta) / (P1 / (1 - tau_in) ** NUM_in)
                            PA3 = 0.5 * tau_ap * (1 - tau_ont) * (1 - tau_sta) / (P1 / (1 - tau_in) ** NUM_in)
                            PS3 = tau_sta * (1 - tau_ont) * (1 - tau_ap) / (P1 / (1 - tau_in) ** NUM_in)
                            Ts2 = PO3 * Ts_ont + 2 * PA3 * Ts_ap + PS3 * Ts_sta
                            DATA_down = PA3 * DATA_ont
                            list_tau_in[z] = tau_in
                            list_f[z] = f
                            list_P1[z] = P1
                            list_Ts2[z] = Ts2
                            list_DATA_down[z] = DATA_down
                            flag[z] = 1
                        else:
                            tau_in = list_tau_in[z] 
                            f = list_f[z]
                            P1 = list_P1[z]
                            Ts2 = list_Ts2[z]
                            DATA_down = list_DATA_down[z]
                        
                    P2 = NUM_in * tau_in * f / (1 - tau_in)
                    P3 = 1 - P2 - P1 - f
                    mu = 0
                    for n in range(0, 100):
                        for m in range(0, n + 1):
                            for r in range(0, n - m + 1):
                                C1 = math.factorial(n) / (math.factorial(m) * math.factorial(n - m))
                                C2 = math.factorial(n - m) / (math.factorial(r) * math.factorial(n - m - r))
                                mu += (C1 * (f ** m) * C2 * (P3 ** r) * (P1 ** (n - m - r)) * P2) * (m * slot + r * Tc + (n - m - r) * Ts1 + Ts2)
                    TP_down = DATA_down / mu
                    list0.append(TP_down)
                    list1.append(mu)
            
            # 计算平稳分布
            s = int((W + 1) * (W / 2 + 1))
            P = np.zeros(shape = (s, s))
            P[0, 1] = 1
            #r1 = row_mat r2 = row_st
            for r1 in range(1, W):
                P[r1, r1 - 1] = PA1
                P[r1, r1 + 1] = PO1
                P[r1, r1 + W] = PA1
            P[W, W - 1] = 0.5
            P[W, 2 * W] = 0.5
            for r2 in range(1, W + 1):
                for r1 in range(list_beg[r2], list_end[r2] + 1):
                    if list_beg[r2] == list_end[r2]:
                        P[list_beg[r2], list_end[r2] - 1] = 1
                    elif r1 == list_beg[r2]:
                        P[r1, list_beg[r2 - 1] + 1] = PS2
                        P[r1, r1 + 1] = PO2
                    elif r1 == list_end[r2]:
                        P[r1, list_end[r2 - 1]] = PS1
                        P[r1, r1 - 1] = PA2
                        P[r1, list_end[r2 + 1]] = PA2
                    else:
                        P[r1, r1 - 1] = PA3
                        P[r1, r1 + 1] = PO3
                        P[r1, r1 + W - r2] = PA3
                        P[r1, r1 - W + r2 -1] = PS3
            pi = ste(P)
            
            TP_final = 0
            sum_mu_pi = 0
            for t in range(0, len(pi)):
                sum_mu_pi += pi[t] * list1[t]   
            for k in range(0, len(pi)):
                TP_final += (list0[k] * list1[k] * pi[k] / sum_mu_pi)
            list2.append([CW_ont, CW_ap, CW_sta, TP_final])

f0 = open('E:/EDCA_Model/TCP/AP.txt', 'w')
for i in list2:
    f0.write(str(i) + '\n')
