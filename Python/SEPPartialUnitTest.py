

def centralgravity(mu, position):
    return -mu * position / (np.linalg.norm(position)**Decimal(3.0))

def SEPcorrection(mu, r_1, r_2, dr):
    r_1_c = r_1 + dr

    r_12 = r_2 - r_1
    r_12_c = r_2 - r_1_c

    a = centralgravity(mu, r_12)
    a_c = centralgravity(mu, r_12_c)

    # print(a, a_c, a_c-a)

    return a_c - a

def AnalyticalPartial(mu, r_1, r_2, dr):
    r_1_c = r_1 + dr

    r_12 = r_2 - r_1
    r_12_c = r_2 - r_1_c

    d_12 = np.linalg.norm(r_12)
    d_12_c = np.linalg.norm(r_12_c)

    p1 = Decimal(3.0) * np.outer(r_12, r_12_c.T) / d_12_c**5
    p2 = Decimal(-3.0) * np.outer(dr, r_12_c.T) / d_12_c**5
    p3 = Decimal(-3.0) * np.outer(r_12, r_12.T) / d_12**5
    p4 = np.dot(( Decimal(1.0) / d_12**3 - Decimal(1.0) / d_12_c**3 ), np.identity(3, dtype=Decimal))

    print("partial term 1: \n", mu*p1)
    print("partial term 2: \n", mu*p2)
    print("partial term 3: \n", mu*p3)
    print("partial term 4: \n", mu*p4)

    return -mu*(p1+p2+p3+p4)

import numpy as np
from decimal import *

getcontext().prec = 40

mu_test = 1.32712440041939e+20
mu_S = Decimal(1.32712440041939e+20)


r_S = np.asarray(([Decimal(-1058202435.85883), Decimal(-407616171.803058) , Decimal(-143292503.024126) ]))
r_M = np.asarray([Decimal(17776989161.8444), Decimal(-56861189168.2378), Decimal(-32252099174.0247) ])

dr_SEP = np.asarray([Decimal(-0.00068076), Decimal(-0.000483976), Decimal(-0.00019087) ])

da = SEPcorrection(mu_S, r_S, r_M, dr_SEP)

p = Decimal(1000.0)
p_x = np.asarray([p, Decimal(0.0), Decimal(0.0)])
p_y = np.asarray([Decimal(0.0), p, Decimal(0.0)])
p_z = np.asarray([Decimal(0.0), Decimal(0.0), p])

cd_x = ( SEPcorrection(mu_S, r_S+p_x, r_M, dr_SEP) - SEPcorrection(mu_S, r_S-p_x, r_M, dr_SEP) ) / (Decimal(2.0)*p)
cd_y = ( SEPcorrection(mu_S, r_S+p_y, r_M, dr_SEP) - SEPcorrection(mu_S, r_S-p_y, r_M, dr_SEP) ) / (Decimal(2.0)*p)
cd_z = ( SEPcorrection(mu_S, r_S+p_z, r_M, dr_SEP) - SEPcorrection(mu_S, r_S-p_z, r_M, dr_SEP) ) / (Decimal(2.0)*p)

cd = np.vstack([cd_x, cd_y, cd_z])

partial = AnalyticalPartial(mu_S, r_S, r_M, dr_SEP)

print("central difference: \n", cd)
print("partial: \n", partial)
print("partial - central difference: \n", partial-cd)


# code_da = np.asarray([Decimal(-4.01575647734909e-12),  Decimal(1.20189969088358e-12),  Decimal(1.04582835447342e-12)])
# print(code_da - da)






