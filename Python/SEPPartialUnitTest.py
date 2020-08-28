

def CentralGravity(mu, position):
    return -mu * position / (np.linalg.norm(position)**Decimal(3.0))


def SEPcorrection(mu, r_1, r_2, dr):

    r_1_c = r_1 + dr
    r_12 = r_2 - r_1
    r_12_c = r_2 - r_1_c

    a = CentralGravity(mu, r_12)
    a_c = CentralGravity(mu, r_12_c)

    # print("a_c: ", a_c)
    # print("a: ", a)
    # print("dif: ", a_c-a)

    return a_c - a


def PartialWrtPosition(mu, r_1, r_2, dr):

    r_1_c = r_1 + dr
    r_12 = r_2 - r_1
    r_12_c = r_2 - r_1_c
    d_12 = np.linalg.norm(r_12)
    d_12_c = np.linalg.norm(r_12_c)

    p1 = Decimal(3.0) * np.outer(r_12, r_12_c.T) / d_12_c**5
    p2 = Decimal(-3.0) * np.outer(dr, r_12_c.T) / d_12_c**5
    p3 = Decimal(-3.0) * np.outer(r_12, r_12.T) / d_12**5
    p4 = np.dot(( Decimal(1.0) / d_12**3 - Decimal(1.0) / d_12_c**3 ), np.identity(3, dtype=Decimal))

    return -mu*(p1+p2+p3+p4)


def PartialWrtMu(mu, r_1, r_2, dr, acc):

    r_1_c = r_1 + dr
    r_12 = r_2 - r_1
    r_12_c = r_2 - r_1_c
    d_12 = np.linalg.norm(r_12)
    d_12_c = np.linalg.norm(r_12_c)

    p1  = acc/mu
    p2a = np.identity(3, dtype=Decimal) / (d_12_c**3)
    p2b = Decimal(3.0) * np.outer(r_12_c, r_12_c.T) / (d_12_c**5)
    p2  = np.matmul((p2a - p2b), dr)

    # print("p1 \n", p1)
    # print("p2a \n", p2a)
    # print("p2b \n", p2b)
    # print("p2 \n", p2)

    return p1 - p2

    # c = dr / mu
    #
    # p1 = r_12_c / (d_12_c**3)
    # p2 = Decimal(3.0) * np.dot(r_12_c.T, r_12_c) * c / (mu * d_12_c**5)
    # p3 = c / (mu * d_12_c**3)
    # p4 = r_12 / d_12**3

    # return -p1 + p2 - p3 + p4

def PartialWrtEta(mu, r_1, r_2, dr, n):

    c = dr / n

    r_1_c = r_1 + c*n
    r_12 = r_2 - r_1
    r_12_c = r_2 - r_1_c
    d_12 = np.linalg.norm(r_12)
    d_12_c = np.linalg.norm(r_12_c)

    # p1 = -mu * c / (d_12_c**3)
    # p2 = Decimal(3.0) * mu * np.dot(r_12_c.T, r_12_c) * c / (d_12_c**5)

    brackets1 = np.identity(3, dtype=Decimal) / (d_12_c**3)
    brackets2 = Decimal(3.0) * np.outer(r_12_c, r_12_c.T) / (d_12_c**5)
    partial = mu * np.matmul( (brackets1 - brackets2) , dr) / n

    return partial


def CentralDifferenceWrtPos(p, mu, r_1, r_2, dr):

    p_x = np.asarray([p, Decimal(0.0), Decimal(0.0)])
    p_y = np.asarray([Decimal(0.0), p, Decimal(0.0)])
    p_z = np.asarray([Decimal(0.0), Decimal(0.0), p])

    cd_x = (SEPcorrection(mu, r_1 + p_x, r_2, dr) - SEPcorrection(mu_S, r_1 - p_x, r_2, dr)) / (
                Decimal(2.0) * p)
    cd_y = (SEPcorrection(mu, r_1 + p_y, r_2, dr) - SEPcorrection(mu_S, r_1 - p_y, r_2, dr)) / (
                Decimal(2.0) * p)
    cd_z = (SEPcorrection(mu, r_1 + p_z, r_2, dr) - SEPcorrection(mu_S, r_1 - p_z, r_2, dr)) / (
                Decimal(2.0) * p)

    return np.vstack([cd_x, cd_y, cd_z])


def CentralDifferenceWrtMu(p, mu, r_1, r_2, dr):

    dr_up = dr*mu/(mu+p)
    dr_down = dr*mu/(mu-p)

    # print(dr, dr_up, dr_down)

    cd = (SEPcorrection(mu + p, r_1, r_2, dr_up) - SEPcorrection(mu_S - p, r_1, r_2, dr_down)) / (
                Decimal(2.0) * p)
    return cd

def CentralDifferenceWrtEta(p, mu, r_1, r_2, dr, n):

    dr_up = dr*(n+p)/n
    dr_down = dr*(n-p)/n

    cd = (SEPcorrection(mu, r_1, r_2, dr_up) - SEPcorrection(mu_S, r_1, r_2, dr_down)) / (
                Decimal(2.0) * p)
    return cd



import numpy as np
from decimal import *

getcontext().prec = 33

mu_test = 1.32712440041939e+20
mu_S = Decimal(1.32712440041939e+20)

r_S = np.asarray(([Decimal(-1058202435.85883), Decimal(-407616171.803058) , Decimal(-143292503.024126) ]))
r_M = np.asarray([Decimal(17776989161.8444), Decimal(-56861189168.2378), Decimal(-32252099174.0247) ])

dr_SEP_nvfalse = np.asarray([Decimal(2.92724249666266), Decimal(2.33826700072565), Decimal(0.898587531648646) ])
dr_SEP_nvtrue = np.asarray([Decimal(0.390298999555446), Decimal(0.311768933430425), Decimal(0.119811670886616) ])

dr_SEP = dr_SEP_nvfalse

da = SEPcorrection(mu_S, r_S, r_M, dr_SEP)

p_pos = Decimal(1000.0)
cd_pos = CentralDifferenceWrtPos(p_pos, mu_S, r_S, r_M, dr_SEP)

partial_pos = PartialWrtPosition(mu_S, r_S, r_M, dr_SEP)

print("central difference wrt position: \n", cd_pos)
print("partial wrt position: \n", partial_pos)
print("partial - central difference wrt position: \n", partial_pos-cd_pos)


# p_mu = Decimal(1E19)
# cd_mu = CentralDifferenceWrtMu(p_mu, mu_S, r_S, r_M, dr_SEP)
# partial_mu = PartialWrtMu(mu_S, r_S, r_M, dr_SEP, da)
#
# print("central difference wrt mu: \n", cd_mu)
# print("partial wrt mu: \n", partial_mu)
# print("partial - central difference wrt mu: \n", partial_mu-cd_mu)


# p_eta = Decimal(1.0E-4)
# eta = Decimal(1.0E-3)
# cd_eta = CentralDifferenceWrtEta(p_eta, mu_S, r_S, r_M, dr_SEP, eta)
# partial_eta = PartialWrtEta(mu_S, r_S, r_M, dr_SEP, eta)
#
# print("central difference wrt eta: \n", cd_eta)
# print("partial wrt eta: \n", partial_eta)
# print("partial - central difference wrt eta: \n", partial_eta-cd_eta)








