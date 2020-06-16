

import numpy as np

# n = 10000
# p = 12
# c = 2
# P  = np.random.rand(p,p)
# Hx = np.random.rand(n,p)
# W_diag = np.random.rand(n)
# W  = np.diag( W_diag )
# Hc = np.random.rand(n,c)
# C  = np.random.rand(c,c)

dir_output = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/Schettino2015/'

P = np.genfromtxt(dir_output + 'InitialCovarianceMatrix.dat')
C = np.genfromtxt(dir_output + 'ConsiderParameterAprioriCovariance.dat')[6:,6:]
Hx = np.genfromtxt(dir_output + 'test_unnormalizedPartialDerivatives.dat')[::5,:]
Hc = np.genfromtxt(dir_output + 'test_partialDerivativesOfConsiderParameters.dat')[::5,:]
W_diag = np.genfromtxt(dir_output + 'ObservationWeightDiagonal.dat')[::5]
W  = np.diag( W_diag )

n = len(Hx)
p = len(P)
c = len(C)


# regular multiplication
concov1_term1 = np.matmul(np.matmul(P, np.transpose(Hx)),W)
concov1_term2 = np.matmul(np.matmul(Hc, C), np.transpose(Hc))
concov1_term3 = np.transpose(concov1_term1)

concov1 = np.matmul(np.matmul(concov1_term1, concov1_term2), concov1_term3)


# rowwise
concov2 = np.zeros((p,p))

for i in range(n):
    Hxi = Hx[i,:]
    Hci = Hc[i,:]
    Wi = W_diag[i]

    concov2_term1 = np.matmul(P, np.transpose(Hxi)) * Wi
    concov2_term2 = np.matmul(np.matmul(Hci, C), np.transpose(Hci))
    concov2_term3 = np.transpose(concov2_term1)

    concov2i = np.matmul(concov2_term1 * concov2_term2, concov2_term3)
    concov2 = concov2 + concov2i


# reduced amount of observations
step = 2
Hxr = Hx[0::step,:]
Hcr = Hc[0::step,:]
Wr  = W_diag[0::step]
new_n = np.size(Hxr,axis=0)

concov3_term1 = np.matmul(P, np.transpose(Hxr))
for i in range(new_n):
    concov3_term1[:,i] *= Wr[i]

concov3_term2 = np.matmul(np.matmul(Hcr, C), np.transpose(Hcr))
concov3_term3 = np.transpose(concov3_term1)

concov3 = np.matmul(np.matmul(concov3_term1, concov3_term2), concov3_term3)


# try different multiplication order
concov4_term1 = np.matmul(np.matmul(P, np.transpose(Hx)),W)
concov4_term3 = np.transpose(concov4_term1)

concov4_term2a = np.matmul(concov4_term1, Hc)
concov4_term2c = np.matmul(np.transpose(Hc), concov4_term3)

concov4 = np.matmul(np.matmul(concov4_term2a, C), concov4_term2c)



error1 = concov1.diagonal()
error2 = concov2.diagonal()
error3 = concov3.diagonal()
error4 = concov4.diagonal()

print(n, new_n)

print(error1)
# print(error2)
# print(error3)
print(error4)
print(error4/error1)
print((error4-error1)/error1)

# print(error3-error1)
# print(error3/error1)
# print((error3-error1)/error1)

