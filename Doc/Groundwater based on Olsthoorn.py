import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import math
from scipy.linalg import expm, sqrtm
from scipy.sparse.linalg import spsolve

font_size = 14
plt.rcParams.update({'font.size': font_size})

# c, (nLay by nSec) matrix of vertical resistance values of top layer and aquitards
c = np.array([[175,   10,    150,   80,    75,    50],      # --> invloed redelijk: weerstand rijnstrangen van 100 naar 125 -> 3.5891m2/d naar 3.3848m2/d
              [3.2e4, 3.2e4, 3.2e4, 3.2e4, 1.0e4, 5.0e3],   # --> invloed minimaal bij veranderen in "likely range" (3.640m2/d naar 3.566m2/d)
              [1,     5.0e2, 1.0e3, 1.0e3, 1.0e3, 5.0e2]])  # --> invloed nihil

# T, (nLay by nSec) matrix of transmissivity values
T = np.array([[750, 250, 750, 750, 750, 750],
              [250, 500, 100, 50,  50,  250],
              [250, 250, 250, 250, 250, 250]])

# h, (1 by nSec) vector of heads on top of each sections
h = np.array([-1, 1, -0.5, 2.5, -1, -1.5])

# Q, (nNod by nSec) matrix of nodal injections [L2/T]
Q = np.zeros_like(c[:, 0:-1])

# x, (nNod by 1) vector of coordinates of intersection points except +/-inf
x = np.array([-3000, -2500, 0, 1500, 3000])

# X, vector of points where values will be computed
X = np.arange(-4000, 4010, 10)


def nsecn(x=x, kD=T, c=c, h=h, Q=Q, X=X, plot_phi=False, plot_q=True, plot_s=False):
    nLay = len(kD[:, 0])
    nSec = len(kD[0, :])

    # include the outer sections to infinity
    a = np.zeros((nLay, 1))
    Q = np.concatenate((a, Q, a), axis=1)
    x = np.append(x, math.inf)
    x = np.append(-math.inf, x)
    Nx = len(x)
    H = np.ones((nLay, 1)) * h

    # Mid-section points are used to compute relative coordinates within sections
    xMidSec = 0.5 * (x[:-1] + x[1:])
    xMidSec[0] = x[1]
    xMidSec[-1] = x[-2]

    # System matrices for all sections
    A = np.zeros((nLay, nLay, nSec))
    for iSec in range(nSec):
        a = 1 / (kD[:, iSec] * c[:, iSec])
        p = np.append(c[1:nLay, iSec], math.inf)
        b = 1 / (kD[:, iSec] * p)
        A[:, :, iSec] = np.diag(a + b) - np.diag(a[1:nLay], -1) - np.diag(b[0:nLay - 1], 1)

    # Generating and filling the coefficient matrix C
    C = np.zeros((nLay * (2 * (Nx - 2)), nLay * (2 * (Nx - 2) + 2)))  # coefficient matrix
    R = np.zeros((nLay * (2 * (Nx - 2)), 1))  # right hand side vector

    nNod = nSec - 1
    for i in range(nNod):
        # i is section N, also the left node number of the section
        # j is right the node number of the section
        # ii and jj point to the position within the total coefficient matrix
        j = i + 1
        ii = 2 * nLay * (i)
        jj = ii + nLay
        C[ii:jj, ii:jj] = +expm(-(x[j] - xMidSec[i]) * sqrtm(A[:, :, i]))
        C[ii:jj, ii + nLay:jj + nLay] = +expm(+(x[j] - xMidSec[i]) * sqrtm(A[:, :, i]))
        C[ii:jj, ii + 2 * nLay: jj + 2 * nLay] = -expm(-(x[j] - xMidSec[j]) * sqrtm(A[:, :, j]))
        C[ii:jj, ii + 3 * nLay: jj + 3 * nLay] = -expm(+(x[j] - xMidSec[j]) * sqrtm(A[:, :, j]))
        R[ii:jj] = np.vstack(-H[:, i] + H[:, j])

        C[ii + nLay:jj + nLay, ii:jj] = np.matmul(np.matmul(-np.diag(kD[:, i]), sqrtm(A[:, :, i])),
                                                  expm(-(x[j] - xMidSec[i]) * sqrtm(A[:, :, i])))
        C[ii + nLay:jj + nLay, ii + nLay:jj + nLay] = np.matmul(np.matmul(+np.diag(kD[:, i]), sqrtm(A[:, :, i])),
                                                                expm(+(x[j] - xMidSec[i]) * sqrtm(A[:, :, i])))
        C[ii + nLay:jj + nLay, ii + 2 * nLay:jj + 2 * nLay] = np.matmul(
            np.matmul(+np.diag(kD[:, j]), sqrtm(A[:, :, j])), expm(-(x[j] - xMidSec[j]) * sqrtm(A[:, :, j])))
        C[ii + nLay:jj + nLay, ii + 3 * nLay:jj + 3 * nLay] = np.matmul(
            np.matmul(-np.diag(kD[:, j]), sqrtm(A[:, :, j])), expm(+(x[j] - xMidSec[j]) * sqrtm(A[:, :, j])))
        R[ii + nLay:jj + nLay] = np.vstack(Q[:, j])

    # Solve the system, using all layers and leaving out the outer column as they have no freedom, because the sections extend to infinity
    COEF = np.vstack(spsolve((C[:, nLay:-nLay]), R))
    COEF = np.concatenate((np.zeros((nLay, 1)), COEF, np.zeros((nLay, 1))))

    # output:
    # phi [H] = computed heads, a (nLay by length(X)) matrix
    # q [L2/T] = computed flows, a (nLay by length(X)) matrix
    # s [L/T] = downward positive seepage rate through top of each layer, a nLay by length(X) matrix
    phi = np.zeros((nLay, len(X)))
    q = np.zeros((nLay, len(X)))
    s = np.zeros((nLay, len(X)))

    for i in range(len(X)):
        iSec = np.nonzero(np.logical_and(X[i] > x[:-1], X[i] <= x[1:]))
        iSec = iSec[0]
        iSec = iSec[0]
        k = 2 * nLay * (iSec)
        l = k + nLay

        C1 = np.matmul(expm(-(X[i] - xMidSec[iSec]) * sqrtm(A[:, :, iSec])), COEF[k:l])
        C2 = np.matmul(expm(+(X[i] - xMidSec[iSec]) * sqrtm(A[:, :, iSec])), COEF[k + nLay:l + nLay])
        C3 = np.matmul(sqrtm(A[:, :, iSec]), C1)
        C4 = np.matmul(sqrtm(A[:, :, iSec]), C2)

        phi[:, i] = np.hstack(C1) + np.hstack(C2) + (H[:, iSec])
        q[:, i] = np.hstack(np.matmul(np.diag(T[:, iSec]), (C3 - C4)))

        sNet = np.matmul(np.matmul(np.diag(T[:, iSec]), sqrtm(A[:, :, iSec])), (C3 + C4))
        s[nLay - 1, i] = sNet[nLay - 1]
        for iLay in np.arange(nLay - 2, -1, -1):
            s[iLay, i] = sNet[iLay] + s[iLay + 1, i]

    if plot_phi is True:
        plt.figure(figsize=(15, 10))
        plt.plot(X, phi[0, :], label='Head in first aquifer', c='b', linestyle='--')
        plt.plot(X, phi[1, :], label='Head in second aquifer', c='b', linestyle='-.')
        plt.plot(X, phi[2, :], label='Head in third aquifer', c='b', linestyle=':')
        plt.axvline(0, c='grey', linestyle=':')
        plt.axvline(1500, c='grey', linestyle=':')
        plt.axvline(-3000, c='grey', linestyle=':')
        plt.axvline(-2500, c='grey', linestyle=':')
        v = -1.5
        plt.plot([0, 1500], [v - 0.1, v - 0.1], c='grey')
        plt.plot([-3000, -2500], [v - 0.1, v - 0.1], c='grey')
        plt.text(750, v, 'Rijnstrangen', color='grey', horizontalalignment='center')
        plt.text(-2750, v, 'Rhine', color='grey', horizontalalignment='center')
        plt.title('Heads in aquifers', fontsize=22)
        plt.xlabel('Distance along cross-section [m]')
        plt.ylabel('Head [m]')
        plt.legend(loc='upper left', fontsize=14)
        plt.show()

    if plot_q is True:
        qsum = q[0, :] + q[1, :] + q[2, :]
        qleft = qsum[400]
        qright = qsum[550]
        plt.figure(figsize=(15, 10))
        plt.plot(X, q[0, :], label='q in first aquifer', c='b', linestyle='--')
        plt.plot(X, q[1, :], label='q in second aquifer', c='b', linestyle='-.')
        plt.plot(X, q[2, :], label='q in third aquifer', c='b', linestyle=':')
        plt.plot(X, qsum, label='Total q', c='black', linestyle='-')
        plt.plot([0, 1500], [qleft, qright], 'o', c='black')
        plt.text(80, qleft, round(qleft, 3))
        plt.text(1580, qright, round(qright, 3))
        plt.axvline(0, c='grey', linestyle=':')
        plt.axvline(1500, c='grey', linestyle=':')
        plt.axvline(-3000, c='grey', linestyle=':')
        plt.axvline(-2500, c='grey', linestyle=':')
        v = qleft - 0.8
        plt.plot([0, 1500], [v - 0.1, v - 0.1], c='grey')
        plt.plot([-3000, -2500], [v - 0.1, v - 0.1], c='grey')
        plt.text(750, v, 'Rijnstrangen', color='grey', horizontalalignment='center')
        plt.text(-2750, v, 'Rhine', color='grey', horizontalalignment='center')
        plt.suptitle('Groundwater flux', fontsize=22, y=0.95)
        plt.title('positive = northward; negative = southward', fontsize=18)
        plt.xlabel('Distance along cross-section [m]')
        plt.ylabel('Flux [m2/d]')
        plt.legend(loc='upper left', fontsize=14)
        plt.show()

    if plot_s is True:
        plt.figure(figsize=(15, 10))
        plt.plot(X, s[0, :], label='Seepage first aquifer', c='b', linestyle='--')
        plt.plot(X, s[1, :], label='Seepage second aquifer', c='b', linestyle='-.')
        plt.plot(X, s[2, :], label='Seepage third aquifer', c='b', linestyle=':')
        plt.axvline(0, c='grey', linestyle=':')
        plt.axvline(1500, c='grey', linestyle=':')
        plt.axvline(-3000, c='grey', linestyle=':')
        plt.axvline(-2500, c='grey', linestyle=':')
        v = -0.055
        plt.plot([0, 1500], [v - 0.002, v - 0.002], c='grey')
        plt.plot([-3000, -2500], [v - 0.002, v - 0.002], c='grey')
        plt.text(750, v, 'Rijnstrangen', color='grey', horizontalalignment='center')
        plt.text(-2750, v, 'Rhine', color='grey', horizontalalignment='center')
        plt.title('Seepage', fontsize=22)
        plt.xlabel('Distance along cross-section [m]')
        plt.ylabel('Seepage [m/d]')
        plt.legend(loc='upper left', fontsize=14)
        plt.show()

    return


nsecn(plot_phi=True, plot_q=True, plot_s=True)

