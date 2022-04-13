import numpy as np
import matplotlib.pyplot as plt

import help_laba

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

from math import sqrt, log

def import_file(C, L, S):
    file = open('C{}_L{}_S{}_dL0,5.dat'.format(C, L, S))
    file_str = file.read().splitlines()
    lam = []
    trans = []
    for i in range(2, len(file_str)):
        lam.append(float(file_str[i].split()[0].replace(',', '.')))
        trans.append(float(file_str[i].split()[1].replace(',', '.')))

    file.close()
    return lam, trans


def spectr_plot(lam, trans, C, L, S):
    plt.figure(figsize=(10, 6))
    plt.plot(lam, trans, lw=0.7)
    plt.xlabel(r'$\lambda$, nm', size=16)
    plt.ylabel(r'Transmittance, %', size=16)
    plt.grid()
    plt.xticks(np.arange(480, 680, 5))
    plt.yticks(np.arange(45, 81, 1))
    plt.savefig(r'specrt_C{}_L{}_S{}_dL0,5.png'.format(C, L, S))
    plt.close()


def parabola(x, y, x_,  v__, k):
    C, B, A, s_C, s_B, s_A = help_laba.parabola(x, y)
    plt.plot(x_, C * np.power(x_, 2) + B * x_ + A, '-', color='red', lw=1)
    plt.plot(x, y, '.', color='black')
    plt.grid()
    plt.xlabel('$v^{,}+1/2$', size=14)
    plt.ylabel(r'$\nu + {}$,'.format(k) + ', см$^{-1}$', size=14)
    plt.savefig('parabola_v''{}.png'.format(v__))
    plt.close()
    nu_00 = A + 0.5 * B - 0.5*214.5
    s_nu_00 = sqrt(s_A ** 2 + (0.5 * s_B) ** 2)
    print('Параболическая регрессия для v__='+str(v__))
    print('A = ' + str(A) + '+-' + str(s_A))
    print('B = ' + str(B) + '+-' + str(s_B))
    print('C = ' + str(C) + '+-' + str(s_C))
    print('nu_00 = ' + str(nu_00) + '+-' + str(s_nu_00))

    return nu_00, s_nu_00, A, s_A, B, s_B, C, s_C


lam_dt_45, trans_dt_45 = import_file('45', '480_655', '0.1')
lam_dt_55, trans_dt_55 = import_file('55', '480_655', '0.1')
lam_dt_70, trans_dt_70 = import_file('70', '400_700', '0.1')

lam_ds_05, trans_ds_05 = import_file('55', '480_655', '0.5')
lam_ds_1, trans_ds_1 = import_file('55', '480_655', '1')
lam_ds_5, trans_ds_5 = import_file('55', '480_655', '5')

spectr_plot(lam_dt_45, trans_dt_45, '45', '480_655', '0.1')
spectr_plot(lam_dt_55, trans_dt_55, '55', '480_655', '0.1')
spectr_plot(lam_dt_70[401:1280], trans_dt_70[401:1280], '70', '400_700', '0.1')
spectr_plot(lam_ds_05, trans_ds_05, '55', '480_655', '0.1')
spectr_plot(lam_ds_1, trans_ds_1, '55', '480_655', '1')
spectr_plot(lam_ds_5, trans_ds_5, '55', '480_655', '5')

# v''=0 с 20 по 35
l1 = np.array([559.000000000009, 556.200000000008, 553.400000000008, 550.800000000007, 548.400000000006,
               545.800000000006, 543.600000000005, 541.200000000005, 539.000000000004, 537.000000000004,
               534.800000000003, 533.000000000003, 531.000000000003, 529.200000000002, 527.400000000002,
               525.800000000001])
# v''=1 с 1 по 22
l2 = np.array([636.400000000026, 632.400000000026, 627.800000000025, 624.000000000024, 619.600000000023,
               615.400000000022, 611.200000000021, 607.40000000002, 603.400000000019, 599.600000000018,
               595.800000000017, 592.000000000016, 588.600000000016, 585.200000000015, 581.600000000014,
               578.000000000013, 574.800000000012, 571.600000000012, 568.600000000011, 565.60000000001,
               562.80000000001, 560.000000000009])

nu1 = 1/(l1 * 10**(-7))
nu2 = 1/(l2 * 10**(-7))

print('Регрессия v''=0' + str(nu1) + ' см-1')
print('Регрессия v''=1' + str(nu2) + ' см-1')


Delta_G = []
v_ = []
v1 = np.arange(1, 23, 1)
v0 = np.arange(20, 36, 1)

for i in range(len(l2)-1):
    Delta_G.append(nu2[i+1] - nu2[i])
    v_.append(v1[i+1])

for i in range(len(l1)-1):
    Delta_G.append(nu1[i+1] - nu1[i])
    v_.append(v0[i+1])

a_, b_, s_b_, s_a_ = help_laba.LSM(np.array(v_[2:]), np.array(Delta_G[2:]))
x1 = np.arange(1, 36, 1)
y1 = a_ * x1 + b_
plt.plot(x1, y1, color='blue', lw=1)
plt.plot(v_, Delta_G, '.', color='blue')
plt.xlabel('$v^{,}+1$', size=14)
plt.ylabel(r'$\Delta G^{,}_{v^{,}+1/2}$, см$^{-1}$', size=14)
plt.grid()
plt.savefig('lunear.png')
plt.close()

print('Линейная регрессия')
print('omega e = ' + str(b_) + '+-' + str(s_b_) + ' см-1')
print('omega chi e = ' + str(a_/(-2)) + '+-' + str(1/2 * s_a_) + ' см-1')
print('Корр матрица = ' + str(np.corrcoef(v_, Delta_G)))
print('Корр коэф = ' + str(np.corrcoef(v_, Delta_G)[0][1]))


y_0 = nu1 + 1/2 * 214.5
x_0 = v0 + 1/2

y_1 = nu2 + 3/2 * 214.5
x_1 = v1 + 1/2

nu_00_0, s_nu_00_0, A0, s_A0, B0, s_B0, C0, s_C0 = parabola(x_0, y_0, np.arange(20, 36.05, 0.05), 0, 1/2 * 214.5)
nu_00_1, s_nu_00_1, A1, s_A1, B1, s_B1, C1, s_C1 = parabola(x_1, y_1, np.arange(1, 24.05, 0.05), 1, 3/2 * 214.5)

nu_00 = np.mean(np.array([nu_00_0, nu_00_1]))
s_nu_00 = sqrt(np.std(np.array([nu_00_0, nu_00_1])) ** 2 + s_nu_00_0 ** 2)

A = np.mean(np.array([A0, A1]))
s_A = sqrt(np.std(np.array([A0, A1])) ** 2 +s_A0 ** 2)

B = np.mean(np.array([B0, B1]))
s_B = sqrt(np.std(np.array([B0, B1])) ** 2 +s_B0 ** 2)

C = np.mean(np.array([C0, C1]))
s_C = sqrt(np.std(np.array([C0, C1])) ** 2 +s_C1 ** 2)

print('Средние значения для параболической регрессии')
print('A = ' + str(A) + '+-' + str(s_A))
print('B = ' + str(B) + '+-' + str(s_B))
print('C = ' + str(C) + '+-' + str(s_C))
print('nu_00 = ' + str(nu_00) + '+-' + str(s_nu_00))

delta_nu = np.zeros(len(nu1)-1)
for i in range(len(nu1)-1):
    delta_nu[i] = nu1[i+1] - nu1[i]

c, b, a, s_c, s_b, s_a = help_laba.parabola(delta_nu, nu1[0:len(nu1)-1])
X = np.arange(55, 95, 0.05)
plt.plot(X, c * np.power(X, 2) + b * X + a, '-', color='red', lw=1)
plt.plot(delta_nu, nu1[0:len(nu1)-1], '.', color='black')
plt.ylabel(r'$\nu$, см$^{-1}$', size=14)
plt.xlabel(r'$\Delta \nu$, см$^{-1}$', size=14)
plt.grid()
plt.savefig('parabola_nu.png')
plt.close()

print('nu гр = ' + str(a) + '+-' + str(s_a))

print('На основе параболической регресии были взяты omega e и omega e chi')
print('Определение D0')

D0_a = B ** 2/ (4 * (-1 * C))
s_D0_a = D0_a * sqrt((2 * s_B/B) ** 2 + (s_C/C) ** 2)

print('a) = ' + str(D0_a) + '+-' + str(s_D0_a))

D0_b = a - nu_00
s_D0_b = sqrt(s_a ** 2 + s_nu_00 ** 2)

print('b) = ' + str(D0_b) + '+-' + str(s_D0_b))

v_ = np.array(v_) - 1/2
k, n, p, s_k, s_n, s_p = help_laba.parabola(v_, Delta_G)
xx = np.arange(1, 35, 0.05)
plt.plot(xx, k * np.power(xx, 2) + n * xx + p, '-', color='red', lw=1)
plt.plot(v_, Delta_G, '.', color='black')
plt.grid()
plt.xlabel(r'$v^{,} + 1/2$', size=14)
plt.ylabel(r'$\Delta G^{,}_{v^{,}+1/2}$, см$^{-1}$', size=14)
plt.savefig('Berja_Shponer.png')
plt.close()

print('Экстраполяция параболой kx^2+nx+p (Берджа-Шпонер)')
print('k =' + str(k) + '+-' + str(s_k))
print('n =' + str(n) + '+-' + str(s_n))
print('p =' + str(p) + '+-' + str(s_p))

print('Пересечение с нулем в 45.2672578989056, площадь через вольфрам 3412.21 - D0c ')

D0__ = a - 7603

print('Энергия диссоциации D0__ =' + str(D0__) + '+-' + str(s_a))

nu_max = 1/(531.000000000003 * 10**(-7))

print('Nu max (32) = ' + str(nu_max))

U_ = nu_max - nu_00 + b_/2
s_U_ = sqrt(s_nu_00 ** 2 + (1/2 * s_b_) ** 2)

De_ = D0_a + b_/2
s_De_ = sqrt(s_D0_a ** 2 + (1/2 * s_b_) ** 2)

beta_ = 0.12177 * b_ * sqrt(126.90447/(2 * De_))
s_beta = beta_ * sqrt((s_b_/b_)**2 + (1/2 * s_De_/De_) ** 2)

re_ = 2.67 + 1/beta_ * log(1 + sqrt(U_/De_))
s_re_ = re_ * sqrt((s_beta/beta_)**2 + (1/2*s_U_/U_)**2 + (1/2*s_De_/De_)**2)

print('U_ = ' + str(U_) + '+-' + str(s_U_))
print('De_ = ' + str(De_) + '+-' + str(s_De_))
print('beta_ = ' + str(beta_) + '+-' + str(s_beta))
print('re_ = ' + str(re_) + '+-' + str(s_re_))


r = np.arange(2.3, 8, 0.00001)
U_r = De_*np.power((1 - np.exp(-beta_ * (r - re_))), 2)

plt.plot(r, U_r, color='blue', lw=1)
plt.xlabel(r'$r, \AA$', size=14)
plt.ylabel(r'$U(r)$, см$^{-1}$', size=14)
plt.grid()
plt.savefig('Morse.png')
plt.close()

plt.figure(figsize=(10, 6))
plt.plot(lam_dt_45, trans_dt_45, label=r'$45 C^\circ$', lw=0.7)
plt.plot(lam_dt_55, trans_dt_55, label=r'$55 C^\circ$', lw=0.7)
plt.plot(lam_dt_70[401:1280], trans_dt_70[401:1280], label=r'$70 C^\circ$', lw=0.7)
plt.xticks(np.arange(480, 680, 20))
plt.xlabel(r'$\lambda$, nm', size=16)
plt.ylabel(r'Transmittance, %', size=16)
plt.legend(prop={'size': 16})
# plt.grid()
plt.savefig('разные темп.png')
plt.close()

plt.figure(figsize=(10, 6))
plt.plot(lam_dt_55, trans_dt_55, label='0,1 нм', lw=0.7)
plt.plot(lam_ds_05, trans_ds_05, label='0,5 нм', lw=0.7)
plt.plot(lam_ds_1, trans_ds_1, label='1 нм', lw=0.7)
plt.plot(lam_ds_5, trans_ds_5, label='5 нм', lw=0.7)
plt.xticks(np.arange(480, 680, 20))
plt.xlabel(r'$\lambda$, nm', size=16)
plt.ylabel(r'Transmittance, %', size=16)
plt.legend(prop={'size': 16})
plt.savefig('разные щели.png')
plt.close()
