import numpy as np
import matplotlib.pyplot as plt
import os
import help_laba

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

from math import sqrt

def import_file(name):
    file = open(os.path.join('files', name))
    file_str = file.read().splitlines()
    lam = []
    abs = []
    for i in range(8, len(file_str)):
        lam.append(float(file_str[i].split()[0].replace(',', '.')))
        abs.append(float(file_str[i].split()[1].replace(',', '.')))

    file.close()
    return lam, abs

lam0_05, abs0_05 = import_file(r'воздух 0.5 без апод.ascii')
lam0_1, abs0_1 = import_file(r'воздух 1 без апод.ascii')
lam0_2, abs0_2 = import_file(r'воздух 2 без апод.ascii')
lam0_4, abs0_4 = import_file(r'воздух 4 без апод.ascii')

lam0_a05, abs0_a05 = import_file(r'воздух 0.5 апод.ascii')
lam0_a4, abs0_a4 = import_file(r'воздух 4 апод.ascii')


# срез [3117:3679]
plt.figure(figsize=(10, 8))
plt.plot(lam0_4, abs0_4, '-', label='4', lw=0.5)
plt.plot(lam0_2, abs0_2, '-', label='2', lw=0.5)
plt.plot(lam0_1, abs0_1, '-', label='1', lw=0.5)
plt.plot(lam0_05, abs0_05, '-', label='0.5', lw=0.5)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.legend(prop={'size': 14})
plt.savefig('воздух без апод.png')
plt.close()

plt.figure(figsize=(10, 8))
plt.plot(lam0_a4, abs0_a4, '-', label='4', lw=1)
plt.plot(lam0_05, abs0_05, '-', label='0.5', lw=1)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.legend(prop={'size': 14})
plt.savefig('воздух апод.png')
plt.close()

# plt.plot(lam0_a4[3117:3679], abs0_a4[3117:3679], '-', label='4 с апод.', lw=0.5)
# plt.plot(lam0_4[3117:3679], abs0_4[3117:3679], '-', label='4 без апод.', lw=0.5)
plt.plot(lam0_a05[3117:3679], abs0_a05[3117:3679], '-', label='0.5', lw=0.5)
plt.plot(lam0_05[3117:3679], abs0_05[3117:3679], '-', label='0.5 без апод.', lw=0.5)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.legend(prop={'size': 14})
plt.savefig('сравнение апод и без апод 0.5.png')
plt.close()

plt.plot(lam0_05[2289-9:2912-9], abs0_05[2289-9:2912-9], '-', label='0.5 без апод.', lw=1)
plt.xlim(3600, 3900)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('h20 val.png')
plt.close()

plt.plot(lam0_05, abs0_05, '-', label='0.5 без апод.', lw=1)
plt.xlim(1200, 2100)
plt.ylim(0.5 * 10**7)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('h20 def.png')
plt.close()

plt.plot(lam0_05, abs0_05, '-', label='0.5 без апод.', lw=1)
plt.xlim(2260, 2400)
plt.ylim(0.25 * 10**7, 2.5 * 10**7)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('co2 val.png')
plt.close()

plt.plot(lam0_05, abs0_05, '-', label='0.5 без апод.', lw=1)
plt.xlim(640, 700)
plt.ylim(0, 0.7 * 10**7)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('co2 def.png')
plt.close()

lam_ahcl, abs_ahcl = import_file(r'hcl 0.5 апод.ascii')
lam_hcl, abs_hcl = import_file(r'hcl 0.5 без апод.ascii')

plt.plot(lam_hcl, abs_hcl, '-', color='red', lw=1)
plt.xlim(2600, 3150)
plt.ylim(0.6, 1)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
# plt.grid()
plt.savefig('hcl без апод.png')
plt.close()

# P-ветвь
j_p = np.arange(1, 12, 1)
u_p = np.array([2864.75, 2843.05, 2820.86, 2798.68, 2775.53, 2751.89, 2727.3, 2702.7, 2677.62, 2651.1, 2625.05])

# R-ветвь
j_r = np.arange(0, 11)
u_r = np.array([2905.74, 2925.52, 2944.32, 2963.13, 2980.5, 2997.86, 3014.26, 3029.69, 3044.64,
                3059.11, 3064.9])

print('Task2')

print('R')
print('j = ' + str(j_r))
print('u = ' + str(u_r))
print('P')
print('j = ' + str(j_p))
print('u = ' + str(u_p))





lam_31, abs_31 = import_file('нпво 3.1 1 без апод.ascii')
lam_32, abs_32 = import_file('нпво 3.2 1 без апод.ascii')
lam_33, abs_33 = import_file('нпво 3.3 1 без апод.ascii')

plt.figure(figsize=(10, 6))
plt.plot(lam_31, abs_31, '-', color='green', lw=1)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('нпво 3.1 без апод.png')
plt.close()

plt.figure(figsize=(10, 6))
plt.plot(lam_32, abs_32, '-', color='green', lw=1)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('нпво 3.2 без апод.png')
plt.close()

plt.figure(figsize=(10, 6))
plt.plot(lam_33, abs_33, '-', color='green', lw=1)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('нпво 3.3 без апод.png')
plt.close()

del_F_ = u_r[1:] - u_p[0:10]
j_ = j_p[0:10] + 1/2

del_F__ = u_r[0:10] - u_p[1:]
j__ = j_p[1:]-1

a1, b1, sa1, sb1 = help_laba.LSM(j_, del_F_)
a2, b2, sa2, sb2 = help_laba.LSM(j__, del_F__)

x = np.arange(0, 13, 1)

plt.plot(x, a1 * x + b1, '-', color='darkblue', lw=0.7)
plt.plot(j_, del_F_, '.', color='red', markersize=8)
plt.xlabel(r'$j+1/2$', size=14)
plt.ylabel(r'$\Delta_2 F^{,}(j)$, см$^{-1}$', size=14)
plt.grid()
plt.savefig(r'd2F_.png')
plt.close()

plt.plot(x, a2 * x + b2, '-', color='darkblue', lw=0.7)
plt.plot(j__, del_F__, '.', color='red', markersize=8)
plt.xlabel(r'$j+1/2$', size=14)
plt.ylabel(r'$\Delta_2 F^{,,}(j)$, см$^{-1}$', size=14)
plt.grid()
plt.savefig(r'd2F__.png')
plt.close()

print('ax+b')
print('d2F_ 4B1 = '+str(a1) + '+-' + str(sa1) + ', b1 = ' + str(b1) + '+-' + str(sb1))
print('d2F__ 4B0 = '+str(a2) + '+-' + str(sa2) + ', b2 = ' + str(b2) + '+-' + str(sb2))

B1 = a1/4
B0 = a2/4

sB1 = sa1/4
sB0 = sa2/4

print('B1='+str(B1)+'+-'+str(sB1))
print('B0='+str(B0)+'+-'+str(sB0))

nev_ = del_F_ - a1*j_ + b1
nev__ = del_F__ - a2*j__ + b2

plt.plot(j_, nev_, '.', markersize=8)
plt.xlabel(r'$j+1/2$', size=14)
plt.ylabel(r'Невязка $\Delta_2 F^{,}(j)$, см$^{-1}$', size=14)
plt.grid()
plt.savefig(r'd2F_ невязка.png')
plt.close()

plt.plot(j__, nev__, '.', markersize=8)
plt.xlabel(r'$j+1/2$', size=14)
plt.ylabel(r'Невязка $\Delta_2 F^{,,}(j)$, см$^{-1}$', size=14)
plt.grid()
plt.savefig(r'd2F__ невязка.png')
plt.close()

y1 = del_F_/j_
y2 = del_F__/j__

x1 = np.power(j_, 2)
x2 = np.power(j__, 2)

k1, m1, sk1, sm1 = help_laba.LSM(x1, y1)
k2, m2, sk2, sm2 = help_laba.LSM(x2, y2)

X = np.arange(0, 110)

plt.plot(X, k1*X+m1, '-', color='darkblue', lw=0.7)
plt.plot(x1, y1, '.', color='red', markersize=8)
plt.xlabel(r'$(j+1/2)^2$', size=14)
plt.ylabel(r'$\Delta_2 F^{,}(j)/(j+1/2)$, см$^{-1}$', size=14)
plt.grid()
plt.savefig(r'd2F_ sq.png')
plt.close()

plt.plot(X, k2*X+m2, '-', color='darkblue', lw=0.7)
plt.plot(x2, y2, '.', color='red', markersize=8)
plt.xlabel(r'$(j+1/2)^2$', size=14)
plt.ylabel(r'$\Delta_2 F^{,,}(j)/(j+1/2)$, см$^{-1}$', size=14)
plt.grid()
plt.savefig(r'd2F__ sq.png')
plt.close()

print('kx+m')
print('4B1 = ' + str(m1) + '+-' + str(sm1))
print('4B0 = ' + str(m2) + '+-' + str(sm2))
print('-8D1 = ' + str(k1) + '+-' + str(sk1))
print('-8D1 = ' + str(k2) + '+-' + str(sk2))

B1_ = m1/4
sB1_ = sm1/4
B0_ = m2/4
sB0_ = sm2/4

D1 = -k1/8
sD1 = sk1/8
D0 = -k2/8
sD0 = sk2/8

print('B1 = ' + str(B1_) + '+-' + str(sB1_))
print('B0 = ' + str(B0_) + '+-' + str(sB0_))
print('D1 = ' + str(D1) + '+-' + str(sD1))
print('D1 = ' + str(D0) + '+-' + str(sD0))

a_p, b_p, c_p, sa_p, sb_p, sc_p = help_laba.parabola(j_p+1, u_p)
a_r, b_r, c_r, sa_r, sb_r, sc_r = help_laba.parabola(j_r+1, u_r)

x_ = np.arange(0, 14)
y_p = a_p * np.power(x_, 2) + b_p * x_ + c_p
y_r = a_r * np.power(x_, 2) + b_r * x_ + c_r

plt.plot(j_p+1, u_p, '.', color='red', markersize=8)
plt.plot(x_, y_p, '-', color='black', lw=1)
plt.xlabel(r'$(j+1)^2$', size=14)
plt.ylabel(r'$\nu_p$, см$^{-1}$ ', size=14)
plt.grid()
plt.savefig('для omega e p.png')
plt.close()

plt.plot(j_r+1, u_r, '.', color='red', markersize=8)
plt.plot(x_, y_r, '-', color='black', lw=1)
plt.xlabel(r'$(j+1)^2$', size=14)
plt.ylabel(r'$\nu_r$, см$^{-1}$ ', size=14)
plt.grid()
plt.savefig('для omega e r.png')
plt.close()

print('P-ветвь, параметры параболы ax^2+bx+c')
print('a = ' + str(a_p) + '+-' + str(sa_p))
print('b = ' + str(b_p) + '+-' + str(sb_p))
print('c = ' + str(c_p) + '+-' + str(sc_p))

print('R-ветвь, параметры параболы ax^2+bx+c')
print('a = ' + str(a_r) + '+-' + str(sa_r))
print('b = ' + str(b_r) + '+-' + str(sb_r))
print('c = ' + str(c_r) + '+-' + str(sc_r))

omega_e_p = c_p/(1-0.01)
s_omega_e_p = sc_p/(1-0.01)
omega_e_r = c_r/(1-0.01)
s_omega_e_r = sc_r/(1-0.01)

print('P omega e = ' + str(omega_e_p) + '+-' + str(s_omega_e_p))
print('R omega e = ' + str(omega_e_r) + '+-' + str(s_omega_e_r))

mu = 1*35.5/36.5 * 1.66 * 10**(-24)
B = (B1 + B0)/2
sB = sqrt(np.std([B1, B0]) * 2 + sB1 ** 2)
r = sqrt(2.79 * 10**(-39) / (B * mu))
s_r = r * 1/2 * sB/B

print('r = ' + str(r * 10**(8)) + '+-' + str(s_r * 10**8))

omega_e = np.mean([omega_e_r, omega_e_p])
s_omega_e = sqrt(np.std([omega_e_r, omega_e_p])**2 + s_omega_e_r ** 2)

print('среднее omega e = ' + str(omega_e) + '+-' + str(s_omega_e))

k = mu * (omega_e * 2 * 3.14 * 3 * 10**8) ** 2
sk = k * 2 * s_omega_e/omega_e
print('k = ' + str(k) + '+-' + str(sk))

mu_37 = 37/38 * 1.66 * 10**(-24)

omega_e_37 = 1/(2 * 3.14 * 3 * 10**8) * sqrt(k/mu_37)
s_omega_e37 = omega_e_37 * 1/2 * sk/k
print('omega e 37 = ' + str(omega_e_37) + '+-' + str(s_omega_e37))

s_delta_nu = sqrt(s_omega_e ** 2+ s_omega_e37**2)

print('delta nu = ' + str(omega_e - omega_e_37) + '+-' + str(s_delta_nu))

R = omega_e/(omega_e - omega_e_37)
sR = sqrt((s_omega_e/omega_e)**2 + (s_delta_nu/(omega_e - omega_e_37)) ** 2)

print('R = ' + str(R) + '+-' + str(sR))

print('Задание с Ar, Kr')

file = open(os.path.join('files', 'Ar_matrix.dpt'))
file_str = file.read().splitlines()
lam_Ar = []
abs_Ar = []
for i in range(len(file_str)):
    lam_Ar.append(float(file_str[i].split()[0]))
    abs_Ar.append(float(file_str[i].split()[1]))
file.close()


file = open(os.path.join('files', 'Kr_matrix.dpt'))
file_str = file.read().splitlines()
lam_Kr = []
abs_Kr = []
for i in range(len(file_str)):
    lam_Kr.append(float(file_str[i].split()[0]))
    abs_Kr.append(float(file_str[i].split()[1]))
file.close()

plt.figure(figsize=(10, 6))
plt.plot(lam_Ar, abs_Ar, '-', color='green', lw=1)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('Ar.png')
plt.close()

plt.figure(figsize=(10, 6))
plt.plot(lam_Kr, abs_Kr, '-', color='green', lw=1)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('Kr.png')
plt.close()


plt.plot(lam_Ar, abs_Ar, '-', color='green', lw=1)
plt.xlim(1000, 1300)
plt.ylim(-0.05, 0.2)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('Ar_cut.png')
plt.close()


d_ar = 3 / (2 * 1.3 * (1167.65203 - 1113.37585))

print('d_ar = ' + str(d_ar) + 'см')

plt.plot(lam_Kr, abs_Kr, '-', color='green', lw=1)
plt.hlines(0.00109, 1000, 1300)
plt.xlim(1000, 1300)
plt.ylim(-0.008, 0.002)
plt.xlabel(r'$\nu$, см$^{-1}$ ', size=14)
plt.ylabel(r'Поглощение', size=14)
plt.savefig('Kr_cut.png')
plt.close()

d_kr = 2 /(2 * 1.4 * (1290.48189-1114.44216))
print('d_kr' + str(d_kr))