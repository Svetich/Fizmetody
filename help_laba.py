import numpy as np
from scipy.optimize import curve_fit
from math import sqrt

def LSM(x, y):
    function = lambda x, a, b: a*x + b
    popt, pcov = curve_fit(function, xdata=x, ydata=y)

    sigma_a = np.sqrt(pcov[0,0])
    sigma_b = np.sqrt(pcov[1, 1])

    return popt[0], popt[1], sigma_a, sigma_b


def chi_sq(x, y, err):
    function = lambda x, a, b: a * x + b
    popt, pcov = curve_fit(function, xdata=x, ydata=y, sigma=err)

    sigma_a = np.sqrt(pcov[0, 0])
    sigma_b = np.sqrt(pcov[1, 1])

    return popt[0], popt[1], sigma_a, sigma_b


def parabola(x, y):
    function = lambda x, a, b, c: a * np.power(x, 2) + b * x + c
    popt, pcov = curve_fit(function, xdata=x, ydata=y)

    sigma_a = np.sqrt(pcov[0, 0])
    sigma_b = np.sqrt(pcov[1, 1])
    sigma_c = np.sqrt(pcov[2, 2])

    return popt[0], popt[1], popt[2], sigma_a, sigma_b, sigma_c


def Gauss(x, y):
    function = lambda x, a, b, c: a * np.exp(-np.power((x-b), 2)/(2 * c ** 2))
    popt, pcov = curve_fit(function, xdata=x, ydata=y)

    sigma_a = np.sqrt(pcov[0, 0])
    sigma_b = np.sqrt(pcov[1, 1])
    sigma_c = np.sqrt(pcov[2, 2])

    return popt[0], popt[1], popt[2], sigma_a, sigma_b, sigma_c

def Lor(x, y):
    function = lambda x_, a, b, c: a * b/(x_**2 + b) + c
    popt, pcov = curve_fit(function, xdata=x, ydata=y)

    sigma_a = np.sqrt(pcov[0, 0])
    sigma_b = np.sqrt(pcov[1, 1])
    sigma_c = np.sqrt(pcov[2, 2])


    return popt[0], popt[1], popt[2], sigma_a, sigma_b, sigma_c