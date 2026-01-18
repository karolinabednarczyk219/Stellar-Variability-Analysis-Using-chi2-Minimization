#Projekt - Pracownia Komputerowa
#Karolina Bednarczyk
#Nr indeksu: 459196
#Wyznaczenie okresu zmienności gwiazd metodą minimalizacji chi2



import math
import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
import os


#import data

def load_data(file):
    data = np.loadtxt(file)
    t = data[:, 0]
    m = data[:, 1]
    sigma = data[:, 2]
    return t, m, sigma

# weighted mean

def weighted_mean(m, sigma):
    weights = 1 / sigma**2
    A = np.sum(weights * m) / np.sum(weights)
    m_0 = m - A
    return m_0

# a and b fitting using chi2 minimization

def fit_ab_chi2(t, m_0, sigma, P):
    n = len(t)
    S_ss = np.sum(np.sin(2 * np.pi * t / P)**2 / sigma**2)
    S_sc = np.sum(np.sin(2 * np.pi * t / P) * np.cos(2 * np.pi * t / P) / sigma**2)
    S_cc = np.sum(np.cos(2 * np.pi * t / P)**2 / sigma**2)
    S_ms = np.sum(m_0 * np.sin(2 * np.pi * t / P) / sigma**2)
    S_mc = np.sum(m_0 * np.cos(2 * np.pi * t / P) / sigma**2)

    delta = S_ss * S_cc - S_sc**2

    if delta < 1e-12:
        return 0, 0, np.inf

    a = (S_ms * S_cc - S_mc * S_sc) / delta
    b = (S_mc * S_ss - S_ms * S_sc) / delta

    chi2 = np.sum(((m_0 - (a * np.sin(2 * np.pi * t / P)+ b * np.cos(2 * np.pi * t / P)))/ sigma)**2)

    return a, b, chi2
# periods of oscillation

def period_search(t, m_0, sigma, Pmin=0.1, Pmax=100.0, N=1000000):
    fmin = 1.0 / Pmax
    fmax = 1.0 / Pmin
    freqs = np.linspace(fmin, fmax, N)
    periods = 1.0 / freqs

    chi2_values = np.zeros(N)
    a_values = np.zeros(N)
    b_values = np.zeros(N)

    for i, P in enumerate(periods):
        a, b, chi2 = fit_ab_chi2(t, m_0, sigma, P)
        a_values[i] = a
        b_values[i] = b
        chi2_values[i] = chi2

    return periods, chi2_values, a_values, b_values

def plot_star(file, P_best, a, b):
    t, m, sigma = load_data(file)
    m_0 = weighted_mean(m, sigma)
    periods, chi2_values, _, _ = period_search(t, m_0, sigma, N=20000)
    star_id = os.path.basename(file).replace(".dat", "")

    # chi2 plot

    plt.figure(figsize=(7, 4))
    plt.plot(periods, chi2_values, lw=1)
    plt.axvline(P_best, color='r', ls='--', label=f"P = {P_best:.4f} d")
    plt.xscale("log")
    plt.xlabel("Period P [days]")
    plt.ylabel(r"$\chi^2$")
    plt.title(f"Star {star_id}: dependence of $\\chi^2(P)$")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # light curve in phase plot
   
    phi = (t / P_best) % 1.0
    phi_plot = np.concatenate([phi, phi + 1.0])
    m_plot = np.concatenate([m_0, m_0])
    sigma_plot = np.concatenate([sigma, sigma])

    # two period model 
    
    phi_model = np.linspace(0, 2, 1000)
    model = a * np.sin(2*np.pi*phi_model) + b * np.cos(2*np.pi*phi_model)

    plt.figure(figsize=(7, 4))
    plt.errorbar(phi_plot, m_plot, yerr=sigma_plot, fmt='.', ms=4)
    plt.plot(phi_model, model, 'r', lw=2)

    plt.xlabel("Phase")
    plt.ylabel("Δ brightness [mag]")
    plt.title(f"Star {star_id}: light curve (2 periods, P = {P_best:.4f} d)")
    plt.gca().invert_yaxis()
    plt.xlim(0, 2)
    plt.tight_layout()
    plt.show()


# main program

files = sorted(glob.glob('/Users/karolinabednarczyk/Downloads/01/*.dat'))

results = []

for file in files:
    t, m, sigma = load_data(file)
    m_0 = weighted_mean(m, sigma)

    periods, chi2_values, a_values, b_values = period_search(t, m_0, sigma)
    id = np.argmin(chi2_values)

    P_best = periods[id]
    a_best = a_values[id]
    b_best = b_values[id]
    chi2_min = chi2_values[id]

    N = len(t)
    chi2_red = chi2_min / (N - 2)

    chi2_mean = np.mean(chi2_values)
    contrast = (chi2_mean - chi2_min) / chi2_mean

    # variable star criteria

    is_variable = (chi2_red < 3) and (contrast > 0.1)

    # save results
    
    results.append((
        file,
        P_best,
        a_best,
        b_best,
        chi2_min,
        chi2_red,
        contrast,
        is_variable
    ))

    print(f"{file}  P={P_best:.4f}  variable={is_variable}")

# Save results to a text file

with open("results.txt", "w") as f:
    f.write("File\tP_best\ta\tb\tchi2\tchi2_red\tcontrast\tvariable\n")
    for r in results:
        f.write(
            f"{r[0]}\t{r[1]:.6f}\t{r[2]:.6f}\t{r[3]:.6f}\t"
            f"{r[4]:.3f}\t{r[5]:.3f}\t{r[6]:.3f}\t{str(r[7])}\n"
        )
    
# variable stars

star_id = 0

variable_stars = [r for r in results if r[7]==True]
print("Number of variable stars:", len(variable_stars))

for r in variable_stars:
    file = r[0]
    P_best = r[1]
    a = r[2]
    b = r[3]

    plot_star(file, P_best, a, b)
    plot_star(file, P_best*2, a, b)



