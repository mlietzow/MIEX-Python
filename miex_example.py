import miex
import numpy as np
import matplotlib.pyplot as plt

# ====================================================================================================

# EXAMPLE 1
# calculate scattering matrix elements for single particle size

# half number of scattering angles (0...pi/2)
nang = 91
nang2 = 2 * nang - 1
# complex refractive index
m = 1.5 + 0.0j
# size parameter
x = 10.0

# run miex and calculate scattering matrix elements
result = miex.shexqnn2(x, m, nang, doSA=True)
matrix = miex.scattering_matrix_elements(result[7], result[8])

# plot the results
fig, ax = plt.subplots(2, 1, sharex=True, layout='constrained')
theta = np.linspace(0, 180, nang2)

ax[0].plot(theta, matrix[0]/matrix[0][0])
ax[0].set_ylabel('S11')
ax[0].set_yscale('log')

ax[1].plot(theta, -matrix[1]/matrix[0])
ax[1].set_xlabel('scattering angle')
ax[1].set_ylabel('-S12 / S11')

ax[0].set_title(f'scattering matrix elements for x = {x} and m = {m}')
plt.show()

# ====================================================================================================

# EXAMPLE 2
# calculate scattering matrix elements and average over size distribution

# half number of scattering angles (0...pi/2)
nang = 91
nang2 = 2 * nang - 1
# number of size parameters
nsize = 100
# complex refractive index
m = 1.5 + 0.0j

# array of size parameter and size distribution
x_arr = np.geomspace(1e-2, 1e1, nsize)
nx_arr = x_arr**(-3.5)

# run miex and calculate scattering matrix elements
S_11_tmp = np.zeros((nsize, nang2))
S_12_tmp = np.zeros((nsize, nang2))
for i, x in enumerate(x_arr):
    result = miex.shexqnn2(x, m, nang, doSA=True)
    matrix = miex.scattering_matrix_elements(result[7], result[8])
    S_11_tmp[i] = matrix[0]
    S_12_tmp[i] = matrix[1]

S_11 = np.zeros(nang2)
S_12 = np.zeros(nang2)
for i in range(nang2):
    S_11[i] = np.trapz(S_11_tmp[:,i] * nx_arr, x=x_arr)
    S_12[i] = np.trapz(S_12_tmp[:,i] * nx_arr, x=x_arr)

# plot the results
fig, ax = plt.subplots(2, 1, sharex=True, layout='constrained')
theta = np.linspace(0, 180, nang2)

ax[0].plot(theta, S_11/S_11[0])
ax[0].set_ylabel('S11')
ax[0].set_yscale('log')

ax[1].plot(theta, -S_12/S_11)
ax[1].set_xlabel('scattering angle')
ax[1].set_ylabel('-S12 / S11')

ax[0].set_title(f'scattering matrix elements for m = {m}')
plt.show()

# ====================================================================================================

# EXAMPLE 3
# calculate efficiency factors for multiple particle sizes

# complex refractive index
m = 1.5 + 0.01j
# number of size parameters
nsize = 1000

# array of size parameter
x_arr = np.geomspace(1e-1, 1e2, nsize)

# run miex and calculate efficiency factors
Q_ext = np.zeros_like(x_arr)
Q_abs = np.zeros_like(x_arr)
Q_sca = np.zeros_like(x_arr)
for i, x in enumerate(x_arr):
    result = miex.shexqnn2(x, m)
    Q_ext[i] = result[0]
    Q_abs[i] = result[1]
    Q_sca[i] = result[2]

#plot the results
fig, ax = plt.subplots(1, 1, layout='constrained')
ax.plot(x_arr, Q_ext, label='ext')
ax.plot(x_arr, Q_abs, label='abs')
ax.plot(x_arr, Q_sca, label='sca')

ax.set_xlabel('size parameter')
ax.set_xscale('log')

ax.set_ylabel('efficiency factor')
ax.set_yscale('log')

ax.legend()
ax.set_title(f'efficiency factors for m = {m}')
plt.show()

# ====================================================================================================
