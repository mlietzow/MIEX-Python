import miex
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st


st.set_page_config(page_title='MIEX', page_icon=None, layout="wide")
st.title('MIEX')
st.text('MIEX is a Mie scattering code for large grains')

st.divider()

col1, col2 = st.columns(2)
with col1:
    ri_real = st.number_input('Real part of refractive index:', value=1.50, format='%f', step=0.01, min_value=1.0)
with col2:
    ri_imag = st.number_input('Imaginary part of refractive index:', value=0.0, format='%f', step=0.01, min_value=0.0)

wavelength = st.number_input(r'Wavelength $\lambda$ / m:', value=1e-6, format='%e', step=1e-6)

st.divider()

radio_button = st.radio(
    'Single grain size or grain size distribution $n(r) \propto n^q$:',
    options=['single', 'distribution'], horizontal=True)

if radio_button == 'single':
    particle_size = st.number_input(r'Grain size $r$ / m:', value=1e-6, format='%e', step=1e-6)
    size_parameter = 2.0 * np.pi * particle_size / wavelength
    st.write(r'Size parameter $x = 2 \pi r / \lambda = $', size_parameter)
else:
    col1, col2 = st.columns(2)
    with col1:
        particle_size_min = st.number_input(r'Minimum grain size $r_{\rm min}$ / m:', value=1e-7, format='%e', step=1e-6)
        size_parameter_min = 2.0 * np.pi * particle_size_min / wavelength
        st.write(r'Minimum size parameter $x_{\rm min} = 2 \pi r_{\rm min} / \lambda = $', size_parameter_min)
    with col2:
        particle_size_max = st.number_input(r'Maximum grain size $r_{\rm max}$ / m:', value=1e-5, format='%e', step=1e-6)
        size_parameter_max = 2.0 * np.pi * particle_size_max / wavelength
        st.write(r'Maximum size parameter $x_{\rm max} = 2 \pi r_{\rm max} / \lambda = $', size_parameter_max)
    
    exponent = st.number_input(r'Size distribtion Exponent $q$', value=-3.5, format='%f', step=0.1)
    nsize = st.number_input(r'Number of size bins:', value=100, format='%d', step=1, min_value=1)

    size_parameter = np.geomspace(size_parameter_min, size_parameter_max, nsize)
    size_distribution = size_parameter**exponent

st.divider()

doSA = st.checkbox('Calculate scattering matrix elements', value=False)
nang = st.number_input(
    r'Half number of scattering angles in the intervall $0$ to $\pi/2$ (equidistantly distributed):',
    value=91, format='%d', step=1, min_value=1, disabled=not doSA)

st.divider()

if st.button('Run'):
    nang2 = 2 * nang - 1
    ri = complex(ri_real, ri_imag)
    if radio_button == 'single':
        result = miex.shexqnn2(size_parameter, ri, nang, doSA)
        matrix = miex.scattering_matrix_elements(result[7], result[8])
        st.write(r'$Q_{\rm ext} =$', result[0])
        st.write(r'$Q_{\rm abs} =$', result[1])
        st.write(r'$Q_{\rm sca} =$', result[2])
        st.write(r'$Q_{\rm bk} =$', result[3])
        st.write(r'$Q_{\rm pr} =$', result[4])
        st.write(r'$A =$', result[5])
        st.write(r'$g_{\rm sca} =$', result[6])

        if doSA:
            fig, ax = plt.subplots(2, 1, sharex=True)
            theta = np.linspace(0, 180, nang2)

            ax[0].plot(theta, matrix[0]/matrix[0][0])
            ax[0].set_ylabel('S11')
            ax[0].set_yscale('log')

            ax[1].plot(theta, -matrix[1]/matrix[0])
            ax[1].set_xlabel('scattering angle')
            ax[1].set_ylabel('-S12 / S11')

            st.pyplot(fig)
    else:
        Q_ext_tmp = np.zeros(nsize)
        Q_abs_tmp = np.zeros(nsize)
        Q_sca_tmp = np.zeros(nsize)
        Q_bk_tmp = np.zeros(nsize)
        g_sca_tmp = np.zeros(nsize)
        S_11_tmp = np.zeros((nsize, nang2))
        S_12_tmp = np.zeros((nsize, nang2))
        for s, size_param in enumerate(size_parameter):
            result = miex.shexqnn2(size_param, ri, nang, doSA)
            matrix = miex.scattering_matrix_elements(result[7], result[8])
            Q_ext_tmp[s] = result[0]
            Q_abs_tmp[s] = result[1]
            Q_sca_tmp[s] = result[2]
            Q_bk_tmp[s] = result[3]
            g_sca_tmp[s] = result[6]
            S_11_tmp[s] = matrix[0]
            S_12_tmp[s] = matrix[1]
        
        Q_ext = np.trapz(Q_ext_tmp * size_distribution, x=size_parameter)
        Q_abs = np.trapz(Q_abs_tmp * size_distribution, x=size_parameter)
        Q_sca = np.trapz(Q_sca_tmp * size_distribution, x=size_parameter)
        Q_bk = np.trapz(Q_bk_tmp * size_distribution, x=size_parameter)
        g_sca = np.trapz(g_sca_tmp * size_distribution, x=size_parameter)
        Q_pr = Q_ext - g_sca * Q_sca
        Albedo = Q_sca / Q_ext

        S_11 = np.zeros(nang2)
        S_12 = np.zeros(nang2)
        for i in range(nang2):
            S_11[i] = np.trapz(S_11_tmp[:,i] * size_distribution, x=size_parameter)
            S_12[i] = np.trapz(S_12_tmp[:,i] * size_distribution, x=size_parameter)

        st.write(r'$Q_{\rm ext} =$', result[0])
        st.write(r'$Q_{\rm abs} =$', result[1])
        st.write(r'$Q_{\rm sca} =$', result[2])
        st.write(r'$Q_{\rm bk} =$', result[3])
        st.write(r'$Q_{\rm pr} =$', result[4])
        st.write(r'$A =$', result[5])
        st.write(r'$g_{\rm sca} =$', result[6])

        if doSA:
            fig, ax = plt.subplots(2, 1, sharex=True)
            theta = np.linspace(0, 180, nang2)

            ax[0].plot(theta, S_11/S_11[0])
            ax[0].set_ylabel('S11')
            ax[0].set_yscale('log')

            ax[1].plot(theta, -S_12/S_11)
            ax[1].set_xlabel('scattering angle')
            ax[1].set_ylabel('-S12 / S11')

            st.pyplot(fig)
