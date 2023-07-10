import miex
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st


# Converter for lambda/n/k database files
def conv(line):
    return line.replace(b'D', b'e')


st.set_page_config(page_title='MIEX', page_icon=None, layout="wide")
st.title('MIEX')
st.text('MIEX is a Mie scattering code for large grains')

st.divider()

radio_wavelength = st.radio(
    'Single wavelength or upload a dust data files (three columns: wavelength/micron real imag):',
    options=['single', 'upload'], horizontal=True)

if radio_wavelength == 'single':
    col1, col2 = st.columns(2)
    with col1:
        ri_real = np.array([st.number_input('Real part of refractive index:', value=1.50, format='%f', step=0.01, min_value=1.0)])
    with col2:
        ri_imag = np.array([st.number_input('Imaginary part of refractive index:', value=0.0, format='%f', step=0.01, min_value=0.0)])

    wavelength = np.array([st.number_input(r'Wavelength $\lambda$ [micron]:', value=1.0, format='%f', step=0.1)])
    nlam = 1
    abun = np.array([1.0])
else:
    col1, col2 = st.columns(2)
    with col1:
        ncomp = st.number_input('Number of chemical components:', value=1, format='%d', step=1, min_value=1)
    with col2:
        nlam = st.number_input('Number wavelegnths:', value=100, format='%d', step=1, min_value=1)
    fnames = []
    abun = np.ones(ncomp) * 100.0
    for icomp in range(ncomp):
        with col1:
            fnames.append(st.file_uploader(f'Choose {icomp+1}. component:', key=f'file{icomp}'))
        with col2:
            abun[icomp] = st.number_input(f'Relative abundance of the {icomp+1}. component [%]:', value=100.0, format='%f', step=1.0, key=f'abun{icomp}')
    
    abun /= 100.0

st.divider()

radio_grain = st.radio(
    'Single grain size or grain size distribution $n(r) \propto n^q$:',
    options=['single', 'distribution'], horizontal=True)

if radio_grain == 'single':
    radmin = st.number_input(r'Grain radius $r$ [micron]:', value=1.0, format='%f', step=0.1)

    radmax = radmin
    exponent = 0.0
    nrad = 1
else:
    col1, col2 = st.columns(2)
    with col1:
        radmin = st.number_input(r'Minimum grain size $r_{\rm min}$ [micron]:', value=0.01, format='%f', step=0.1)
        exponent = st.number_input(r'Size distribtion Exponent $q$', value=-3.5, format='%f', step=0.1)
    with col2:
        radmax = st.number_input(r'Maximum grain size $r_{\rm max}$ [micron]:', value=1.0, format='%f', step=0.1)
        nrad = st.number_input(r'Number of size bins:', value=100, format='%d', step=1, min_value=1)

st.divider()

col1, col2 = st.columns(2)
with col1:
    doSA = st.checkbox('Calculate scattering matrix elements', value=False)
with col2:
    nang = st.number_input(
        r'Half number of scattering angles in the intervall $0$ to $\pi/2$ (equidistantly distributed):',
        value=91, format='%d', step=1, min_value=1, disabled=not doSA)

st.divider()

if st.button('Run'):
    wavelength = np.zeros(nlam)
    ri_real = np.zeros((ncomp, nlam))
    ri_imag = np.zeros((ncomp, nlam))

    q_ext = np.zeros(nlam)
    q_sca = np.zeros(nlam)
    q_abs = np.zeros(nlam)
    q_bk = np.zeros(nlam)
    # q_pr = np.zeros(nlam)

    c_ext = np.zeros(nlam)
    c_sca = np.zeros(nlam)
    c_abs = np.zeros(nlam)
    c_bk = np.zeros(nlam)
    # c_pr = np.zeros(nlam)

    albedo = np.zeros(nlam)
    g_sca = np.zeros(nlam)

    nang2 = 2 * nang - 1
    S11 = np.zeros((nang2, nlam))
    S12 = np.zeros((nang2, nlam))
    S33 = np.zeros((nang2, nlam))
    S34 = np.zeros((nang2, nlam))

    if np.sum(abun) != 1:
        st.warning('Warning: The sum of the relative abundances is not 100 %')

    # read lambda/n/k database
    for icomp in range(ncomp):
        if fnames[icomp] is None:
            st.error('Error: Dust data file missing!')
            st.stop()
        w, n, k = np.loadtxt(fnames[icomp], unpack=True, comments='#', converters=conv)
        wavelength = w[:nlam]
        ri_real[icomp] = n[:nlam]
        ri_imag[icomp] = k[:nlam]
        # first two lines of file give information about the content of the file; no need to save in a variable

    placeholder = st.empty()
    with placeholder.container():
        progress_text = 'Calculating optical properties ...'
        progress_bar = st.progress(0, text=progress_text)

    radminlog = np.log10(radmin)
    radmaxlog = np.log10(radmax)
    if nrad > 1:
        steplog = (radmaxlog - radminlog) / (nrad - 1.0)
    else:
        steplog = 0.0

    counter = 0
    for ilam in range(nlam):
        weisum = 0.0
        wrad = 0.0
        wqsc = 0.0
        refmed = 1.0

        for icomp in range(ncomp):
            for irad in range(nrad):
                # show progress every 1 per cent
                counter += 1
                if int(counter % ((nlam * ncomp * nrad) / 10)) == 0:
                    progress_bar.progress(int(counter / (nlam * ncomp * nrad) * 100), text=progress_text)

                # current radius / radius interval
                rad = 10.0**(radminlog + irad * steplog)
                rad1 = 10.0**(radminlog + (irad + 1.0) * steplog)
                if nrad > 1:
                    delrad = rad1 - rad
                else:
                    delrad = 1.0

                # size parameter
                x = 2.0 * np.pi * rad * refmed / wavelength[ilam]

                # complex refractive index
                ri = complex(ri_real[icomp,ilam], ri_imag[icomp,ilam]) / refmed

                # derive the scattering parameters
                q_extx, q_absx, q_scax, q_bkx, q_prx, albedox, g_scax, S1x, S2x = miex.shexqnn2(x, ri, nang, doSA)

                # update average values
                weight = abun[icomp] * rad**exponent * delrad
                weisum = weisum + weight

                wradx = np.pi * (rad * 1.0e-6)**2 * weight
                wqscx = wradx * q_scax

                wrad += wradx
                wqsc += wqscx

                c_ext[ilam] += q_extx * wradx
                c_sca[ilam] += q_scax * wradx
                c_bk[ilam] += q_bkx * wradx
                c_abs[ilam] += q_absx * wradx

                q_ext[ilam] += q_extx * wradx
                q_sca[ilam] += q_scax * wradx
                q_bk[ilam] += q_bkx * wradx
                q_abs[ilam] += q_absx * wradx

                g_sca[ilam] += g_scax * wqscx

                S11x, S12x, S33x, S34x = miex.scattering_matrix_elements(S1x, S2x)

                S11[:,ilam] += S11x * weight
                S12[:,ilam] += S12x * weight
                S33[:,ilam] += S33x * weight
                S34[:,ilam] += S34x * weight

        c_ext[ilam] /= weisum
        c_sca[ilam] /= weisum
        c_bk[ilam] /= weisum
        c_abs[ilam] /= weisum

        q_ext[ilam] /= wrad
        q_sca[ilam] /= wrad
        q_bk[ilam] /= wrad
        q_abs[ilam] /= wrad

        S11[:,ilam] /= weisum
        S12[:,ilam] /= weisum
        S33[:,ilam] /= weisum
        S34[:,ilam] /= weisum

        albedo[ilam] = c_sca[ilam] / c_ext[ilam]
        g_sca[ilam] /= wqsc

    placeholder.info('Plotting results ...')

    data_dict = {
        'wavelength': wavelength,
        'Qext': q_ext,
        'Qabs': q_abs,
        'Qsca': q_sca,
        'Qbk': q_bk,
        'Qpr': q_ext - g_sca * q_sca,
        'A': albedo,
        'gsca': g_sca,
    }

    if nlam == 1:
        df = pd.DataFrame(data_dict)
        st.dataframe(
            df, use_container_width=True, hide_index=True,
            column_config={key: st.column_config.NumberColumn(format='%e') for key in data_dict},
            )
    else:
        fig, ax = plt.subplots(2, 1, sharex=True, layout='constrained')

        ax[0].plot(wavelength, q_ext, label='extinction')
        ax[0].plot(wavelength, q_abs, label='absorption')
        ax[0].plot(wavelength, q_sca, label='scattering')
        ax[0].plot(wavelength, q_bk, label='backscattering')
        ax[0].plot(wavelength, q_ext - g_sca * q_sca, label='radiation pressure')

        ax[0].set_ylabel('efficiency factor')
        ax[0].set_yscale('log')
        ax[0].legend()

        ax[1].plot(wavelength, albedo, label='single scattering albedo')
        ax[1].plot(wavelength, q_sca, label='scattering assymetry factor')
        
        # ax[1].set_yscale('log')
        ax[1].set_xlabel('wavelength [micron]')
        ax[1].set_xscale('log')
        ax[1].legend()

        st.pyplot(fig, use_container_width=True)

    if doSA:
        if nlam == 1:
            fig, ax = plt.subplots(2, 2, sharex=True, layout='constrained')
            theta = np.linspace(0, 180, nang2)

            fig.suptitle(f'wavelength: {wavelength[0]} [micron]')

            ax[0,0].plot(theta, S11[:,0] / S11[0,0])
            ax[0,0].set_ylabel('S11 / S11(0 deg)')
            ax[0,0].set_yscale('log')

            ax[0,1].plot(theta, -S12[:,0] / S11[:,0])
            ax[0,1].set_ylabel('-S12 / S11')
            ax[0,1].yaxis.set_label_position("right")
            ax[0,1].yaxis.tick_right()

            ax[1,0].plot(theta, S33[:,0] / S11[:,0])
            ax[1,0].set_xlabel('scattering angle [deg]')
            ax[1,0].set_ylabel('S33 / S11')

            ax[1,1].plot(theta, S34[:,0] / S11[:,0])
            ax[1,1].set_xlabel('scattering angle [deg]')
            ax[1,1].set_ylabel('S34 / S11')
            ax[1,1].yaxis.set_label_position("right")
            ax[1,1].yaxis.tick_right()

        else:
            fig, ax = plt.subplots(2, 2, sharex=True, sharey=True, layout='constrained')

            theta = np.linspace(0, 180, nang2)

            im = ax[0,0].pcolormesh(theta, wavelength, np.transpose(S11 / S11[0,:]), vmax=1, vmin=0, shading='nearest')
            ax[0,0].set_title('S11 / S11(0 deg)')
            ax[0,0].set_ylabel('wavelength [micron]')
            ax[0,0].set_yscale('log')
            fig.colorbar(im, ax=ax[0,0])

            im = ax[0,1].pcolormesh(theta, wavelength, np.transpose(-S12 / S11), vmin=-1, vmax=1, cmap='bwr')
            ax[0,1].set_title('-S12 / S11')
            fig.colorbar(im, ax=ax[0,1])

            im = ax[1,0].pcolormesh(theta, wavelength, np.transpose(S33 / S11), vmin=-1, vmax=1, cmap='bwr')
            ax[1,0].set_title('S33 / S11')
            ax[1,0].set_xlabel('scattering angle [deg]')
            ax[1,0].set_ylabel('wavelength [micron]')
            fig.colorbar(im, ax=ax[1,0])

            im = ax[1,1].pcolormesh(theta, wavelength, np.transpose(S34 / S11), vmin=-1, vmax=1, cmap='bwr')
            ax[1,1].set_title('S34 / S11')
            ax[1,1].set_xlabel('scattering angle [deg]')
            fig.colorbar(im, ax=ax[1,1])

        st.pyplot(fig, use_container_width=True)

    placeholder.success('Done!')
