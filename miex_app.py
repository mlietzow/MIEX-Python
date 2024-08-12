from src.miex import miex
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st


# Converter for lambda/n/k database files
def conv(line):
    return line.replace("D", "e")


st.set_page_config(page_title="MIEX", page_icon=None)
st.title("MIEX App")
st.write("This app is a Mie scattering code for large grains written in Python and based on [MIEX](https://ui.adsabs.harvard.edu/abs/2018ascl.soft10019W) by [Wolf & Voshchinnikov (2004)](https://ui.adsabs.harvard.edu/abs/2004CoPhC.162..113W).")
st.write("""
    The following quantities for

    1. single grain sizes / chemical components and
    2. mixtures of chemically different grains with a size distribution

    can be calculated:

    - Scattering matrix elements $S_{11}$, $S_{12}$, $S_{33}$, and $S_{34}$,
    - Extinction efficiency factor ($Q_\\mathrm{ext}$) and Extinction cross-section ($C_\\mathrm{ext}$),
    - Scattering efficiency factor ($Q_\\mathrm{sca}$) and Scattering cross-section ($C_\\mathrm{sca}$),
    - Absorption efficiency factor ($Q_\\mathrm{abs}$) and Absorption cross-section ($C_\\mathrm{abs}$),
    - Backscattering efficiency factor ($Q_\\mathrm{bk}$) and Backscattering cross-section ($C_\\mathrm{bk}$),
    - Radiation pressure efficiency factor ($Q_\\mathrm{pr}$),
    - Albedo,
    - Scattering asymmetry factor ($g$).
""")

st.divider()

radio_wavelength = st.radio(
    "Single wavelength or upload dust data files (three columns: wavelength/micron real imag):",
    options=["single", "upload"],
    horizontal=True,
)

col1, col2 = st.columns(2)
if radio_wavelength == "single":
    col1, col2 = st.columns(2)
    with col1:
        input_ri_real = np.array(
            [
                st.number_input(
                    "Real part of refractive index:",
                    value=1.50,
                    format="%e",
                    step=0.01,
                    min_value=0.0,
                )
            ]
        )
        input_wavelength = np.array(
            [
                st.number_input(
                    "Wavelength $\\lambda$ [micron]:", value=1.0, format="%e", step=0.1, min_value=0.0
                )
            ]
        )
    with col2:
        input_ri_imag = np.array(
            [
                st.number_input(
                    "Imaginary part of refractive index:",
                    value=0.0,
                    format="%e",
                    step=0.01,
                    min_value=0.0,
                )
            ]
        )

    nlam = 1
    ncomp = 1
    abun = np.array([1.0])
else:
    st.info(
        "all data files have to contain the refractive index and the same wavelength distribution"
    )
    with col1:
        ncomp = st.number_input(
            "Number of chemical components:", value=1, format="%d", step=1, min_value=1
        )
    with col2:
        nlam = st.number_input(
            "Number of wavelengths:", value=100, format="%d", step=1, min_value=1
        )
    fnames = []
    abun = np.ones(ncomp) * 100.0
    for icomp in range(ncomp):
        with col1:
            fnames.append(
                st.file_uploader(f"Choose {icomp+1}. component:", key=f"file{icomp}")
            )
        with col2:
            abun[icomp] = st.number_input(
                f"Relative abundance of the {icomp+1}. component [%]:",
                value=100.0,
                format="%f",
                step=1.0,
                key=f"abun{icomp}",
                min_value=0.0,
                max_value=100.0,
            )

    abun /= 100.0

st.divider()

radio_grain = st.radio(
    "Single grain size or grain size distribution:",
    options=["single", "distribution"],
    horizontal=True,
)

col1, col2 = st.columns(2)
if radio_grain == "single":
    with col1:
        radmin = st.number_input(
            "Grain radius $r$ [micron]:", value=1.0, format="%e", step=0.1, min_value=0.0
        )

    radmax = radmin
    exponent = 0.0
    nrad = 1
    dist_type = ""
else:
    dist_type = st.radio(
        "Distribution type",
        options=["Power law", "Power law with exponential decay"],
        horizontal=True,
    )
    st.caption("Power law: $n(r) \\propto r^q$")
    st.caption(
        "Power law with exponential decay: $n(r) \\propto r^q \\times \\exp(-r / p)$"
    )
    with col1:
        radmin = st.number_input(
            "Minimum grain size $r_{\\rm min}$ [micron]:",
            value=0.01,
            format="%e",
            step=0.1,
        )
        exponent = st.number_input(
            "Size distribtion exponent $q$", value=-3.5, format="%f", step=0.1, max_value=0.0
        )
        if "exponential" in dist_type:
            parameter2 = st.number_input(
                "Exponential decay parameter $p$", value=1.0, format="%e", step=0.1, min_value=0.0
            )
    with col2:
        radmax = st.number_input(
            "Maximum grain size $r_{\\rm max}$ [micron]:",
            value=1.0,
            format="%e",
            step=0.1,
        )
        nrad = st.number_input(
            "Number of size bins:", value=100, format="%d", step=1, min_value=1
        )

st.divider()

col1, col2 = st.columns(2)
with col1:
    doSA = st.checkbox("Calculate scattering matrix elements", value=False)
with col2:
    nang2 = st.number_input(
        "Number of scattering angles in the intervall $0$ to $\\pi$ (equidistantly distributed, must be odd):",
        value=91,
        format="%d",
        step=2,
        min_value=1,
        disabled=not doSA,
    )

st.divider()

col1, col2 = st.columns(2)
with col1:
    run_miex = st.button("Run")
with col2:
    plot_res = st.checkbox("Plot results", value=False)

st.divider()
st.write("The original source code written in FORTRAN90 is distributed under the [CPC license](https://www.elsevier.com/about/policies/open-access-licenses/elsevier-user-license/cpc-license). Modified and ported to Python with permission of S. Wolf.")

if run_miex:
    # ---------------------------------------------------------------------------------------------------
    # 2. Read data files & Prepare the calculations
    # ---------------------------------------------------------------------------------------------------
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

    S11 = np.zeros((nang2, nlam))
    S12 = np.zeros((nang2, nlam))
    S33 = np.zeros((nang2, nlam))
    S34 = np.zeros((nang2, nlam))

    if np.sum(abun) != 1:
        st.warning("Warning: The sum of the relative abundances is not 100 %")

    if nang2 % 2 == 1:
        nang = int(0.5 * (nang2 - 1) + 1)
    else:
        st.error("Error: Number of scattering angles must be odd!")
        st.stop()

    if radio_wavelength == "single":
        wavelength = input_wavelength
        ri_real[0] = input_ri_real
        ri_imag[0] = input_ri_imag
    else:
        # read lambda/n/k database
        for icomp in range(ncomp):
            if fnames[icomp] is None:
                st.error("Error: Dust data file missing!")
                st.stop()
            w, n, k = np.loadtxt(
                fnames[icomp], unpack=True, comments="#", converters=conv
            )

            if isinstance(w, float):
                if nlam == 1:
                    wavelength = np.array([w])
                    ri_real[icomp] = np.array([n])
                    ri_imag[icomp] = np.array([k])
                else:
                    st.error(
                        "Error: Number of defined wavelengths is larger than the number of wavelengths given in the file"
                    )
                    st.stop()
            else:
                if nlam <= len(w):
                    wavelength = w[:nlam]
                    ri_real[icomp] = n[:nlam]
                    ri_imag[icomp] = k[:nlam]
                else:
                    st.error(
                        "Error: Number of defined wavelengths is larger than the number of wavelengths given in the file"
                    )
                    st.stop()

    # define radial step width
    radminlog = np.log10(radmin)
    radmaxlog = np.log10(radmax)
    if nrad > 1:
        steplog = (radmaxlog - radminlog) / (nrad - 1.0)
    else:
        steplog = 0.0

    # ---------------------------------------------------------------------------------------------------
    # 3. Run the Mie scattering routines
    # ---------------------------------------------------------------------------------------------------
    counter = 0
    placeholder = st.empty()
    with placeholder.container():
        progress_text = "Calculating optical properties ..."
        progress_bar = st.progress(0, text=progress_text)
    for ilam in range(nlam):
        weisum = 0.0
        wrad = 0.0
        wqsc = 0.0
        refmed = 1.0

        for icomp in range(ncomp):
            for irad in range(nrad):
                # show progress every 10 per cent
                counter += 1
                if int(counter % ((nlam * ncomp * nrad) / 10)) == 0:
                    progress_bar.progress(
                        int(counter / (nlam * ncomp * nrad) * 100), text=progress_text
                    )

                # current radius / radius interval
                rad = 10.0 ** (radminlog + irad * steplog)
                rad1 = 10.0 ** (radminlog + (irad + 1.0) * steplog)
                if nrad > 1:
                    delrad = rad1 - rad
                else:
                    delrad = 1.0

                # size parameter
                x = 2.0 * np.pi * rad * refmed / wavelength[ilam]

                # complex refractive index
                ri = complex(ri_real[icomp, ilam], ri_imag[icomp, ilam]) / refmed

                try:
                    # derive the scattering parameters
                    q_extx, q_absx, q_scax, q_bkx, q_prx, albedox, g_scax, S1x, S2x = (
                        miex.shexqnn2(x, ri, nang, doSA)
                    )
                except Exception as e:
                    st.error(e)
                    st.stop()

                # update average values
                weight = abun[icomp] * rad**exponent * delrad
                if "exponential" in dist_type:
                    weight *= np.exp(-rad / parameter2)
                weisum = weisum + weight

                wradx = np.pi * (rad * 1.0e-6) ** 2 * weight
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

                S11[:, ilam] += S11x * weight
                S12[:, ilam] += S12x * weight
                S33[:, ilam] += S33x * weight
                S34[:, ilam] += S34x * weight

        c_ext[ilam] /= weisum
        c_sca[ilam] /= weisum
        c_bk[ilam] /= weisum
        c_abs[ilam] /= weisum

        q_ext[ilam] /= wrad
        q_sca[ilam] /= wrad
        q_bk[ilam] /= wrad
        q_abs[ilam] /= wrad

        S11[:, ilam] /= weisum
        S12[:, ilam] /= weisum
        S33[:, ilam] /= weisum
        S34[:, ilam] /= weisum

        albedo[ilam] = c_sca[ilam] / c_ext[ilam]
        g_sca[ilam] /= wqsc

    # ---------------------------------------------------------------------------------------------------
    # 4. Save the results
    # ---------------------------------------------------------------------------------------------------
    output_file = "# *** PROJECT PARAMETERS ***\n"
    output_file += "#\n"
    output_file += f"# Number of wavelengths            : {nlam}\n"
    output_file += f"# Number of chemical components    : {ncomp}\n"
    if radio_wavelength == "single":
        output_file += f"# Real part of refractive index    : {ri_real[0,0]}\n"
        output_file += f"# Imag part of refractive index    : {ri_imag[0,0]}\n"
        output_file += f"# Wavelength [micron]              : {wavelength[0]}\n"
    else:
        output_file += "# Relative abundances [%]          :\n"
        for icomp in range(ncomp):
            output_file += f"# {icomp+1}. component: {abun[icomp]*100}\n"
        output_file += "# Name(s) of the dust data file(s) :\n"
        for icomp in range(ncomp):
            output_file += f"# {icomp+1}. component: {fnames[icomp].name}\n"
    if radio_grain == "single":
        output_file += f"# Grain radius [micron]            : {radmin}\n"
    else:
        output_file += f"# Minimum grain radius [micron]    : {radmin}\n"
        output_file += f"# Maximum grain radius [micron]    : {radmax}\n"
        output_file += f"# Size distribution exponent       : {exponent}\n"
        if "exponential" in dist_type:
            output_file += f"# Exponential decay parameter      : {parameter2}\n"
        output_file += f"# Number of size bins              : {nrad}\n"
    if doSA:
        output_file += f"# Number of scattering angles      : {nang2}\n"
    output_file += "#\n"
    output_file += "#\n"
    # -----------------------------------------------------------------------------------------------
    output_file += "# *** RESULTS ***\n"
    output_file += "#\n"
    output_file += (
        "# 1. Wavelength [micron]         Q_ext                    C_ext [m^2]\n"
    )
    for ilam in range(nlam):
        output_file += f"{wavelength[ilam]} {q_ext[ilam]} {c_ext[ilam]}\n"
    output_file += (
        "# =============================================================== #\n"
    )

    output_file += (
        "# 2. Wavelength [micron]         Q_sca                    C_sca [m^2]\n"
    )
    for ilam in range(nlam):
        output_file += f"{wavelength[ilam]} {q_sca[ilam]} {c_sca[ilam]}\n"
    output_file += (
        "# =============================================================== #\n"
    )

    output_file += (
        "# 3. Wavelength [micron]         Q_bk                      C_bk [m^2]\n"
    )
    for ilam in range(nlam):
        output_file += f"{wavelength[ilam]} {q_bk[ilam]} {c_bk[ilam]}\n"
    output_file += (
        "# =============================================================== #\n"
    )

    output_file += (
        "# 4. Wavelength [micron]         Q_abs                    C_abs [m^2]\n"
    )
    for ilam in range(nlam):
        output_file += f"{wavelength[ilam]} {q_abs[ilam]} {c_abs[ilam]}\n"
    output_file += (
        "# =============================================================== #\n"
    )

    output_file += "# 5. Wavelength [micron], Albedo\n"
    for ilam in range(nlam):
        output_file += f"{wavelength[ilam]} {albedo[ilam]}\n"
    output_file += (
        "# =============================================================== #\n"
    )

    output_file += "# 6. Wavelength [micron], Scattering asymmetry factor g\n"
    for ilam in range(nlam):
        output_file += f"{wavelength[ilam]} {g_sca[ilam]}\n"
    output_file += (
        "# =============================================================== #\n"
    )

    output_file += "# 7. Wavelength [micron]         Q_pr\n"
    for ilam in range(nlam):
        output_file += f"{wavelength[ilam]} {q_ext[ilam] - g_sca[ilam] * q_sca[ilam]}\n"
    output_file += (
        "# =============================================================== #\n"
    )

    output_file += (
        "# 8. Wavelength [micron]     theta [degree]         S11-S12-S33-S34\n"
    )
    if doSA:
        for ilam in range(nlam):
            for iang in range(nang2):
                # angx: scattering angle [degree]
                angx = iang * 180.0 / (nang2 - 1.0)
                output_file += f"{wavelength[ilam]} {angx} {S11[iang,ilam]} {S12[iang,ilam]} {S33[iang,ilam]} {S34[iang,ilam]}\n"
    else:
        output_file += "# not calculated.\n"

    # ---------------------------------------------------------------------------------------------------
    # 5. Plot the results
    # ---------------------------------------------------------------------------------------------------

    if plot_res:
        placeholder.info("Plotting results ...")

        if nlam == 1:
            data_dict = {
                "wavelength": wavelength,
                "Q_ext": q_ext,
                "C_ext [m^2]": c_ext,
                "Q_abs": q_abs,
                "C_abs [m^2]": c_abs,
                "Q_sca": q_sca,
                "C_sca [m^2]": c_sca,
                "Q_bk": q_bk,
                "C_bk [m^2]": c_bk,
                "Qpr": q_ext - g_sca * q_sca,
                "A": albedo,
                "gsca": g_sca,
            }

            df = pd.DataFrame(data_dict)
            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    key: st.column_config.NumberColumn(format="%e") for key in data_dict
                },
            )

        else:
            fig, ax = plt.subplots(
                3, 1, sharex=True, figsize=(6.4, 6.4), layout="constrained"
            )

            ax[0].plot(wavelength, c_ext, label="extinction")
            ax[0].plot(wavelength, c_abs, label="absorption")
            ax[0].plot(wavelength, c_sca, label="scattering")
            ax[0].plot(wavelength, c_bk, label="backscattering")

            ax[0].set_ylabel(r"cross section [m$^2$]")
            ax[0].set_yscale("log")
            ax[0].legend()

            ax[1].plot(wavelength, q_ext, label="extinction")
            ax[1].plot(wavelength, q_abs, label="absorption")
            ax[1].plot(wavelength, q_sca, label="scattering")
            ax[1].plot(wavelength, q_bk, label="backscattering")
            ax[1].plot(wavelength, q_ext - g_sca * q_sca, label="radiation pressure")

            ax[1].set_ylabel("efficiency factor")
            ax[1].set_yscale("log")
            ax[1].legend()

            ax[2].plot(wavelength, albedo, label="single scattering albedo")
            ax[2].plot(wavelength, g_sca, label="scattering assymetry factor")

            # ax[2].set_yscale("log")
            ax[2].set_xlabel("wavelength [micron]")
            ax[2].set_xscale("log")
            ax[2].legend()

            st.pyplot(fig, use_container_width=True)

        if doSA:
            if nlam == 1:
                fig, ax = plt.subplots(
                    2, 2, sharex=True, figsize=(6.4, 4.8), layout="constrained"
                )
                theta = np.linspace(0, 180, nang2)

                fig.suptitle(f"wavelength: {wavelength[0]} [micron]")

                ax[0, 0].plot(theta, S11[:, 0] / S11[0, 0])
                ax[0, 0].set_ylabel("S11 / S11(0 deg)")
                ax[0, 0].set_yscale("log")

                ax[0, 1].plot(theta, -S12[:, 0] / S11[:, 0])
                ax[0, 1].set_ylabel("-S12 / S11")
                ax[0, 1].yaxis.set_label_position("right")
                ax[0, 1].yaxis.tick_right()

                ax[1, 0].plot(theta, S33[:, 0] / S11[:, 0])
                ax[1, 0].set_xlabel("scattering angle [deg]")
                ax[1, 0].set_ylabel("S33 / S11")

                ax[1, 1].plot(theta, S34[:, 0] / S11[:, 0])
                ax[1, 1].set_xlabel("scattering angle [deg]")
                ax[1, 1].set_ylabel("S34 / S11")
                ax[1, 1].yaxis.set_label_position("right")
                ax[1, 1].yaxis.tick_right()

            else:
                fig, ax = plt.subplots(
                    2,
                    2,
                    sharex=True,
                    sharey=True,
                    figsize=(6.4, 4.8),
                    layout="constrained",
                )

                theta = np.linspace(0, 180, nang2)

                im = ax[0, 0].pcolormesh(
                    theta,
                    wavelength,
                    np.transpose(S11 / S11[0, :]),
                    vmax=1,
                    vmin=0,
                    shading="nearest",
                )
                ax[0, 0].set_title("S11 / S11(0 deg)")
                ax[0, 0].set_ylabel("wavelength [micron]")
                ax[0, 0].set_yscale("log")
                fig.colorbar(im, ax=ax[0, 0])

                im = ax[0, 1].pcolormesh(
                    theta,
                    wavelength,
                    np.transpose(-S12 / S11),
                    vmin=-1,
                    vmax=1,
                    cmap="bwr",
                )
                ax[0, 1].set_title("-S12 / S11")
                fig.colorbar(im, ax=ax[0, 1])

                im = ax[1, 0].pcolormesh(
                    theta,
                    wavelength,
                    np.transpose(S33 / S11),
                    vmin=-1,
                    vmax=1,
                    cmap="bwr",
                )
                ax[1, 0].set_title("S33 / S11")
                ax[1, 0].set_xlabel("scattering angle [deg]")
                ax[1, 0].set_ylabel("wavelength [micron]")
                fig.colorbar(im, ax=ax[1, 0])

                im = ax[1, 1].pcolormesh(
                    theta,
                    wavelength,
                    np.transpose(S34 / S11),
                    vmin=-1,
                    vmax=1,
                    cmap="bwr",
                )
                ax[1, 1].set_title("S34 / S11")
                ax[1, 1].set_xlabel("scattering angle [deg]")
                fig.colorbar(im, ax=ax[1, 1])

            st.pyplot(fig, use_container_width=True)

    with placeholder.container():
        st.success("Done!")
        st.download_button("Download results", output_file, file_name="results.txt")
