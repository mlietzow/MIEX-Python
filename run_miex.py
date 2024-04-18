import miex
import numpy as np
import sys
import os
import time

# ****************************************************************************************************
# MIEX: MIE SCATTERING CODE FOR LARGE GRAINS (https://doi.org/10.1016/j.cpc.2004.06.070)
#                    _________________________________________________________________________________
#                    Contact information:   wolf@astrophysik.uni-kiel.de (Sebastian Wolf)
#                                           rbrunngraeber@astrophysik.uni-kiel.de (Robert Brunngräber)
#                                           mlietzow@astrophysik.uni-kiel.de (Moritz Lietzow)
# ====================================================================================================
# GENERAL CODE DESCRIPTION
# ------------------------
#
# Based on Mie scattering theorie, the following quantities for
#
#       (a) single grain sizes / chemical components
#   and (b) mixtures of chemically different grains with a size distributions
#
# can be calculated:
#
#   * Scattering matrix elements S11, S12, S33, and S34
#   * Extinction effiency factor         (Q_ext) & Extinction cross section     (C_ext)
#   * Scattering effiency factor         (Q_sca) & Scattering cross section     (C_sca)
#   * Absorption effiency factor         (Q_abs) & Absorption cross section     (C_abs)
#   * Backscattering effiency factor     (Q_bk)  & Backscattering cross section (C_bk)
#   * Radiation pressure effiency factor (Q_pr)
#   * Albedo
#   * Scattering assymetry factor (g).
#
# ____________________________________________________________________________________________________
# The optical data of the grains have to be provided in files with the following tabular form
#   * first  row: wavelength [micron]
#   * second row: n (=real[ri])
#   * second row: k (=imag[ri]).
#
# Rem.: For astrophysical applications, such tables can be found, e.g., at
#       http://www.astro.uni-jena.de/Users/database/entry.html
#       ('Jena-Petersburg Database of Optical Constants'). See, also,
#       N.V.Voshchinnikov: 'Optics of Cosmic Dust', Astrophysics and Space Physics Review 12,  1 (2002)
#       for further references.
#
# ====================================================================================================


def conv(line):
    # Converter for lambda/n/k database files
    return line.replace(b"D", b"e")


def main(input_filename):
    # ---------------------------------------------------------------------------------------------------
    # 0. General settings
    # ---------------------------------------------------------------------------------------------------
    with open(input_filename) as input_file:
        # Real refractive index of the surrouding medium
        refmed = float(input_file.readline())
        print(f"Real refractive index of surrouding medium    : {refmed}")

        # ---------------------------------------------------------------------------------------------------
        # 1. Get main parameters
        # ---------------------------------------------------------------------------------------------------
        nlam = int(input_file.readline())
        ncomp = int(input_file.readline())
        print(f"Number of wavelengths                         : {nlam}")
        print(f"Number of chemical components                 : {ncomp}")

        fnames = []
        abun = np.ones(ncomp) * 100.0

        # ---------------------------------------------------------------------------------------------------
        print("Name of the dust data files (lambda/n/k data)")
        print("    [all data files have to contain the refractive]")
        print("    [index and the same wavelength distribution   ]")
        for icomp in range(ncomp):
            fname = input_file.readline().rstrip()
            print(f"    {icomp+1}. component : {fname}")
            fnames.append(fname)

        if ncomp > 1:
            print("Relative abundances of the different components [%]")
            for icomp in range(ncomp):
                ab = float(input_file.readline())
                print(f"    {icomp+1}. component : {ab}")
                abun[icomp] = ab

        abun /= 100.0

        ask1 = int(input_file.readline())
        print("-1- Single grain size")
        print(f"-2- Grain size distribution                   : {ask1}")
        if ask1 == 1:
            radmin = float(input_file.readline())
            radmax = radmin
            nrad = 1
            alpha = 0.0
            print(f"    Grain radius [micron]                     : {radmin}")
        else:
            radmin = float(input_file.readline())
            radmax = float(input_file.readline())
            alpha = float(input_file.readline())
            nrad = int(input_file.readline())
            print(f"    Minimum grain size [micron]               : {radmin}")
            print(f"    Maximum grain size [micron]               : {radmax}")
            print(f"    Size distribution exponent                : {alpha}")
            print(f"    Number of size bins                       : {nrad}")

        ask2 = int(input_file.readline())
        print(f"Calculate scattering matrix elements (0=n/1=y): {ask2}")
        if ask2 == 1:
            nang2 = int(input_file.readline())
            print("Number of scattering angles in the interval")
            print(" [0 deg,180 deg];  odd number!")
            print(f" [example: '181' -> step width = 1 deg]       : {nang2}")
            if nang2 % 2 == 1:
                nang = int((nang2 - 1) / 2 + 1)
                doSA = True
            else:
                raise Exception("Number of scattering angles must be odd! Aborting.")
        else:
            nang = 1
            nang2 = 1
            doSA = False

        fresult = input_file.readline().rstrip()
        print(f"Project name                                  : {fresult}")

        ask3 = int(input_file.readline())
        print(f"Save results in separate files (0=n/1=y)      : {ask3}")
        if ask3 == 0:
            svsep = False
        else:
            svsep = True

    # ---------------------------------------------------------------------------------------------------
    # 2. Read data files & Prepare the calculations
    # ---------------------------------------------------------------------------------------------------
    wavelength = np.zeros(nlam)
    n_real = np.zeros((ncomp, nlam))
    k_imag = np.zeros((ncomp, nlam))

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

    print(" >> Calculation started ...")

    # read lambda/n/k database
    for icomp in range(ncomp):
        w, n, k = np.loadtxt(
            "ri-data/" + fnames[icomp], comments="#", unpack=True, converters=conv
        )
        wavelength = w[:nlam]
        n_real[icomp] = n[:nlam]
        k_imag[icomp] = k[:nlam]
        # first two lines of file give information about the content of the file; no need to save in a variable

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
    print("0 %", end="\r")
    for ilam in range(nlam):
        weisum = 0.0
        wrad = 0.0
        wqsc = 0.0

        for icomp in range(ncomp):
            for irad in range(nrad):
                # show progress every 10 per cent
                counter += 1
                if int(counter % ((nlam * ncomp * nrad) / 10)) == 0:
                    print(int(counter / (nlam * ncomp * nrad) * 100), "%", end="\r")

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
                ri = complex(n_real[icomp, ilam], k_imag[icomp, ilam]) / refmed

                # derive the scattering parameters
                q_extx, q_absx, q_scax, q_bkx, q_prx, albedox, g_scax, S1x, S2x = (
                    miex.shexqnn2(x, ri, nang, doSA)
                )

                # update average values
                weight = abun[icomp] * rad**alpha * delrad
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
    os.makedirs(os.path.dirname("results/" + fresult), exist_ok=True)
    with open("results/" + fresult, "w") as output_file:
        output_file.write("# *** PROJECT PARAMETERS ***\n")
        output_file.write("#\n")
        output_file.write(f"# Number of wavelengths            : {nlam}\n")
        output_file.write(f"# Number of chemical components    : {ncomp}\n")
        output_file.write("# Relative abundances [%]          :\n")
        for icomp in range(ncomp):
            output_file.write(f"# {icomp+1}. component: {abun[icomp]*100}\n")
        output_file.write("# Name(s) of the dust data file(s) :\n")
        for icomp in range(ncomp):
            output_file.write(f"# {icomp+1}. component: {fnames[icomp]}\n")
        output_file.write(f"# Minimum grain radius [micron]    : {radmin}\n")
        output_file.write(f"# Maximum grain radius [micron]    : {radmax}\n")
        output_file.write(f"# Size distribution exponent       : {alpha}\n")
        output_file.write(f"# Number of size bins              : {nrad}\n")
        if doSA:
            output_file.write(f"# Number of scattering angles      : {nang2}\n")
        output_file.write("#\n")
        output_file.write("#\n")
        # -----------------------------------------------------------------------------------------------
        output_file.write("# *** RESULTS ***\n")
        output_file.write("#\n")
        output_file.write(
            "# 1. Wavelength [micron]         Q_ext                    C_ext [m^2]\n"
        )
        for ilam in range(nlam):
            output_file.write(f"{wavelength[ilam]} {q_ext[ilam]} {c_ext[ilam]}\n")
        output_file.write(
            "# =============================================================== #\n"
        )

        output_file.write(
            "# 2. Wavelength [micron]         Q_sca                    C_sca [m^2]\n"
        )
        for ilam in range(nlam):
            output_file.write(f"{wavelength[ilam]} {q_sca[ilam]} {c_sca[ilam]}\n")
        output_file.write(
            "# =============================================================== #\n"
        )

        output_file.write(
            "# 3. Wavelength [micron]         Q_bk                      C_bk [m^2]\n"
        )
        for ilam in range(nlam):
            output_file.write(f"{wavelength[ilam]} {q_bk[ilam]} {c_bk[ilam]}\n")
        output_file.write(
            "# =============================================================== #\n"
        )

        output_file.write(
            "# 4. Wavelength [micron]         Q_abs                    C_abs [m^2]\n"
        )
        for ilam in range(nlam):
            output_file.write(f"{wavelength[ilam]} {q_abs[ilam]} {c_abs[ilam]}\n")
        output_file.write(
            "# =============================================================== #\n"
        )

        output_file.write("# 5. Wavelength [micron], Albedo\n")
        for ilam in range(nlam):
            output_file.write(f"{wavelength[ilam]} {albedo[ilam]}\n")
        output_file.write(
            "# =============================================================== #\n"
        )

        output_file.write("# 6. Wavelength [micron], Scattering asymmetry factor g\n")
        for ilam in range(nlam):
            output_file.write(f"{wavelength[ilam]} {g_sca[ilam]}\n")
        output_file.write(
            "# =============================================================== #\n"
        )

        output_file.write("# 7. Wavelength [micron]         Q_pr\n")
        for ilam in range(nlam):
            output_file.write(
                f"{wavelength[ilam]} {q_ext[ilam] - g_sca[ilam] * q_sca[ilam]}\n"
            )
        output_file.write(
            "# =============================================================== #\n"
        )

        output_file.write(
            "# 8. Wavelength [micron]     theta [degree]         S11-S12-S33-S34\n"
        )
        if doSA:
            for ilam in range(nlam):
                for iang in range(nang2):
                    # angx: scattering angle [degree]
                    angx = iang * 180.0 / (nang2 - 1.0)
                    output_file.write(
                        f"{wavelength[ilam]} {angx} {S11[iang,ilam]} {S12[iang,ilam]} {S33[iang,ilam]} {S34[iang,ilam]}\n"
                    )
        else:
            output_file.write("# not calculated.\n")

    # ---------------------------------------------------------------------------------------------------
    if svsep:
        # 1. Extinction efficiency factor
        with open("results/" + fresult + ".q_ext", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {q_ext[ilam]}\n")

        # 2. Extinction cross section [m^2]
        with open("results/" + fresult + ".c_ext", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {c_ext[ilam]}\n")

        # 3. Scattering efficiency factor
        with open("results/" + fresult + ".q_sca", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {q_sca[ilam]}\n")

        # 4. Scattering cross section [m^2]
        with open("results/" + fresult + ".c_sca", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {c_sca[ilam]}\n")

        # 5. Backscattering effiency factor
        with open("results/" + fresult + ".q_bk", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {q_bk[ilam]}\n")

        # 6. Backscattering cross section [m^2]
        with open("results/" + fresult + ".c_bk", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {c_bk[ilam]}\n")

        # 7. Absorption efficieny factor
        with open("results/" + fresult + ".q_abs", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {q_abs[ilam]}\n")

        # 8. Absorption cross section [m^2]
        with open("results/" + fresult + ".c_abs", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {c_abs[ilam]}\n")

        # 9. Albedo
        with open("results/" + fresult + ".alb", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {albedo[ilam]}\n")

        # 10. Scattering asymmetry factor g
        with open("results/" + fresult + ".g_sca", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(f"{wavelength[ilam]} {g_sca[ilam]}\n")

        # 11. Q_pr
        with open("results/" + fresult + ".q_pr", "w") as output_file:
            for ilam in range(nlam):
                output_file.write(
                    f"{wavelength[ilam]} {q_ext[ilam] - g_sca[ilam] * q_sca[ilam]}\n"
                )

        # 12. Scattering Matrix elements S11, S12, S33, and S34
        if doSA:
            # S11
            with open("results/" + fresult + ".S11", "w") as output_file:
                for ilam in range(nlam):
                    output_file.write(f"{wavelength[ilam]}\n")
                    for iang in range(nang2):
                        # angx: scattering angle [degree]
                        angx = iang * 180.0 / (nang2 - 1.0)
                        output_file.write(f"{angx} {S11[iang,ilam]}\n")

            # S12
            with open("results/" + fresult + ".S12", "w") as output_file:
                for ilam in range(nlam):
                    output_file.write(f"{wavelength[ilam]}\n")
                    for iang in range(nang2):
                        # angx: scattering angle [degree]
                        angx = iang * 180.0 / (nang2 - 1.0)
                        output_file.write(f"{angx} {S11[iang,ilam]}\n")

            # S33
            with open("results/" + fresult + ".S33", "w") as output_file:
                for ilam in range(nlam):
                    output_file.write(f"{wavelength[ilam]}\n")
                    for iang in range(nang2):
                        # angx: scattering angle [degree]
                        angx = iang * 180.0 / (nang2 - 1.0)
                        output_file.write(f"{angx} {S11[iang,ilam]}\n")

            # S34
            with open("results/" + fresult + ".S34", "w") as output_file:
                for ilam in range(nlam):
                    output_file.write(f"{wavelength[ilam]}\n")
                    for iang in range(nang2):
                        # angx: scattering angle [degree]
                        angx = iang * 180.0 / (nang2 - 1.0)
                        output_file.write(f"{angx} {S11[iang,ilam]}\n")

    print("    ... done.")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise Exception(f"Input file missing")
    else:
        start_time = time.time()
        main(sys.argv[1])
        print(f"--- {time.time() - start_time:.3f} seconds ---")
