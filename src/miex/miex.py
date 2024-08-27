import numpy as np
from numba import njit


# ****************************************************************************************************
# MIEX: MIE SCATTERING CODE FOR LARGE GRAINS (https://doi.org/10.1016/j.cpc.2004.06.070)
#                    _________________________________________________________________________________
#                    Contact information:   wolf@astrophysik.uni-kiel.de (Sebastian Wolf)
#                                           rbrunngraeber@astrophysik.uni-kiel.de (Robert Brunngräber)
#                                           mlietzow@astrophysik.uni-kiel.de (Moritz Lietzow-Sinjen)
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
#   * Scattering assymetry factor (g)
#
# ____________________________________________________________________________________________________
# The optical data of the grains have to be provided in files with the following tabular form
#   * first  column: wavelength [micron]
#   * second column: n (=real[ri])
#   * third  column: k (=imag[ri]).
#
# Rem.: For astrophysical applications, such tables can be found, e.g., at
#       http://www.astro.uni-jena.de/Users/database/entry.html
#       ('Jena-Petersburg Database of Optical Constants'). See, also,
#       N.V.Voshchinnikov: 'Optics of Cosmic Dust', Astrophysics and Space Physics Review 12, 1 (2002)
#       for further references.
#
# ====================================================================================================


@njit(cache=True)
def aa2(ax, ri, num, ru):
    """Calculations of the ratio of derivative to the function for Bessel functions of half order with complex argument: J'(n)/J(n).
    The calculations are given by the recursive expression 'from top to bottom' beginning from n=num.
    This routine is based on the routine 'aa' published by
    N.V.Voshchinnikov: 'Optics of Cosmic Dust', Astrophysics and Space Physics Review 12, 1 (2002)

    Parameters
    ----------
    ax : float
        inverse size_parameter

    ri : complex float
        complex refractive index

    num : int
        number for subroutine

    ru : complex array like
        ru-array of results
    """

    s = ax / ri
    ru[num - 1] = (num + 1.0) * s

    for i in range(num - 1, 0, -1):
        ru[i - 1] = (i + 1.0) * s - 1.0 / (ru[i] + (i + 1.0) * s)


@njit(cache=True)
def shexqnn2(x, ri, nang=2, doSA=False, nterm=2e7, eps=1.0e-20, xmin=1.0e-06):
    """Derive quantities for a single size parameter and chemical composition.
    This routine is based on the routine 'shexqnn' published by
    N.V.Voshchinnikov: 'Optics of Cosmic Dust', Astrophysics and Space Physics Review 12, 1 (2002)

    Parameters
    ----------
    x : float
        size parameter = 2 * pi * radius / wavelength

    ri : complex float
        complex refractive index

    nang : int, optional, default = 2
        half number of scattering angles theta in the intervall 0...pi/2 (equidistantly distributed)

    doSA : bool, optional, default = False
        calculation of the scattering amplitudes

    nterm : int, optional, default = 2e7
        Maximum number of terms to be considered

    eps : float, optional, default = 1.0e-20
        Accuracy to be achieved

    xmin : float, optional, default = 1.0e-06
        Minimum size parameter

    Returns
    -------
    0:  Q_ext : float
        extinction efficiency

    1:  Q_abs : float
        absorption efficiency

    2:  Q_sca : float
        scattering efficiency

    3:  Q_bk : float
        backscattering efficiency

    4:  Q_pr : float
        radiation pressure efficiency

    5:  albedo : float
        single scattering albedo

    6:  g_sca : float
        scattering asymmetry factor

    7:  SA_1 : ndarray
        scattering amplitude function. The length of the array is 2*nang-1

    8:  SA_2 : ndarray
        scattering amplitude function. The length of the array is 2*nang-1

    Raises
    ------
    ValueError
        if Mie scattering limit is exceeded or if required terms exceed maximum number of terms
    ValueError
        if somehow calculations result in NaN
    """

    fact = np.array([1.0, 1.0e250])
    factor = 1.0e250

    if x <= xmin:
        raise ValueError(
            "Mie scattering limit",
            xmin,
            "exceeded, current size parameter:",
            x,
            "decrease default value of the argument 'xmin'",
        )

    ax = 1.0 / x
    b = 2.0 * ax**2
    ss = complex(0.0, 0.0)
    s3 = complex(0.0, -1.0)
    an = 3.0

    # Define the number for subroutine aa2 [Loskutov (1971)]
    y = np.abs(ri) * x
    num = int(1.25 * y + 15.5)
    if y < 1.0:
        num = int(7.5 * y + 9.0)
    elif y > 100.0 and y < 50000.0:
        num = int(1.0625 * y + 28.5)
    elif y >= 50000.0:
        num = int(1.005 * y + 50.5)

    if num > nterm:
        raise ValueError(
            "Maximum number of terms:",
            nterm,
            "number of terms required: ",
            num,
            "increase default value of the argument 'nterm'",
        )

    # Logarithmic derivative to Bessel function (complex argument)
    ru = np.zeros(num, dtype=np.complex128)
    aa2(ax, ri, num, ru)

    # ------------------------------------------------------------------------------------------
    # FIRST TERM
    # ------------------------------------------------------------------------------------------

    # Bessel functions
    ass = 1.0 / np.sqrt(0.5 * np.pi * ax)
    w1 = 2.0 / np.pi * ax
    Si = np.sin(x) * ax
    Co = np.cos(x) * ax

    # n = 0
    besJ0 = Si * ass
    besY0 = -Co * ass
    iu0 = 0

    # n = 1
    besJ1 = (Si * ax - Co) * ass
    besY1 = (-Co * ax - Si) * ass
    iu1 = 0
    iu2 = 0

    # Mie coefficients
    s = ru[0] / ri + ax
    s1 = s * besJ1 - besJ0
    s2 = s * besY1 - besY0
    ra0 = s1 / (s1 - s3 * s2)  # coefficient a_1

    s = ru[0] * ri + ax
    s1 = s * besJ1 - besJ0
    s2 = s * besY1 - besY0
    rb0 = s1 / (s1 - s3 * s2)  # coefficient b_1

    # efficiency factors
    r = -1.5 * (ra0 - rb0)
    Q_ext = an * np.real(ra0 + rb0)
    Q_sca = an * (np.abs(ra0) ** 2 + np.abs(rb0) ** 2)

    # scattering amplitude functions
    nang2 = 2 * nang - 1
    SA_1 = np.zeros(nang2, dtype=np.complex128)
    SA_2 = np.zeros(nang2, dtype=np.complex128)
    if doSA:
        mu = np.cos(np.linspace(0, np.pi, nang2))
        fpi0 = np.zeros(nang2)
        fpi1 = np.ones(nang2)

        fac = 1.5
        ftau = mu * fpi1 - 2.0 * fpi0

        SA_1 += fac * (ra0 * fpi1 + rb0 * ftau)
        SA_2 += fac * (ra0 * ftau + rb0 * fpi1)

        fpi1_tmp = fpi1
        fpi1 = fpi1 * mu * 3.0
        fpi1 = fpi1 - fpi0 * 2.0
        fpi0 = fpi1_tmp

    # ------------------------------------------------------------------------------------------
    # 2., 3., ... num
    # ------------------------------------------------------------------------------------------
    z = -1.0

    for iterm in range(2, num + 1):
        an = an + 2.0
        an2 = an - 2.0

        # Bessel functions
        if iu1 == iu0:
            besY2 = an2 * ax * besY1 - besY0
        else:
            besY2 = an2 * ax * besY1 - besY0 / factor

        if np.abs(besY2) > 1.0e300:
            besY2 = besY2 / factor
            iu2 = iu1 + 1

        # rbrunngraeber 10/14: Changed from besJ2 = (w1 + besY2 * besJ1) / besY1,
        # because (besY2 * besJ1) could become very large (1e+300) for large grain sizes,
        # (besY2 / besY1) is about 1, suggested by fkirchschlager
        besJ2 = besY2 / besY1
        besJ2 = w1 / besY1 + besJ2 * besJ1

        # Mie coefficients
        s = ru[iterm - 1] / ri + iterm * ax

        s1 = s * besJ2 / fact[iu2] - besJ1 / fact[iu1]
        s2 = s * besY2 * fact[iu2] - besY1 * fact[iu1]
        # coefficient a_n, (n=iterm)
        ra1 = s1 / (s1 - s3 * s2)

        s = ru[iterm - 1] * ri + iterm * ax
        s1 = s * besJ2 / fact[iu2] - besJ1 / fact[iu1]
        s2 = s * besY2 * fact[iu2] - besY1 * fact[iu1]
        # coefficient b_n, (n=iterm)
        rb1 = s1 / (s1 - s3 * s2)

        # efficiency factors
        z = -z
        rr = z * (iterm + 0.5) * (ra1 - rb1)
        r += rr
        ss += (iterm - 1.0) * (iterm + 1.0) / iterm * (
            ra0 * np.conj(ra1) + rb0 * np.conj(rb1)
        ) + an2 / iterm / (iterm - 1.0) * (ra0 * np.conj(rb0))
        qq = an * np.real(ra1 + rb1)
        Q_ext += qq
        Q_sca += an * (np.abs(ra1) ** 2 + np.abs(rb1) ** 2)

        # leaving-the-loop with error criterion
        if np.isnan(Q_ext):
            raise ValueError("Q_ext is not a number")

        # leaving-the-loop criterion
        if np.abs(qq / Q_ext) < eps:
            break

        # Bessel functions
        besJ0 = besJ1
        besJ1 = besJ2
        besY0 = besY1
        besY1 = besY2
        iu0 = iu1
        iu1 = iu2
        ra0 = ra1
        rb0 = rb1

        # scattering amplitude functions
        if doSA:
            fac = (2.0 * iterm + 1.0) / (iterm * (iterm + 1.0))
            ftau = iterm * mu * fpi1 - (iterm + 1.0) * fpi0

            SA_1 += fac * (ra0 * fpi1 + rb0 * ftau)
            SA_2 += fac * (ra0 * ftau + rb0 * fpi1)

            fpi1_tmp = fpi1
            fpi1 = fpi1 * mu * (2.0 + 1.0 / iterm)
            fpi1 = fpi1 - fpi0 * (1.0 + 1.0 / iterm)
            fpi0 = fpi1_tmp

    # efficiency factors (final calculations)
    Q_ext *= b
    Q_sca *= b
    Q_bk = 2.0 * b * np.abs(r) ** 2
    Q_pr = Q_ext - 2.0 * b * np.real(ss)
    Q_abs = Q_ext - Q_sca
    albedo = Q_sca / Q_ext
    g_sca = (Q_ext - Q_pr) / Q_sca

    return Q_ext, Q_abs, Q_sca, Q_bk, Q_pr, albedo, g_sca, SA_1, SA_2


def scattering_matrix_elements(SA_1, SA_2):
    """Calculations of the scattering matrix elements

    Parameters
    ----------
    SA_1 : array_like
        scattering amplitude function

    SA_2 : array_like
        scattering amplitude function

    Returns
    -------
    S_11, S_12, S_33, S_34 : array_like
        scattering matrix elements
    """

    S_11 = 0.5 * (np.abs(SA_2) ** 2 + np.abs(SA_1) ** 2)
    S_12 = 0.5 * (np.abs(SA_2) ** 2 - np.abs(SA_1) ** 2)
    S_33 = np.real(SA_2 * np.conj(SA_1))
    S_34 = np.imag(SA_2 * np.conj(SA_1))

    # If scattering angle is 0: SA_1 = SA_2 -> S_12 = S_34 = 0
    # Numerically, this is not always true -> set to 0 by hand
    S_12[0] = 0.0
    S_34[0] = 0.0

    # If scattering angle is pi: SA_1 = -SA_2 -> S_12 = S_34 = 0
    # Numerically, this is not always true -> set to 0 by hand
    S_12[-1] = 0.0
    S_34[-1] = 0.0

    return S_11, S_12, S_33, S_34
