import miex.miex
import unittest
import numpy as np


class test_miex(unittest.TestCase):
    """compare results by
    C.F.Bohren & D.R.Huffman 'Absorption and Scattering of Light by Small Particles' (1998)
    https://doi.org/10.1002/9783527618156
    (see Appendix A, page 482)
    and by
    W.J.Wiscombe: 'Mie Scattering Calculations: Advances in Technique and Fast, Vector-speed Computer Codes' (1979)
    https://doi.org/10.5065/D6ZP4414
    (see Appendix, A22)
    """

    def test_bh(self):
        ri = complex(1.55, 0)
        lambda0 = 0.6328
        radius = 0.525
        x = 2.0 * np.pi * radius / lambda0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=11, doSA=True)

        self.assertAlmostEqual(result["Q_ext"], 3.10543, msg="incorrect Q_ext", delta=0.00001)
        self.assertAlmostEqual(result["Q_abs"], 0.00000, msg="incorrect Q_abs", delta=0.00001)
        self.assertAlmostEqual(result["Q_sca"], 3.10543, msg="incorrect Q_sca", delta=0.00001)
        self.assertAlmostEqual(result["Q_bk"], 2.92534, msg="incorrect Q_bk", delta=0.00001)

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = -scat_mat["S_12"] / S_11
        S_33 = scat_mat["S_33"] / S_11
        S_34 = scat_mat["S_34"] / S_11
        S_11 /= S_11[0]

        bh_S_11 = np.array(
            [
                0.100000e01,
                0.785390e00,
                0.356897e00,
                0.766119e-01,
                0.355355e-01,
                0.701845e-01,
                0.574313e-01,
                0.219660e-01,
                0.125959e-01,
                0.173750e-01,
                0.124601e-01,
                0.679093e-02,
                0.954239e-02,
                0.863419e-02,
                0.227421e-02,
                0.543998e-02,
                0.160243e-01,
                0.188852e-01,
                0.195254e-01,
                0.301676e-01,
                0.383189e-01,
            ]
        )
        bh_S_12 = np.array(
            [
                0.000000e00,
                -0.459811e-02,
                -0.458541e-01,
                -0.364744e00,
                -0.534997e00,
                0.959953e-02,
                0.477927e-01,
                -0.440604e00,
                -0.831996e00,
                0.341670e-01,
                0.230462e00,
                -0.713472e00,
                -0.756255e00,
                -0.281215e00,
                -0.239612e00,
                -0.850804e00,
                -0.706334e00,
                -0.891081e00,
                -0.783319e00,
                -0.196194e00,
                0.000000e00,
            ]
        )
        bh_S_33 = np.array(
            [
                0.100000e01,
                0.999400e00,
                0.986022e00,
                0.843603e00,
                0.686967e00,
                0.959825e00,
                0.985371e00,
                0.648043e00,
                0.203255e00,
                0.795354e00,
                0.937497e00,
                -0.717397e-02,
                -0.394748e-01,
                0.536251e00,
                0.967602e00,
                0.187531e00,
                0.495254e00,
                0.453277e00,
                -0.391613e00,
                -0.962069e00,
                -0.100000e01,
            ]
        )
        bh_S_34 = np.array(
            [
                0.000000e00,
                0.343261e-01,
                0.160184e00,
                0.394076e00,
                -0.491787e00,
                -0.280434e00,
                0.163584e00,
                0.621216e00,
                -0.516208e00,
                -0.605182e00,
                0.260742e00,
                0.700647e00,
                -0.653085e00,
                -0.795835e00,
                0.795798e-01,
                -0.490882e00,
                -0.505781e00,
                -0.226817e-01,
                0.482752e00,
                0.189556e00,
                0.000000e00,
            ]
        )

        for i, j in zip(S_11, bh_S_11):
            self.assertAlmostEqual(i, j, msg="incorrect S_11", delta=0.000001)
        for i, j in zip(S_12, bh_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.000001)
        for i, j in zip(S_33, bh_S_33):
            self.assertAlmostEqual(i, j, msg="incorrect S_33", delta=0.000001)
        for i, j in zip(S_34, bh_S_34):
            self.assertAlmostEqual(i, j, msg="incorrect S_34", delta=0.000001)

    def test_w_1(self):
        ri = complex(1.5, 0)
        x = 10.0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=19, doSA=True)

        self.assertAlmostEqual(
            result["Q_ext"], 2.881999, msg="incorrect Q_ext", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_abs"], 0.000000, msg="incorrect Q_abs", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_sca"], 2.881999, msg="incorrect Q_sca", delta=0.000001
        )
        self.assertAlmostEqual(
            result["g_sca"], 0.742913, msg="incorrect g_sca", delta=0.000001
        )

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = scat_mat["S_12"]
        S_12 /= S_11

        w_S_11 = np.array(
            [
                0.520856e04,
                0.428132e04,
                0.234775e04,
                0.818121e03,
                0.149547e03,
                0.139195e02,
                0.768071e02,
                0.128393e03,
                0.789487e02,
                0.293697e02,
                0.518546e02,
                0.631908e02,
                0.341567e02,
                0.256714e02,
                0.320167e02,
                0.235919e02,
                0.194168e02,
                0.184988e02,
                0.917521e01,
                0.963035e01,
                0.136986e02,
                0.485859e01,
                0.280440e01,
                0.853248e01,
                0.439826e01,
                0.570995e01,
                0.113872e02,
                0.304142e01,
                0.675713e01,
                0.242055e02,
                0.159589e02,
                0.161562e02,
                0.582167e02,
                0.679231e02,
                0.286866e02,
                0.249122e02,
                0.423766e02,
            ]
        )
        w_S_12 = np.array(
            [
                0.0000,
                -0.0239,
                -0.0963,
                -0.2080,
                -0.2887,
                0.4685,
                0.0005,
                -0.0892,
                -0.0767,
                0.7251,
                0.4834,
                -0.0592,
                -0.0163,
                0.9793,
                0.3242,
                -0.3951,
                0.5014,
                0.9514,
                -0.0269,
                0.1127,
                0.7525,
                0.7257,
                0.8993,
                0.5987,
                -0.4844,
                0.3364,
                0.9396,
                0.6460,
                0.6282,
                0.9286,
                0.7664,
                -0.4108,
                0.2385,
                0.5276,
                0.6986,
                0.0830,
                0.0000,
            ]
        )

        for i, j in zip(S_11, w_S_11):
            self.assertAlmostEqual(
                i / j, 1.0, msg=f"incorrect S_11, {i} != {j}", delta=0.00001
            )
        for i, j in zip(S_12, w_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.0001)

    def test_w_2(self):
        ri = complex(1.5, 0)
        x = 100.0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=19, doSA=True)

        self.assertAlmostEqual(
            result["Q_ext"], 2.094388, msg="incorrect Q_ext", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_abs"], 0.000000, msg="incorrect Q_abs", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_sca"], 2.094388, msg="incorrect Q_sca", delta=0.000001
        )
        self.assertAlmostEqual(
            result["g_sca"], 0.818246, msg="incorrect g_sca", delta=0.000001
        )

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = scat_mat["S_12"]
        S_12 /= S_11

        w_S_11 = np.array(
            [
                0.275508e08,
                0.106440e06,
                0.586911e05,
                0.258367e05,
                0.202131e05,
                0.192755e05,
                0.159327e05,
                0.538269e04,
                0.866710e04,
                0.540267e04,
                0.305275e04,
                0.510959e04,
                0.167021e04,
                0.202680e04,
                0.960828e03,
                0.122361e04,
                0.394465e03,
                0.254063e03,
                0.505356e03,
                0.507283e03,
                0.131614e03,
                0.222694e03,
                0.133922e03,
                0.705029e02,
                0.401041e02,
                0.100503e03,
                0.167811e03,
                0.232785e03,
                0.243533e03,
                0.235398e03,
                0.186086e03,
                0.121783e04,
                0.353283e04,
                0.239888e04,
                0.272592e04,
                0.622951e03,
                0.434048e04,
            ]
        )
        w_S_12 = np.array(
            [
                0.0000,
                -0.0107,
                0.0098,
                0.0712,
                -0.0248,
                -0.0730,
                -0.0540,
                0.2001,
                0.0919,
                -0.1301,
                0.4545,
                0.1026,
                0.4679,
                -0.2130,
                0.7357,
                -0.0054,
                0.8365,
                -0.2535,
                -0.0016,
                -0.4214,
                -0.5649,
                -0.3412,
                0.3830,
                -0.0700,
                -0.8726,
                0.1822,
                0.5235,
                0.4015,
                0.3292,
                0.0578,
                0.5225,
                -0.8534,
                -0.8664,
                0.6181,
                0.8582,
                0.5069,
                0.0000,
            ]
        )

        for i, j in zip(S_11, w_S_11):
            self.assertAlmostEqual(
                i / j, 1.0, msg=f"incorrect S_11, {i} != {j}", delta=0.00001
            )
        for i, j in zip(S_12, w_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.0001)

    def test_w_3(self):
        ri = complex(1.5, 0)
        x = 1000.0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=19, doSA=True)

        self.assertAlmostEqual(
            result["Q_ext"], 2.013945, msg="incorrect Q_ext", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_abs"], 0.000000, msg="incorrect Q_abs", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_sca"], 2.013945, msg="incorrect Q_sca", delta=0.000001
        )
        self.assertAlmostEqual(
            result["g_sca"], 0.827882, msg="incorrect g_sca", delta=0.000001
        )

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = scat_mat["S_12"]
        S_12 /= S_11

        w_S_11 = np.array(
            [
                0.253568e12,
                0.482785e07,
                0.883654e06,
                0.163897e07,
                0.136169e07,
                0.113470e07,
                0.103882e07,
                0.929666e06,
                0.863859e06,
                0.766220e06,
                0.266763e06,
                0.324223e06,
                0.282935e06,
                0.237449e06,
                0.160593e06,
                0.525112e05,
                0.789951e05,
                0.306843e05,
                0.248427e05,
                0.214530e05,
                0.106283e05,
                0.152392e05,
                0.197453e05,
                0.912219e04,
                0.714740e04,
                0.791914e04,
                0.952024e04,
                0.109979e05,
                0.666722e04,
                0.168603e05,
                0.211737e05,
                0.559288e04,
                0.208814e06,
                0.267605e06,
                0.392870e06,
                0.224420e06,
                0.257577e07,
            ]
        )
        w_S_12 = np.array(
            [
                0.0000,
                -0.0455,
                -0.1300,
                -0.0275,
                0.0909,
                0.0469,
                0.0701,
                0.0406,
                -0.0202,
                -0.0694,
                0.5051,
                0.1317,
                -0.0010,
                -0.1018,
                -0.0523,
                0.3454,
                -0.0192,
                -0.3916,
                0.0127,
                -0.8087,
                -0.8588,
                -0.6411,
                -0.6625,
                -0.5131,
                0.3474,
                -0.4821,
                -0.2850,
                -0.2910,
                0.3604,
                -0.3089,
                -0.4512,
                0.4230,
                -0.0666,
                0.0290,
                -0.8299,
                -0.0267,
                0.0000,
            ]
        )

        for i, j in zip(S_11, w_S_11):
            self.assertAlmostEqual(
                i / j, 1.0, msg=f"incorrect S_11, {i} != {j}", delta=0.00001
            )
        for i, j in zip(S_12, w_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.0001)

    def test_w_4(self):
        ri = complex(1.5, 0)
        x = 5000.0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=19, doSA=True)

        self.assertAlmostEqual(
            result["Q_ext"], 2.008650, msg="incorrect Q_ext", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_abs"], 0.000000, msg="incorrect Q_abs", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_sca"], 2.008650, msg="incorrect Q_sca", delta=0.000001
        )
        self.assertAlmostEqual(
            result["g_sca"], 0.829592, msg="incorrect g_sca", delta=0.000001
        )

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = scat_mat["S_12"]
        S_12 /= S_11

        w_S_11 = np.array(
            [
                0.157609e15,
                0.394653e08,
                0.229931e08,
                0.339037e08,
                0.262577e08,
                0.236482e08,
                0.155814e08,
                0.323654e08,
                0.234225e08,
                0.178268e08,
                0.673414e07,
                0.586469e07,
                0.709207e07,
                0.262123e07,
                0.310573e07,
                0.158540e07,
                0.219096e07,
                0.101618e07,
                0.567675e06,
                0.323205e06,
                0.241279e06,
                0.249034e06,
                0.143024e06,
                0.417106e06,
                0.253470e06,
                0.199198e06,
                0.273150e06,
                0.254040e06,
                0.185353e06,
                0.363299e06,
                0.272294e06,
                0.256060e06,
                0.130619e08,
                0.372204e07,
                0.453092e07,
                0.205242e07,
                0.237786e09,
            ]
        )
        w_S_12 = np.array(
            [
                0.0000,
                -0.0223,
                0.1179,
                0.1247,
                0.0279,
                0.1168,
                0.2597,
                -0.0712,
                -0.0175,
                -0.0512,
                0.5594,
                0.3873,
                0.0148,
                0.8293,
                0.0368,
                0.0580,
                -0.1909,
                -0.3772,
                -0.4089,
                -0.7399,
                -0.7719,
                -0.6412,
                0.9228,
                -0.8142,
                -0.1234,
                0.0435,
                -0.3815,
                -0.1249,
                -0.1602,
                -0.6980,
                -0.7087,
                -0.0846,
                -0.8202,
                -0.1764,
                0.2883,
                0.6997,
                0.0000,
            ]
        )

        for i, j in zip(S_11, w_S_11):
            self.assertAlmostEqual(
                i / j, 1.0, msg=f"incorrect S_11, {i} != {j}", delta=0.00001
            )
        for i, j in zip(S_12, w_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.0001)

    def test_w_5(self):
        ri = complex(1.5, 0.1)
        x = 10.0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=19, doSA=True)

        self.assertAlmostEqual(
            result["Q_ext"], 2.459791, msg="incorrect Q_ext", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_abs"], 1.224646, msg="incorrect Q_abs", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_sca"], 1.235144, msg="incorrect Q_sca", delta=0.000001
        )
        self.assertAlmostEqual(
            result["g_sca"], 0.922350, msg="incorrect g_sca", delta=0.000001
        )

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = scat_mat["S_12"]
        S_12 /= S_11

        w_S_11 = np.array(
            [
                0.379171e04,
                0.300320e04,
                0.141624e04,
                0.313014e03,
                0.124235e02,
                0.296988e02,
                0.273164e02,
                0.112113e02,
                0.109517e02,
                0.607843e01,
                0.220902e01,
                0.632075e01,
                0.646946e01,
                0.225394e01,
                0.215826e01,
                0.392848e01,
                0.299433e01,
                0.163623e01,
                0.183556e01,
                0.209544e01,
                0.166228e01,
                0.137914e01,
                0.153058e01,
                0.147431e01,
                0.116521e01,
                0.135300e01,
                0.174359e01,
                0.136826e01,
                0.798073e00,
                0.974236e00,
                0.133396e01,
                0.141816e01,
                0.148012e01,
                0.126487e01,
                0.106733e01,
                0.172292e01,
                0.231818e01,
            ]
        )
        w_S_12 = np.array(
            [
                0.0000,
                -0.0014,
                -0.0068,
                -0.0301,
                -0.6183,
                -0.2518,
                -0.2817,
                -0.7359,
                -0.6378,
                -0.5132,
                -0.9235,
                -0.7194,
                -0.6077,
                -0.1744,
                -0.2426,
                -0.7454,
                -0.6373,
                0.3019,
                -0.0893,
                -0.8614,
                -0.6653,
                0.2706,
                0.0790,
                -0.7132,
                -0.8966,
                0.1033,
                0.3819,
                -0.0370,
                -0.6271,
                0.0599,
                0.3753,
                0.1218,
                -0.2643,
                -0.6463,
                -0.8175,
                -0.2177,
                0.0000,
            ]
        )

        for i, j in zip(S_11, w_S_11):
            self.assertAlmostEqual(
                i / j, 1.0, msg=f"incorrect S_11, {i} != {j}", delta=0.00001
            )
        for i, j in zip(S_12, w_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.0001)

    def test_w_6(self):
        ri = complex(1.5, 0.1)
        x = 100.0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=19, doSA=True)

        self.assertAlmostEqual(
            result["Q_ext"], 2.089822, msg="incorrect Q_ext", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_abs"], 0.957688, msg="incorrect Q_abs", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_sca"], 1.132134, msg="incorrect Q_sca", delta=0.000001
        )
        self.assertAlmostEqual(
            result["g_sca"], 0.950392, msg="incorrect g_sca", delta=0.000001
        )

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = scat_mat["S_12"]
        S_12 /= S_11

        w_S_11 = np.array(
            [
                0.273645e08,
                0.911574e05,
                0.130172e05,
                0.314304e04,
                0.121824e04,
                0.911319e03,
                0.801673e03,
                0.629347e03,
                0.465786e03,
                0.370932e03,
                0.317391e03,
                0.269848e03,
                0.230451e03,
                0.202758e03,
                0.180401e03,
                0.162813e03,
                0.149203e03,
                0.138475e03,
                0.130113e03,
                0.123599e03,
                0.118559e03,
                0.114659e03,
                0.111687e03,
                0.109437e03,
                0.107749e03,
                0.106503e03,
                0.105600e03,
                0.104961e03,
                0.104521e03,
                0.104229e03,
                0.104044e03,
                0.103935e03,
                0.103877e03,
                0.103850e03,
                0.103840e03,
                0.103837e03,
                0.103837e03,
            ]
        )
        w_S_12 = np.array(
            [
                0.0000,
                -0.0247,
                -0.0848,
                -0.2375,
                -0.4798,
                -0.5484,
                -0.5537,
                -0.6268,
                -0.7490,
                -0.8418,
                -0.8905,
                -0.9395,
                -0.9797,
                -0.9960,
                -0.9944,
                -0.9750,
                -0.9395,
                -0.8885,
                -0.8273,
                -0.7576,
                -0.6831,
                -0.6072,
                -0.5316,
                -0.4586,
                -0.3897,
                -0.3257,
                -0.2673,
                -0.2148,
                -0.1684,
                -0.1278,
                -0.0932,
                -0.0643,
                -0.0409,
                -0.0229,
                -0.0101,
                -0.0025,
                0.0000,
            ]
        )

        for i, j in zip(S_11, w_S_11):
            self.assertAlmostEqual(
                i / j, 1.0, msg=f"incorrect S_11, {i} != {j}", delta=0.00001
            )
        for i, j in zip(S_12, w_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.0001)

    def test_w_7(self):
        ri = complex(1.5, 0.1)
        x = 1000.0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=19, doSA=True)

        self.assertAlmostEqual(
            result["Q_ext"], 2.019703, msg="incorrect Q_ext", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_abs"], 0.912770, msg="incorrect Q_abs", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_sca"], 1.106932, msg="incorrect Q_sca", delta=0.000001
        )
        self.assertAlmostEqual(
            result["g_sca"], 0.950880, msg="incorrect g_sca", delta=0.000001
        )

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = scat_mat["S_12"]
        S_12 /= S_11

        w_S_11 = np.array(
            [
                0.255002e12,
                0.103489e07,
                0.270449e06,
                0.145265e06,
                0.101955e06,
                0.803929e05,
                0.646676e05,
                0.528211e05,
                0.436516e05,
                0.364909e05,
                0.308618e05,
                0.264252e05,
                0.229245e05,
                0.201609e05,
                0.179799e05,
                0.162601e05,
                0.149061e05,
                0.138426e05,
                0.130097e05,
                0.123601e05,
                0.118559e05,
                0.114671e05,
                0.111696e05,
                0.109441e05,
                0.107752e05,
                0.106505e05,
                0.105601e05,
                0.104961e05,
                0.104520e05,
                0.104227e05,
                0.104042e05,
                0.103933e05,
                0.103874e05,
                0.103846e05,
                0.103836e05,
                0.103834e05,
                0.103834e05,
            ]
        )
        w_S_12 = np.array(
            [
                0.0000,
                -0.0682,
                -0.1681,
                -0.2830,
                -0.3931,
                -0.4861,
                -0.5789,
                -0.6682,
                -0.7517,
                -0.8267,
                -0.8909,
                -0.9417,
                -0.9765,
                -0.9939,
                -0.9928,
                -0.9737,
                -0.9379,
                -0.8878,
                -0.8265,
                -0.7571,
                -0.6829,
                -0.6069,
                -0.5315,
                -0.4586,
                -0.3897,
                -0.3258,
                -0.2674,
                -0.2149,
                -0.1684,
                -0.1279,
                -0.0932,
                -0.0643,
                -0.0409,
                -0.0229,
                -0.0101,
                -0.0025,
                0.0000,
            ]
        )

        for i, j in zip(S_11, w_S_11):
            self.assertAlmostEqual(
                i / j, 1.0, msg=f"incorrect S_11, {i} != {j}", delta=0.00001
            )
        for i, j in zip(S_12, w_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.0001)

    def test_w_8(self):
        ri = complex(1.5, 0.1)
        x = 5000.0

        result = miex.get_mie_coefficients(x=x, ri=ri, nang=19, doSA=True)

        self.assertAlmostEqual(
            result["Q_ext"], 2.006775, msg="incorrect Q_ext", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_abs"], 0.907582, msg="incorrect Q_abs", delta=0.000001
        )
        self.assertAlmostEqual(
            result["Q_sca"], 1.099193, msg="incorrect Q_sca", delta=0.000001
        )
        self.assertAlmostEqual(
            result["g_sca"], 0.950650, msg="incorrect g_sca", delta=0.000001
        )

        scat_mat = miex.get_scattering_matrix_elements(result["SA_1"], result["SA_2"])
        S_11 = scat_mat["S_11"]
        S_12 = scat_mat["S_12"]
        S_12 /= S_11

        w_S_11 = np.array(
            [
                0.157315e15,
                0.772728e07,
                0.417917e07,
                0.311291e07,
                0.245545e07,
                0.197572e07,
                0.160555e07,
                0.131668e07,
                0.109011e07,
                0.911909e06,
                0.771474e06,
                0.660667e06,
                0.573177e06,
                0.504089e06,
                0.449555e06,
                0.406550e06,
                0.372691e06,
                0.346095e06,
                0.325266e06,
                0.309020e06,
                0.296412e06,
                0.286689e06,
                0.279248e06,
                0.273608e06,
                0.269384e06,
                0.266266e06,
                0.264005e06,
                0.262403e06,
                0.261300e06,
                0.260568e06,
                0.260105e06,
                0.259832e06,
                0.259684e06,
                0.259616e06,
                0.259591e06,
                0.259585e06,
                0.259585e06,
            ]
        )
        w_S_12 = np.array(
            [
                0.0000,
                -0.1103,
                -0.1975,
                -0.2927,
                -0.3902,
                -0.4856,
                -0.5788,
                -0.6680,
                -0.7513,
                -0.8264,
                -0.8906,
                -0.9414,
                -0.9763,
                -0.9936,
                -0.9926,
                -0.9735,
                -0.9378,
                -0.8878,
                -0.8264,
                -0.7570,
                -0.6829,
                -0.6069,
                -0.5315,
                -0.4586,
                -0.3897,
                -0.3258,
                -0.2674,
                -0.2149,
                -0.1684,
                -0.1279,
                -0.0932,
                -0.0643,
                -0.0409,
                -0.0229,
                -0.0101,
                -0.0025,
                0.0000,
            ]
        )

        for i, j in zip(S_11, w_S_11):
            self.assertAlmostEqual(
                i / j, 1.0, msg=f"incorrect S_11, {i} != {j}", delta=0.00001
            )
        for i, j in zip(S_12, w_S_12):
            self.assertAlmostEqual(i, j, msg="incorrect S_12", delta=0.0001)


if __name__ == "__main__":
    unittest.main(verbosity=2)
