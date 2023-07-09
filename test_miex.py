import miex
import numpy as np
import unittest

class test_miex(unittest.TestCase):
    '''
        compare results by
        C.F.Bohren & D.R.Huffman 'Absorption and Scattering of Light by Small Particles' (1998)
        https://doi.org/10.1002/9783527618156
        (see Appendix A, page 482)
        and by
        W.J.Wiscombe: 'Mie Scattering Calculations: Advances in Technique and Fast, Vector-speed Computer Codes' (1979)
        https://doi.org/10.5065/D6ZP4414
        (see Appendix, A22)
    '''

    def test_bh(self):
        ri = complex(1.55, 0)
        lambda0 = 0.6328
        radius = 0.525
        x = 2.0 * np.pi * radius / lambda0

        result = miex.shexqnn2(x=x, ri=ri, nang=11, doSA=True)

        self.assertAlmostEqual(result[0], 3.10543, msg='incorrect Q_ext', delta=0.00001)
        self.assertAlmostEqual(result[1], 0.00000, msg='incorrect Q_abs', delta=0.00001)
        self.assertAlmostEqual(result[2], 3.10543, msg='incorrect Q_sca', delta=0.00001)
        self.assertAlmostEqual(result[3], 2.92534, msg='incorrect Q_bk', delta=0.00001)

        S_11, S_12, S_33, S_34 = miex.scattering_matrix_elements(result[7], result[8])
        S_12 /= -S_11
        S_33 /= S_11
        S_34 /= S_11
        S_11 /= S_11[0]
        
        bh_S_11 = np.array([
            0.100000e+01, 0.785390e+00, 0.356897e+00, 0.766119e-01, 0.355355e-01, 0.701845e-01, 0.574313e-01, 0.219660e-01, 0.125959e-01,
            0.173750e-01, 0.124601e-01, 0.679093e-02, 0.954239e-02, 0.863419e-02, 0.227421e-02, 0.543998e-02, 0.160243e-01, 0.188852e-01,
            0.195254e-01, 0.301676e-01, 0.383189e-01])
        bh_S_12 = np.array([
            0.000000e+00,-0.459811e-02,-0.458541e-01,-0.364744e+00,-0.534997e+00, 0.959953e-02, 0.477927e-01,-0.440604e+00,-0.831996e+00,
            0.341670e-01, 0.230462e+00,-0.713472e+00,-0.756255e+00,-0.281215e+00,-0.239612e+00,-0.850804e+00,-0.706334e+00,-0.891081e+00,
           -0.783319e+00,-0.196194e+00, 0.000000e+00])
        bh_S_33 = np.array([
            0.100000e+01, 0.999400e+00, 0.986022e+00, 0.843603e+00, 0.686967e+00, 0.959825e+00, 0.985371e+00, 0.648043e+00, 0.203255e+00,
            0.795354e+00, 0.937497e+00,-0.717397e-02,-0.394748e-01, 0.536251e+00, 0.967602e+00, 0.187531e+00, 0.495254e+00, 0.453277e+00,
           -0.391613e+00,-0.962069e+00,-0.100000e+01])
        bh_S_34 = np.array([
            0.000000e+00, 0.343261e-01, 0.160184e+00, 0.394076e+00,-0.491787e+00,-0.280434e+00, 0.163584e+00, 0.621216e+00,-0.516208e+00,
           -0.605182e+00, 0.260742e+00, 0.700647e+00,-0.653085e+00,-0.795835e+00, 0.795798e-01,-0.490882e+00,-0.505781e+00,-0.226817e-01,
            0.482752e+00, 0.189556e+00, 0.000000e+00])
        
        for i, j in zip(S_11, bh_S_11):
            self.assertAlmostEqual(i, j, msg='incorrect S_11', delta=0.000001)
        for i, j in zip(S_12, bh_S_12):
            self.assertAlmostEqual(i, j, msg='incorrect S_12', delta=0.000001)
        for i, j in zip(S_33, bh_S_33):
            self.assertAlmostEqual(i, j, msg='incorrect S_33', delta=0.000001)
        for i, j in zip(S_34, bh_S_34):
            self.assertAlmostEqual(i, j, msg='incorrect S_34', delta=0.000001)

    def test_w_1(self):
        ri = complex(1.5, 0)
        x = 10.0

        result = miex.shexqnn2(x=x, ri=ri)

        self.assertAlmostEqual(result[0], 2.881999, msg='incorrect Q_ext', delta=0.000001)
        self.assertAlmostEqual(result[1], 0.000000, msg='incorrect Q_abs', delta=0.000001)
        self.assertAlmostEqual(result[2], 2.881999, msg='incorrect Q_sca', delta=0.000001)
        self.assertAlmostEqual(result[6], 0.742913, msg='incorrect g_sca', delta=0.000001)

        # TODO: compare intensity and degree of polarization

    def test_w_2(self):
        ri = complex(1.5, 0)
        x = 100.0

        result = miex.shexqnn2(x=x, ri=ri)

        self.assertAlmostEqual(result[0], 2.094388, msg='incorrect Q_ext', delta=0.000001)
        self.assertAlmostEqual(result[1], 0.000000, msg='incorrect Q_abs', delta=0.000001)
        self.assertAlmostEqual(result[2], 2.094388, msg='incorrect Q_sca', delta=0.000001)
        self.assertAlmostEqual(result[6], 0.818246, msg='incorrect g_sca', delta=0.000001)

        # TODO: compare intensity and degree of polarization

    def test_w_3(self):
        ri = complex(1.5, 0)
        x = 1000.0

        result = miex.shexqnn2(x=x, ri=ri)

        self.assertAlmostEqual(result[0], 2.013945, msg='incorrect Q_ext', delta=0.000001)
        self.assertAlmostEqual(result[1], 0.000000, msg='incorrect Q_abs', delta=0.000001)
        self.assertAlmostEqual(result[2], 2.013945, msg='incorrect Q_sca', delta=0.000001)
        self.assertAlmostEqual(result[6], 0.827882, msg='incorrect g_sca', delta=0.000001)

        # TODO: compare intensity and degree of polarization

    def test_w_4(self):
        ri = complex(1.5, 0)
        x = 5000.0

        result = miex.shexqnn2(x=x, ri=ri)

        self.assertAlmostEqual(result[0], 2.008650, msg='incorrect Q_ext', delta=0.000001)
        self.assertAlmostEqual(result[1], 0.000000, msg='incorrect Q_abs', delta=0.000001)
        self.assertAlmostEqual(result[2], 2.008650, msg='incorrect Q_sca', delta=0.000001)
        self.assertAlmostEqual(result[6], 0.829592, msg='incorrect g_sca', delta=0.000001)

        # TODO: compare intensity and degree of polarization
    
    def test_w_5(self):
        ri = complex(1.5, 0.1)
        x = 10.0

        result = miex.shexqnn2(x=x, ri=ri)

        self.assertAlmostEqual(result[0], 2.459791, msg='incorrect Q_ext', delta=0.000001)
        self.assertAlmostEqual(result[1], 1.224646, msg='incorrect Q_abs', delta=0.000001)
        self.assertAlmostEqual(result[2], 1.235144, msg='incorrect Q_sca', delta=0.000001)
        self.assertAlmostEqual(result[6], 0.922350, msg='incorrect g_sca', delta=0.000001)

        # TODO: compare intensity and degree of polarization

    def test_w_6(self):
        ri = complex(1.5, 0.1)
        x = 100.0

        result = miex.shexqnn2(x=x, ri=ri)

        self.assertAlmostEqual(result[0], 2.089822, msg='incorrect Q_ext', delta=0.000001)
        self.assertAlmostEqual(result[1], 0.957688, msg='incorrect Q_abs', delta=0.000001)
        self.assertAlmostEqual(result[2], 1.132134, msg='incorrect Q_sca', delta=0.000001)
        self.assertAlmostEqual(result[6], 0.950392, msg='incorrect g_sca', delta=0.000001)

        # TODO: compare intensity and degree of polarization

    def test_w_7(self):
        ri = complex(1.5, 0.1)
        x = 1000.0

        result = miex.shexqnn2(x=x, ri=ri)

        self.assertAlmostEqual(result[0], 2.019703, msg='incorrect Q_ext', delta=0.000001)
        self.assertAlmostEqual(result[1], 0.912770, msg='incorrect Q_abs', delta=0.000001)
        self.assertAlmostEqual(result[2], 1.106932, msg='incorrect Q_sca', delta=0.000001)
        self.assertAlmostEqual(result[6], 0.950880, msg='incorrect g_sca', delta=0.000001)

        # TODO: compare intensity and degree of polarization
    
    def test_w_8(self):
        ri = complex(1.5, 0.1)
        x = 5000.0

        result = miex.shexqnn2(x=x, ri=ri)

        self.assertAlmostEqual(result[0], 2.006775, msg='incorrect Q_ext', delta=0.000001)
        self.assertAlmostEqual(result[1], 0.907582, msg='incorrect Q_abs', delta=0.000001)
        self.assertAlmostEqual(result[2], 1.099193, msg='incorrect Q_sca', delta=0.000001)
        self.assertAlmostEqual(result[6], 0.950650, msg='incorrect g_sca', delta=0.000001)

        # TODO: compare intensity and degree of polarization

if __name__ == '__main__':
    unittest.main(verbosity=2)
