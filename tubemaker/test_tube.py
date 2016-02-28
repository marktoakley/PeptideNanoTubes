'''
@author:  Mark Oakley
'''
import unittest
import numpy as np

from tubemaker.nanotube import orient_coords, build_ring, build_tube

class RandomTest(unittest.TestCase):
    
    def setUp(self):
        self.res_coords = np.array([[  3.32577000e+00,   1.54790900e+00,  -1.60720400e-06],
                                    [  3.90940700e+00,   7.23611000e-01,  -2.73988200e-06],
                                    [  3.97004800e+00,   2.84579500e+00,  -1.31116300e-07],
                                    [  3.67166300e+00,   3.40012900e+00,  -8.89820000e-01],
                                    [  3.57696500e+00,   3.65383800e+00,   1.23214300e+00],
                                    [  3.87748400e+00,   3.11579500e+00,   2.13119700e+00],
                                    [  4.07505900e+00,   4.62301700e+00,   1.20578600e+00],
                                    [  2.49699500e+00,   3.80107500e+00,   1.24137900e+00],
                                    [  5.48554100e+00,   2.70520700e+00,  -4.39875500e-06],
                                    [  6.00882400e+00,   1.59317500e+00,  -8.44976800e-06]])
        self.res_coords = orient_coords(self.res_coords)
        
    def test_ring(self):
        res_coords = self.res_coords.copy()
        coords = build_ring(8,res_coords)
        self.assertEquals(80,len(coords))
        
    def test_tube(self):
        res_coords = self.res_coords.copy()
        coords = build_tube(4, 8, res_coords)
        self.assertEquals(320,len(coords))    
        
if __name__ == "__main__":
    unittest.main()