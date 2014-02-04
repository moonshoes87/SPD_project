#!/usr/bin/env python
# test_test.py: -*- Python -*-  DESCRIPTIVE TEXT.
# 
#  Copyright (c) 2014 Lori Jonestrask
#  Author: Lori Jonestrask (mintblue87@gmail.com) .

print __name__
import unittest
from ..lib import lib_arai_plot_statistics


class FakeTest(unittest.TestCase):

    def test_me(self):
        self.assertEqual(1,1)



if __name__ == "__main__":
    unittest.main()
