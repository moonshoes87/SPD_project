#!/usr/bin/env python


# test_fake.py: -*- Python -*-  DESCRIPTIVE TEXT.
# 
#  Copyright (c) 2014 Lori Jonestrask
#  Author: Lori Jonestrask (mintblue87@gmail.com) .

import os
import unittest

print os.getcwd()

class DoThings(unittest.TestCase):
    
    def test_for_true(self):
        self.assertEqual(True, True)

    def test_for_false(self):
        self.assertFalse(False)
