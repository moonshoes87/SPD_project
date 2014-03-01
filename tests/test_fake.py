#!/usr/bin/env python


# test_fake.py: -*- Python -*-  DESCRIPTIVE TEXT.
# 
#  Copyright (c) 2014 Lori Jonestrask
#  Author: Lori Jonestrask (mintblue87@gmail.com) .

import os
import unittest

print "monkeying around"
print os.getcwd()
print "end of monkey time"

class DoThings(unittest.TestCase):
    
    def test_for_true(self):
        self.assertEqual(True, True)

    def test_for_false(self):
        self.assertFalse(False)
