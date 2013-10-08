#!/usr/bin/env python

import unittest
import spd


#class ToRomanBadInput(unittest.TestCase):                            
#    def testTooLarge(self):                                          
#        """toRoman should fail with large input"""                   
#        self.assertRaises(roman.OutOfRangeError, roman.toRoman, 4000)


print spd.thing.s

class TestSpecimenValues(unittest.TestCase):

    def test_name_length(self):
        self.assertEqual(len(spd.thing.s), 12) 
        
    def test_thing_is_pie(self):
        self.assertRaises(ValueError, spd.YorkRegression, spd.thing)



if __name__ == "__main__":
    unittest.main()
