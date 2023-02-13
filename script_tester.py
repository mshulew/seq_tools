#!/usr/bin/env python3
# coding: utf-8

"""
opens a file and allows testing for script development
"""

import sys

        
if __name__ == "__main__":

# parse arguments
    input_filename = sys.argv[1]
    
    with open(input_filename, 'r') as input_file:
        for line in input_file:
            print(line.rstrip())
        

  
    
