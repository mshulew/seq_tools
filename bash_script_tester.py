#!/usr/bin/env python3
# coding: utf-8

"""
opens a file and allows testing for script development
"""

import subprocess

        
if __name__ == "__main__":

        output_working_dir = '/home/ubuntu/data/test_output/complete_rnaseq_output_2/Lib9_r1'

        cmd = 'samtools view {} -F 272 > {}'.format(output_working_dir + '/Aligned.sortedByCoord.out.bam',output_working_dir + '/temp.sam')
        process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        process.wait()
        result = subprocess.run(['wc','-l',output_working_dir + '/temp.sam'], stdout=subprocess.PIPE)
        print(result)
        report_data = int(str(result.stdout).split()[0].strip("b'"))
        print(report_data)
        cmd = 'rm {}'.format(output_working_dir + '/temp.sam')
        process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        process.wait()
  
    
