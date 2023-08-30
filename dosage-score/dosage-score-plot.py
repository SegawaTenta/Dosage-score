#!/usr/bin/env python3

import os
import shutil
import sys
import argparse
from utils import time_stamp, clean_cmd, call_log
from dosage_plot import Dosage_plot

class Dosage_score_plot(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description='dosage-score pipeline 2023/8/21')
        
        parser.add_argument('-r1','--ref_fa1',
                            required=True,
                            type=str,
                            help='Referance fasta file.')

        parser.add_argument('-b','--bam_info',
                            required=True,
                            type=str,
                            help='BAM info file.')

        parser.add_argument('-o','--output_dir',
                            required=False,
                            type=str,
                            default="dosage_score",
                            help='Output directory. [dosage_score]')

        parser.add_argument('-chr_position','--chr_position',
                            required=False,
                            type=str,
                            default="",
                            help='Selecting genomic region.')

        parser.add_argument('-window_size','--window_size',
                            required=False,
                            type=int,
                            default=2000000,
                            help='Minumum plot in window size. [2000000]')
        
        parser.add_argument('-c','--color_info',
                            required=False,
                            type=str,
                            default="",
                            help='color info file.')
        
        parser.add_argument('-min_rate','--min_rate',
                            required=False,
                            type=float,
                            default=0.05,
                            help='Minimum ratio window.')

        args = parser.parse_args()
        self.ref_fa1=args.ref_fa1
        self.bam_info=args.bam_info
        self.output_dir=args.output_dir
        self.chr_position=args.chr_position
        self.window_size=args.window_size
        self.color_info=args.color_info
        self.min_rate=args.min_rate

    def step1(self):
        dosage_pot = Dosage_plot(self.output_dir,
                                self.ref_fa1,
                                self.bam_info,
                                self.chr_position,
                                self.window_size,
                                self.color_info,
                                self.min_rate
                                )
        dosage_pot.run()

    def run(self):
        os.makedirs('{}'.format(self.output_dir), exist_ok=True)
        os.makedirs('{}/4_plot'.format(self.output_dir), exist_ok=True)
        self.step1()

if __name__ == '__main__':
    print(time_stamp(), 'start to run dosage-score plot.', flush=True)
    Dosage_score_plot().run()
    print(time_stamp(), 'dosage-score plot successfully finished.\n', flush=True)

