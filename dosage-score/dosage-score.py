#!/usr/bin/env python3

import os
import shutil
import sys
import argparse
from utils import time_stamp, clean_cmd, call_log
from linked_regions import Linked_regions
from sliding_window_calculation import Sliding_window_calculation
from dosage_score_calculation import Dosage_score_calculation


class Dosage_score(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description='dosage-score pipeline 2023/8/21')
        
        parser.add_argument('-r1','--ref_fa1',
                            required=True,
                            type=str,
                            help='Referance fasta file.')
        
        parser.add_argument('-r2','--ref_fa2',
                            required=False,
                            type=str,
                            default="",
                            help='Masked referance fasta file.')
        
        parser.add_argument('-l','--link_file',
                            required=False,
                            type=str,
                            help='Link file compared between 2 genome fasta.')

        parser.add_argument('-g','--genome_info',
                            required=True,
                            type=str,
                            help='Genome info file.')

        parser.add_argument('-b','--bam_info',
                            required=True,
                            type=str,
                            help='BAM info file.')

        parser.add_argument('-o','--output_dir',
                            required=False,
                            type=str,
                            default="dosage_score",
                            help='Output directory. [dosage_score]')

        parser.add_argument('-t','--thread',
                            required=False,
                            type=int,
                            default=1,
                            help='Thread. [1]')

        parser.add_argument('-minMQ','--minMQ',
                            required=False,
                            type=int,
                            default=60,
                            help='Minimum mapping quality for an alignment to be used. [60]')


        parser.add_argument('-minBQ','--minBQ',
                            required=False,
                            type=int,
                            default=13,
                            help='Minimum base quality for a base to be considered. [13]')

        """
        parser.add_argument('-chr_position','--chr_position',
                            required=False,
                            type=str,
                            default="",
                            help='Selecting genomic region.')
        """

        parser.add_argument('-window_size','--window_size',
                            required=False,
                            type=int,
                            default=2000000,
                            help='Minumum plot in window size. [2000000]')
        
        parser.add_argument('-step_size','--step_size',
                            required=False,
                            type=int,
                            default=500000,
                            help='Minumum plot in window size. [500000]')

        args = parser.parse_args()
        self.ref_fa1=args.ref_fa1
        self.ref_fa2=args.ref_fa2
        self.genome_info=args.genome_info
        self.link_file=args.link_file
        self.bam_info=args.bam_info
        self.output_dir=args.output_dir
        self.thread=args.thread
        self.minMQ=args.minMQ
        self.minBQ=args.minBQ
        # self.chr_position=args.chr_position
        self.chr_position=""
        self.window_size=args.window_size
        self.step_size=args.step_size

        if not self.ref_fa2:
            self.ref_fa2=self.ref_fa1

    def step1(self):
        # link ファイルの重複削除とフォーマット変換
        os.makedirs('{}/1_linked_regions'.format(self.output_dir), exist_ok=True)
        linked_regions = Linked_regions(self.output_dir,
                                        self.ref_fa1,
                                        self.link_file,
                                        self.chr_position)
        linked_regions.run()

    def step2(self):
        # mpileupでdepthを算出してwindowごとに計算
        os.makedirs('{}/2_sliding_window_average'.format(self.output_dir), exist_ok=True)
        sliding_window_calculation = Sliding_window_calculation(self.output_dir,
                                                                self.ref_fa1,
                                                                self.ref_fa2,
                                                                self.genome_info,
                                                                self.bam_info,
                                                                self.thread,
                                                                self.minMQ,
                                                                self.minBQ,
                                                                self.window_size,
                                                                self.step_size
                                                                )
        sliding_window_calculation.run()

    def step3(self):
        # mpileupでdepthを算出してwindowごとに計算
        os.makedirs('{}/3_dosage_score'.format(self.output_dir), exist_ok=True)
        dosage_score_calculation = Dosage_score_calculation(self.output_dir,
                                                            self.ref_fa1,
                                                            self.genome_info,
                                                            self.bam_info
                                                            )
        dosage_score_calculation.run()

    def run(self):
        os.makedirs('{}'.format(self.output_dir), exist_ok=True)
        os.makedirs('{}/log'.format(self.output_dir), exist_ok=True)
        log_txt = open("{0}/log/info.txt".format(self.output_dir), "w")
        log_txt.write("ref_fa1\t{0}\n".format(self.ref_fa1))
        log_txt.write("ref_fa2\t{0}\n".format(self.ref_fa2))
        log_txt.write("genome_info\t{0}\n".format(self.genome_info))
        log_txt.write("link_file\t{0}\n".format(self.link_file))
        log_txt.write("bam_info\t{0}\n".format(self.bam_info))
        log_txt.write("output_dir\t{0}\n".format(self.output_dir))
        log_txt.write("thread\t{0}\n".format(self.thread))
        log_txt.write("minMQ\t{0}\n".format(self.minMQ))
        log_txt.write("minBQ\t{0}\n".format(self.minBQ))
        log_txt.write("chr_position\t{0}\n".format(self.chr_position))
        log_txt.write("window_size\t{0}\n".format(self.window_size))
        log_txt.write("step_size\t{0}\n".format(self.step_size))

        self.step1()
        self.step2()
        self.step3()

if __name__ == '__main__':
    print(time_stamp(), 'start to run dosage-score.', flush=True)
    Dosage_score().run()
    print(time_stamp(), 'dosage-score successfully finished.\n', flush=True)

