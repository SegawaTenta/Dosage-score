import sys
import os
import re
import subprocess as sbp
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from dosage_score.utils import time_stamp, clean_cmd, call_log

class Dosage_score_calculation(object):
    def __init__(self,output_dir,ref_fa1,genome_info,bam_info):
        self.output_dir = output_dir
        self.ref_fa1 = ref_fa1
        self.genome_info = genome_info
        self.bam_info = bam_info

    def run(self):

        print(time_stamp(),
              'start to dosage-score calculation.',
              flush=True)
        

        genome_info_open = open(self.genome_info, "r")
        hash_chr={}
        for genome_info_line in genome_info_open :
            genome_info_line=genome_info_line.replace('\n', '')
            colom = genome_info_line.split()
            mark = str(colom[0])
            chrs = str(colom[1])
            hash_chr[mark]=chrs

        bam_info_open = open(self.bam_info, "r")
        hash_bam={}
        line_count=0
        for bam_info_line in bam_info_open :
            bam_info_line=bam_info_line.replace('\n', '')
            line_count=line_count+1

            colom = bam_info_line.split()
            expected_genome = str(colom[3])

            colom2 = expected_genome.split(",")
            for i in colom2:
                if line_count in hash_bam:
                    hash_bam[line_count]=hash_bam[line_count]+","+hash_chr[i]
                else:
                    hash_bam[line_count]=hash_chr[i]
                


        
        
        
        out_txt = open("{0}/3_dosage_score/dosage_score.txt".format(self.output_dir), "w")

        df = pd.read_csv("{0}/2_sliding_window_average/sliging_window_average.txt".format(self.output_dir), sep='\t', header=None)
        medium_list=[]
        for i in range(1,line_count+1):

            cont_chr = hash_bam[i].split(",")
            filtered_df = df[df[0].isin(cont_chr)]
            median_value = filtered_df[2+i].median()
            medium_list.append(median_value)


        sliding_window_txt_open = open("{0}/2_sliding_window_average/sliging_window_average.txt".format(self.output_dir), "r")

        for txt_line in sliding_window_txt_open :
            txt_line=txt_line.replace('\n', '')
            colom = txt_line.split()
            sample_num=len(colom)-3
            chr = str(colom[0])
            posi = str(colom[1])
            a_posi = str(colom[2])

            out_txt.write("{0}\t{1}\t{2}".format(chr,posi,a_posi))

            for i in range(1,sample_num+1):
                depth=float(colom[2+i])
                dosage_score=(depth/medium_list[i-1])*2
                dosage_score=int(dosage_score*1000)/1000
                out_txt.write("\t{0}".format(dosage_score))

            out_txt.write("\n")
        
        print(time_stamp(),
              'dosage-score calculation successfully finished.',
              flush=True)