import sys
import os
import re
import subprocess as sbp
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from utils import time_stamp, clean_cmd, call_log

class Sliding_window_calculation(object):
    def __init__(self,output_dir,ref_fa1,ref_fa2,genome_info,bam_info,thread,minMQ,minBQ,window_size,step_size):
        self.output_dir = output_dir
        self.ref_fa1 = ref_fa1
        self.ref_fa2 = ref_fa2
        self.genome_info = genome_info
        self.bam_info = bam_info
        self.thread = thread
        self.minMQ = minMQ
        self.minBQ = minBQ
        self.window_size = window_size
        self.step_size = step_size

    def run(self):

        print(time_stamp(),
              'start to sliding window calculation.',
              flush=True)

        os.makedirs('{0}/2_sliding_window_average/region_txt'.format(self.output_dir), exist_ok=True)

        region_txt=[]
        fai = open(self.ref_fa1 + '.fai', "r")
        for fai_line in fai :
            fai_line=fai_line.replace('\n', '')
            colom = fai_line.split()
            chr = str(colom[0])
            e_posi = int(colom[1])
            s_posi = 0
            
            while s_posi <= e_posi + self.step_size:
                # print(chr,s_posi)
                # print(s_posi-(self.window_size)/2,s_posi+(self.window_size)/2)
                txt = open(self.output_dir + '/1_linked_regions/analysis_region.txt', "r")
                s_window=int(s_posi-(self.window_size)/2)
                e_window=int(s_posi+(self.window_size)/2)
                if s_window<1:
                        s_window=1
                output_txt = open("{0}/2_sliding_window_average/region_txt/{1}.{2}.{3}".format(self.output_dir,
                                                                                                    chr,
                                                                                                    s_window,
                                                                                                    e_window
                                                                                                    ), "w")
                out_count=0

                for txt_line in txt :
                    txt_line=txt_line.replace('\n', '')
                    colom2 = txt_line.split()
                    target_chr = str(colom2[0])
                    target_s_posi = int(colom2[1])
                    target_e_posi = int(colom2[2])

                    if target_chr == chr and not (target_e_posi < s_window or target_s_posi > e_window):
                        if target_s_posi < s_window:
                            target_s_posi = s_window

                        if target_e_posi > e_window:
                            target_e_posi = e_window

                        if target_s_posi<1:
                            target_s_posi=1

                        # print(chr,s_posi,target_chr,target_s_posi,target_e_posi)
                        output_txt.write("{0}\t{1}\t{2}\n".format(target_chr,target_s_posi,target_e_posi))
                        out_count=out_count+1

                output_txt.close()
                if out_count==0:
                    os.remove("{0}/2_sliding_window_average/region_txt/{1}.{2}.{3}".format(self.output_dir,
                                                                                            chr,
                                                                                            s_window,
                                                                                            e_window
                                                                                            ))
                else:
                    region_txt.append("{0}/2_sliding_window_average/region_txt/{1}.{2}.{3}".format(self.output_dir,
                                                                                            chr,
                                                                                            s_window,
                                                                                            e_window
                                                                                            ))
                    
                s_posi = s_posi + self.step_size
        







        def mpileup(a_region_txt):
            genome_info_open = open(self.genome_info, "r")
            hash_chr={}
            for genome_info_line in genome_info_open :
                genome_info_line=genome_info_line.replace('\n', '')
                colom = genome_info_line.split()
                mark = str(colom[0])
                chrs = str(colom[1])

                if "," in chrs:
                    colom2 = chrs.split(",")
                    for chr in colom2:
                        hash_chr[chr]=mark
                else:
                    hash_chr[chrs]=mark

            bam_info_open = open(self.bam_info, "r")
            hash_bam={}
            bams=""
            line_count=0
            for bam_info_line in bam_info_open :
                bam_info_line=bam_info_line.replace('\n', '')
                line_count=line_count+1

                colom = bam_info_line.split()
                sample_name = str(colom[0])
                bam_path = str(colom[1])
                stat = str(colom[2])
                expected_genome = str(colom[3])

                bams=bams+" "+bam_path
                if stat == "Control":
                    hash_bam[line_count]=expected_genome

            region_txt_open = open(a_region_txt, "r")
            file_name = os.path.basename(a_region_txt)
            output_txt = open("{0}/2_sliding_window_average/sliging_window_total/{1}".format(self.output_dir,
                                                                                                file_name
                                                                                                ), "w")

            posi_count=0
            sum_count=[]
            for txt_line in region_txt_open :
                txt_line=txt_line.replace('\n', '')
                colom = txt_line.split()
                chr = str(colom[0])
                s_posi = int(colom[1])
                e_posi = int(colom[2])

                cmd="samtools mpileup -q {0} -Q {1} -r {2}:{3}-{4} -f {5} {6}".format(self.minMQ,
                                                                                      self.minBQ,
                                                                                      chr,
                                                                                      s_posi,
                                                                                      e_posi,
                                                                                      self.ref_fa2,
                                                                                      bams)
                
                cmd_processsbp=sbp.Popen(cmd,stdout=sbp.PIPE,stderr=sbp.DEVNULL,shell=True)
                for line in cmd_processsbp.stdout:
                    line=line.decode().strip()
                    array=line.split()
                    chr=array[0]
                    posi=int(array[1])
                    base=array[2]
                    sample_num=int((len(array)-3)/3)
                    
                    depth_list=[]
                    out_witch=0
                    if not base=="N":
                        for i in range(1,sample_num+1):
                            sample_depth=int(array[i*3])
                            if i in hash_bam:
                                if hash_chr[chr]==hash_bam[i]:
                                    if sample_depth==0:
                                        out_witch=out_witch+1
                                else:
                                    if sample_depth>0:
                                        out_witch=out_witch+1
                            
                            depth_list.append(sample_depth)

                        if out_witch==0:
                            posi_count=posi_count+1
                            if len(sum_count) == 0:
                                sum_count=depth_list
                            else:
                                new_list=[]
                                for index, i in enumerate(sum_count):
                                    sum_depth=i+depth_list[index]
                                    new_list.append(sum_depth)
                                sum_count=new_list

            depth="\t".join(map(str, sum_count))
            output_txt.write("{0}\t{1}\n".format(posi_count,depth))

                    
        os.makedirs('{0}/2_sliding_window_average/sliging_window_total'.format(self.output_dir), exist_ok=True)
        
        with ThreadPoolExecutor(max_workers=self.thread) as executor:
            executor.map(mpileup, region_txt) 
        
        







        output_txt = open("{0}/2_sliding_window_average/sliging_window_average.txt".format(self.output_dir), "w")
        
        for i in region_txt:
            file_name = os.path.basename(i)

            parts = file_name.split('.')
            e_posi = int(parts[-1])
            s_posi = int(parts[-2])
            chr = '.'.join(parts[:-2])

            txt_open = open("{0}/2_sliding_window_average/sliging_window_total/{1}".format(self.output_dir,
                                                                                            file_name), "r")
            
            for txt_line in txt_open:
                txt_line=txt_line.replace('\n', '')
                colom = txt_line.split()
                length=len(colom)
                position_num=int(colom[0])

                if not position_num==0:
                    middole_posi=int(e_posi-self.window_size/2)
                    output_txt.write("{0}\t{1}\t{2}".format(chr,middole_posi,position_num))
                    
                    for ii in range(1,length):
                        total_depth=int(colom[ii])
                        ave=int((total_depth/position_num)*1000)/1000
                        output_txt.write("\t{0}".format(ave))

                    output_txt.write("\n")

                

            



        
        print(time_stamp(),
              'sliding window calculation successfully finished.',
              flush=True)