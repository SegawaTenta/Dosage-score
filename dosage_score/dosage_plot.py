import sys
import os
import re
import subprocess as sbp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from concurrent.futures import ThreadPoolExecutor
from dosage_score.utils import time_stamp, clean_cmd, call_log

class Dosage_plot(object):
    def __init__(self,output_dir,ref_fa1,bam_info,chr_position,window_size,color_info,min_rate):
        self.output_dir=output_dir
        self.ref_fa1=ref_fa1
        self.bam_info=bam_info
        self.chr_position=chr_position
        self.window_size=window_size
        self.color_info=color_info
        self.min_rate=min_rate

    def run(self):

        hash_max_posi={}
        fai_info_open = open("{0}.fai".format(self.ref_fa1), "r")
        for fai_info_line in fai_info_open :
            fai_info_line=fai_info_line.replace('\n', '')
            colom = fai_info_line.split()
            chr = str(colom[0])
            max_posi = int(colom[1])
            hash_max_posi[chr]=max_posi

        score_info_open = open("{0}/3_dosage_score/dosage_score.txt".format(self.output_dir), "r")
        out_txt = open("{0}/4_plot/filtered_dosage_score.txt".format(self.output_dir), "w")

        for score_info_line in score_info_open :
            score_info_line=score_info_line.replace('\n', '')
            colom = score_info_line.split()
            chr = str(colom[0])
            posi = int(colom[1])
            analysis_posi = int(colom[2])

            s_posi=posi-self.window_size/2
            e_posi=posi+self.window_size/2

            if s_posi<1:
                s_posi=1

            if e_posi>hash_max_posi[chr]:
                e_posi=hash_max_posi[chr]

            possible_posi=e_posi-s_posi+1
            if analysis_posi/possible_posi>=self.min_rate:
                out_txt.write("{}\n".format(score_info_line))

        out_txt.close()

        name_list=[]
        bam_info_open = open(self.bam_info, "r")
        for bam_info_line in bam_info_open :
            bam_info_line=bam_info_line.replace('\n', '')
            colom = bam_info_line.split()
            name = str(colom[0])
            name_list.append(name)


        hash_color={}
        if self.color_info:
            color_info_open = open(self.color_info, "r")
            for color_info_line in color_info_open :
                color_info_line=color_info_line.replace('\n', '')
                colom = color_info_line.split()
                chr = str(colom[0])
                color = str(colom[1])
                hash_color[chr]=color


        ######plot#########
        if not self.chr_position:
            sample_num=len(name_list)
            df=pd.read_table("{0}/4_plot/filtered_dosage_score.txt".format(self.output_dir),header=None)

            chr_item = df[0].unique()
            chr_num=len(chr_item)
            fig_col=5
            fig_raw=math.ceil(chr_num/5)
            maxposi=max(df[1])/1000000

            for i in range(sample_num):
                count1=0
                fig = plt.figure(figsize=(20,8),dpi=300)
                plt.subplots_adjust(wspace=0.4, hspace=0.6)

                for chr in chr_item:
                    count1=count1+1
                    ax = fig.add_subplot(fig_raw, fig_col, count1)
                    plt.title(chr,fontsize=9)
                    plt.xlim(0,maxposi)
                    plt.ylim(0,4.5)
                    plt.xticks(np.arange(0, maxposi, 20))
                    plt.xlabel("Chr position(Mb)",fontsize=9)
                    plt.ylabel("Dosage-score",fontsize=9)
                    plt.gca().spines['right'].set_visible(False)
                    plt.gca().spines['top'].set_visible(False)

                    df11=df[df[0]==chr]

                    plt.hlines(1, 0,maxposi, color='black', linestyle='dotted',linewidth=0.5)
                    plt.hlines(2, 0,maxposi, color='black', linestyle='dotted',linewidth=0.5)
                    plt.hlines(3, 0,maxposi, color='black', linestyle='dotted',linewidth=0.5)
                    plt.hlines(4, 0,maxposi, color='black', linestyle='dotted',linewidth=0.5)

                    if chr in hash_color:
                        ax.plot(df11[1]/1000000, df11[i+3], color=hash_color[chr],linewidth=3)
                    else:
                        ax.plot(df11[1]/1000000, df11[i+3], color="#000000",linewidth=3)


                fig.savefig("{0}/4_plot/{1}.png".format(self.output_dir,name_list[i]))
                plt.clf()

        else:
            sample_num=len(name_list)
            colom=self.chr_position.split(":")
            
            df=pd.read_table("{0}/4_plot/filtered_dosage_score.txt".format(self.output_dir),header=None)
            
            df11=df[df[0]==colom[0]]
            filtered_df = df11[(df11[1] >= int(colom[1])) & (df11[1] <= int(colom[2]))]

            for i in range(sample_num):
                fig = plt.figure(figsize=(8,4),dpi=300)
                plt.subplots_adjust(wspace=0.4, hspace=0.6)
                plt.title(colom[0],fontsize=9)
                plt.xlim(int(colom[1])/1000000,int(colom[2])/1000000)
                plt.ylim(0,4.5)
                plt.yticks([0, 1, 2, 3, 4])
                plt.xlabel("Chr position (Mb)",fontsize=9)
                plt.ylabel("Dosage-score",fontsize=9)
                plt.gca().spines['right'].set_visible(False)
                plt.gca().spines['top'].set_visible(False)
                plt.hlines(1, 0,int(colom[2])/1000000, color='black', linestyle='dotted',linewidth=0.5)
                plt.hlines(2, 0,int(colom[2])/1000000, color='black', linestyle='dotted',linewidth=0.5)
                plt.hlines(3, 0,int(colom[2])/1000000, color='black', linestyle='dotted',linewidth=0.5)
                plt.hlines(4, 0,int(colom[2])/1000000, color='black', linestyle='dotted',linewidth=0.5)

                if chr in hash_color:
                    plt.plot(filtered_df[1]/1000000, filtered_df[i+3], color=hash_color[colom[0]],linewidth=3)
                else:
                    plt.plot(filtered_df[1]/1000000, filtered_df[i+3], color="#000000",linewidth=3)

                fig.savefig("{0}/4_plot/{1}_{2}_{3}_{4}.png".format(self.output_dir,name_list[i],str(colom[0]),str(colom[1]),str(colom[2])))
                plt.clf()


            








