import sys
import os
import re
import subprocess as sbp
from dosage_score.utils import time_stamp, clean_cmd, call_log

class Linked_regions(object):

    def __init__(self,output_dir,ref_fa,link_file,chr_position):
        self.output_dir = output_dir
        self.ref_fa = ref_fa
        self.link_file = link_file
        self.chr_position = chr_position

    def make_fai(self):
        cmd = 'samtools faidx {} \
            >> {}/log/samtools.log \
            2>&1'.format(self.ref_fa, self.output_dir)
        
        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.output_dir, 'samtools', cmd)
            sys.exit(1)

    def run(self):

        print(time_stamp(),
              'start to linked region selected.',
              flush=True)

        if not os.path.exists(self.ref_fa + '.fai'):
            self.make_fai()

        if not self.link_file:
            fai = open(self.ref_fa + '.fai', "r")
            output_txt = open("{0}/1_linked_regions/analysis_region.txt".format(self.output_dir), "w")
            for fai_line in fai :
                fai_line=fai_line.replace('\n', '')
                colom = fai_line.split()
                chr = str(colom[0])
                end_posi = str(colom[1])

                if self.chr_position:
                    colom_c = self.chr_position.split(":")
                    target_chr=colom_c[0]
                    target_s_posi=colom_c[1]
                    target_e_posi=colom_c[2]
                    if target_chr==chr:
                        if target_s_posi<1:
                            target_s_posi=1
                        if target_e_posi>end_posi:
                            target_e_posi=end_posi
                        output_txt.write("{0}\t{1}\t{2}\n".format(target_chr,target_s_posi,target_e_posi))
                else:
                    output_txt.write("{0}\t1\t{1}\n".format(chr,end_posi))

        else:

            hash_delite={}
            hashA={}
            hashB={}
            line_count=0

            txt = open(self.link_file, "r")
            output_txt = open("{0}/1_linked_regions/analysis_region.txt".format(self.output_dir), "w")
            for txt_line in txt :
                txt_line=txt_line.replace('\n', '')
                line_count=line_count+1
                colom = txt_line.split()
                chrA = str(colom[0])
                s_posiA = int(colom[1])
                e_posiA = int(colom[2])
                chrB = str(colom[3])
                s_posiB = int(colom[4])
                e_posiB = int(colom[5])

                if e_posiA < s_posiA:
                    s_posiA, e_posiA = e_posiA, s_posiA

                if e_posiB < s_posiB:
                    s_posiB, e_posiB = e_posiB, s_posiB

                for i in range(1,line_count):
                    colomA = hashA[i].split(":")
                    keep_chrA=colomA[0]
                    keep_s_posiA=int(colomA[1])
                    keep_e_posiA=int(colomA[2])

                    if keep_chrA == chrA:
                        if keep_e_posiA>=s_posiA and keep_s_posiA<=e_posiA:
                            hash_delite[i]=1
                            hash_delite[line_count]=1

                    colomB = hashB[i].split(":")
                    keep_chrB=colomB[0]
                    keep_s_posiB=int(colomB[1])
                    keep_e_posiB=int(colomB[2])

                    if keep_chrB == chrB:
                        if keep_e_posiB>=s_posiB and keep_s_posiB<=e_posiB:
                            hash_delite[i]=1
                            hash_delite[line_count]=1

                hashA[line_count]=chrA + ":" + str(s_posiA) + ":" + str(e_posiA)
                hashB[line_count]=chrB + ":" + str(s_posiB) + ":" + str(e_posiB)

            txt.close()

            line_count=0
            txt = open(self.link_file, "r")
            for txt_line in txt :
                txt_line=txt_line.replace('\n', '')
                line_count=line_count+1
                if not line_count in hash_delite:
                    colom = txt_line.split()
                    chrA = str(colom[0])
                    s_posiA = int(colom[1])
                    e_posiA = int(colom[2])
                    chrB = str(colom[3])
                    s_posiB = int(colom[4])
                    e_posiB = int(colom[5])

                    if e_posiA < s_posiA:
                        s_posiA, e_posiA = e_posiA, s_posiA

                    if e_posiB < s_posiB:
                        s_posiB, e_posiB = e_posiB, s_posiB

                    if self.chr_position:
                        colom_c = self.chr_position.split(":")
                        target_chr=colom_c[0]
                        target_s_posi=int(colom_c[1])
                        target_e_posi=int(colom_c[2])

                        if target_chr==chrA:
                            if target_s_posi>e_posiA or target_e_posi>s_posiA:
                                pass
                            else:
                                if target_s_posi>s_posiA:
                                    s_posiA=target_s_posi
                                if target_e_posi<e_posiA:
                                    e_posiA=target_e_posi
                                
                                output_txt.write("{0}\t{1}\t{2}\n".format(target_chr,str(s_posiA),str(e_posiA)))

                        if target_chr==chrB:
                            if target_s_posi>e_posiB or target_e_posi>s_posiB:
                                pass
                            else:
                                if target_s_posi>s_posiB:
                                    s_posiB=target_s_posi
                                if target_e_posi<e_posiB:
                                    e_posiB=target_e_posi
                                
                                output_txt.write("{0}\t{1}\t{2}\n".format(target_chr,str(s_posiB),str(e_posiB)))

                    else:
                        output_txt.write("{0}\t{1}\t{2}\n{3}\t{4}\t{5}\n".format(chrA,str(s_posiA),str(e_posiA),chrB,str(s_posiB),str(e_posiB)))

        print(time_stamp(),
              'linked region selected successfully finished.',
              flush=True)