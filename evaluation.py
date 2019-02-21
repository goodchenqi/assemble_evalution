# -*- coding:UTF-8 -*-
#!/thinker/storage/software/Genome/anaconda2/bin/python
from __future__ import unicode_literals
# -*- coding: utf-8 -*-

############packages usage#############
import argparse
import re
import os
import gzip
import time
import math
import sys 
import logging
import glob
import random
from collections import defaultdict,OrderedDict 
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont
except:
    print('bad')
import pysam
from config import *
#####Description####
usage = '''
Author : chenqi
Email  : chenqi@gooalgene.com
Date   : 2018-9-3
Version: v3.1
Description:
    
Example: 
    python %s 

''' % (__file__[__file__.rfind(os.sep) + 1:])

#####HelpFormat#####
class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def run_cmd(cmd):
    #may be here a problem but it is ok ,youjust need to change it 
    logging.info(cmd)
    flag = os.system(cmd)
    if flag == 0 or flag == 256:
        return 0
    sys.exit(0)

def contiue_run(tag,workpath):
    trace = workpath + '/trace_tmp/trace.txt'
    status = []
    key = []
    code = []
    filePath = []
    exit_flag = False
    if os.path.exists(trace):
        with open(trace) as f:
            f.readline()
            for line in f:
                key.append(line.split('\t')[3])
                status.append(line.split('\t')[4])
                code.append(line.split('\t')[5])
                filePath.append(line.split('\t')[1])
        count = 0
        for i in key:
            if tag in i and code[count] == '0' and status[count] == 'COMPLETED':
                exit_flag = True
                tmp_path = workpath+'/work/%s*/*'%filePath[count]
                tmp_path_list = glob.glob(r'%s'%tmp_path)
                for file in tmp_path_list:
                    os.system('ln -s %s'%file)
                break
            count += 1
    if exit_flag == False:
        return 1
    elif exit_flag == True:
        return 0

def fa_or_fq(file):
    for i in os.popen('less -S %s |head -1'%file):
        if i.startswith('>'):
            return 'fa'
        elif i.startswith('@'):
            return 'fq'
        else:
            print('please check your input file:%s\n\n\t you must input fastq or fasta format\n')
            sys.exit(0)

def stat_N50(file):
    #base_info
    total_contig,total_scaffold,base_info,sequence=[],[],[0,0,0,0,0],[]
    for line in open(file):
        if line.startswith('>'):
            if sequence:
                count=0
                scaffold=''.join(sequence)
                for i in ['A','T','C','G','N']:
                    base_info[count]+=scaffold.count(i)
                    count+=1
                total_scaffold.append(len(scaffold))
                for contig in re.split(r'[n|N]+',scaffold):
                    total_contig.append(len(contig))
            sequence=[]
            continue
        sequence.append(line.strip().upper())
    scaffold=''.join(sequence)
    count=0
    for i in ['A','T','C','G','N']:
        base_info[count]+=scaffold.count(i)
        count+=1
    total_scaffold.append(len(scaffold))
    for contig in re.split(r'[n|N]+',scaffold):
        total_contig.append(len(contig))
    sequence=[]
    #total_contig:[13201252, 791743, 316486, 34071, 167849, 7000, 28949, 136, 7265, 21549, 37522]
    #total_scaffold:[30427671, 19698289, 23459830, 18585056, 26975502, 366924, 154478]
    total_contig=sorted(total_contig,reverse=True)
    total_scaffold=sorted(total_scaffold,reverse=True)
    #contig_max_length,contig_total_length,scaffold_max_length,scaffold_total_length
    contig_max_length,contig_total_length,scaffold_max_length,scaffold_total_length=total_contig[0],sum(total_contig),total_scaffold[0],sum(total_scaffold)
    total_base,total_base_contig=sum(base_info),sum(base_info)-base_info[-1]
    base_info.append(base_info[2]+base_info[3])
    base_info.append(total_base)
    N_info_scaffold=[total_base*0.9,total_base*0.8,total_base*0.7,total_base*0.6,total_base*0.5]
    N_info_contig=[total_base_contig*i for i in [0.9,0.8,0.7,0.6,0.5]]
    count,check,contig_sum,scaffold_sum=0,1,0,0
    num_100_contig,num_2k_contig,num_100_scaffold,num_2k_scaffold=0,0,0,0
    contig_N_seq,scaffold_N_seq,contig_i,scaffold_i=[0,0,0,0,0],[0,0,0,0,0],1,1
    while 1:
        if count>=len(total_contig):
            break
        elif count>=len(total_scaffold) and check:
            check=0
        else:
            #count_contig:num of total_contig,13201252, 791743, 316486, 34071
            count_contig=total_contig[count]
            contig_sum+=total_contig[count]
            if count_contig>=100:num_100_contig+=1
            if count_contig>=2000:num_2k_contig+=1
            if check:
                #count_scaffold:num of total_scaffold,30427671, 19698289, 23459830
                count_scaffold=total_scaffold[count]
                scaffold_sum+=total_scaffold[count]
                if count_scaffold>=100:num_100_scaffold+=1
                if count_scaffold>=2000:num_2k_scaffold+=1
                if N_info_scaffold and scaffold_sum>=N_info_scaffold[-1]:
                    scaffold_N_seq[-scaffold_i]=[total_scaffold[count],count+1]
                    N_info_scaffold.pop()
                    scaffold_i+=1
            if N_info_contig and contig_sum>=N_info_contig[-1]:
                contig_N_seq[-contig_i]=[total_contig[count],count+1]
                N_info_contig.pop()
                contig_i+=1
            count+=1
    return [[num_100_contig,num_2k_contig,num_100_scaffold,num_2k_scaffold,len(total_contig),len(total_scaffold)],contig_N_seq,scaffold_N_seq,[[contig_max_length,scaffold_max_length],[contig_total_length,scaffold_total_length]],base_info]

def mapped_result(file,samtools):
    #first you need to check the software (samtools)
    if not os.path.exists(file):
        print('sorry!may be the depth %s didn\'t exists,please check it!!'%file)
        sys.exit(0)
    mapped_rate=[]
    for i in os.popen('%s flagstat -@ 6 %s'%(samtools,file)): 
        if re.findall(r'mapped',i):
            mapped_rate.append(re.split(r'[(|)]',i.strip())[1].split(':')[0].strip())
        elif re.findall(r'properly paired',i):
            info=re.split(r'[(|)]',i.strip())[1].split(':')[0].strip()
            if info=='N/A' or info=='-nan%':
                return mapped_rate
            else:
                mapped_rate.append(info)
                return mapped_rate

###############################draw gc&&depth picture###################################
def total_point_num_get(file):
    total_base=0
    with open(file,'r') as f:
        for line in f:
            if line.startswith('>'):continue
            total_base+=len(line.strip())
    str_num=str(int(total_base/2000))
    return int(str_num[0])*10**(len(str_num)-1)

def ref_deal(ref,point_lenth,depth_file,outdir,draw_soft):
    total_gc,gc_len,count_len,list_base,chr_name={},0,0,[],''
    with open(ref,'r') as f:
        for line in f:
            if line.startswith('>'):
                if chr_name:
                    count=0
                    for i in list_base:
                        count+=len(i)
                        if count>point_lenth:
                            count=count-point_lenth
                            gc_len+=(i[:-count].upper().count('G')+i[:-count].upper().count('C'))
                            total_gc[chr_name].append([float(gc_len)/float(point_lenth)])
                            gc_len=0
                        elif count==point_lenth:
                            count=0
                            gc_len+=(i.upper().count('G')+i.upper().count('C'))
                            total_gc[chr_name].append([float(gc_len)/float(point_lenth)])
                            gc_len=0
                        elif count<point_lenth:
                            gc_len+=(i.upper().count('G')+i.upper().count('C'))
                    if count:
                        total_gc[chr_name].append([float(gc_len)/float(point_lenth)])

                chr_name=re.split(r'\s+',line.strip().lstrip('>'))[0]
                if  chr_name not in total_gc:total_gc[chr_name]=[]
                list_base=[]
                continue
            list_base.append(line.strip())
    count=0
    for i in list_base:
        count+=len(i)
        if count>=point_lenth:
            count=count-point_lenth
            gc_len+=(i[:-count].upper().count('G')+i[:-count].upper().count('C'))
            total_gc[chr_name].append([float(gc_len)/float(point_lenth)])
            gc_len=0
        elif count==point_lenth:
            count=0
            total_gc[chr_name].append([float(gc_len)/float(point_lenth)])
            gc_len=0
        elif count<point_lenth:
            gc_len+=(i.upper().count('G')+i.upper().count('C'))
    list_base=[]
    if count:
        total_gc[chr_name].append([float(gc_len)/float(point_lenth)])
    
    ##################depth deal####################total_gc:key:chr   value:[[GC,depth].....]
    dict_chr={}
    with open(depth_file,'r') as f:
        for line in f:
            info=line.strip().split('\t')
            if info[0] not in dict_chr:
                dict_chr[info[0]]=0
                count_point,count,repeat_point_lenth=0,0,point_lenth

            if int(info[1])==repeat_point_lenth:
                repeat_point_lenth+=point_lenth
                total_gc[info[0]][count_point].append(int(count/point_lenth))
                count=0
                count_point+=1
            elif int(info[1])>repeat_point_lenth:
                repeat_point_lenth+=point_lenth
                total_gc[info[0]][count_point].append(int(count/point_lenth))
                count=int(info[2])
                count_point+=1
            count+=int(info[2])
    w=open('%s/result.txt'%outdir,'w')
    for key,value in total_gc.items():
        for i in value:
            try:
                w.write('1\t%d\t%.6f\n'%(i[1],i[0]))
            except:
                pass
    w.close()
    data_file('%s/result.txt'%outdir)
    cmd='Rscript %s %s/result.txt %s/result.png'%(draw_soft,outdir,outdir)
    run_cmd(cmd)

def data_file(infile):
    v1,v2=[],[]
    with open(infile,'r') as f:
        for line in f:
            info=line.strip().split('\t')
            v1.append(int(info[1]))
            v2.append(float(info[2]))
    v1=np.array(v1)
    v2=np.array(v2)
    v1_std,v1_mean=np.std(v1),np.mean(v1)
    v2_std,v2_mean=np.std(v2),np.mean(v2)
    v1_del=np.where(v1>(v1_mean+3*v1_std))
    v1_del_low=np.where(v1<(v1_mean-3*v1_std))
    v2_del=np.where(v2>(v2_mean+3*v2_std))
    v2_del_low=np.where(v2<(v2_mean-3*v2_std))
    v_total=[]
    for s in range(len(v1)):
        if (s in v1_del_low[0]) or (s in v2_del_low[0]) or (s in v1_del[0]) or (s in v2_del[0]):
            print(v1[s],v2[s])
            continue
        v_total.append([v1[s],v2[s]])
    w=open(infile,'w')
    for i in v_total:
        w.write('1\t%d\t%f\n'%(i[0],i[1]))
    w.close()

def depth_result(file):
    total_line,cov_line,total_depth,sta,end,depth_all=0,[0,0,0,0],{},0,0,0
    for line in open(file,'r'):
        info=line.strip().split('\t')
        total_line+=1
        depth_all+=int(info[2])
        if int(info[2])>0:
            cov_line[0]+=1
        if int(info[2])>=4:
            cov_line[1]+=1
        if int(info[2])>=10:
            cov_line[2]+=1
        if int(info[2])>=20:
            cov_line[3]+=1
        if info[0] not in total_depth:
            if total_depth:
                depth_average=int(depth/(end-sta))
                total_depth[key]=depth_average
            key=info[0]
            total_depth[key]=0
            sta=int(info[1])
            depth=0
        end=int(info[1])
        depth+=int(info[2])
    average_depth=float(depth_all)/float(total_line)
    coverage_rate=['%.6f'%(float(info)/float(total_line)) for info in cov_line]
    return total_depth,coverage_rate,average_depth

####################################################################################

############here need to change###########
def result_file(info,file_name):
    info_save,check=[],3
    cog_id=info.split('\t')[0]
    if os.path.exists(file_name):
        for line in open(file_name,'r'):
            if line.startswith(cog_id):
                if cog_id=='<blast_nt' or cog_id=='<busco_evo':
                    if line.strip()==info.strip():
                        check=0
                        break
                    else:
                        check=2
                        continue
                elif cog_id=='<3_mapped':
                    if line.strip()==info.strip():
                        check=0
                        break
                    if info.split('\t')[1]==line.split('\t')[1]:
                        info_save.append(info)
                        check=2
                        continue
                    elif check!=2:
                        check=1  
                elif cog_id=='<2_mapped':
                    if line.strip()==info.strip():
                        check=0
                        break
                    if info.split('\t')[1]==line.split('\t')[1]:
                        info_save.append(info)
                        check=2
                        continue
                    elif check!=2:
                        check=1
                elif cog_id=='<mrna_mapped':
                    if line.strip()==info.strip():
                        check=0
                        break
                    if info.split('\t')[1]==line.split('\t')[1]:
                        info_save.append(info)
                        check=2
                        continue
                    elif check!=2:
                        check=1
            info_save.append(line)
        if check==1:
            w=open(file_name,'a')
            w.write(info)
            w.close()
        elif check==2:
            w=open(file_name,'w')
            for info in info_save:
                w.write(info)
            w.close()
        elif check==3:
            w=open(file_name,'a')
            w.write(info)
            w.close()
        info_save=[]
    else:
        w=open(file_name,'w')
        w.write(info)
        w.close()

################main function################
############report html#############
class Report():
    def getopt(self):
        parser = argparse.ArgumentParser(
            formatter_class=HelpFormatter, description=usage)
        parser.add_argument(
            'func',choices=['ReportHtml'])
        parser.add_argument(
            '-i','--input',help='input file[infofile]',type=str)
        parser.add_argument(
            '-r','--ref',help='input ref file',type=str,dest='ref')
        parser.add_argument(
            '-o','--outdir',help='output dir',type=str,dest='outdir')
        # parser.add_argument('--block',help='circos window',type=int,dest='block')
        parser.add_argument(
            '--species',help='input species name',type=str,dest='species',nargs='+')
        parser.add_argument('-tag',dest='tag',type=str)
        parser.add_argument('-workpath',dest='workpath',type=str)
        args = parser.parse_args()
        if not args.input:
            print('input file must be given!!!')
            sys.exit(0)
        elif not args.outdir:
            print('outdir must be given!!!')
            sys.exit(0)
        elif not args.ref:
            print('ref file must be given!!!')
            sys.exit(0)
        elif not args.species:
            print('species name must be given!!!')
            sys.exit(0)
        return args

    def draw_pie(self,data):
        fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
        
        
        total=sum(data)
        info_name=['A','T','C','G','N']
        recipe = ['{0}:{1: .2f}%'.format(info_name[i],float(data[i]*100)/float(total)) for i in range(5) ]
        
        wedges, texts = ax.pie(data)
        
        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="-"),
                  bbox=bbox_props, zorder=0, va="center")
        
        for i, p in enumerate(wedges):
            ang = (p.theta2 - p.theta1)/2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))
            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = "angle,angleA=0,angleB={}".format(ang)
            kw["arrowprops"].update({"connectionstyle": connectionstyle})
            ax.annotate(recipe[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                         horizontalalignment=horizontalalignment, **kw)
        
        
        plt.show()
        ####save pie picture######
        plt.savefig('result_pie.png')

    def report_result(self,result_file,outdir,ref_file,species_name):
        ##################get each step info###################
        stat_result=stat_N50(ref_file)
        stat_info,base_info,count='','',1
        for i in range(90,40,-10):
            stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">N%d</td><td style="text-align: center;">%d</td><td style="text-align: center;">%d</td><td style="text-align: center;">%d</td><td style="text-align: center;">%d</td></tr>'''%(i,stat_result[1][-count][0],stat_result[1][-count][1],stat_result[2][-count][0],stat_result[2][-count][1])
            count+=1
            
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Total length</td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td></tr>'''%(stat_result[3][1][0],stat_result[3][1][1])
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Number(>=100bp)</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td></tr>'''%(stat_result[0][0],stat_result[0][2])
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Number(>=2kb)</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td></tr>'''%(stat_result[0][1],stat_result[0][3])
    #     stat_info+='''
    #         <tr style="height: 5px;">
    # <td style="text-align: center;">Total number</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td></tr>'''%(stat_result[0][4],stat_result[0][5])
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Max length</td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td></tr>'''%(stat_result[3][0][0],stat_result[3][0][1])
        #picture change  
    #     base_total=['A','T','C','G','N','GC','Total']
    #     for i in range(len(base_total)):
    #         base_info+='''
    # <tr style="height: 5px;">
    # <td style="text-align: center;">{0}</td><td style="text-align: center;">{1}</td><td style="text-align: center;">{2: .2f}%</td></tr>'''.format(base_total[i],stat_result[-1][i],(float(stat_result[-1][i])/float(stat_result[-1][-1]))*100)
        base_info="<img src=\"result_pie.png\">"
        self.draw_pie(stat_result[-1][:5])

        info_2_mapp,info_2_hm,info_3,info_nt,png_2,png_3,info_mrna,info_busco='','','','','','','',''
        png2_check,png3_check=0,0
        dick_check={'<3_mapped':0,'<2_mapped':0,'<blast_nt':0,'<mrna_mapped':0,'<busco_evo':0}
        dict_func={'<3_mapped':'Reads3Stat','<2_mapped':'Reads2Stat','<blast_nt':'BlastNtMapp','<mrna_mapped':'MrnaMapped','<busco_evo':'BuscoEvaluation'}
        for line in open(result_file,'r'):
            info=line.strip().split('\t')
            dick_check[info[0]]+=1
        ################2 mapp info#############
            if line.startswith('<2_mapped'):
                info_2_mapp+='''
                <tr style="height: 5px;">
    <td width="300px" style="text-align: center;">{0}</td><td width="700px" style="text-align: center;">{1}</td><td width="700px" style="text-align: center;">{2}</td><td width="700px" style="text-align: center;">{3: .2f}</td><td width="700px" style="text-align: center;">{4}</td><td width="700px" style="text-align: center;">{5}</td><td width="700px" style="text-align: center;">{6}</td><td width="700px" style="text-align: center;">{7}</td></tr>'''.format(info[1],info[2].split(':')[0],info[2].split(':')[1],float(info[3].split(':')[1]),info[4].split(':')[1],info[4].split(':')[2],info[4].split(':')[3],info[4].split(':')[4])
                if info[1]=='Total':
                    info_2_hm+='''
            <tr style="height: 5px;">
        <td style="text-align: center;">{0}</td><td style="text-align: center;">{1: .2f}%</td><td style="text-align: center;">{2: .2f}%</td><td style="text-align: center;">{3: .2f}%</td><td style="text-align: center;">{4: .2f}%</td></tr>'''.format(info[1],float(info[5].split(':')[1])*100/float(stat_result[-1][-1]),float(info[7].split(':')[1])*100/float(stat_result[-1][-1]),float(info[4].split(':')[1])*100/float(stat_result[-1][-1]),float(info[6].split(':')[1])*100/float(stat_result[-1][-1]))
                    png_2+="<img src=\"assemble2/result.png\">"

        #############3 mapp info#############
            elif line.startswith('<3_mapped'):
                info_3+='''
            <tr style="height: 5px;">
    <td width="300px" style="text-align: center;">{0}</td><td width="700px" style="text-align: center;">{1}</td><td width="700px" style="text-align: center;">{2: .2f}</td><td width="700px" style="text-align: center;">{3}</td><td width="700px" style="text-align: center;">{4}</td><td width="700px" style="text-align: center;">{5}</td><td width="700px" style="text-align: center;">{6}</td></tr>'''.format(info[1],info[2].split(':')[1],float(info[3].split(':')[1]),info[4].split(':')[1],info[4].split(':')[2],info[4].split(':')[3],info[4].split(':')[4])

                png_3+="<img src=\"assemble3/result.png\">"

        #############mrna mapp info##########
            elif line.startswith('<mrna_mapped'):
                info_mrna+='''
            <tr style="height: 5px;">
        <td width="300px" style="text-align: center;">{0}</td><td width="700px" style="text-align: center;">{1}</td><td width="700px" style="text-align: center;">{2}</td></tr>'''.format(info[1],info[2].split(':')[1],info[2].split(':')[2])

        #############nt mapp info############ total_base = total p_ctg.fa base num
            elif line.startswith('<blast_nt'):
                for nt_line in info[3:]:
                    nt_info=re.split(r'[:|,]',nt_line)
                    if '' in nt_info:
                        nt_info.remove('')
                    info_nt+='''
            <tr style="height: 5px;">
        <td style="text-align: center;">{0}</td><td style="text-align: center;">{1}</td><td style="text-align: center;">{2}</td><td style="text-align: center;">{3: .3f}%</td></tr>'''.format(nt_info[0],format(int(nt_info[-2]),','),format(int(info[2].split(':')[1]),','),float(nt_info[-1])*100)

        ###########busco evoluation info#####
            elif line.startswith('<busco_evo'):
                for busco_line in info[1:]:
                    busco_info=busco_line.split(',')
                    info_busco+='''
                <tr style="height: 5px;">
        <td style="text-align: center;">{0}</td><td style="text-align: center;">{1}</td><td style="text-align: center;">{2}</td></tr>'''.format(busco_info[0],format(int(busco_info[1]),','),busco_info[2])
        
        for key,value in dick_check.items():
            if not value:
                print('please check the \"%s\",you can restart function \"\"!!!!'%(key,dict_func[key]))
                sys.exit(0)
        w=open('%s/assemble_report.html'%outdir,'w')
        info='''<!DOCTYPE html>
<html>
<head>
    <TITLE>基因组组装评估报告</TITLE>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />   
    
    <link href="src/css/bootstrap.min.css" type="text/css" rel="stylesheet" />
    <link rel="stylesheet" type="text/css" href="src/css/css.css" />
    <link href="src/css/nav.css" type="text/css" rel="stylesheet" />
    <link href="src/js/fancybox/jquery.fancybox.css" rel="stylesheet" type="text/css" /> 
    <link href="src/css/jscrollpane_gooalgene.css" rel="stylesheet" type="text/css" />
    <script src="src/js/jquery-1.11.3.min.js" type="text/javascript"></script>
    <script src="src/js/nav.js" type="text/javascript"></script>
    <script src="src/js/bootstrap.min.js" type="text/javascript"></script>
    <script src="src/js/fancybox/jquery.fancybox.pack.js" type="text/javascript"></script>
    <script src="src/js/albumSider.js" type="text/javascript"></script>
    <script src="src/js/report.js" type="text/javascript"></script>
    <script src="src/js/js.js" type="text/javascript"></script>
    <script src="src/js/jquery.jscrollpane.min.js" type="text/javascript"></script>
    <script src="src/js/jquery.mousewheel.js" type="text/javascript"></script>
    <script src="src/js/jquery.treeview.js" type="text/javascript"></script>
</head>
<!DOCTYPE html>
<body>
    <div style="width: 97%;margin: 0 auto;clear: both;">
        <a href="http://www.gooalgene.com" target="blank"><img style="display: inline;width:16%;" src="src/images/logo.png" title="logo" alt="logo.png" /></a>
        <h1><i>{0}</i> 基因组组装评估报告</h1>
    </div>
    <HR style="FILTER: alpha(opacity=100,finishopacity=0,style=3)" color=#987cb9 SIZE=3>
    <div class="container" style="width: 100%;">
        <div class="row" >
            <div role="complementary" class="col-md-3" style="page-break-after:always;">
                <nav class="bs-docs-sidebar">
                <p style="font:1.3em Microsoft YaHei;text-indent:1em;font-weight:bold;line-height: 1.5;">目录栏</p>
                    <ul class="nav bs-docs-sidenav">
                        <li><a href="#a1">1.项目背景</a></li>
                        <li><a href="#b2">2.组装评估分析流程</a></li>
                        <li><a href="#c3">3.组装结果统计</a></li>
                        <li><a href="#d4">4.数据正确性评估</a></li>
                        <li><a href="#e5">5.序列一致性评估</a></li>
                        <li><a href="#f6">6.组装完整性评估</a> 
                            <ul class="nav">
                                <li><a href="#f61" class="js-active-bg">5.1 保守基因完整性评估</a></li>
                            </ul>
                            <ul class="nav">
                                <li><a href="#f62" class="js-active-bg">5.1 基因区完整性评估</a></li>
                            </ul>
                        </li>
                        <li><a href="#g7">7.准确性及杂合度评估</a>
                        </li>
                        <li><a href="#h8">8.基因组圈图</a></li>
                    </ul>
                </nav>
            </div>  
    <div role="main" class="col-md-9" >     
                <div class="tMainBody ">                    
                    <h2 id="a1">1.项目背景</h2>
                    <div class="studyUnit">
                        <p>物种名称:{1}</p>
                        <p class="paragraph">分析内容：将该基因组与NT数据库、其二三代数据、转录组数据进行比对，并通过BUSCO等方法，整体地评估{2}基因组的组装效果。</p>
                    </div>
                </div>
                <div class="tMainBody ">
                    <h2 id="b2">2.组装评估分析流程</h2>
                    <img src='基因组组装评估流程图.png'>
                    <p class="name">图1.组装评估分析流程图
                </div>
                <div class="tMainBody ">
                    <h2 id="c3">3. 组装结果统计</h2>
                    <p class="paragraph">组装得到基因组序列的基本情况统计见下表:</p>
                    <div class="studyUnit">
                        <p class="name">表1. 组装结果序列统计</p>
<table class="table table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
    <td width="300px" style="text-align: center;"></td><td width="600px" style="text-align: center;">Contig</td><td width="600px" style="text-align: center;">Scaffold</td></tr></table>
<table class="table table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
    <td width="300px" style="text-align: center;"></td><td style="text-align: center;">Length(bp)</td><td style="text-align: center;">Number</td><td style="text-align: center;">Size(bp)</td><td style="text-align: center;">Number</td></tr>{3}</table></div>
        <p class="premark"><b>注：</b>
            <br/>&emsp;&emsp;(1) N90 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 90%的时候，那一条序列的长度（bp）。
            <br/>&emsp;&emsp;(2) N80 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 80%的时候，那一条序列的长度（bp）。
            <br/>&emsp;&emsp;(3) N70 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 70%的时候，那一条序列的长度（bp）。
            <br/>&emsp;&emsp;(4) N60 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 60%的时候，那一条序列的长度（bp）。
            <br/>&emsp;&emsp;(5) N50 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 50%的时候，那一条序列的长度（bp）。
            <br/>&emsp;&emsp;(6) Max length：组装得到的最长的序列的长度（bp）。
            <br/>&emsp;&emsp;(7) Total  length：组装得到的序列总长（bp）。
            <br/>&emsp;&emsp;(8) Number(>=100bp)：组装得到长度大于100bp的序列条数。
            <br/>&emsp;&emsp;(9) Number(>=2kb)：组装得到长度大于2kb的序列条数。
            <br/>&emsp;&emsp;(10) Total number：组装得到的序列条数。
        </p>
                    <p class="paragraph">对基因组序列进行A、T、C、G的含量统计分析，统计结果见下图:</p>
                        <div class="studyUnit">{4}
                        <p class="name">图1 .碱基组成分布情况统计</p>
</div>
    </div>
            <div class="tMainBody ">
                <h2 id="d4">4 数据正确性评估</h2>
                    <p class="paragraph">为了确认组装结果属于目标物种，以1000 bp为步长将基因组序列打断，把打断后的序列通过软件Blast比对到NCBI核苷酸数据库（NT库），表3展示了比对次数排名前五的物种。</p>
                    <div class="studyUnit">
                        <p class="name">表3.打断后序列与NT库比对情况部分展示</p>
<table class="table table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
     <td style="text-align: center;">Species</td><td style="text-align: center;">Blast number</td><td style="text-align: center;">Total blast number</td><td style="text-align: center;">Total (%)</td></tr>{5}</table>
            </div>   
        </div>       
            <div class="tMainBody ">
                <h2 id="e5">5 序列一致性评估</h2>
                    <p class="paragraph">分别使用软件bwa和blasr将二、三代数据比对到XXX基因组上，统计比对率以及覆盖率，这两个值越高则认为组装结果与reads的一致性越高，组装的效果也相应的更好。</p>
                    <div class="studyUnit">
                        <p class="name">表4.二代reads比对结果统计</p>
<table class="table left table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;"><td width="500px" style="text-align: center;">Sample name</td><td width="700px" style="text-align: center;">Mapping rate(%)</td><td width="700px" style="text-align: center;">Paired mapping rate(%)</td><td width="700px" style="text-align: center;">Average sequencing depth</td><td width="700px" style="text-align: center;">Coverage(%)</td><td width="700px" style="text-align: center;">Coverage at least 4X(%)</td><td width="700px" style="text-align: center;">Coverage at least 10X(%)</td><td width="700px" style="text-align: center;">Coverage at least 20X(%)</td></tr>{6}</table>
        <p class="premark"><b>注：</b>
            <br/>&emsp;&emsp;(1) Mapping rate：能够比对到基因组的二代数据比例；
            <br/>&emsp;&emsp;(2) Paired mapping rate：能够成对地比对到基因组的二代数据比例；
            <br/>&emsp;&emsp;(3) Average sequence depth：二代数据的平均测序深度；
            <br/>&emsp;&emsp;(4) Coverage：基因组上被二代数据覆盖区域占整个基因组的比例；
            <br/>&emsp;&emsp;(5) Coverage at least 4X：基因组上至少被4X二代数据覆盖区域占整个基因组的比例；
            <br/>&emsp;&emsp;(6) Coverage at least 10X：基因组上至少被10X二代数据覆盖区域占整个基因组的比例；
            <br/>&emsp;&emsp;(7) Coverage at least 20X：基因组上至少被20X二代数据覆盖区域占整个基因组的比例。
        </p>
                    </div>
                    <p class="paragraph"></p>
                    <div class="studyUnit">
                        <p class="name">表5.三代reads比对结果统计</p>
<table class="table left table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;"><td width="500px" style="text-align: center;">Sample name</td><td width="700px" style="text-align: center;">Mapping rate(%)</td><td width="700px" style="text-align: center;">Average sequencing depth</td><td width="700px" style="text-align: center;">Coverage(%)</td><td width="700px" style="text-align: center;">Coverage at least 4X(%)</td><td width="700px" style="text-align: center;">Coverage at least 10X(%)</td><td width="700px" style="text-align: center;">Coverage at least 20X(%)</td></tr>{7}</table>
        <p class="premark"><b>注：</b>
            <br/>&emsp;&emsp;(1) Mapping rate：能够比对到基因组的三代数据比例；
            <br/>&emsp;&emsp;(2) Paired mapping rate：能够成对地比对到基因组的三代数据比例；
            <br/>&emsp;&emsp;(3) Average sequence depth：三代数据的平均测序深度；
            <br/>&emsp;&emsp;(4) Coverage：基因组上被三代数据覆盖区域占整个基因组的比例；
            <br/>&emsp;&emsp;(5) Coverage at least 4X：基因组上至少被4X三代数据覆盖区域占整个基因组的比例；
            <br/>&emsp;&emsp;(6) Coverage at least 10X：基因组上至少被10X三代数据覆盖区域占整个基因组的比例；
            <br/>&emsp;&emsp;(7) Coverage at least 20X：基因组上至少被20X三代数据覆盖区域占整个基因组的比例。
        </p>
                    </div>
            <p class="paragraph"></p>      
            <p class="paragraph">同时绘制GC与深度分布关联分析统计图，展示、评估测序的均匀性，图中横坐标表示GC含量，纵坐标表示覆盖深度，右侧展示了contig覆盖深度分布，上方展示了GC含量分布，中间大图为根据contigs的GC分布及覆盖深度信息绘制的散点图，其中颜色深浅用于反映散点图中点的密度。</p>
                        <div class="studyUnit">
                                {8}{9}{10}                               
                    <p class="name">图1. GC含量与覆盖深度（Depth）关联分析统计图（二代reads比对)</p>
                    </div>
                    <div class="studyUnit">
                        {11}{12}{13}
                    <p class="name">图2. GC含量与覆盖深度（Depth）关联分析统计图（三代reads比对）</p>
            </div>
            <div class="tMainBody ">
                <h2 id="f6">6 组装完整性评估</h2>
                    <p class="paragraph">根据基因组与转录组数据及保守基因集的比对，综合评估组装基因组的完整性。</p>
                    <div class="tMainBody ">
                        <h3 id="f61">6.1.保守基因完整性评估</h3>
                        <p class="paragraph">BUSCO（Benchmarking Universal Single-Copy Orthologs）是利用直系同源数据库OrthoDB对基因组组装的完整性进行定量评估的软件。 BUSCO抽样了数百个基因组，从中选择单拷贝直系同源＞90%的基因构建了六种主要的系统进化分枝的基因集。</p>
                        <p class="paragraph">本次评估采用的BUSCO基因集为{14}。BUSCO评估结果见表6。</p>
                        <div class="studyUnit">
                        <p class="name">表6.BUSCO评估</p>
<table class="table left table-hover table-striped table-bordered" style="width:100%;">
                        <tr style="height: 5px;">
    <td style="text-align: center;">Type</td><td style="text-align: center;">Proteins</td><td style="text-align: center;">Percentage(%)</td></tr>
        <tr style="height: 5px;">{15}</table>
        <p class="premark"><b>注：</b>
            <br/>&emsp;&emsp;(1) Complete BUSCOs：能完整比对BUSCO的基因；
            <br/>&emsp;&emsp;(2) Complete and single-copy BUSCOs：一个BUSCO 能完整比对上一个基因；
            <br/>&emsp;&emsp;(3) Complete and duplicated BUSCOs：一个BUSCO 能完整比对上多个基因；
            <br/>&emsp;&emsp;(4) Fragmented BUSCOs：只有部分序列能比对上BUSCO profile的基因；
            <br/>&emsp;&emsp;(5) Missing BUSCOs：没有能比对上BUSCO profile的基因；
            <br/>&emsp;&emsp;(6) Total BUSCO groups searched：BUSCO groups的基因总数。
        </p>
                    </div>
                </div>
                    <div class="tMainBody ">
                        <h3 id="f62">6.2.基因区完整性评估</h3>
                        <p class="paragraph">通过软件hisat2将转录组数据比对到基因组上，统计比对率，该值越高则表示基因组相对来说组装得越完整。</p>
                        <div class="studyUnit">
                        <p class="name">表7.转录组比对率统计</p>
                        <table class="table left table-hover table-striped table-bordered" style="width:100%;">
                        <tr style="height: 5px;"><td width="500px" style="text-align: center;">Sample name</td><td width="700px" style="text-align: center;">Mapping rate(%)</td><td width="700px" style="text-align: center;">Paired mapping rate(%)</td></tr>{16}</table>
                        <p class="premark"><b>注：</b>
                <br/>&emsp;&emsp;(1) Mapping rate：能够比对到基因组的转录组数据比例；
                <br/>&emsp;&emsp;(2) Paired Mapping rate：能够成对地比对到基因组的转录组数据比例。</p>
                </div>
            </div>
        </div>
            <div class="tMainBody ">
                <h2 id="g7">7 准确性及杂合度评估</h2>
                <p class="paragraph">将二代reads比对到基因组后，通过软件samtools、picard以及GATK识别突变。分别统计SNP、InDel的纯合率和杂合率，纯合率越低表明基因组准确性越高，杂合率越高说明基因组的杂合率越高。</p>
                <div class="studyUnit">
                    <p class="name">表8.纯合率及杂合率统计</p>
                    <table class="table left table-hover table-striped table-bordered" style="width:100%;">
                        <tr style="height: 5px;"><td width="500px" style="text-align: center;">Sample name</td><td style="text-align: center;">Homozygous SNP rate(%)</td><td style="text-align: center;">Homozygous InDel rate(%)</td><td style="text-align: center;">Heterozygous SNP rate(%)</td><td style="text-align: center;">Heterozygous InDel rate(%)</td></tr>{17}</table>
                        <p class="premark"><b>注：</b>
                            <br/>&emsp;&emsp;(1) Homozygous SNP rate：纯合SNP位点占基因组比例；
                            <br/>&emsp;&emsp;(2) Homozygous InDel rate：纯合InDel位点占基因组比例；
                            <br/>&emsp;&emsp;(3) Heterozygous SNP rate：杂合SNP位点占基因组比例；
                            <br/>&emsp;&emsp;(4) Heterozygous InDel rate：杂合InDel位点占基因组比例。
                        </p>
                </div>
            </div>
            <div class="tMainBody ">
                <h2 id="h8">8 基因组圈图</h2>
                <p class="paragraph">绘制基因组圈图整体展示评估结果。</p>
                <img src="circos/circos.png">
                <p class="name">图3. 基因组圈图</p>
                <p class="premark"><b>注： window为{18}k</b>
                    <br/>&emsp;&emsp;A.基因组信息；
                    <br/>&emsp;&emsp;B.GC含量分布；
                    <br/>&emsp;&emsp;C.二代reads深度分布；
                    <br/>&emsp;&emsp;D.三代reads深度分布；
                    <br/>&emsp;&emsp;E.外圈为纯合SNP分布，内圈为杂合SNP分布；
                    <br/>&emsp;&emsp;F.外圈为纯合InDel分布，内圈为杂合InDel分布；
                    <br/>&emsp;&emsp;G. 完整比对BUSCO的基因在基因组上分布情况，蓝色为single-copy BUSCO，红色为duplicated BUSCO。
                </p>
                </div>
            </div>
        </div>
    </div>    
</div>  
</body>
</html>
    '''.format(species_name,species_name,species_name,stat_info,base_info,info_nt,info_2_mapp,info_3,'','',png_2,'','',png_3,'busco_dataset',info_busco,info_mrna,info_2_hm,format(float(total_point_num_get(ref_file))/float(1000),','))
        w.write(info.encode('utf-8'))
        w.close()

    def main(self):
        args=self.getopt()
        if contiue_run(args.tag,args.workpath):
            config=Config()
            check_dir(args.outdir)
            self.report_result(args.input,args.outdir,args.ref,' '.join(args.species))

############3 read mapp#############
class Assemble3_stat():
    def getopt(self):
        parser = argparse.ArgumentParser(
            formatter_class=HelpFormatter, description=usage)
        parser.add_argument(
            'func',choices=['Reads3Stat'])
        parser.add_argument(
            '-r','--ref',help='assemble file that you need to evalue',dest='ref',type=str)
        parser.add_argument(
            '-i','--input',help='input fofn file or read 3 file',type=str,nargs='+')
        parser.add_argument(
            '-o','--outdir',help='output dir',type=str,dest='outdir')
        parser.add_argument(
            '-n','--thread',help='set num of thread',type=int,dest='thread',default=10)
        parser.add_argument(
            '--vf',help='set memory of each thread[unit:G,default:2]',dest='vf',default=2,type=int)
        parser.add_argument(
            '--infofile',help='the each step info file\nuse for report',dest='infofile',type=str)
        parser.add_argument('-tag',dest='tag',type=str)
        parser.add_argument('-workpath',dest='workpath',type=str)
        args = parser.parse_args()
        if not args.ref:
            print('assemble file must be given!!!')
            sys.exit(0)
        elif not args.input:
            print('input fofn file must be given!!!')
            sys.exit()
        elif not args.outdir:
            print('outdir must be given!!!')
            sys.exit(0)
        return args

    def assemble3_stat(self,ref,fofnfile,blasr,samtools,r_gc_depth,outdir,nproc):
        args=self.getopt()
        outdir=outdir+'/assemble3'
        check_dir(outdir)
        #################3 mapped###########
        # cmd='%s %s %s --nproc %s --bam --out %s/result.bam'%(blasr,fofnfile,ref,nproc,outdir)    
        # run_cmd(cmd)
        # # #################3 sort#############
        # cmd='%s sort -@ %s -m %sG -o %s/result.sort.bam %s/result.bam'%(samtools,nproc,args.vf,outdir,outdir)       
        # run_cmd(cmd)
        # # #################3 depth###########
        # cmd='%s depth -aa %s/result.sort.bam>%s/result.depth'%(samtools,outdir,outdir)
        # run_cmd(cmd)

        total_point_num=total_point_num_get(ref)
        ref_deal(ref,total_point_num,'%s/result.depth'%outdir,outdir,r_gc_depth)
        mapped_rate=mapped_result('%s/result.sort.bam'%(outdir),samtools)
        total_depth,coverage_rate,average_depth=depth_result('%s/result.depth'%(outdir))
        if os.path.exists('%s/result.bam'%(outdir)):os.remove('%s/result.bam'%(outdir))
        
        return coverage_rate,mapped_rate,'%s/result.png'%outdir,average_depth

    def main(self):
        args=self.getopt()
        if contiue_run(args.tag,args.workpath):
            config=Config()
            blasr = check_software('blasr') if check_software('blasr') else config.blasr
            samtools = check_software('samtools') if check_software('samtools') else config.samtools
            r_gc_depth=check_software('gc-coverage.r') if check_software('gc-coverage.r') else config.r_gc_depth
            check_dir(args.outdir)
            
            coverage_rate3,mapped_rate3,png_3,average_depth=self.assemble3_stat(args.ref,' '.join(args.input),blasr,samtools,r_gc_depth,args.outdir,args.thread)
            sample_name=args.input[0].split('/')[-1].split('_')[0] if len(args.input)>1 else ''
            info='<3_mapped\tTotal\tmapped_rate:%s\taverage_depth:%.3f\tcoverage_rate:%s\tpng_3:%s\n'%(':'.join(mapped_rate3),average_depth,':'.join(coverage_rate3),png_3)
            file_info_get='%s/result.txt'%args.outdir if not args.infofile else args.infofile
            result_file(info,file_info_get)

class Assemble2_stat():
    def getopt(self):
        parser = argparse.ArgumentParser(
            formatter_class=HelpFormatter, description=usage)
        parser.add_argument(
            'func',choices='Reads2Stat')
        parser.add_argument(
            '-r','--ref',help='assemble file that you need to evalue',dest='ref',type=str)
        parser.add_argument(
            '-i','--input',help='input file[2_read_1,2_read_2]',nargs=2,type=str)
        parser.add_argument(
            '-o','--outdir',help='output dir',type=str,dest='outdir')
        parser.add_argument(
            '-n','--thread',help='set num of thread',type=int,dest='thread',default=10)
        parser.add_argument(
            '--bams',help='input all bams file',type=str,dest='bams',nargs='+')
        parser.add_argument(
            '--prefix',help='input sample name',type=str,dest='prefix')
        parser.add_argument(
            '--vf',help='set memory of each thread[unit:G,default:2]',dest='vf',default=2,type=int)
        parser.add_argument(
            '--infofile',help='the each step info file\nuse for report',dest='infofile',type=str)
        parser.add_argument('-tag',dest='tag',type=str)
        parser.add_argument('-workpath',dest='workpath',type=str)
        args = parser.parse_args()
        if not args.ref:
            print('assemble file must be given!!!')
            sys.exit(0)
        elif not args.input and not args.bams:
            print('input file must be given!!!')
            sys.exit(0)
        elif not args.outdir:
            print('outdir must be given!!!')
            sys.exit(0)
        return args

    def bam_stat(self,ref,bams,samtools,bwa,r_gc_depth,thread_num,outdir):
        args=self.getopt()
        outdir=outdir+'/assemble2'
        check_dir(outdir)
        config=Config()
        gatk=check_software('GenomeAnalysisTK.jar') if check_software('GenomeAnalysisTK.jar') else config.gatk
        picard=check_software('picard.jar') if check_software('picard.jar') else config.picard
        samtools=check_software('samtools') if check_software('samtools') else config.samtools
        cmd="%s merge %s/total_result.bam %s"%(samtools,outdir,' '.join(bams))
        if len(bams)>1:
            run_cmd(cmd)
        else:
            cmd='mv %s %s/total_result.bam'%(bams[0],outdir)
            run_cmd(cmd)

        cmd="%s sort -m %sG -@ %s -o %s/result.sorted.bam %s/total_result.bam"%(samtools,args.vf,thread_num,outdir,outdir)
        run_cmd(cmd)

        ##################depth######################
        cmd='%s depth -aa %s/result.sorted.bam>%s/result.depth'%(samtools,outdir,outdir)
        run_cmd(cmd)
        #######################picard sort##############
        cmd="java -Djava.io.tmpdir=%s/tmp -jar %s SortSam VALIDATION_STRINGENCY=LENIENT INPUT=%s/result.sorted.bam OUTPUT=%s/result.sort.bam SORT_ORDER=coordinate"%(outdir,picard,outdir,outdir)
        run_cmd(cmd)

        ######################add header################
        cmd="java -jar %s AddOrReplaceReadGroups I=%s/result.sort.bam O=%s/result.bam RGID=4 RGLB=result RGPL=illumina RGPU=result RGSM=20"%(picard,outdir,outdir)
        run_cmd(cmd)

        #######################samtools index###########
        cmd="%s index %s/result.bam"%(samtools,outdir)
        run_cmd(cmd)

        ######################calling##################
        cmd="java -Xmx4G -jar %s -T HaplotypeCaller -R %s -I %s/result.bam --genotyping_mode DISCOVERY -o %s/result.gvcf -nct %s"%(gatk,ref,outdir,outdir,thread_num)
        run_cmd(cmd)

        ######################select snp && indel and stat#######
        cmd='java -Xmx4G -jar %s -T SelectVariants -o %s/result.snp.vcf --variant %s/result.gvcf -selectType SNP -R %s'%(gatk,outdir,outdir,ref)
        run_cmd(cmd)
        homog,mix,chroms,total_site,rate_total=[0,0],[0,0],[],[0,0],[0,0,0,0]
        for line in open('%s/result.snp.vcf'%(outdir),'r'):
            if line.startswith('#'):continue
            chrs=line.split('\t')[0]
            if chrs not in chroms:chroms.append(chrs)
            types=line.strip().split('\t')[-1].split(':')[0]
            if types=='1/1':
                homog[0]+=1
            elif types=='0/1':
                mix[0]+=1
            total_site[0]+=1
        rate_total[1]=homog[0]
        rate_total[0]=mix[0]

        cmd='java -Xmx2G -jar %s -T SelectVariants -o %s/result.indel.vcf --variant %s/result.gvcf -selectType INDEL -R %s'%(gatk,outdir,outdir,ref)
        run_cmd(cmd)
        for line in open('%s/result.indel.vcf'%(outdir),'r'):
            types=line.strip().split('\t')[-1].split(':')[0]
            if types=='1/1':
                homog[1]+=1
            elif types=='0/1':
                mix[1]+=1
            total_site[1]+=1
        rate_total[3]=homog[1]
        rate_total[2]=mix[1]

        ###################filter#######################
        ###snp filt
        cmd='java -Xmx4g -jar %s -R %s -T VariantFiltration --variant %s/result.snp.vcf --clusterSize 4 --clusterWindowSize 10 --maskName aroundIndel --mask %s/result.indel.vcf -maskExtend 3 --filterName \"lowMQRankSum\" --filterExpression \"MQRankSum < -12.5\" --filterName \"highFS\" --filterExpression \"FS > 60.0\" --filterName \"lowReadPosRankSum\" --filterExpression \"ReadPosRankSum < -8.0\" --filterName \"lowMQ\" --filterExpression \"MQ < 40.0\" --filterName \"lowQD\" --filterExpression \"QD < 2.0\" --out %s/snp.filt.vcf --genotypeFilterName \"lowDP\" --genotypeFilterExpression \"DP < 8.0\"'%(gatk,ref,outdir,outdir,outdir)
        run_cmd(cmd)
        ###indel filt
        cmd='java -Xmx4g -jar %s -T VariantFiltration -R %s -V %s/result.indel.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum <-20.0\" --filterName \"filter\" -o %s/indel.filt.vcf'%(gatk,ref,outdir,outdir)
        run_cmd(cmd)
        
        total_point_num=total_point_num_get(ref)
        ref_deal(ref,total_point_num,'%s/result.depth'%(outdir),outdir,r_gc_depth)
        mapped_rate=mapped_result('%s/result.bam'%(outdir),samtools)
        total_depth,coverage_rate,average_depth=depth_result('%s/result.depth'%(outdir))
        # png_2=gc_depth(total_depth,total_gc,r_gc_depth,'result',outdir)
        if os.path.exists('%s/result.sam'%(outdir)):os.remove('%s/result.sam'%(outdir))
        if os.path.exists('%s/result.sorted.bam'%(outdir)):os.remove('%s/result.sorted.bam'%(outdir))
        if os.path.exists('%s/result.sort.bam'%(outdir)):os.remove('%s/result.sort.bam'%(outdir))
        
        return coverage_rate,mapped_rate,rate_total,'%s/result.png'%outdir,average_depth


    def assemble2_stat(self,ref,tartget1,tartget2,samtools,bwa,thread_num,outdir):
        args=self.getopt()
        outdir=outdir+'/assemble2'
        check_dir(outdir)
        config=Config()
        #####################software###############
        gatk=check_software('GenomeAnalysisTK.jar') if check_software('GenomeAnalysisTK.jar') else config.gatk
        picard=check_software('picard.jar') if check_software('picard.jar') else config.picard
        ##################bwa index###############
        cmd='%s index %s'%(bwa,ref)
        run_cmd(cmd)
        ###################samtools faidx##########
        cmd='%s faidx %s'%(samtools,ref)
        run_cmd(cmd)
        ##################picard index############
        cmd='java -jar %s CreateSequenceDictionary R=%s O=%s.dict'%(picard,ref,'.'.join(ref.split('.')[:-1]))
        run_cmd(cmd)
        ##################bwa mem##################
        cmd = "%s mem -t %s %s %s %s > %s/result.sam"%(bwa,thread_num,ref,tartget1,tartget2,outdir)
        run_cmd(cmd)
        cmd="%s view -@ %s -bSh %s/result.sam -o %s/%s.bam"%(samtools,thread_num,outdir,outdir,args.prefix)
        run_cmd(cmd)
        cmd="%s sort -m %sG -@ %s -o %s/result.sorted.bam %s/%s.bam"%(samtools,args.vf,thread_num,outdir,outdir,args.prefix)
        run_cmd(cmd)

        ##################depth######################
        cmd='%s depth -aa %s/result.sorted.bam>%s/result.depth'%(samtools,outdir,outdir)
        run_cmd(cmd)
        
        mapped_rate=mapped_result('%s/%s.bam'%(outdir,args.prefix),samtools)
        total_depth,coverage_rate,average_depth=depth_result('%s/result.depth'%(outdir))
        if os.path.exists('%s/result.sam'%(outdir)):os.remove('%s/result.sam'%(outdir))
        if os.path.exists('%s/result.sorted.bam'%(outdir)):os.remove('%s/result.sorted.bam'%(outdir))
        if os.path.exists('%s/result.sort.bam'%(outdir)):os.remove('%s/result.sort.bam'%(outdir))
        
        return coverage_rate,mapped_rate,average_depth
    
    def main(self):
        args=self.getopt()
        if contiue_run(args.tag,args.workpath):
            config=Config()
            check_dir(args.outdir)
            samtools = check_software('samtools') if check_software('samtools') else config.samtools
            r_gc_depth=check_software('gc-coverage.r') if check_software('gc-coverage.r') else config.r_gc_depth
            bwa=check_software('bwa') if check_software('bwa') else config.bwa
            if args.input:
                coverage_rate2,mapped_rate2,average_depth=self.assemble2_stat(args.ref,args.input[0],args.input[1],samtools,bwa,args.thread,args.outdir)
                info='<2_mapped\t%s\tmapped_rate:%s\taverage_depth:%.3f\tcoverage_rate:%s\n'%(args.prefix,':'.join(mapped_rate2),average_depth,':'.join(coverage_rate2))
            elif args.bams:
                coverage_rate2,mapped_rate2,total_rate,png_2,average_depth=self.bam_stat(args.ref,args.bams,samtools,bwa,r_gc_depth,args.thread,args.outdir)
                info='<2_mapped\tTotal\tmapped_rate:%s\taverage_depth:%.3f\tcoverage_rate:%s\tsnp_mix:%d\tsnp_homog:%d\tindel_mix:%d\tindel_homog:%d\tpng_2:%s\n'%(':'.join(mapped_rate2),average_depth,':'.join(coverage_rate2),total_rate[0],total_rate[1],total_rate[2],total_rate[3],png_2)
            file_info_get='%s/result.txt'%args.outdir if not args.infofile else args.infofile
            result_file(info,file_info_get)

class Assemble_mrna():
    def getopt(self):
        parser = argparse.ArgumentParser(
            formatter_class=HelpFormatter, description=usage)
        parser.add_argument(
            'func',choices=['MrnaMapped'])
        parser.add_argument(
            '-r','--ref',help='assemble file that you need to evalue',dest='ref',type=str,required=True)
        parser.add_argument(
            '-i','--input',help='input file[mrna_read_1,mrna_read_2]',dest='input',nargs=2,type=str,required=True)
        parser.add_argument(
            '-o','--outdir',help='output dir',type=str,dest='outdir',required=True)
        parser.add_argument(
            '-n','--thread',help='set num of thread',type=int,dest='thread',default=10)
        parser.add_argument(
            '--prefix',help='input sample name',type=str,dest='prefix')
        parser.add_argument(
            '--infofile',help='the each step info file\nuse for report',dest='infofile',type=str)
        parser.add_argument('-tag',dest='tag',type=str)
        parser.add_argument('-workpath',dest='workpath',type=str)
        args=parser.parse_args()
        return args

    def assemble_mrna(self,ref,hisat2_build,hisat2,mrna_1,mrna_2,thr_num,samtools,outdir):
        outdir=outdir+'/assemble_mrna'
        check_dir(outdir)
        #################hisat2 index###################
        cmd='%s -p %s %s %s'%(hisat2_build,thr_num,ref,'.'.join(ref.split('.')[:-1]))
        run_cmd(cmd)
        ##from here
        cmd='%s -p %s -x %s -1 %s -2 %s -S %s/result.sam'%(hisat2,thr_num,'.'.join(ref.split('.')[:-1]),mrna_1,mrna_2,outdir)
        run_cmd(cmd)
        cmd='%s view -bS %s/result.sam>%s/result.bam'%(samtools,outdir,outdir)
        run_cmd(cmd)
        mapped_rate=mapped_result('%s/result.bam'%(outdir),samtools)
        if os.path.exists('%s/result.sam'%(outdir)):os.remove('%s/result.sam'%(outdir))
        return mapped_rate

    def main(self):
        args=self.getopt()
        if contiue_run(args.tag,args.workpath):
            config=Config()
            check_dir(args.outdir)
            samtools = check_software('samtools') if check_software('samtools') else config.samtools
            hisat2=check_software('hisat2') if check_software('hisat2') else config.hisat2
            hisat2_build=check_software('hisat2-build') if check_software('hisat2-build') else config.hisat2_build
            
            mapped_ratem=self.assemble_mrna(args.ref,hisat2_build,hisat2,args.input[0],args.input[1],args.thread,samtools,args.outdir)
            sample_name=args.input[0].split('/')[-1].split('_')[0]
            info='<mrna_mapped\t%s\tmapped_rate:%s\n'%(args.prefix,':'.join(mapped_ratem))
            file_info_get='%s/result.txt'%args.outdir if not args.infofile else args.infofile
            result_file(info,file_info_get)

class Blast_nt():
    def getopt(self):
        parser = argparse.ArgumentParser(
            formatter_class=HelpFormatter, description=usage)
        parser.add_argument(
            'func',choices=['BlastNtMapp'])
        parser.add_argument(
            '-r','--ref',help='assemble file that you need to evalue',dest='ref',type=str,required=True)
        parser.add_argument(
            '-i','--input',help='input file nt_mapped:[nt_database]',dest='input',type=str,default='/mnt/beegfs/data/database/01.nt_nr/FASTA/nt')
        parser.add_argument(
            '-o','--outdir',help='output dir',type=str,dest='outdir',required=True)
        parser.add_argument(
            '-n','--thread',help='set num of thread',type=int,dest='thread',default=10)
        parser.add_argument(
            '--infofile',help='the each step info file\nuse for report',dest='infofile',type=str)
        parser.add_argument('-tag',dest='tag',type=str)
        parser.add_argument('-workpath',dest='workpath',type=str)
        args=parser.parse_args()
        return args

    def blast_nt(self,blastn,nt_database,assemble_file,thr_num,outdir):
        outdir=outdir+'/blast_nt'
        output_prefix=assemble_file.split('/')[-1]
        check_dir(outdir)
        read_200,total_read={},0
        if fa_or_fq(assemble_file)=='fa':
            info=''
            w=open('%s/new_%s'%(outdir,output_prefix),'w')
            assemble_file_info=open(assemble_file,'r') if 'gz' not in assemble_file.split('.') else gzip.open(assemble_file,'rb').decode('utf-8')
            for line in assemble_file_info:
                if line.startswith('>'):
                    if info:
                        count_info=0
                        if len(info)<1000:
                            count_info+=1
                            continue
                        sta,count=0,0
                        for i in range(1000,len(info),1000):
                            w.write(key+'_%s\n'%count)
                            w.write(info[sta:i]+'\n')
                            sta=i
                            count+=1
                            
                        total_read+=count
                        if len(info)-i>200:read_200[key]=info[i:]
                    key=re.split(r'\s+',line.strip())[0]
                    info=''
                    continue
                info+=line.strip()               
            count=0
            if total_read<20000:
                for key,value in read_200.items():
                    
                    w.write(key+'\n'+value+'\n')
                    if count+total_read>=20000:
                        break
                    count+=1
            w.close()
        output_prefix=assemble_file.split('/')[-1].split('.')[0]
        cmd='%s -db %s -query %s/new_%s -out %s/%s.out -evalue 0.00001 -outfmt \'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle\'  -num_threads %d'%(blastn,nt_database,outdir,assemble_file.split('/')[-1],outdir,output_prefix,thr_num)
        run_cmd(cmd)
        result_info,total_mapp_read={},{}
        for line in open('%s/%s.out'%(outdir,output_prefix),'r'):
            info=line.strip().split('\t')
            if info[0] not in total_mapp_read:
                total_mapp_read[info[0]]=0
            elif info[0] in total_mapp_read:
                continue
            species_mapp=re.split(r'\s+',info[-1])[1] if 'PREDICTED:' not in info[-1] else re.split(r'\s+',info[-1])[2]
            if species_mapp not in result_info:result_info[species_mapp]={}
            if info[0] not in result_info[species_mapp]:result_info[species_mapp][info[0]]=0
        total_mapp_num=len(total_mapp_read)
        total_mapp_read.clear()
        for key,value in result_info.items():
            result_info[key]=[float(len(value))/float(total_mapp_num),len(value)]
        total_read_num=total_read if total_read>=20000 else total_read+count
        return result_info,total_mapp_num,total_read_num

    def main(self):
        args=self.getopt()
        if contiue_run(args.tag,args.workpath):
            config=Config()
            check_dir(args.outdir)
            #here to start
            blastn=check_software('blastn') if check_software('blastn') else config.blastn
            check_dir(args.outdir)
            dict_nt,total_mapp_num,total_read_num=self.blast_nt(blastn,args.input,args.ref,args.thread,args.outdir)
            info='<blast_nt\ttotal_read_num:%s\ttotal_mapp_num:%s'%(total_read_num,total_mapp_num)
            result_nt=sorted(dict_nt.items(),key=lambda item:item[1][1],reverse=True)
            count,stat=0,0
            for key,value in result_nt:
                stat+=value[0]
                if count<5:
                    info+='\t%s:%d,%f'%(key,value[1],value[0])
                else:
                    if stat>0.8:
                        break
                    elif value[0]<0.05:
                        break
                    else:
                        info+='\t%s:%d,%f'%(key,value[1],value[0])
                count+=1
            info+='\n'
            file_info_get='%s/result.txt'%args.outdir if not args.infofile else args.infofile
            result_file(info,file_info_get)


class Busco_evo():
    def getopt(self):
        parser = argparse.ArgumentParser(formatter_class=HelpFormatter, description=usage)
        parser.add_argument('func',choices=['BuscoEvaluation'])
        parser.add_argument('-o', '--output', help='output dir', dest='output', type=str)
        parser.add_argument('-g','--genome', help='input genome file', dest='genome', type=str)
        parser.add_argument('--species', help='please input species name', dest='species',type=str,nargs='*')
        parser.add_argument('--augustus_dataset',help='please input augustus dataset',dest='augustus_dataset',type=str)
        parser.add_argument('--busco_dataset',help='please input busco dataset',dest='busco_dataset',type=str)
        parser.add_argument('--cpu', help='set cpu num', dest='cpu', default=40,type=int)
        parser.add_argument('--vf', help='set Virtual Memory', dest='vf', default='2g',type=str)
        parser.add_argument('--infofile',help='the each step info file use for report',dest='infofile',type=str)
        parser.add_argument('-tag',dest='tag',type=str)
        parser.add_argument('-workpath',dest='workpath',type=str)
        args = parser.parse_args()
        if not args.output:
            print('config must be given!!!!!![-o]')
            sys.exit(1)
        elif not args.genome:
            print("genome file must be given!!!!![-g]")
            sys.exit(1)
        elif not args.species and not (args.busco_dataset and args.augustus_dataset):
            print('species[--species] or (augustus_dataset[--augustus_dataset] and busco_dataset[--busco_dataset]) must be given!!!')
            sys.exit(1)

        return args

    def run_gene_mode(self,genome,output,busco,cpu,august_dataset,busco_dataset,env):
        out_prefix=genome.split('/')[-1].split('.')[0]
        cmd='gzip -dc %s >%s && '%(genome,os.path.basename(genome).rstrip('.gz'))
        if '.gz' in genome:run_cmd(cmd)
        check_dir(output)
        outdir=output+'/1.busco_geno'
        check_dir(outdir)
        genome=os.path.abspath(genome)
        cmd='source %s && cd %s && python %s --in %s --cpu %s --out %s --force --mode genome --lineage_path %s -sp %s --blast_single_core && cd -'%(env,outdir,busco,genome,cpu,out_prefix,busco_dataset,august_dataset)
        print(cmd)
        run_cmd(cmd)
        
        count,info_busco=0,'<busco_evo'
        for line in open('%s/run_%s/short_summary_%s.txt'%(outdir,out_prefix,out_prefix),'r'):
            if line.startswith('#') or line.startswith('\n') or line.startswith(' '):
                continue
            info = line.strip().split('\t')
            if len(info)==1:
                rate_all=re.split(r'[[|\]|,]',info[0])
                rate_all.remove('')
            elif len(info)==2:
                rate_write=rate_all[count].split(':')[1] if 'n' not in rate_all[count] else '100%'
                info_busco+='\t%s,%s,%s'%(info[1],info[0],rate_write)
                count+=1
        info_busco+='\n'
        return info_busco

    def busco_august_find(self,species):
        config=Config()
        name_file,nodes_file,augustus_file=config.name_dmp,config.node_dmp,config.augustus_info
        busco_dict=config.busco_dict
        augustus_dataset=config.augustus_dataset
        dict_augustus_level={}
        for line in open(augustus_file,'r'):
            info=line.strip().split(';')
            dict_augustus_level[info[0]]=OrderedDict()
            for i in info[1].split(',')[:-1]:
                dict_augustus_level[info[0]][i.split(":")[0]]=i.split(':')[1]

        augustus,busco,id_name='','',''
        for line in open(name_file,'r'):
            info=line.strip().split('|')
            name=info[1].strip().lower()
            if name == species.lower():
                id_name=info[0].strip()
                if id_name in augustus_dataset:
                    augustus=augustus_dataset[id_name]
                elif id_name in busco_dict:
                    busco=busco_dict[id_name]
                break
        if id_name:
            dict_nodes_high={}
            for line in open(nodes_file,'r'):
                info=line.strip().split('|')
                query_id=info[0].strip()
                advanced_id=info[1].strip()
                if query_id not in dict_nodes_high:dict_nodes_high[query_id]=advanced_id
            if busco and augustus:
                return augustus,busco
            else:
                count=0
                if augustus:
                    while count<20000:
                        count+=1
                        id_name=dict_nodes_high[id_name]
                        if id_name in busco_dict:
                            busco=busco_dict[id_name]
                            break
                else:
                    while count<20000:
                        count+=1
                        #######call augustus dataset#####
                        if not augustus:
                            for key,value in dict_augustus_level.items():
                                if id_name in value:
                                    augustus=key
                        id_name=dict_nodes_high[id_name]
                        if id_name in busco_dict:
                            busco=busco_dict[id_name]
                        if augustus and busco:
                            break
                if count>=20000:
                    print('sorry!we get no ideal,may be you can find the %s busco and augustus dataset by yourself!'%species)
        else:
            print('%s is wrong Latin name,please check your input name'%species)
            sys.exit(0)
        return config.augustus_dataset[augustus],busco

    def main(self):
        config=Config()
        args=self.getopt()
        if contiue_run(args.tag,args.workpath):
            if args.species:
                augustus,busco=self.busco_august_find(' '.join(args.species))
                info_busco=self.run_gene_mode(args.genome,args.output,config.busco,args.cpu,augustus,busco,config.env)
            elif args.busco_dataset and args.augustus_dataset:
                info_busco=self.run_gene_mode(args.genome,args.output,config.busco,args.cpu,args.augustus_dataset,args.busco_dataset,config.env)
            file_info_get='%s/result.txt'%args.output if not args.infofile else args.infofile
            result_file(info_busco,file_info_get)

class Circos():
    def get_yeast(self,ref_len, select_chrom, out_put):
        w=open(out_put,'w')
        col_lst=[]
        with open(ref_len,'r') as f:
            for line in f:
                info=line.strip().split('\t')
                if info[0] not in select_chrom:continue
                w.write('chr\t-\t%s\t%s\t0\t%s\t%s,%s,%s\n'%(info[0],info[0],info[1],random.randint(0,255),random.randint(0,255),random.randint(0,255)))
        w.close()

    def variation_circos(self,vcf_file,color,chroms,outdir,type_out,block,fa_len):
        w_homo=open('%s/%s_homo.circos'%(outdir,type_out),'w')
        w_mix=open('%s/%s_mix.circos'%(outdir,type_out),'w')
        chrs_all,count=defaultdict(dict),0
        for line in open(vcf_file,'r'):
            if line.startswith("#"): continue
            if "PASS" not in line: continue
            info=line.strip().split('\t')
            chrom , loc, GT=info[0], int(info[1]), info[-1].split(':')[0]
            if chrom not in chroms:continue
            if chrom not in chrs_all:
                if chrs_all:
                    for i in range(len(lst_homo)):
                        w_homo.write('{0}\t{1}\t{2}\t{3}\n'.format(chrs,block*i,block*(i+1),lst_homo[i]))
                    for i in range(len(lst_mix)):
                        w_mix.write('{0}\t{1}\t{2}\t{3}\n'.format(chrs,block*i,block*(i+1),lst_mix[i]))
                chrs=chrom
                block_num = int(fa_len[chrs]/block)
                lst_homo = [0] * (block_num + 1)
                lst_mix = [0] * (block_num + 1)
                chrs_all[chrom]=0
            if GT=='0/0':continue
            if GT=='1/1':
                pos=int(loc/block)
                lst_homo[pos]+=1
            else:
                pos=int(loc/block)
                lst_mix[pos]+=1
        if chrs in chroms:
            for i in range(len(lst_homo)):
                w_homo.write('{0}\t{1}\t{2}\t{3}\n'.format(chrs,block*i,block*(i+1),lst_homo[i]))
            for i in range(len(lst_mix)):
                w_mix.write('{0}\t{1}\t{2}\t{3}\n'.format(chrs,block*i,block*(i+1),lst_mix[i]))
        w_mix.close()
        w_homo.close()

    def get_depth_dic(self,depth_f, color,chroms,block,fa_len,outdir,type_out):
        w=open('%s/depth_%s.circos'%(outdir,type_out),'w')
        chrs_all=defaultdict(dict)
        max_depth=0
        for line in open(depth_f,'r'):
            info=line.strip().split("\t")
            chrom, loc, dep = info[0],int(info[1]),int(info[2])
            if chrom not in chroms:continue
            if chrom not in chrs_all:
                if chrs_all:
                    sort_lst=sorted(lst)
                    end_lst=int(len(sort_lst)*0.9)
                    value_max9=sort_lst[end_lst]
                    max_depth=value_max9 if value_max9>max_depth else max_depth
                    for i in range(len(lst)):
                        w.write('{0}\t{1}\t{2}\t{3}\tfill_color={4}\n'.format(chrs,block*i+1,block*(i+1),lst[i],color))
                    lst=[]
                chrs=chrom
                chrs_all[chrom]=0
                block_num = int(fa_len[chrs]/block)
                lst = [0] * (block_num + 1)
            pos=int(loc/block)
            lst[pos]+=dep
        sort_lst=sorted(lst)
        end_lst=int(len(sort_lst)*0.9)
        value_max9=sort_lst[end_lst]
        max_depth=value_max9 if value_max9>max_depth else max_depth
        for i in range(len(lst)):
            w.write('{0}\t{1}\t{2}\t{3}\tfill_color={4}\n'.format(chrs,block*i+1,block*(i+1),lst[i],color))
        w.close()
        return max_depth

    def depth_conf_result(self,max_depth_2,max_depth_3,outdir):
        w=open('%s/depth_2.conf'%outdir,'w')
        w.write('''type        = histogram
file        = ./circos/depth_2.circos
r1          = 0.83r
r0          = 0.73r
min         = 0 
max         = %s
extend_bin  = no
fill_color  = 179,30,68
color       = 179,30,68
thickness   = 0 
orientation = out 
            '''%max_depth_2)
        w.close()
        w=open('%s/depth_3.conf'%outdir,'w')
        w.write('''type        = histogram
file        = ./circos/depth_3.circos
r1          = 0.71r
r0          = 0.61r
min         = 0 
max         = %s
extend_bin  = no
fill_color  = 239,59,96
color       = 239,59,96
thickness   = 0 
orientation = out 
            '''%max_depth_3)
        w.close()

    def GC_density(self,infile,chroms,block,fa_len,outdir,color):
        out_dic,lis_base,chrs={},[],''
        max_GC,min_GC=0,100000000000000000000000
        result={}
        for line in open(infile):
            if line.startswith('>'):
                if chrs and  lis_base:
                    if chrs in chroms:
                        lenth=len(lis_base[0])
                        block_num = int(fa_len[chrs]/block)
                        result[chrs]=[]
                        for i in range(block_num+1):
                            sta_index=block*i/lenth
                            end_index=block*(i+1)/lenth+1
                            sta_pos=0
                            end_pos=block
                            sequence=''.join(lis_base[sta_index:end_index+1])[sta_pos:end_pos]
                            result[chrs].append(float((sequence.count('G')+sequence.count('C')))/float(block))
                        num_max,num_min=max(result[chrs]),min(result[chrs])
                        max_GC=num_max if max_GC<num_max else max_GC
                        min_GC=num_min if max_GC>num_min else min_GC
                        lis_base=[]
                    else:
                        lis_base=[]
                else:
                    lis_base=[]
                chrs=re.split(r'\s+',line.lstrip('>'))[0]
                continue
            lis_base.append(line.strip())
        if chrs in chroms:
            lenth=len(lis_base[0])
            block_num = int(fa_len[chrs]/block)
            result[chrs]=[]
            for i in range(block_num+1):
                sta_index=block*i/lenth
                end_index=block*(i+1)/lenth+1
                sta_pos=0
                end_pos=block
                sequence=''.join(lis_base[sta_index:end_index+1])[sta_pos:end_pos]
                result[chrs].append(float((sequence.count('G')+sequence.count('C')))/float(block))
            num_max,num_min=max(result[chrs]),min(result[chrs])
            max_GC=num_max if max_GC<num_max else max_GC
            min_GC=num_min if max_GC>num_min else min_GC
        lis_base=[]
        output=open('%s/GC.circos'%outdir,'w')
        for k, v in result.items():
            count = 0
            for i in v:
                output.write("{0}\t{1}\t{2}\t{3}\tfill_color={4}\n".format(k, count*block+1, (count+1)*block, i,color))
                count += 1
        output.close()
        return max_GC,min_GC

    def busco_circos(self,infile,outdir):
        w_single=open('%s/single.circos'%outdir,'w')
        w_duplic=open('%s/duplic.circos'%outdir,'w')
        for line in open(infile):
            if line.startswith('#'):continue
            info=line.strip().split('\t')
            if info[1]=='Complete':
                w_single.write('\t'.join(info[2:5])+'\tfill_color=blue\n')
            elif info[1]=='Duplicated':
                w_duplic.write('\t'.join(info[2:5])+'\tfill_color=red\n')
        w_single.close()
        w_duplic.close()

    def config_file(self,chroms,outdir,block_num):
        w=open('%s/chromosomes.conf'%outdir,'w')
        w.write('''chromosomes_reverse = /genome|mito/
chromosomes_units           = %d
chromosomes_display_default = yes
chromosomes_order  = %s'''%(block_num,','.join(chroms)))
        w.close()
        w=open('%s/ideogram.conf'%outdir,'w')
        w.write('''<ideogram>
#<spacing>
#default = 5u
#</spacing>
<spacing>

default = 0.0025r
break   = 0.5r
<pairwise %s;%s> 
        spacing = 8r
</pairwise>
</spacing>
<<include ideogram.position.conf>>
<<include ideogram.label.conf>>
<<include bands.conf>>

</ideogram>'''%(chroms[0],chroms[-1]))
        w.close()

    def conf_GC(self,outdir,max_GC,min_GC):
        w=open('%s/GC.conf'%outdir,'w')
        w.write('''type        = histogram
file        = ./circos/GC.circos
r1          = 0.95r
r0          = 0.85r
min         = %s 
max         = %s
extend_bin  = no
fill_color  = 122,63,97
color       = 122,63,97
thickness   = 0 
orientation = out 
            '''%(min_GC,max_GC))
        w.close()

    def picture_draw_info(self,png_file):
        img = Image.open(png_file)

        draw = ImageDraw.Draw(img)
        ttfront = ImageFont.truetype('/usr/share/fonts/wqy-zenhei/wqy-zenhei.ttc',90)
        draw.text((1480,280),'A',fill=(25,25,25),font=ttfront)
        draw.text((1480,420),'B',fill=(25,25,25),font=ttfront)
        draw.text((1480,570),'C',fill=(25,25,25),font=ttfront)
        draw.text((1480,700),'D',fill=(25,25,25),font=ttfront)
        draw.text((1480,830),'E',fill=(25,25,25),font=ttfront)
        draw.text((1480,970),'F',fill=(25,25,25),font=ttfront)
        draw.text((1480,1110),'G',fill=(25,25,25),font=ttfront)
        img.show()
        img.save(png_file)

    def circos(self):
        args=self.getopt()
        config=Config()
        outdir=args.outdir+'/circos'
        check_dir(outdir)   
        color={'snp_vcf': '240,139,61','indel_vcf': '207,49,113','depth_2': '179,30,68','depth_3':'239,59,96','GC':'122,63,97'}
        fa = pysam.FastaFile(args.references)
        fa_len = dict(zip(fa.references, fa.lengths))
        chroms=self.gene_chr(fa_len)
        args.block=total_point_num_get(args.references)
        self.config_file(list(chroms.keys()),outdir,args.block)
        self.get_yeast(args.references+'.fai', chroms.keys(),'%s/yeast.txt'%outdir)
        max_GC,min_GC=self.GC_density(args.references,chroms,args.block,fa_len,outdir,color['GC'])
        self.conf_GC(outdir,max_GC,min_GC)
        self.variation_circos(args.snp_vcf,color['snp_vcf'],chroms,outdir,'snp',args.block,fa_len)
        self.variation_circos(args.indel_vcf,color['indel_vcf'],chroms,outdir,'indel',args.block,fa_len)
        max_depth_2=self.get_depth_dic(args.depth_2, color['depth_2'],chroms,args.block,fa_len,outdir,'2')
        max_depth_3=self.get_depth_dic(args.depth_3, color['depth_3'],chroms,args.block,fa_len,outdir,'3')
        self.depth_conf_result(max_depth_2,max_depth_3,outdir)
        self.busco_circos(args.busco_tsv,outdir)
        os.system('cp -r %s/* %s'%(config.circos_conf,outdir))
        os.system('%s -conf %s/circos.conf'%(config.circos_soft,outdir))
        self.picture_draw_info('circos.png')


    def getopt(self):
        parser = argparse.ArgumentParser(
            formatter_class=HelpFormatter, description=usage)
        parser.add_argument(
            'func',choices=['Circos'])
        parser.add_argument('-r','--references',dest='references',help='input ref gene',type=str,required=True)
        parser.add_argument('-o','--outdir',dest='outdir',help='outdir',type=str,required=True)
        parser.add_argument('--block',dest='block',help='input the size of block',type=int)
        parser.add_argument('--snp_vcf',dest='snp_vcf',help='input snp_vcf file',type=str,required=True)
        parser.add_argument('--indel_vcf',dest='indel_vcf',help='input indel_vcf',type=str,required=True)
        parser.add_argument('--depth_2',dest='depth_2',help='input depth_2 file',type=str,required=True)
        parser.add_argument('--depth_3',dest='depth_3',help='input depth_3 file',type=str,required=True)
        parser.add_argument('--busco_tsv',dest='busco_tsv',help='input busco_tsv file',type=str,required=True)
        parser.add_argument('-tag',dest='tag',type=str)
        parser.add_argument('-workpath',dest='workpath',type=str)
        args = parser.parse_args()
        return args

    def gene_chr(self,fa_len):
        chr_len_sorted=sorted(fa_len.items(),key=lambda item:item[1],reverse=True)
        total_len=sum(list(fa_len.values()))
        
        judge_chr=[total_len/100,total_len/200,total_len/2]
        count_len,chr_name_all,check=0,OrderedDict(),0
        ###########scaffold###########
        for key,value in chr_len_sorted:
            if value>judge_chr[0]:
                chr_name_all[key]=0
                count_len+=value
        if count_len<judge_chr[2]:
            check=1
        if check:
            #######contig###########
            count_len,chr_name_all,check=0,OrderedDict(),0
            if chr_len_sorted[0][1]<judge_chr[1]:
                for key,value in chr_len_sorted:
                    if count_len>9:break
                    chr_name_all[key]=0
                    count_len+=1
            else:
                for key,value in chr_len_sorted:
                    if value>judge_chr[1]:
                        chr_name_all[key]=0
                    else:
                        break
        return chr_name_all

    def main(self):
        args=self.getopt()
        if contiue_run(args.tag,args.workpath):
            self.circos()

usage = '''
Author : chenqi
Email  : chenqi@gooalgene.com
Date   : 2018-11-4
Link   : 
Version: v1.0
Description:
    Assemble Evaluation!
Example:
    python %s <functions> [options]... 
<functions>
    Reads3Stat                     Stat 3 mapped result
    Reads2Stat                     Stat 2 mapped result
    MrnaMapped                     Stat mrna mapped result
    BlastNtMapp                    Stat nt mapped result
    BuscoEvaluation                Busco evaluetion
    Circos                         Draw circos picture
    ReportHtml                     Result Report
''' % (__file__[__file__.rfind(os.sep) + 1:])

def main():
    if len(sys.argv)==1 or sys.argv[1]=='-h' or sys.argv[1]=='--help':
        print(usage)
        sys.exit(0)
    if sys.argv[1]=='Reads3Stat':
        process=Assemble3_stat()
        process.main()
    elif sys.argv[1]=='Reads2Stat':
        process=Assemble2_stat()
        process.main()
    elif sys.argv[1]=='MrnaMapped':
        process=Assemble_mrna()
        process.main()
    elif sys.argv[1]=='BlastNtMapp':
        process=Blast_nt()
        process.main()
    elif sys.argv[1]=='BuscoEvaluation':
        process=Busco_evo()
        process.main()
    elif sys.argv[1]=='ReportHtml':
        process=Report()
        process.main()
    elif sys.argv[1]=='Circos':
        process=Circos()
        process.main()
    else:
        print('''
sorry! we didn\'t create the function now,please email to my mailbow and may be i will create it

please print \'-h\' or \'--help\' for help,chose function fitst!!!
    
    thank you so much for using it

email:chenqi@gooalgene.com
            ''')
        sys.exit(0)

if __name__=='__main__':
    try:
        t1 = time.time()
        time1 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t1))
        print('Start at : ' + time1)
        main()
        t2 = time.time()
        time2 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t2))
        print('End at : ' + time2)
        t3 = t2 - t1
        print('Spend time: ' + fmt_time(t3))
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! >_< See you!\n")
        sys.exit(0)

