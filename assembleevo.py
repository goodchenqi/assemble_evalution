#!/thinker/storage/software/Genome/anaconda2/bin/python
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
try:
    import pysam
    import fire
    import pandas as pd
except:
    print('import error: fire and panda can\'t import ,maybe you can run commend \"export PYTHONPATH=/thinker/storage/software/lib/python2/lib/python2.7/site-packages:/thinker/storage/software/lib/python2/lib64/python2.7/site-packages/networkx/networkx:/thinker/storage/software/lib/python2/lib64/python2.7/site-packages:$PYTHONPATH && python /thinker/storage/software/lib/python2/lib64/python2.7/site-packages/networkx/setup.py build_ext --inplace --force\"')
    sys.exit(0)
from collections import defaultdict
from config import *
#####Description####
usage = '''
Author : chenqi
Email  : chenqi@gooalgene.com
Date   : 2018-9-3
Version: v1.0
Description:
    
Example: 
    python %s 

''' % (__file__[__file__.rfind(os.sep) + 1:])

#####HelpFormat#####
class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

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
    total_line,total_contig,total_scaffold,base_info,sequence=0,[],[],[0,0,0,0,0],[]
    for line in open(file):
        if line.startswith('>'):
            if sequence:
                count=0
                scaffold=''.join(sequence)
                for i in ['A','T','C','G','N']:
                    base_info[count]+=scaffold.count(i)
                total_scaffold.append(len(scaffold))
                for contig in re.split(r'[n|N]+',scaffold):
                    total_contig.append(len(contig))
            sequence=[]
            continue
        sequence.append(line.strip().upper())
    count=0
    scaffold=''.join(sequence)
    for i in ['A','T','C','G','N']:
        base_info[count]+=scaffold.count(i)
    total_scaffold.append(len(scaffold))
    for contig in re.split(r'[n|N]+',scaffold):
        total_contig.append(len(contig))
    sequence=[]

    total_contig=sorted(total_contig,reverse=True)
    total_scaffold=sorted(total_scaffold,reverse=True)
    contig_max_length,contig_total_length,scaffold_max_length,scaffold_total_length=total_contig[0],sum(total_contig),total_scaffold[0],sum(total_scaffold)
    total_base,total_base_contig=sum(base_info),sum(base_info)-base_info[-1]
    N_info_scaffold=[total_base*0.9,total_base*0.8,total_base*0.7,total_base*0.6,total_base*0.5]
    N_info_contig=[total_base_contig*i for i in [0.9,0.8,0.7,0.6,0.5]]
    count,check,count_contig,count_scaffold=0,1,0,0
    num_100_contig,num_2k_contig,num_100_scaffold,num_2k_scaffold=0,0,0,0
    contig_N_seq,scaffold_N_seq,contig_i,scaffold_i=[0,0,0,0,0],[0,0,0,0,0],1,1
    while 1:
        if count>=len(total_contig):
            break
        elif count>=len(total_scaffold):
            check=0
        else:
            count_contig+=total_contig[count]
            if count_contig>=100:num_100_contig+=1
            if count_contig>=2000:num_2k_contig+=1
            if check:
                count_scaffold+=total_scaffold[count]
                if count_scaffold>=100:num_100_scaffold+=1
                if count_scaffold>=2000:num_2k_scaffold+=1
                if N_info_scaffold and count_scaffold>=N_info_scaffold[-1]:
                    scaffold_N_seq[-scaffold_i]=[total_scaffold[count],count+1]
                    N_info_scaffold.pop()
                    scaffold_i+=1
            if N_info_contig and count_contig>=N_info_contig[-1]:
                contig_N_seq[-contig_i]=[total_contig[count],count+1]
                N_info_contig.pop()
                contig_i+=1
            count+=1

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

def mapped_result(file,samtools):
    #first you need to check the software (samtools)
    if not os.path.exists(file):
        print('sorry!may be the depth %s didn\'t exists,please check it!!'%file)
        sys.exit(0)
    mapped_rate=[]
    for i in os.popen('%s flagstat %s'%(samtools,file)): 
        if re.findall(r'mapped',i):
            mapped_rate.append(re.split(r'[(|)]',i.strip())[1].split(':')[0].strip())
        elif re.findall(r'properly paired',i):
            info=re.split(r'[(|)]',i.strip())[1].split(':')[0].strip()
            if info=='N/A' or info=='-nan%':
                return mapped_rate
            else:
                mapped_rate.append(info)
                return mapped_rate

def ref_deal(ref):
    input_type=fa_or_fq(ref)
    total_gc,total_len,gc_len={},0,0
    if input_type=='fa':
        for line in open(ref,'r'):
            if line.strip().startswith('>'):
                if total_len!=0:
                    gc_rate=float(gc_len)/float(total_len)
                    total_gc[key]=gc_rate
                key=re.split(r'\s+',line)[0][1:]
                if key not in total_gc:total_gc[key]=0
                continue
            total_len+=len(line.strip())
            gc_len+=(line.upper().count('G')+line.upper().count('C'))
    return total_gc

def gc_depth(depth,gc,draw_soft,output_prefix,outdir):
    total_result={}
    w=open('%s/%s.txt'%(outdir,output_prefix),'w')
    for key,value in depth.items():
        kes=str(value)+'\t%.6f'%(gc[key])
        if kes not in total_result:total_result[kes]=0
        total_result[kes]+=1
    depth.clear()
    for key,value in total_result.items():
        w.write(str(value)+'\t'+key+'\n')
    w.close()
    cmd='Rscript %s %s/%s.txt %s/%s.png'%(draw_soft,outdir,output_prefix,outdir,output_prefix)
    run_cmd(cmd)
    return '%s/%s.png'%(outdir,output_prefix)

def resume_info(file):
    count,resume_list=0,[]
    if os.path.exists(file):
        for line in open(file,'r'):
            resume_list.append(int(line.strip()))
        try:
            resume=max(resume_list)
            return resume
        except:
            return 0
    else:
        return 0

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
        parser.add_argument(
            '--species',help='input species name',type=str,dest='species')
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

    def report_result(self,result_file,outdir,ref_file,species_name):
        ##################get each step info###################
        stat_result=stat_N50(ref_file)
        stat_info,base_info,count='','',0
        for i in range(90,40,-10):
            stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">N%d</td><td style="text-align: center;">%d</td><td style="text-align: center;">%d</td><td style="text-align: center;">%d</td><td style="text-align: center;">%d</td></tr>'''%(i,stat_result[2][count][0],stat_result[2][count][1],stat_result[1][count][0],stat_result[1][count][1])
            count+=1
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Max length</td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td></tr>'''%(stat_result[3][0][0],stat_result[3][0][1])
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Total length</td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td></tr>'''%(stat_result[3][1][0],stat_result[3][1][1])
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Number(>=100bp)</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td></tr>'''%(stat_result[0][0],stat_result[0][2])
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Number(>=2kb)</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td></tr>'''%(stat_result[0][1],stat_result[0][3])
        stat_info+='''
            <tr style="height: 5px;">
    <td style="text-align: center;">Total number</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td><td style="text-align: center;"></td><td style="text-align: center;">%d</td></tr>'''%(stat_result[0][4],stat_result[0][5])
        base_total=['A','T','C','G','N','GC','Total']
        for i in range(len(base_total)):
            base_info+='''
    <tr style="height: 5px;">
    <td style="text-align: center;">{0}</td><td style="text-align: center;">{1}</td><td style="text-align: center;">{2: .2f}%</td></tr>'''.format(base_total[i],stat_result[-1][i],(float(stat_result[-1][i])/float(stat_result[-1][-1]))*100)
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
                info_2_hm+='''
            <tr style="height: 5px;">
        <td style="text-align: center;">{0}</td><td style="text-align: center;">{1: .2f}%</td><td style="text-align: center;">{2: .2f}%</td><td style="text-align: center;">{3: .2f}%</td><td style="text-align: center;">{4: .2f}%</td></tr>'''.format(info[1],float(info[5].split(':')[1])*100,float(info[7].split(':')[1])*100,float(info[4].split(':')[1])*100,float(info[6].split(':')[1])*100)

                png_2+="<li href=\"{0}\"><img src=\"{1}\" /><span>{2}</span></li>".format(info[-1].split(':')[1],info[-1].split(':')[1],info[1])
                if not png2_check:
                    png_2_random=info[-1].split(':')[1]
                png2_check=1

        #############3 mapp info#############
            elif line.startswith('<3_mapped'):
                info_3+='''
            <tr style="height: 5px;">
    <td width="300px" style="text-align: center;">{0}</td><td width="700px" style="text-align: center;">{1}</td><td width="700px" style="text-align: center;">{2: .2f}</td><td width="700px" style="text-align: center;">{3}</td><td width="700px" style="text-align: center;">{4}</td><td width="700px" style="text-align: center;">{5}</td><td width="700px" style="text-align: center;">{6}</td></tr>'''.format(info[1],info[2].split(':')[1],float(info[3].split(':')[1]),info[4].split(':')[1],info[4].split(':')[2],info[4].split(':')[3],info[4].split(':')[4])

                png_3+="<li href=\"{0}\"><img src=\"{1}\" /><span>{2}</span></li>".format(info[-1].split(':')[1],info[-1].split(':')[1],info[1])

                if not png3_check:
                    png_3_random=info[-1].split(':')[1]
                png3_check=1

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
                        <li><a href="#b2">2.组装结果统计</a></li>
                        <li><a href="#c3">3.数据正确性评估</a></li>
                        <li><a href="#d4">4.序列一致性评估</a></li>
                        <li><a href="#e5">5.组装完整性评估</a> 
                            <ul class="nav">
                                <li><a href="#e51" class="js-active-bg">5.1 保守基因完整性评估</a></li>
                            </ul>
                            <ul class="nav">
                                <li><a href="#e52" class="js-active-bg">5.1 基因区完整性评估</a></li>
                            </ul>
                        </li>
                        <li><a href="#f6">6.准确性及杂合度评估</a>
                        </li>
                        <li><a href="#g7">7.基因组圈图</a></li>
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
                    <h2 id="b2">2. 组装结果统计</h2>
                    <p class="paragraph">组装得到基因组序列的基本情况统计见下表:</p>
                    <div class="studyUnit">
                        <p class="name">表1. 组装结果序列统计</p>
<table class="table table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
    <td style="text-align: center;"></td><td style="text-align: center;">contig</td><td style="text-align: center;">scaffold</td></tr></table>
<table class="table table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
    <td style="text-align: center;"></td><td style="text-align: center;">Length(bp)</td><td style="text-align: center;">Number</td><td style="text-align: center;">Size(bp)</td><td style="text-align: center;">Number</td></tr>{3}</table></div>
                    <p class="paragraph">注：</p>
                    <p class="paragraph">(1) N90 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 90%的时候，那一条序列的长度（bp）。</p>
                    <p class="paragraph">(2) N80 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 80%的时候，那一条序列的长度（bp）。</p>
                    <p class="paragraph">(3) N70 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 70%的时候，那一条序列的长度（bp）。</p>
                    <p class="paragraph">(4) N60 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 60%的时候，那一条序列的长度（bp）。</p>
                    <p class="paragraph">(5) N50 length：将组装好的序列按照从长到短进行排列并进行累加，当累加长度达到总长的 50%的时候，那一条序列的长度（bp）。</p>
                    <p class="paragraph">(6) Max length：组装得到的最长的序列的长度（bp）。</p>
                    <p class="paragraph">(7) Total  length：组装得到的序列总长（bp）。</p>
                    <p class="paragraph">(8) Number(>=100bp)：组装得到长度大于100bp的序列条数。</p>
                    <p class="paragraph">(9) Number(>=2kb)：组装得到长度大于2kb的序列条数。</p>
                    <p class="paragraph">(10) Total number：组装得到的序列条数。</p>
                    <p class="paragraph">对基因组序列进行A、T、C、G的含量统计分析，统计结果见下表:</p>
                        <div class="studyUnit">
                        <p class="name">表2.碱基组成分布情况统计</p>
<table class="table table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
    <td style="text-align: center;">Base</td><td style="text-align: center;">Number</td><td style="text-align: center;">% of genome</td></tr>{4}</table></div>
    </div>
            <div class="tMainBody ">
                <h2 id="c3">3 数据正确性评估</h2>
                    <p class="paragraph">为了确认组装结果属于目标物种，以1000 bp为步长将基因组序列打断，把打断后的序列通过软件Blast比对到NCBI核苷酸数据库（NT库），表3展示了比对次数排名前五的物种。</p>
                    <div class="studyUnit">
                        <p class="name">表3.打断后序列与NT库比对情况部分展示</p>
<table class="table table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
    <td style="text-align: center;">Species</td><td style="text-align: center;">Blast number</td><td style="text-align: center;">Total blast number</td><td style="text-align: center;">Total (%)</td></tr>{5}</table>
            </div>   
        </div>       
            <div class="tMainBody ">
                <h2 id="d4">4 序列一致性评估</h2>
                    <p class="paragraph">分别使用软件bwa和blasr将二、三代数据比对到XXX基因组上，统计比对率以及覆盖率，这两个值越高则认为组装结果与reads的一致性越高，组装的效果也相应的更好。</p>
                    <div class="studyUnit">
                        <p class="name">表4.二代reads比对结果统计</p>
<table class="table left table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
    <td width="300px" style="text-align: center;">species</td><td width="700px" style="text-align: center;">Mapping rate(%)</td><td width="700px" style="text-align: center;">Paired mapping rate(%)</td><td width="700px" style="text-align: center;">Average sequencing depth</td><td width="700px" style="text-align: center;">Coverage(%)</td><td width="700px" style="text-align: center;">Coverage at least 4X(%)</td><td width="700px" style="text-align: center;">Coverage at least 10X(%)</td><td width="700px" style="text-align: center;">Coverage at least 20X(%)</td></tr>{6}</table>
                    </div>
                    <p class="paragraph"></p>
                    <div class="studyUnit">
                        <p class="name">表5.三代reads比对结果比对结果统计</p>
<table class="table left table-hover table-striped table-bordered" style="width: 100%">
        <tr style="height: 5px;">
    <td width="300px" style="text-align: center;">species</td><td width="700px" style="text-align: center;">Mapping rate(%)</td><td width="700px" style="text-align: center;">Average sequencing depth</td><td width="700px" style="text-align: center;">Coverage(%)</td><td width="700px" style="text-align: center;">Coverage at least 4X(%)</td><td width="700px" style="text-align: center;">Coverage at least 10X(%)</td><td width="700px" style="text-align: center;">Coverage at least 20X(%)</td></tr>{7}</table>
                    </div>
            <p class="paragraph"></p>      
            <p class="paragraph">同时绘制GC与深度分布关联分析统计图，展示、评估测序的均匀性，图中横坐标表示GC含量，纵坐标表示覆盖深度，右侧展示了contig覆盖深度分布，上方展示了GC含量分布，中间大图为根据contigs的GC分布及覆盖深度信息绘制的散点图，其中颜色深浅用于反映散点图中点的密度。</p>
            <div class="albumSlider lane">
                                <div class="selectImg">
                                    <span class="img_name">S60_per_base_sequence_quality_1.png</span>
                                    <div class="fright">
                                        <input type="text">
                                        <input value="搜索" class="confirm" type="button">
                                        <input value="复位" class="showall" type="button">
                                    </div>
                                </div>
                                <div class="fullview"><a href="{8}"><img src="{9}" /></a></div>
                                <div class="slider">
                                    <div class="button movebackward" title="向上"></div>
                                    <div class="imglistwrap">
                                        <ul class="imglist">          
                                            {10}
                                        </ul>
                                    </div>
                                    <div class="button moveforward" title="向下"></div>
                                </div>
                    </div>
                    <p class="name">图1. GC含量与覆盖深度（Depth）关联分析统计图（二代reads比对)</p>
                    <p class="paragraph"></p>
                    <div class="studyUnit">
                        <div class="studyUnit">
                        <p class="paragraph"> </p>
                        <div class="albumSlider lane">
                                <div class="selectImg">
                                    <span class="img_name">S60_per_base_sequence_quality_1.png</span>
                                    <div class="fright">
                                        <input type="text">
                                        <input value="搜索" class="confirm" type="button">
                                        <input value="复位" class="showall" type="button">
                                    </div>
                                </div>
                                <div class="fullview"><a href="{11}"><img src="{12}" /></a></div>
                                <div class="slider">
                                    <div class="button movebackward" title="向上"></div>
                                    <div class="imglistwrap">
                                        <ul class="imglist">          
                                            {13}
                                        </ul>
                                    </div>
                                    <div class="button moveforward" title="向下"></div>
                                </div>
                    </div>
                    <p class="name">图2. GC含量与覆盖深度（Depth）关联分析统计图（三代reads比对）</p>
                    </div>
            </div>
            <div class="tMainBody ">
                <h2 id="e5">5 组装完整性评估</h2>
                    <p class="paragraph">根据基因组与转录组数据及保守基因集的比对，综合评估组装基因组的完整性。</p>
                    <div class="tMainBody ">
                        <h3 id="e51">5.1.保守基因完整性评估</h3>
                        <p class="paragraph">BUSCO（Benchmarking Universal Single-Copy Orthologs）是利用直系同源数据库OrthoDB对基因组组装的完整性进行定量评估的软件。 BUSCO抽样了数百个基因组，从中选择单拷贝直系同源＞90%的基因构建了六种主要的系统进化分枝的基因集。</p>
                        <p class="paragraph">本次评估采用的BUSCO基因集为{14}。BUSCO评估结果见表6。</p>
                        <div class="studyUnit">
                        <p class="name">表6.BUSCO评估</p>
<table class="table left table-hover table-striped table-bordered" style="width:100%;">
                        <tr style="height: 5px;">
    <td style="text-align: center;">Type</td><td style="text-align: center;">Proteins</td><td style="text-align: center;">Percentage(%)</td></tr>
        <tr style="height: 5px;">{15}</table>
            <p class="paragraph">注：(1) Complete BUSCOs：能完整比对BUSCO的基因； </p>
            <p class="paragraph">  (2) Complete and single-copy BUSCOs：一个BUSCO 能完整比对上一个基因； </p>
            <p class="paragraph">  (3) Complete and duplicated BUSCOs：一个BUSCO 能完整比对上多个基因； </p>
            <p class="paragraph">  (4) Fragmented BUSCOs：只有部分序列能比对上BUSCO profile的基因； </p>
            <p class="paragraph">  (5) Missing BUSCOs：没有能比对上BUSCO profile的基因； </p>
            <p class="paragraph">  (6) Total BUSCO groups searched：BUSCO groups的基因总数。</p>
                    </div>
                </div>
                    <div class="tMainBody ">
                        <h3 id="e52">5.2.基因区完整性评估</h3>
                        <p class="paragraph">通过软件hisat2将转录组数据比对到基因组上，统计比对率，该值越高则表示基因组相对来说组装得越完整。</p>
                        <div class="studyUnit">
                        <p class="name">表7.转录组比对率统计</p>
                        <table class="table left table-hover table-striped table-bordered" style="width:100%;">
                        <tr style="height: 5px;">
                        <td style="text-align: center;">Sample ID</td><td style="text-align: center;">Mapping rate(%)</td><td style="text-align: center;">Paired mapping rate(%)</td></tr>{16}</table>
                </div>
            </div>
        </div>
            <div class="tMainBody ">
                <h2 id="f6">6 准确性及杂合度评估</h2>
                <p class="paragraph">将二代reads比对到基因组后，通过软件samtools、picard以及GATK识别突变。分别统计SNP、InDel的纯合率和杂合率，纯合率越低表明基因组准确性越高，杂合率越高说明基因组的杂合率越高。</p>
                <div class="studyUnit">
                    <p class="name">表8.纯合率及杂合率统计</p>
                    <table class="table left table-hover table-striped table-bordered" style="width:100%;">
                        <tr style="height: 5px;">
    <td style="text-align: center;">Sample ID</td><td style="text-align: center;">Homozygous SNP rate(%)</td><td style="text-align: center;">Homozygous InDel rate(%)</td><td style="text-align: center;">Heterozygous SNP rate(%)</td><td style="text-align: center;">Heterozygous InDel rate(%)</td></tr>{17}</table>
                </div>
            </div>
            <div class="tMainBody ">
                <h2 id="g7">7 基因组圈图</h2>
                    </table>
                </div>
            </div>
        </div>
    </div>    
</div>  
</body>
</html>
    '''.format(species_name,species_name,species_name,stat_info,base_info,info_nt,info_2_mapp,info_3,png_2_random,png_2_random,png_2,png_3_random,png_3_random,png_3,'busco_dataset',info_busco,info_mrna,info_2_hm)
        w.write(info)
        w.close()

    def main(self):
        args=self.getopt()
        config=Config()
        check_dir(args.outdir)
        self.report_result(args.input,args.outdir,args.ref,args.species)

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
            '--resume',help='please input resume file',dest='resume',type=str)
        parser.add_argument(
            '--infofile',help='the each step info file\nuse for report',dest='infofile',type=str)
        parser.add_argument(
            '--bamfile',help='input sort bam file',dest='bamfile',type=str)
        args = parser.parse_args()
        if not args.ref:
            print('assemble file must be given!!!')
            sys.exit(0)
        elif not args.bamfile and not args.input:
            print('bamfile or input fofn file must be given!!!')
            sys.exit()
        elif not args.outdir:
            print('outdir must be given!!!')
            sys.exit(0)
        return args

    def assemble3_stat(self,ref,fofnfile,blasr,samtools,r_gc_depth,total_gc,outdir,w,resum,nproc):
        outdir=outdir+'/assemble3'
        check_dir(outdir)
        count=0
        output_prefix=fofnfile.split('/')[-1].split('_')[0]
        # #################3 mapped###########
        cmd='%s %s %s --nproc %s --bam --out %s/%s.bam'%(blasr,fofnfile,ref,nproc,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        # #################3 sort#############
        cmd='%s sort -o %s/%s.sort.bam %s/%s.bam'%(samtools,outdir,output_prefix,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        # #################3 depth###########
        cmd='%s depth %s/%s.sort.bam>%s/%s.depth'%(samtools,outdir,output_prefix,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        mapped_rate=mapped_result('%s/%s.sort.bam'%(outdir,output_prefix),samtools)
        total_depth,coverage_rate,average_depth=depth_result('%s/%s.depth'%(outdir,output_prefix))
        png_3=gc_depth(total_depth,total_gc,r_gc_depth,output_prefix,outdir)
        if os.path.exists('%s/%s.bam'%(outdir,output_prefix)):os.remove('%s/%s.bam'%(outdir,output_prefix))
        count+=1
        return coverage_rate,mapped_rate,png_3,average_depth

    def bam_stat(self,bamfile,samtools,r_gc_depth,total_gc,outdir,w,resum):
        outdir=outdir+'/assemble3'
        check_dir(outdir)
        count=0
        output_prefix=target1.split('/')[-1].split('_')[0]
        cmd='%s depth %s>%s/%s.depth'%(samtools,bamfile,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        mapped_rate=mapped_result('%s/%s.sort.bam'%(outdir,output_prefix),samtools)
        total_depth,coverage_rate,average_depth=depth_result('%s/%s.depth'%(outdir,output_prefix))
        png_3=gc_depth(total_depth,total_gc,r_gc_depth,output_prefix,outdir)
        if os.path.exists('%s/%s.bam'%(outdir,output_prefix)):os.remove('%s/%s.bam'%(outdir,output_prefix))
        count+=1
        return coverage_rate,mapped_rate,png_3,average_depth

    def main(self):
        args=self.getopt()
        config=Config()
        blasr = check_software('blasr') if check_software('blasr') else config.blasr
        samtools = check_software('samtools') if check_software('samtools') else config.samtools
        r_gc_depth=check_software('gc-coverage.r') if check_software('gc-coverage.r') else config.r_gc_depth
        check_dir(args.outdir)
        total_gc= ref_deal(args.ref)
        resume=0 if not args.resume else resume_info(args.resume)
        w=open('%s/3_resume.txt'%args.outdir,'w')
        if args.bamfile:
            coverage_rate3,mapped_rate3,png_3,average_depth=self.bam_stat(args.bamfile,samtools,r_gc_depth,total_gc,args.outdir,w,resume)
        elif args.input:
            coverage_rate3,mapped_rate3,png_3,average_depth=self.assemble3_stat(args.ref,' '.join(args.input),blasr,samtools,r_gc_depth,total_gc,args.outdir,w,resume,args.thread)
        w.close()
        sample_name=args.input[0].split('/')[-1].split('_')[0] if len(args.input)>1 else ''
        info='<3_mapped\t%s\tmapped_rate:%s\taverage_depth:%.3f\tcoverage_rate:%s\tpng_3:%s\n'%(sample_name,':'.join(mapped_rate3),average_depth,':'.join(coverage_rate3),png_3)
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
        parser.add_argument("-v", dest="var_type", choices=["snp", "indel","sv", "cnv", "depth"],help='variant type [required].(2_mapped)',default='snp',type=str)
        parser.add_argument(
            '-n','--thread',help='set num of thread',type=int,dest='thread',default=10)
        parser.add_argument("-b", dest='block', default=50000, type=int,
            help='block length [default:50000].(2_mapped)')
        parser.add_argument(
            '--resume',help='please input resume file',dest='resume',type=str)
        parser.add_argument(
            '--infofile',help='the each step info file\nuse for report',dest='infofile',type=str)
        parser.add_argument(
            '--circos',help='chose if you want draw circos(only chromosomes)[y|n]',choices=['y','n'],default='n',type=str)
        args = parser.parse_args()
        if not args.ref:
            print('assemble file must be given!!!')
            sys.exit(0)
        elif not args.input:
            print('input file must be given!!!')
            sys.exit(0)
        elif not args.outdir:
            print('outdir must be given!!!')
            sys.exit(0)
        return args

    def get_yeast(self,ref_len, select_chrom, out_put):
        ''' refernece file: ref_id{Tab}ref_len '''
        # col_lst = ['106,61,154','128,177,211','139,58,98','177,89,40',
        #     '179,222,105','70,186,218','202,178,214','227,26,28',
        #     '251,128,114','145,51,154','86,69,46','179,222,105',
        #     '118,167,201','128,177,211','252,205,229','102,186,218',
        #     '139,58,98','251,128,114','253,180,98']
        for line in open(ref_len,'r'):
            col_lst.append('%d,%d,%d'%([random.randint(0,255),random.randint(0,255),random.randint(0,255)]))    
        chroms = pd.Series(select_chrom)
        intable = pd.read_table(ref_len, header=None)
        intable = intable.iloc[:,0:2]
        good_table = intable[intable[0].isin(chroms)]

        colu = pd.Series(col_lst[:len(good_table)], name="Col")
        good_table = pd.concat([good_table, colu], axis=1)
        info = lambda lst: pd.Series(["contig", "-", lst[0], lst[0], '0', lst[1], lst['Col']])
        result = good_table.apply(info, axis=1)
        result.to_csv(out_put, sep='\t', header=None, index=False)

    def variants_to_circosconf(self,chromf,outdir,output_prefix):
        w=open('%s/circos.conf'%outdir,'w')
        w.write('''<colors>
<<include colors.conf>>
<<include etc/brewer.conf>>
</colors>

<fonts>
<<include etc/fonts.conf>>
</fonts>


<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

chromosomes_reverse = /genome|mito/
chromosomes_units           = 100
chromosomes_display_default = yes
chromosomes_order  = {0}

### single genomes
# specify the karyotype file here - try other karyotypes in data/karyotype
karyotype = ./yeast.txt

<plots>

<plot>
type        = histogram
file        = {1}/snp.circos
r1          = 0.92r
r0          = 0.84r
min         = 0
#max         = 0.2
extend_bin  = no
fill_color  = dred
color       = dred
thickness   = 0
orientation = out
</plot>

<plot>
type        = histogram
file        = {2}/insert.circos
r1          = 0.76r
r0          = 0.68r
min         = 0
#max         = 0.2
extend_bin  = no
fill_color  = green
color       = green
thickness   = 0
orientation = out
</plot>

<plot>
type        = histogram
file        = {3}/delete.circos
r1          = 0.60r
r0          = 0.52r
min         = 0
#max         = 0.2
extend_bin  = no
fill_color  = blue
color       = blue
thickness   = 0
orientation = out
</plot>

<plot>
type        = histogram
file        = {4}/sv.circos
r1          = 0.44r
r0          = 0.36r
min         = 0
#max         = 0.2
extend_bin  = no
fill_color  = black
color       = black
thickness   = 0
orientation = out
</plot>

<plot>
type        = histogram
file        = {5}/cnv.circos
r1          = 0.28r
r0          = 0.20r
min         = 0
#max         = 0.2
extend_bin  = no
fill_color  = purple
color       = purple
thickness   = 0
orientation = out
</plot>


</plots>

<<include etc/housekeeping.conf>>
'''.format(chromf,outdir,outdir,outdir,outdir,outdir))

    def snp_sv_density(self,vcf, chroms):
        out_dic = {}
        with open(vcf) as inf:
            for line in inf:
                
                if line.startswith("#"): continue
                if "PASS" not in line: continue
                llst = line.split("\t")
                chrom, loc = llst[0], int(llst[1])
                if chrom not in chroms: continue
                out_dic.setdefault(chrom, []).append(loc)
        return out_dic

    def sv_density(self,vcf, chroms):
        out_dic = {}
        inf = gzip.open(vcf) if vcf.endswith("gz") else open(vcf)

        for line in inf:
            if line.startswith("#"): continue
            if "SVTYPE=BND" in line: continue
            llst = line.split("\t")
            chrom, loc = llst[0], int(llst[1])
            if chrom not in chroms: continue
            out_dic.setdefault(chrom, []).append(loc)
        inf.close()
        return out_dic

    def indel_density(self,vcf, chroms):
        ins_dic, del_dic = {}, {}
        with open(vcf) as inf:
            for line in inf:
                if line.startswith("#") or "PASS" not in line: continue     
                llst = line.split("\t")
                chrom, loc, ref, alt = llst[0], int(llst[1]), llst[2], llst[3]
                if chrom not in chroms: continue
                if len(ref) < len(alt):
                    ins_dic.setdefault(chrom, []).append(loc)
                else:
                    del_dic.setdefault(chrom, []).append(loc)
        return ins_dic, del_dic

    def get_depth_dic(self,depth_f, chroms):
        outdic = defaultdict(dict)
        with open(depth_f) as depf:
            for line in depf:
                chrom, loc, dep = line.strip().split()
                if chrom not in chroms: continue
                outdic[chrom].update({int(loc):int(dep)})
        return outdic

    def output(self,var_dic, block, var, color, fa_len,outdir):
        outname = outdir+'/'+var + '.circos'
        output = open(outname, 'w')
        result = {}
        lst = []
        if var == 'depth':
            for k, v in var_dic.items():
                block_num = int(fa_len[k]/block)
                lst = [0.01] * (block_num + 1)
                for loc, dep in v.items():
                    pos = int(loc /block)
                    lst[pos] += dep
                result[k] = lst
        else:
            for k, v in var_dic.items():
                block_num = int(fa_len[k]/block)
                lst = [0.01] * (block_num + 1)
                for i in v:
                    pos = int(i / block)
                    lst[pos] += 1
                result[k] = lst
            #print result

        #print result
        for k, v in result.items():
            count = 0
            for i in v:
                output.write("{0}\t{1}\t{2}\t{3}\tfill_color={4}\n".format(k, count*block+1, (count+1)*block, i, color))
                count += 1
                
    def main_result(self,ref, infile, chromf, var_type, block ,outdir):
        chr = [chromf]    
        fa = pysam.FastaFile(ref)
        fa_len = dict(zip(fa.references, fa.lengths))

        color = {"snp":"red", "insert":"green", "delete":"blue", 
            "sv":"black", 'cnv':"purple", "depth":"gray"}
        if var_type == "snp":
            snp_dic = self.snp_sv_density(infile, chr)
            self.output(snp_dic, block, "snp", color["snp"], fa_len,outdir)
        elif var_type == "indel":
            ins_dic, del_dic = self.indel_density(infile, chr)
            self.output(ins_dic, block, "insert", color["insert"], fa_len,outdir)
            self.output(del_dic, block, 'delete', color["delete"], fa_len,outdir)
        elif var_type == 'sv':
            sv_dic = self.sv_density(infile, chr) # SV density is same as SNP density
            self.output(sv_dic, block, "sv", color["sv"], fa_len,outdir)
        elif var_type == "cnv":
            cnv_dic = self.snp_sv_density(infile, chr) # CNV density issame as SNP density as well
            self.output(cnv_dic, block, "cnv", color["cnv"], fa_len,outdir)
        elif var_type == 'depth':
            dep_dic = self.get_depth_dic(infile, chr)
            self.output(dep_dic, block, "depth", color["depth"], fa_len,outdir)

    def assemble2_stat(self,ref,tartget1,tartget2,samtools,bwa,total_gc,r_gc_depth,thread_num,outdir,w,resum,circos,var_type,block):
        outdir=outdir+'/assemble2'
        check_dir(outdir)
        count=0
        output_prefix=tartget1.split('/')[-1].split('_')[0]
        config=Config()
        #####################software###############
        gatk=check_software('GenomeAnalysisTK.jar') if check_software('GenomeAnalysisTK.jar') else config.gatk
        picard=check_software('picard.jar') if check_software('picard.jar') else config.picard
        #bgzip=check_software('bgzip') if check_software('bgzip') else config.bgzip
        #tabix=check_software('tabix') if check_software('tabix') else config.tabix
        ##################bwa index###############
        cmd='%s index %s'%(bwa,ref)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        ###################samtools faidx##########
        cmd='%s faidx %s'%(samtools,ref)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        ##################picard index############
        cmd='java -jar  %s CreateSequenceDictionary R=%s O=%s.dict'%(picard,ref,'.'.join(ref.split('.')[:-1]))
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        ##################bwa mem##################
        ###########way 1######################
        # cmd = "%s mem -t %s %s %s %s  | %s view -@ %s -bS -t %s.fai - |%s sort -@ %s -m 2000000000 - -o %s/%s.sorted.bam"%(bwa,thread_num,ref,tartget1,tartget2,samtools,thread_num,ref,samtools,thread_num,outdir,output_prefix)
        # if count>resum:
        #     if run_cmd(cmd):sys.exit(1)
        # all_bam=[]
        # for i in glob.glob('%s/%s.sorted.bam*'%(outdir,output_prefix)):
        #     all_bam.append(i.strip())
        # bams=' '.join(all_bam)
        # cmd='%s merge %s/%s.sorted.bam %s'%(samtools,outdir,output_prefix,bams)
        # if count>resum:
        #     if run_cmd(cmd):sys.exit(1)

        ##########way 2#######################
        cmd = "%s mem -t %s -o %s/%s.sam %s %s %s"%(bwa,thread_num,outdir,output_prefix,ref,tartget1,tartget2)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        cmd="%s view -@ %s -bSh %s/%s.sam -o %s/%s.bam"%(samtools,thread_num,outdir,output_prefix,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        cmd="%s sort -m 2G -@ %s -o %s/%s.sorted.bam %s/%s.bam"%(samtools,thread_num,outdir,output_prefix,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)

        ##################depth######################
        cmd='%s depth %s/%s.sorted.bam>%s/%s.depth'%(samtools,outdir,output_prefix,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)

        #######################picard sort##############
        cmd="java -Djava.io.tmpdir=%s/tmp -jar %s SortSam VALIDATION_STRINGENCY=LENIENT INPUT=%s/%s.sorted.bam OUTPUT=%s/%s.sort.bam SORT_ORDER=coordinate"%(outdir,picard,outdir,output_prefix,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)

        ######################add header################
        cmd="java -jar %s AddOrReplaceReadGroups I=%s/%s.sort.bam O=%s/%s.bam RGID=4 RGLB=%s RGPL=illumina RGPU=%s RGSM=20"%(picard,outdir,output_prefix,outdir,output_prefix,output_prefix,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)

        #######################samtools index###########
        cmd="%s index %s/%s.bam"%(samtools,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)

        ######################calling##################
        cmd="java -Xmx2G -jar %s -T HaplotypeCaller -R %s -I %s/%s.bam --genotyping_mode DISCOVERY -o %s/%s.gvcf -nct %s"%(gatk,ref,outdir,output_prefix,outdir,output_prefix,thread_num)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            print('good')
            print(cmd)
            if run_cmd(cmd):sys.exit(1)

        ######################select snp && indel and stat#######
        cmd='java -Xmx2G -jar %s -T SelectVariants -o %s/%s.snp.vcf --variant %s/%s.gvcf -selectType SNP -R %s'%(gatk,outdir,output_prefix,outdir,output_prefix,ref)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        homog,mix,chroms,total_site,rate_total=[0,0],[0,0],[],[0,0],[0,0,0,0]
        for line in open('%s/%s.snp.vcf'%(outdir,output_prefix),'r'):
            if line.startswith('#'):continue
            chrs=line.split('\t')[0]
            if chrs not in chroms:chroms.append(chrs)
            types=line.strip().split('\t')[-1].split(':')[0]
            if types=='1/1':
                homog[0]+=1
            elif types=='0/1':
                mix[0]+=1
            total_site[0]+=1
        rate_total[1]=float(homog[0])/float(total_site[0])
        rate_total[0]=float(mix[0])/float(total_site[0])
        # homog_rate_snp=float(homog[0])/float(total_site[0])
        # mix_rate_snp=float(mix[0])/float(total_site[0])

        cmd='java -Xmx2G -jar %s -T SelectVariants -o %s/%s.indel.vcf --variant %s/%s.gvcf -selectType INDEL -R %s'%(gatk,outdir,output_prefix,outdir,output_prefix,ref)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        for line in open('%s/%s.indel.vcf'%(outdir,output_prefix),'r'):
            types=line.strip().split('\t')[-1].split(':')[0]
            if types=='1/1':
                homog[1]+=1
            elif types=='0/1':
                mix[1]+=1
            total_site[1]+=1
        rate_total[3]=float(homog[1])/float(total_site[1])
        rate_total[2]=float(mix[1])/float(total_site[1])
        # homog_rate_indel=float(homog[1])/float(total_site[1])
        # mix_rate_indel=float(mix[1])/float(total_site[1])

        ###################filter#######################
        cmd='java -Xmx4g -jar %s -R %s -T VariantFiltration --variant %s/%s.snp.vcf --clusterSize 4 --clusterWindowSize 10 --maskName aroundIndel --mask %s/%s.indel.vcf -maskExtend 3 --filterName \"lowMQRankSum\" --filterExpression \"MQRankSum < -12.5\" --filterName \"highFS\" --filterExpression \"FS > 60.0\" --filterName \"lowReadPosRankSum\" --filterExpression \"ReadPosRankSum < -8.0\" --filterName \"lowMQ\" --filterExpression \"MQ < 40.0\" --filterName \"lowQD\" --filterExpression \"QD < 2.0\" --out %s/%s.snp.filt.vcf --genotypeFilterName \"lowDP\" --genotypeFilterExpression \"DP < 8.0\"'%(gatk,ref,outdir,output_prefix,outdir,output_prefix,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        ###############circos if contig don't neet to do this#########################
        if circos=='y':
            output_circos='%s/%s'%(outdir,output_prefix)
            check_dir(output_circos)
            self.get_yeast(ref+'.fai',chroms,'%s/yeast.txt'%output_circos)
            self.main_result(ref, '%s/%s.snp.filt.vcf'%(outdir,output_prefix), chroms, var_type,block,output_circos)
            self.variants_to_circosconf(chroms,output_circos,output_prefix)
            os.system('cp -r %s %s'%(config.circos_conf,output_circos))
            os.system('cd %s && %s && cd -'%(output_circos,config.circos_soft))

        mapped_rate=mapped_result('%s/%s.bam'%(outdir,output_prefix),samtools)
        total_depth,coverage_rate,average_depth=depth_result('%s/%s.depth'%(outdir,output_prefix))
        png_2=gc_depth(total_depth,total_gc,r_gc_depth,output_prefix,outdir)
        if os.path.exists('%s/%s.sam'%(outdir,output_prefix)):os.remove('%s/%s.sam'%(outdir,output_prefix))
        if os.path.exists('%s/%s.sorted.bam'%(outdir,output_prefix)):os.remove('%s/%s.sorted.bam'%(outdir,output_prefix))
        if os.path.exists('%s/%s.sort.bam'%(outdir,output_prefix)):os.remove('%s/%s.sort.bam'%(outdir,output_prefix))
        count+=1
        return coverage_rate,mapped_rate,rate_total,png_2,average_depth
    
    def main(self):
        args=self.getopt()
        config=Config()
        check_dir(args.outdir)
        samtools = check_software('samtools') if check_software('samtools') else config.samtools
        r_gc_depth=check_software('gc-coverage.r') if check_software('gc-coverage.r') else config.r_gc_depth
        bwa=check_software('bwa') if check_software('bwa') else config.bwa
        total_gc= ref_deal(args.ref)
        resume=0 if not args.resume else resume_info(args.resume)
        w=open('%s/2_resume.txt'%args.outdir,'w')
        coverage_rate2,mapped_rate2,total_rate,png_2,average_depth=self.assemble2_stat(args.ref,args.input[0],args.input[1],samtools,bwa,total_gc,r_gc_depth,args.thread,args.outdir,w,resume,args.circos,args.var_type,args.block)
        w.close()

        sample_name=args.input[0].split('/')[-1].split('_')[0]
        info='<2_mapped\t%s\tmapped_rate:%s\taverage_depth:%.3f\tcoverage_rate:%s\tsnp_mix:%.3f\tsnp_homog:%.3f\tindel_mix:%.3f\tindel_homog:%.3f\tpng_2:%s\n'%(sample_name,':'.join(mapped_rate2),average_depth,':'.join(coverage_rate2),total_rate[0],total_rate[1],total_rate[2],total_rate[3],png_2)
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
            '--resume',help='please input resume file',dest='resume',type=str)
        parser.add_argument(
            '--infofile',help='the each step info file\nuse for report',dest='infofile',type=str)
        args=parser.parse_args()
        return args

    def assemble_mrna(self,ref,hisat2_build,hisat2,mrna_1,mrna_2,thr_num,samtools,outdir,w,resum):
        outdir=outdir+'/assemble_mrna'
        check_dir(outdir)
        count=0
        output_prefix=mrna_1.split('/')[-1].split('.')[0].split('_')[0]
        #################hisat2 index###################
        cmd='%s -p %s %s %s'%(hisat2_build,thr_num,ref,ref.split('.')[0])
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        ##from here
        cmd='%s -p %s -x %s -1 %s -2 %s -S %s/%s.sam'%(hisat2,thr_num,ref.split('.')[0],mrna_1,mrna_2,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        cmd='%s view -bS %s/%s.sam>%s/%s.bam'%(samtools,outdir,output_prefix,outdir,output_prefix)
        count+=1
        w.write('%d\n'%count)
        if count>resum:
            if run_cmd(cmd):sys.exit(1)
        mapped_rate=mapped_result('%s/%s.bam'%(outdir,output_prefix),samtools)
        if os.path.exists('%s/%s.sam'%(outdir,output_prefix)):os.remove('%s/%s.sam'%(outdir,output_prefix))
        count+=1
        w.write('%d\n'%count)
        return mapped_rate

    def main(self):
        args=self.getopt()
        config=Config()
        check_dir(args.outdir)
        samtools = check_software('samtools') if check_software('samtools') else config.samtools
        hisat2=check_software('hisat2') if check_software('hisat2') else config.hisat2
        hisat2_build=check_software('hisat2-build') if check_software('hisat2-build') else config.hisat2_build
        resume=0 if not args.resume else resume_info(args.resume)
        w=open('%s/mrna_resume.txt'%args.outdir,'w')
        mapped_ratem=self.assemble_mrna(args.ref,hisat2_build,hisat2,args.input[0],args.input[1],args.thread,samtools,args.outdir,w,resume)
        w.close()
        sample_name=args.input[0].split('/')[-1].split('_')[0]
        info='<mrna_mapped\t%s\tmapped_rate:%s\n'%(sample_name,':'.join(mapped_ratem))
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
            '-i','--input',help='input file nt_mapped:[nt_database]',dest='input',type=str,default='/thinker/storage/data/database/01.nt_nr/FASTA/nt')
        parser.add_argument(
            '-o','--outdir',help='output dir',type=str,dest='outdir',required=True)
        parser.add_argument(
            '-n','--thread',help='set num of thread',type=int,dest='thread',default=10)
        # parser.add_argument(
        #     '--resume',help='please input resume file',dest='resume',type=str)
        parser.add_argument(
            '--infofile',help='the each step info file\nuse for report',dest='infofile',type=str)
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
                        if len(info)-i>200:
                            read_200[key]=info[i:]
                    key=re.split(r'\s+',line.strip())[0]
                    info=''
                    continue
                info+=line.strip()               
            count=0
            if total_read<20000:
                for key,value in read_200.items():
                    count+=1
                    w.write(key+'\n'+value+'\n')
                    if count+total_read>=20000:
                        break
        w.close()
        output_prefix=assemble_file.split('/')[-1].split('.')[0]
        cmd='%s -db %s -query %s/new_%s -out %s/%s.out -evalue 0.00001 -outfmt \'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle\'  -num_threads %d'%(blastn,nt_database,outdir,assemble_file.split('/')[-1],outdir,output_prefix,thr_num)
        run_cmd(cmd)
        result_info,total_mapp_read={},{}
        for line in open('%s/%s.out'%(outdir,output_prefix),'r'):
            info=line.strip().split('\t')
            if info[0] not in total_mapp_read:total_mapp_read[info[0]]=0
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
        config=Config()
        check_dir(args.outdir)
        #here to start
        blastn=check_software('blastn') if check_software('blastn') else config.blastn
        check_dir(args.outdir)
        dict_nt,total_mapp_num,total_read_num=self.blast_nt(blastn,args.input,args.ref,args.thread,args.outdir)
        info='<blast_nt\ttotal_read_num:%s\ttotal_mapp_num:%s'%(total_read_num,total_mapp_num)
        for key,value in dict_nt.items():
            info+='\t%s:%d,%f'%(key,value[1],value[0])
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
        parser.add_argument('--species_all', help='get the all species name[y|n]', dest='species_all',default='n',type=str)
        parser.add_argument('--infofile',help='the each step info file use for report',dest='infofile',type=str)
        args = parser.parse_args()
        return args

    def run_gene_mode(self,genome,output,busco,cpu,august_dataset,busco_dataset,env):
        out_prefix=genome.split('/')[-1].split('.')[0]
        cmd='gzip -dc %s >%s && '%(genome,os.path.basename(genome).rstrip('.gz'))
        if '.gz' in genome:run_cmd(cmd)
        check_dir(output)
        outdir=output+'/1.busco_geno'
        check_dir(outdir)
        genome=os.path.abspath(genome)
        cmd='source %s && cd %s && python %s --in %s --cpu %s --out %s --force --mode genome --lineage_path %s %s --blast_single_core && cd -'%(env,outdir,busco,genome,cpu,out_prefix,busco_dataset,august_dataset)
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

    def main(self):
        busco_dict={'arthropoda': '/thinker/storage/data/database/08.busco/arthropoda_odb9', 'rhizobiales': '/thinker/storage/data/database/08.busco/rhizobiales_odb9', 'microsporidia': '/thinker/storage/data/database/08.busco/microsporidia_odb9', 'gammaproteobacteria': '/thinker/storage/data/database/08.busco/gammaproteobacteria_odb9', 'deltaepsilonsub': '/thinker/storage/data/database/08.busco/deltaepsilonsub_odb9', 'bacteroidetes': '/thinker/storage/data/database/08.busco/bacteroidetes_odb9', 'embryophyta': '/thinker/storage/data/database/08.busco/embryophyta_odb9', 'tenericutes': '/thinker/storage/data/database/08.busco/tenericutes_odb9', 'enterobacteriales': '/thinker/storage/data/database/08.busco/enterobacteriales_odb9', 'metazoa': '/thinker/storage/data/database/08.busco/metazoa_odb9', 'mammalia': '/thinker/storage/data/database/08.busco/mammalia_odb9', 'insecta': '/thinker/storage/data/database/08.busco/insecta_odb9', 'sordariomyceta': '/thinker/storage/data/database/08.busco/sordariomyceta_odb9', 'bacteria': '/thinker/storage/data/database/08.busco/bacteria_odb9', 'eukaryota': '/thinker/storage/data/database/08.busco/eukaryota_odb9', 'protists': '/thinker/storage/data/database/08.busco/protists_ensembl', 'nematoda': '/thinker/storage/data/database/08.busco/nematoda_odb9', 'ascomycota': '/thinker/storage/data/database/08.busco/ascomycota_odb9', 'actinobacteria': '/thinker/storage/data/database/08.busco/actinobacteria_odb9', 'dikarya': '/thinker/storage/data/database/08.busco/dikarya_odb9', 'aves': '/thinker/storage/data/database/08.busco/aves_odb9', 'saccharomyceta': '/thinker/storage/data/database/08.busco/saccharomyceta_odb9', 'spirochaetes': '/thinker/storage/data/database/08.busco/spirochaetes_odb9', 'betaproteobacteria': '/thinker/storage/data/database/08.busco/betaproteobacteria_odb9', 'actinopterygii': '/thinker/storage/data/database/08.busco/actinopterygii_odb9', 'alveolata_stramenophiles': '/thinker/storage/data/database/08.busco/alveolata_stramenophiles_ensembl', 'hymenoptera': '/thinker/storage/data/database/08.busco/hymenoptera_odb9', 'saccharomycetales': '/thinker/storage/data/database/08.busco/saccharomycetales_odb9', 'bacillales': '/thinker/storage/data/database/08.busco/bacillales_odb9', 'proteobacteria': '/thinker/storage/data/database/08.busco/proteobacteria_odb9', 'tetrapoda': '/thinker/storage/data/database/08.busco/tetrapoda_odb9', 'lactobacillales': '/thinker/storage/data/database/08.busco/lactobacillales_odb9', 'cyanobacteria': '/thinker/storage/data/database/08.busco/cyanobacteria_odb9', 'euarchontoglires': '/thinker/storage/data/database/08.busco/euarchontoglires_odb9', 'firmicutes': '/thinker/storage/data/database/08.busco/firmicutes_odb9', 'fungi': '/thinker/storage/data/database/08.busco/fungi_odb9', 'vertebrata': '/thinker/storage/data/database/08.busco/vertebrata_odb9', 'basidiomycota': '/thinker/storage/data/database/08.busco/basidiomycota_odb9', 'endopterygota': '/thinker/storage/data/database/08.busco/endopterygota_odb9', 'pezizomycotina': '/thinker/storage/data/database/08.busco/pezizomycotina_odb9', 'clostridia': '/thinker/storage/data/database/08.busco/clostridia_odb9', 'laurasiatheria': '/thinker/storage/data/database/08.busco/laurasiatheria_odb9', 'diptera': '/thinker/storage/data/database/08.busco/diptera_odb9', 'eurotiomycetes': '/thinker/storage/data/database/08.busco/eurotiomycetes_odb9'}
        augustus_dataset={'leishmania tarentolae': 'leishmania_tarentolae', 'bombus terrestris': 'bombus_terrestris2', 'neurospora': 'neurospora', 'amphimedon': 'amphimedon', 'verticillium albo-atrum': 'verticillium_albo_atrum1', 'conidiobolus coronatus': 'Conidiobolus_coronatus', 'lethenteron camtschaticum': 'japaneselamprey', 'apis dorsata': 'adorsata', 'burkholderia pseudomallei': 'b_pseudomallei', 'verticillium longisporum': 'verticillium_longisporum1', 'rice': 'rice', 'histoplasma capsulatum': 'histoplasma_capsulatum', 'eremothecium gossypii': 'eremothecium_gossypii', 'cryptococcus': 'cryptococcus', 'caenorhabditis': 'caenorhabditis', 'saccharomyces cerevisiae rm11-1a': 'saccharomyces_cerevisiae_rm11-1a_1', 'fusarium graminearum': 'fusarium_graminearum', 'neurospora crassa': 'neurospora_crassa', 'coccidioides immitis': 'coccidioides_immitis', 'xiphophorus maculatus': 'Xipophorus_maculatus', 'coprinus': 'coprinus', 'apis mellifera': 'honeybee1', 'chlamys': 'chlamy2011', 'coprinus cinereus': 'coprinus_cinereus', 'trichinella': 'trichinella', 'rhizopus oryzae': 'rhizopus_oryzae', 'cryptococcus neoformans var. neoformans jec21': 'cryptococcus_neoformans_neoformans_JEC21', 'heliconius melpomene': 'heliconius_melpomene1', 'nicotiana attenuata': 'coyote_tobacco', 'laccaria bicolor': 'laccaria_bicolor', 'aspergillus nidulans': 'aspergillus_nidulans', 'fusarium': 'fusarium', 'chlorella': 'chlorella', 'sulfolobus solfataricus': 'sulfolobus_solfataricus', 'pea aphid': 'pea_aphid', 'petromyzon marinus': 'sealamprey', 'volvox': 'volvox', 'saccharomyces cerevisiae s288c': 'saccharomyces_cerevisiae_S288C', 'triticum aestivum': 'wheat', 'botrytis cinerea': 'botrytis_cinerea', 'human': 'human', 'debaryomyces hansenii': 'debaryomyces_hansenii', 'nasonia': 'nasonia', 'cryptococcus neoformans var. neoformans b-3501a': 'cryptococcus_neoformans_neoformans_B', 'pichia stipitis': 'pichia_stipitis', 'bombus impatiens': 'bombus_impatiens1', 'elephant shark': 'elephant_shark', 'caenorhabditis elegans': 'c_elegans_trsk', 'chicken': 'chicken', 'trichomonas': 'Trichomonas', 'galdieria': 'galdieria', 'chaetomium globosum': 'chaetomium_globosum', 'aspergillus oryzae': 'aspergillus_oryzae', 'saccharomyces': 'saccharomyces', 'parasteatoda': 'parasteatoda', 'maize': 'maize', 'aedes': 'aedes', 'staphylococcus aureus': 's_aureus', 'danio rerio': 'zebrafish', 'sillago': 'Sillago', 'aspergillus fumigatus': 'aspergillus_fumigatus', 'lentinula edodes': 'Lentinula_edodes_ref', 'chlamydomonas': 'chlamydomonas', 'ancylostoma ceylanicum': 'ancylostoma_ceylanicum', 'drosophila melanogaster': 'fly', 'lodderomyces elongisporus': 'lodderomyces_elongisporus', 'rhodnius': 'rhodnius', 'tetrahymena': 'tetrahymena', 'ustilago maydis': 'ustilago_maydis', 'spirodela polyrhiza': 'Spirodela_polyrhiza_new', 'toxoplasma': 'toxoplasma', 'glycine max': 'Gmax_Wm82', 'pneumocystis': 'pneumocystis', 'streptococcus pneumoniae': 's_pneumoniae', 'yarrowia lipolytica': 'yarrowia_lipolytica', 'glycine max x glycine soja': 'Gsoja', 'aspergillus terreus': 'aspergillus_terreus', 'schizosaccharomyces pombe': 'schizosaccharomyces_pombe', 'camponotus floridanus': 'camponotus_floridanus', 'ustilago': 'ustilago', 'cacao': 'cacao', 'phanerochaete chrysosporium': 'phanerochaete_chrysosporium', 'pelteobagrus fulvidraco': 'Pelteobagrus_fulvidraco', 'schistosoma': 'schistosoma2', 'encephalitozoon cuniculi gb-m1': 'encephalitozoon_cuniculi_GB', 'histoplasma': 'histoplasma', 'candida albicans': 'candida_albicans', 'pomacea canaliculata': 'Pomacea_canaliculata', 'cryptococcus gattii vgi': 'cryptococcus_neoformans_gattii', 'bacteria': 'template_prokaryotic', 'brugia': 'brugia', 'candida tropicalis': 'candida_tropicalis', 'arabidopsis': 'arabidopsis', 'escherichia coli k-12': 'E_coli_K12', 'candida guilliermondii': 'candida_guilliermondii', 'plasmodium falciparum': 'pfalciparum', 'thermoanaerobacter tengcongensis': 'thermoanaerobacter_tengcongensis', 'culex': 'culex', 'magnaporthe grisea': 'magnaporthe_grisea', 'tomato': 'tomato', 'tribolium': 'tribolium2012', 'kluyveromyces lactis': 'kluyveromyces_lactis', 'zea mays': 'maize5'}
        config=Config()
        args=self.getopt()
        august_dataset=''
        if args.species_all=='y':
            for line in open(config.species,'r'):
                if line.startswith('#') or line.startswith(' ') or line.startswith('\n'):
                    continue
                info=re.split(r'\s+',line.strip())
                print('\t'.join(info[:2]))
            sys.exit(0)
        if not args.output:
            print('config must be given!!!!!![-o]')
            sys.exit(1)
        elif not args.genome:
            print("genome file must be given!!!!![-g]")
            sys.exit(1)
        elif not args.species and not args.busco_dataset:
            print('species[--species] or (augustus_dataset[--augustus_dataset] and busco_dataset[--busco_dataset]) must be given!!!')
            sys.exit(1)
        if args.busco_dataset:
            if args.busco_dataset not in busco_dict:
                print('%s is wrong busco_dataset,you can chose the ture busco_dataset!!!\n')
                for key in busco_dict.keys():
                    print(key)
                sys.exit(1)
            if args.augustus_dataset not in augustus_dataset and args.augustus_dataset:
                chose_set=input('%s is wrong augustus_dataset,you can chose not to use augustus_dataset[y|n]?')
                if chose_set=='n':
                    for key in augustus_dataset.keys():
                        print(key)
                    sys.exit(1)
            august_dataset='' if not args.augustus_dataset else '--species %s --long'%augustus_dataset[args.augustus_dataset.lower()]
            busco_dataset=busco_dict[args.busco_dataset.lower()]
            info_busco=self.run_gene_mode(args.genome,args.output,config.busco,args.cpu,august_dataset,busco_dataset,config.env)
        else:        
            ################check if species is right###########
            dict_name,dict_node,check,august_check={},{},0,0
            match_name=' '.join(args.species)
            check_species=0
            if '_' not in match_name:
                pass
            else:
                for line_all in open(config.species,'r'):
                    if line_all.startswith('#') or line_all.startswith(' ') or line_all.startswith('\n'):
                        continue
                    info_all=re.split(r'\s+',line_all.strip())
                    if match_name==info_all[0]:
                        match_name=re.sub(r'_',' ',' '.join(info_all[-2:]).rstrip('<').rstrip(' ')) if len(info_all)>2 else re.sub(r'_',' ',info_all[0].rstrip('/'))

            for line in open(config.name_dmp,'r'):
                info=line.strip().split('|')
                if info[0].strip() not in dict_name:dict_name[info[0].strip()]=[]
                dict_name[info[0].strip()].append([info[1].strip(),info[3].strip()])
                if match_name.lower() == info[1].strip().lower():
                    key_id=info[0].strip()
                    check_species=1
            if not check_species:
                print('%s is wrong Latin name,please input right Latin name or you can use parameter\"--species_all y\" to get right input species'%match_name)
                sys.exit(1)
            ################chose dataset#############       
            if match_name.lower() in augustus_dataset:
                august_dataset=augustus_dataset[match_name.lower()]
                august_check=1
            for line in open(config.node_dmp,'r'):
                info=re.sub(r'\s+','',line.strip()).split('|')
                if info[0] not in dict_node:dict_node[info[0]]=[info[1],info[2]]
                if info[0]==key_id:
                    for i in dict_name[key_id]:
                        if i[1]=='scientific name' and i[0].lower() in busco_dict:
                            busco_dataset=busco_dict[i[0].lower()]
                            check=1
                            break
                if check:
                   break
            while 1:
                if key_id=='1' or check==1:
                    break
                else:
                    for i in dict_name[key_id]:
                        if not august_check:
                            if i[1]=='scientific name' and (i[0].lower() in augustus_dataset):
                                august_check=1
                                august_dataset=augustus_dataset[i[0].lower()]

                        if i[1]=='scientific name' and i[0].lower() in busco_dict:
                            busco_dataset=busco_dict[i[0].lower()]
                            check=1
                            break
                    key_id=dict_node[key_id][0]
            if check:
                august_dataset='' if not august_dataset else '--species %s --long'%augustus_dataset[august_dataset]
                info_busco=self.run_gene_mode(args.genome,args.output,config.busco,args.cpu,august_dataset,busco_dataset,config.env)
            else:
                print('%s is wrong Latin name,please input right Latin name or you can use parameter\"--species_all y\" to get right input species'%match_name)
                sys.exit(1)
        file_info_get='%s/result.txt'%args.output if not args.infofile else args.infofile
        result_file(info_busco,file_info_get)

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

