#!/usr/bin/env python
# -*- coding: utf-8 -*-

#####Import Module#####
import logging
import sys
import os
import math
import time
import argparse
import glob

class Config():
    def __init__(self):
        #############busco  config############
        self.env='/thinker/storage/udata/chenqi/assemble_evlu/busco_test/bin/env'
        self.busco='/thinker/storage/software/Genome/busco/scripts/run_BUSCO.py'
        self.species='/thinker/storage/udata/chenqi/assemble_evlu/busco_test/all.species'
        self.name_dmp='/thinker/storage/data/database/01.nt_nr/taxonomy/taxdump/names.dmp'
        self.node_dmp='/thinker/storage/data/database/01.nt_nr/taxonomy/taxdump/nodes.dmp'

        self.samtools='/thinker/storage/bioinfo/02Software/samtools/run/1.3.1/samtools'
        self.blasr='/thinker/storage/software/common/pacbio/smrtlink/install/smrtlink-release_5.1.0.26412/bundles/smrttools/install/smrttools-release_5.1.0.26366/smrtcmds/bin/blasr'
        self.bwa='/thinker/storage/bioinfo/02Software/bwa/run/bwa'
        self.r_gc_depth='/data/tmp_chenqi/gene_assemble/gc-coverage.r*'
        self.hisat2='/thinker/storage/software/Transcriptom/hisat2-2.1.0/hisat2'
        self.hisat2_build='/thinker/storage/software/Transcriptom/hisat2-2.1.0/hisat2-build'
        self.N50_stat='/data/tmp_chenqi/gene_assemble/stat_N50.pl*'
        self.blastn='/thinker/storage/bioinfo/02Software/blast/run/blast-2.4.0/blastn'
        self.gatk='/thinker/storage/bioinfo/02Software/GATK/run/3.8.0/GenomeAnalysisTK.jar'
        self.picard='/thinker/storage/bioinfo/02Software/picard/run/2.2.4/picard.jar'

def fmt_time(spend_time):
    spend_time = int(spend_time)
    day = 24 * 60 * 60
    hour = 60 * 60
    min = 60
    if spend_time < 60:
        return "%ds" % math.ceil(spend_time)
    elif spend_time > day:
        days = divmod(spend_time, day)
        return "%dd%s" % (int(days[0]), fmt_time(days[1]))
    elif spend_time > hour:
        hours = divmod(spend_time, hour)
        return '%dh%s' % (int(hours[0]), fmt_time(hours[1]))
    else:
        mins = divmod(spend_time, min)
        return "%dm%ds" % (int(mins[0]), math.ceil(mins[1]))

def check_software(software_path):
    if os.path.exists(software_path):
        logging.debug("Choose software:" + software_path + "!")
    else:
        output = os.popen('which ' + software_path)
        software_temp = output.read().strip()
        if os.path.exists(software_temp):
            software_path = software_temp
            logging.debug("Choose software:" + software_path + "!")
        else:
            logging.error("Can't locate the " + software_path + "!")
            return 0
    return software_path

def run_cmd(cmd):
    #may be here a problem but it is ok ,youjust need to change it 
    logging.info(cmd)
    flag = os.system(cmd)
    if flag == 0 or flag == 256:
        return 0
    return 1

def check_dir(dirs):
    if not os.path.exists(dirs):os.mkdir(dirs)
