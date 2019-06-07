#!/usr/bin/env python
import re
import os
import natsort
import collections
from ruamel import yaml  # conda install -c conda-forge -n rnaseq ruamel.yaml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', dest='fas', help='path to fastq files', metavar='fq', type=str,  nargs='+')
parser.add_argument('-o', dest='out', default='.')
args = parser.parse_args()

if not os.path.exists('config.yml'):
    f = open('config.yml', 'w')

config = collections.OrderedDict()

config['transcripts_fa'] = '/Users/mdascenzo/workspace/data/Homo_sapiens.GRCh38.rel83.cdna.all.fa'
config['tx2gene_fp'] = '/Users/mdascenzo/workspace/dev/workflows/src/rnaseq/R/tx2gene.EnsDb.Hsapiens.v86.csv'
config['out'] = '/Users/mdascenzo/workspace/dev/rnaseq-sim/analysis'

config['samples'] = collections.OrderedDict()

# currently assumes paired end reads
it = iter(natsort.natsorted(arg.fas))
for x in it:
    paired_reads = (x, next(it))
    sample_id = re.sub('_1$', '', os.path.splitext(os.path.splitext(os.path.basename(paired_reads[0]))[0])[0])

    config['samples'][sample_id] = collections.OrderedDict()
    config['samples'][sample_id]['read1'] = paired_reads[0]
    config['samples'][sample_id]['read2'] = paired_reads[1]

config_file = yaml.dump(config, Dumper=yaml.RoundTripDumper, default_flow_style=False).replace('!!omap', '').replace('- ', '')

fout = open(os.path.join(args.out, 'config.yml'), 'w')
fout.write(config_file)
