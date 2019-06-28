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
parser.add_argument('-s', dest='sim', action='store_true')
parser.add_argument('-n', dest='configfile_name', default='config.yml')

args = parser.parse_args()

# todo: verify input or exit

config = collections.OrderedDict()

config['analysis_name'] = 'rnaseq_analysis'
config['transcripts_fa'] = '/Users/mdascenzo/workspace/data/Homo_sapiens.GRCh38.rel83.cdna.all.fa'
config['tx2gene_fp'] = '/Users/mdascenzo/workspace/data/tx2gene.EnsDb.Hsapiens.v86.csv'

out = args.out
if not os.path.isabs(args.out):
    out = os.path.join(os.getcwd(), args.out)
config['out'] = out

# make paths absolute if not
fas = [os.path.abspath(f) for f in args.fas]

config['options'] = collections.OrderedDict()
if args.sim:
    config['options']['salmon-quant-l'] = 'IU'

config['samples'] = collections.OrderedDict()

# currently assumes paired end reads
it = iter(natsort.natsorted(fas))
for x in it:
    paired_reads = (x, next(it))
    sample_id = re.sub('_1$', '', os.path.splitext(os.path.splitext(os.path.basename(paired_reads[0]))[0])[0])

    config['samples'][sample_id] = collections.OrderedDict()
    config['samples'][sample_id]['read1'] = paired_reads[0]
    config['samples'][sample_id]['read2'] = paired_reads[1]

config_file = yaml.dump(config, Dumper=yaml.RoundTripDumper, default_flow_style=False).replace('!!omap', '').replace('- ', '')

# write config file
config_file_name = os.path.splitext(args.configfile_name)[0] + '.yml'
if not os.path.exists(config_file_name):
    fout = open(config_file_name, 'w')
    fout.write(config_file)
else:
    print('warning: config file exists.')
