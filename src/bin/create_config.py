#!/usr/bin/env python
import re
import os
import natsort
import collections
import argparse
import sys
from ruamel import yaml


__version__ = "0.0.1"

usage = """
create_config.py -f <fasta_files> [options]
version: {__version__}

options:
    -o  path for directory name where analysis output will be placed
    -n  name of configuration file

example usage:
create_config.py -f seq/*.fa.gz -o analysis
""".format(**locals())

parser = argparse.ArgumentParser()
parser.add_argument('-f', dest='fas', help='path to fastq files', metavar='fq', type=str,  nargs='+')
parser.add_argument('-o', dest='out', default='.')
parser.add_argument('-s', dest='sim', action='store_true')
parser.add_argument('-n', dest='configfile_name', default='config.yml')

args = parser.parse_args()

# validate input
if not args.fas:
    print(usage)
    sys.exit()

# set default analysis options
config = collections.OrderedDict()
config['analysis_name'] = 'rnaseq_analysis'
config['star'] = 'no'
config['salmon'] = 'yes'
config['trim'] = 'yes'
config['resources_dir'] = '/Volumes/Precyte1/stage/resources'
config['build'] = 'hg38'
config['genome_uid'] = 'hg38dev'
config['tx_uid'] = 'ensembl_rel83'
config['annotation_gtf'] = 'gencode.v25.primary_assembly.annotation.wERCC92.gtf'
# todo: update with known location
config['tx2gene_fp'] =\
    '/Volumes/Precyte1/stage/resources/transcriptomes/hg38/ensembl_rel86/annotation/tx2gene/tx2gene.EnsDb.Hsapiens.v86.csv'

# set default star options
config['star_sj_db_overhang'] = 74

out = args.out
if not os.path.isabs(args.out):
    out = os.path.join(os.getcwd(), args.out)
config['out'] = out

# make paths absolute if not
fas = [os.path.abspath(f) for f in args.fas]

# set options parameter and default values
config['options'] = collections.OrderedDict()
if args.sim:
    config['options']['salmon-quant-l'] = 'IU'

# configure samples
config['samples'] = collections.OrderedDict()

# set reads for each sample, currently assumes paired end reads
it = iter(natsort.natsorted(fas))
for x in it:
    paired_reads = (x, next(it))
    sample_id = re.sub('_1$', '', os.path.splitext(os.path.splitext(os.path.basename(paired_reads[0]))[0])[0])

    config['samples'][sample_id] = collections.OrderedDict()
    config['samples'][sample_id]['read1'] = paired_reads[0]
    config['samples'][sample_id]['read2'] = paired_reads[1]

# create configuration document
config_file = yaml.dump(config, Dumper=yaml.RoundTripDumper, default_flow_style=False).replace('!!omap', '').replace('- ', '')

# write config file
config_file_name = os.path.splitext(args.configfile_name)[0] + '.yml'
if not os.path.exists(config_file_name):
    fout = open(config_file_name, 'w')
    fout.write(config_file)
else:
    print('warning: config file exists, nothing written.')

# todo: consider copying workflows to local dir
# workflow = 'rnaseq'
# os.path.join(sys.prefix, 'envs', os.environ['CONDA_DEFAULT_ENV'], 'workflows', workflow + '.smk')
