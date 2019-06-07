#!/usr/bin/env python

# conda install -c conda-forge -n rnaseq ruamel.yaml

import re
import io
import os
import natsort
import glob
import collections
from ruamel import yaml
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', '--fas', dest='fas', help='path to fastq files', metavar='fq', type=str)
parser.add_option('-o', dest='out', default='.')

(options, args) = parser.parse_args()

if not os.path.exists('config.yml'):
    f = io.file('config.yml', 'w')

fas = [f for f in glob.glob(options.fas, recursive=False)]

config = collections.OrderedDict()

config['transcripts_fa'] = '/Users/mdascenzo/workspace/data/Homo_sapiens.GRCh38.rel83.cdna.all.fa'
config['tx2gene_fp'] = '/Users/mdascenzo/workspace/dev/pipelines/src/rnaseq/R/tx2gene.EnsDb.Hsapiens.v86.csv'
config['out'] = '/Users/mdascenzo/workspace/dev/rnaseq-sim/analysis'

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

fout = open(os.path.join(options.out, 'config.yml'), 'w')
fout.write(config_file)
