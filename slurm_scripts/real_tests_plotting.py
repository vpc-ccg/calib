#!/usr/bin/env python3
import sys
import os
import plotly
import plotly.graph_objs as go
import math

results_path=sys.argv[1].rstrip('/')
variants_path=sys.argv[2].rstrip('/')
snps_path=sys.argv[3].rstrip('/')

snps = set()
for line in open(variants_path):
    line = line.rstrip().split('\t')
    chrom = line[0]
    pos = int(line[1])
    ref = line[2].upper()
    alt = line[3].upper()
    snps.add(
        (
            chrom,
            pos,
            ref,
            alt
        )
    )
for line in open(snps_path):
    line = line.rstrip().split('\t')
    chrom = line[0]
    pos = int(line[1])
    ref = line[2].upper()
    alt = line[3].upper()
    snps.add(
        (
            chrom,
            pos,
            ref,
            alt
        )
    )


tools=['calib', 'starcode', 'umitools', 'raw']

names=dict(
    calib='Calib',
    rainbow='Rainbow',
    umitools='UMI-tools',
    starcode='starcode-umi',
    raw='No clustering',
)
legend_names=dict(
    calib='Calib',
    rainbow='Rainbow',
    umitools='UMI-tools',
    starcode='starcode-umi',
    raw='No clustering',
)
seq_errs = set()
for tool in tools:
    for f in os.listdir('{}/{}.sinvict/'.format(results_path, tool)):
        seq_errs.add(f.strip('/'))
seq_errs = list(reversed(sorted(seq_errs)))

mutation_calls = dict()
for seq_err in seq_errs:
    mutation_calls[seq_err] = dict()
    for tool in tools:
        sinvict_path = '{}/{}.sinvict/{}/calls_level1.sinvict'.format(results_path, tool, seq_err)
        mutation_calls[seq_err][tool] = dict()
        for line in open(sinvict_path):
            line = line.rstrip().split('\t')
            chrom = line[0]
            pos = int(line[1])
            ref = line[3].upper()
            alt = line[5].upper()
            mutation_calls[seq_err][tool][(chrom,pos,ref,alt,)]=line

colors_1 = [
    '#e41a1c',
    '#377eb8',
    '#4daf4a',
    '#984ea3',
    '#ff7f00',
]

colors_2 = [
    '#f4a3a4',
    '#afcbe3',
    '#b8dfb7',
    '#d6b8dA',
    '#ffcc99',
]
data=list()
height=1200
width=1200
for idx, tool in enumerate(tools):
    all_counts = list()
    validated_counts = list()
    for seq_err in seq_errs:
        mutation_calls_set = set(mutation_calls[seq_err][tool].keys())
        all_counts.append(len(mutation_calls_set))
        validated_counts.append(len(mutation_calls_set.intersection(snps)))
    trace = go.Bar(
        orientation='h',
        y = seq_errs,
        x = validated_counts,
        text = ['{}'.format(count) for count in validated_counts],
        textposition = 'inside',
        name = '{} verified or reported         '.format(names[tool]),
        marker=dict(
            color=colors_1[idx]
        ),
        legendgroup=tool,
    )
    data.append(trace)
    trace = go.Bar(
        orientation='h',
        y = seq_errs,
        x = all_counts,
        text = ['{}'.format(count) for count in all_counts],
        textposition = 'inside',
        name = '{} all calls'.format(names[tool]),
        marker=dict(
            color=colors_2[idx],
        ),
        legendgroup=tool,
    )
    data.append(trace)
layout = go.Layout(
    barmode='group',
    xaxis=dict(
        type='log',
    ),
    yaxis=dict(
        type='category',
    ),
    font=dict(
        size=24,
    ),
    legend=dict(
        traceorder='grouped+reversed',
        #x=0.50,
        #y=0.05,
        tracegroupgap=20,
    ),
)
fig = go.Figure(data=data, layout=layout)
plotly.offline.plot(fig, filename='{}/plots.html'.format(results_path), auto_open=False, image='svg', image_width=width, image_height=height, image_filename='sinvict')
