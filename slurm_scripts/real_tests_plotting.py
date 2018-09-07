#!/usr/bin/env python3
import sys
import os
import plotly
import plotly.graph_objs as go
import math

results_path=sys.argv[1].rstrip('/')
snps_path=sys.argv[2].rstrip('/')
variants_path=sys.argv[3].rstrip('/')

snps = set()
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

tools = set()
for f in os.listdir(results_path):
    if 'sinvict' in f:
        tool = f.split('.')[0]
        tools.add(tool)
tools = sorted(tools)

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
annotations=list()
for idx, tool in enumerate(tools):
    all_counts = list()
    validated_counts = list()
    for seq_err in seq_errs:
        mutation_calls_set = set(mutation_calls[seq_err][tool].keys())
        all_counts.append(len(mutation_calls_set))
        validated_counts.append(len(mutation_calls_set.intersection(snps)))
    trace = go.Bar(
        x = seq_errs,
        y = all_counts,
        # text = ['<b>{}</b>'.format(count) for count in all_counts],
        text = ['{}'.format(count) for count in all_counts],
        textposition = 'auto',
        name = '{} all calls'.format(tool),
        marker=dict(
            color=colors_2[idx],
        ),
    )
    data.append(trace)
    trace = go.Bar(
        x = seq_errs,
        y = validated_counts,
        # text = ['<b>{}</b>'.format(count) for count in validated_counts],
        text = ['{}'.format(count) for count in validated_counts],
        textposition = 'auto',
        name = '{} verified or reported'.format(tool),
        marker=dict(
            color=colors_1[idx]
        ),
    )
    data.append(trace)
layout = go.Layout(
    title='SiNVICT calls with for different tools and sequencing error rates',
    barmode='group',
    yaxis=dict(
            type='log',
    ),
    xaxis=dict(
            type='category',
    ),
    font=dict(
        size=14,
        # family='Times New Roman',
    ),
    legend=dict(
        x=0.1,
        y=0.9,
    ),
)

fig = go.Figure(data=data, layout=layout)
plotly.offline.plot(fig, filename='{}/plots.html'.format(results_path), auto_open=False)
