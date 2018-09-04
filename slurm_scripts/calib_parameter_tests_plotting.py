#!/usr/bin/env python3
import sys
import numpy as np
from plotly.offline import plot
import plotly.graph_objs as go

tsv_path = sys.argv[1]
tsv_lines = open(tsv_path).readlines()

field_to_idx = dict()
for idx, field in enumerate(tsv_lines[0].rstrip().split('\t')):
    field_to_idx[field] = idx

parameters_to_performance = dict()
for line in tsv_lines[1:]:
    line = line.rstrip().split('\t')

    e = int(line[field_to_idx['barcode_error_tolerance']])
    k = int(line[field_to_idx['kmer_size']])
    m = int(line[field_to_idx['minimizers_num']])
    t = int(line[field_to_idx['minimizers_threshold']])
    key = (e, k, m, t)
    if not key in parameters_to_performance:
        parameters_to_performance[key] = dict(ari=list(), time=list(), mem=list())

    ari  = float(line[field_to_idx['ARI']])
    time = float(line[field_to_idx['user_time']])
    mem  = float(line[field_to_idx['mem']])
    parameters_to_performance[key]['ari'].append(ari)
    parameters_to_performance[key]['time'].append(time)
    parameters_to_performance[key]['mem'].append(mem)

m_to_color = {
    3 : '#e41a1c',
    4 : '#377eb8',
    5 : '#4daf4a',
    6 : '#984ea3',
    7 : '#ff7f00',
}

e_to_dash = {
    1 : 'dash',
    2 : 'dot' ,
}

k_to_marker_symbol = {
    4 : 'circle' ,
    8 : 'square',
}

data = list()
for key, val in parameters_to_performance.items():
    e = key[0]
    k = key[1]
    m = key[2]
    t = key[3]

    ari_list = val['ari']
    time_list = val['time']
    mem_list = val['mem']
    min_ari = np.min(ari_list)*100
    max_ari = np.max(ari_list)*100
    avg_ari = np.average(ari_list)*100
    min_time = np.min(time_list)/60
    max_time = np.max(time_list)/60
    avg_time = np.average(time_list)/60
    min_mem = np.min(mem_list)/1024/1024
    max_mem = np.max(mem_list)/1024/1024
    avg_mem = np.average(mem_list)/1024/1024
    if min_ari < 99:
        continue

    text='e = {}; k = {}; m = {}; t = {}<br>Mem (GB): {:4.2f}; {:4.2f}; {:4.2f}<br>User time (min): {:4.2f}; {:4.2f}; {:4.2f}<br>ARI (%): {:4.2f}; {:4.2f}; {:4.2f}'.format(
            e, k, m, t,
            min_mem, max_mem, avg_mem,
            min_time, max_time, avg_time,
            min_ari, max_ari, avg_ari)
    trace = go.Scatter(
        name='e = {}; k = {}; m = {}; t = {}'.format(e, k, m, t),
        x=[max_time],
        y=[min_ari],
        text=[text],
        mode='markers',
        marker=dict(
            symbol= k_to_marker_symbol[k],
            color = m_to_color[m],
            size  = int(avg_mem)
        ),
    )
    data.append(trace)
layout = go.Layout(
    title=tsv_path.split('/')[-1],
    hovermode='closest',
    font=dict(
        size=10,
    ),
    xaxis=dict(
        title='Time in min',
        type='log',
    ),
    yaxis=dict(
        title='ARI score',
    ),
)
fig = go.Figure(data=data, layout=layout)
plot(fig, filename='{}.html'.format(tsv_path), auto_open=False)
