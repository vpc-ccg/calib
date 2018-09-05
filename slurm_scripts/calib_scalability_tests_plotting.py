#!/usr/bin/env python3
import sys
import numpy as np
from dateutil import parser
from plotly.offline import plot
import plotly.graph_objs as go


key_to_symbol = {
    (100000 , 'wall_time'): 'circle',
    (100000 , 'mem')      : 'circle-open',
    (1000000, 'wall_time'): 'square',
    (1000000, 'mem')      : 'square-open',
    (2000000, 'wall_time'): 'diamond',
    (2000000, 'mem')      : 'diamond-open',
}

data = list()
mems = list()
wall_times = list()

output_html = sys.argv[1]
for tsv_path in sys.argv[2:]:
    num_molecules = int(tsv_path.split('molNum')[1].split('/')[0])
    tsv_lines = open(tsv_path).readlines()
    field_to_idx = dict()
    for idx, field in enumerate(tsv_lines[0].rstrip().split('\t')):
        field_to_idx[field] = idx
    thread_count_to_performance = dict()
    for line in tsv_lines[1:]:
        line = line.rstrip().split('\t')
        c = int(line[field_to_idx['log_comment']].split('_')[1])
        mem = float(line[field_to_idx['mem']])/1024/1024
        mem = round(mem, 2)
        wall_time = line[field_to_idx['wall_time']].split(':')
        wall_time = list(reversed(wall_time))
        secs = 0
        for idx, t in enumerate(wall_time):
            secs += float(t)*(60**idx)
        mins = secs/60
        print(wall_time, mins)

        thread_count_to_performance[c] = dict(
            wall_time=mins,
            mem=mem,
        )
        mems.append(mem)
        wall_times.append(mins)
    x_list = list()
    y1_list = list()
    y2_list = list()
    for c in sorted(thread_count_to_performance.keys()):
        x_list.append(c)
        y1_list.append(thread_count_to_performance[c]['wall_time'])
        y2_list.append(thread_count_to_performance[c]['mem'])
    trace = go.Scatter(
        name='Time for molNum = {}'.format(num_molecules),
        x=x_list,
        y=y1_list,
        mode='markers+lines',
        marker=dict(
            symbol=key_to_symbol[(num_molecules,'wall_time')],
            size=10,
        ),
        line=dict(
            color='black',
            dash='solid',
        )
    )
    data.append(trace)
    trace = go.Scatter(
        name='RAM for molNum = {}'.format(num_molecules),
        x=x_list,
        y=y2_list,
        yaxis='y2',
        mode='markers+lines',
        marker=dict(
            symbol=key_to_symbol[(num_molecules,'mem')],
            size=10,
        ),
        line=dict(
            color='black',
            dash='longdash',
        )
    )
    print(key_to_symbol[(num_molecules,'mem')],)
    data.append(trace)
mems.sort()
mems = mems[2:]
wall_times.sort()
wall_times = wall_times[2:]
layout = go.Layout(
    legend=dict(
        x=0.85,
        y=0.80,
        borderwidth=2,
    ),
    hovermode='closest',
    font=dict(
        family='Times New Roman',
    ),
    xaxis=dict(
        title='Number of threads',
        rangemode='nonnegative',
        type='linear',
    ),
    yaxis=dict(
        title='Wall clock time',
        rangemode='nonnegative',
        tickmode='array',
        tickvals=wall_times,
        ticktext=['{:4.2f} min'.format(min) for min in wall_times],
        # tickangle=45,
    ),
    yaxis2=dict(
        title='RAM',
        rangemode='nonnegative',
        overlaying='y',
        side='right',
        anchor='y',
        tickmode='array',
        tickvals=mems,
        ticktext=['{:4.2f} GB'.format(mem) for mem in mems],
        # tickangle=45,
    ),
)
fig = go.Figure(data=data, layout=layout)
plot(fig, filename='{}.html'.format('scalibility_plots'), auto_open=False, image='svg', image_filename='scalibility_plots',image_height=900, image_width=1600)
