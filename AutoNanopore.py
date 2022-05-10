# -*- coding: utf-8 -*-
# @Author  : Xinlong Liu
# @Time    : 2021/10/23
# @Function: Event detection in solid-state nanopore sequencing signals
import pyabf
import numpy as np
import csv
from collections import defaultdict
import pandas as pd
import argparse
import time



def event_detection(filePath,args):
    print('Processing started')
    start_time = time.time()
    save_path = args.output_path
    window_size = args.window_size
    theta = args.theta
    direc_flag = args.signal_direction  # 0 for positive-going signals and 1 for negative-going ones
    split_time = window_size / 1000
    file_index = filePath.split('/')[-1].split('.')[0]
    abf = pyabf.ABF(filePath)
    sampling_rate = abf.dataRate
    print('sampling rate is %s'%sampling_rate)
    x = abf.sweepX
    y = abf.sweepY
    freq = len(x)  # total number of data points
    print('total number of data point is %s'%freq)
    all_time = freq / sampling_rate  # total length in time
    print('total length of the file is %ss'%all_time)
    point_time = all_time / freq  # sampling interval
    split_points = int(freq * split_time / all_time)
    print('%s slices generated' %int(all_time / split_time))
    if direc_flag == 1:
        y = y * (-1)

    positions = []
    for i in range(0, len(x), split_points):
        positions.append(i)
    positions.append(len(x))
    possible_points = defaultdict(list)
    possible_amplitude = []
    possible_peak = []
    num = 0

    for position_index, position in enumerate(positions[:-1]):
        num += 1
        split_area = y[position:positions[position_index + 1]]  # range for the segment
        peak_value = max(split_area)  # highest peak in the segment
        peak_index = np.where(split_area == peak_value)[0][0]  # index of the highest peak in the segment
        peak_index = position + peak_index  # overall index of the highest peak
        flag = 0
        testtt = 0
        for ikk in range(25, 100):
            temprange = y[peak_index - ikk:peak_index]

            tempcoef = (peak_value - np.mean(temprange)) / abs(np.std(temprange))
            if tempcoef > 3:
                end_index = peak_index - ikk
                flag = 1
                break
            if flag == 0:
                end_index = max(0, peak_index - 125)
                testtt += 1
        start_index = max(0, end_index - 125)
        base_line = np.mean(y[start_index:end_index])

        rise = 0.0
        decay = 0.0
        plot_start_index = 0
        plot_end_index = 0
        for index_ in range(peak_index, start_index, -1):
            if y[index_] < base_line + 0.1 * (y[peak_index] - base_line):
                rise = x[peak_index] - x[index_ + 1]
                plot_start_index = index_ + 1
                break
        for index_ in range(peak_index, positions[position_index + 1]):
            if y[index_] <= base_line + 0.1 * (y[peak_index] - base_line):
                decay = x[index_] - x[peak_index]
                plot_end_index = index_
                break
        if rise > 0.0 and decay > 0.0:
            if direc_flag == 1:
                possible_points[x[peak_index]] += [x[peak_index] * 1000, -1 * peak_value - (-1) * base_line,
                                                   -1 * base_line, point_time * 1000 * plot_start_index,
                                                   point_time * 1000 * plot_end_index,
                                                   rise * 1000 + decay * 1000, rise * 1000, decay * 1000,
                                                   plot_start_index, plot_end_index, peak_index, start_index, end_index]
                possible_amplitude.append(base_line - y[peak_index])
                possible_peak.append(x[peak_index] * 1000)
            else:
                possible_points[x[peak_index]] += [x[peak_index] * 1000, peak_value - base_line, base_line,
                                                   point_time * 1000 * plot_start_index,
                                                   point_time * 1000 * plot_end_index,
                                                   rise * 1000 + decay * 1000, rise * 1000, decay * 1000,
                                                   plot_start_index, plot_end_index, peak_index, start_index, end_index]
                possible_amplitude.append(y[peak_index] - base_line)
                possible_peak.append(x[peak_index] * 1000)

    amplitude_Q = np.percentile(possible_amplitude, [25, 50, 75])
    amplitude_q_1 = amplitude_Q[0]
    amplitude_q_3 = amplitude_Q[2]
    amp_IQR = amplitude_q_3 - amplitude_q_1

    print('%s valid slices selected' %len(possible_points))
    f = open(save_path + '/' + str(file_index) + '.csv', 'w', newline='', )
    f_writer = csv.writer(f)
    f_writer.writerow(
        ['Event index', 'Peak Time (ms)', 'Amplitude (nA)', 'Baseline', 'Start time (ms)', 'End time (ms)',
         'Duration (ms)', 'Rise (ms)', 'Decay (ms)', 'Event start index',
         "Event end index", "peak index", "Baseline start index", "Baseline end index"])
    true_index = 1
    for peak_time, values in possible_points.items():
        if direc_flag == 1:
            if values[1] <= amplitude_q_1 - theta * amp_IQR:
                f_writer.writerow([true_index] + values)
                true_index += 1

        else:
            if values[1] >= amplitude_q_3 + theta * amp_IQR:
                f_writer.writerow([true_index] + values)
                true_index += 1

    f.close()
    
    print('Optimizing ...')
    confindex = []
    an_results = pd.read_csv(save_path + '/' + str(file_index) + '.csv')

    dur_Q = np.percentile(an_results['Duration (ms)'], [25, 50, 75])
    dur_q_1 = dur_Q[0]
    dur_q_3 = dur_Q[2]
    dur_IQR = dur_q_3 - dur_q_1

    an_adjusted = an_results[an_results['Duration (ms)'] <= dur_q_3 + 1.5 * dur_IQR]

    if direc_flag == 1:
        an_adjusted['Amplitude (nA)'] = -1*an_adjusted['Amplitude (nA)']
    an_adjusted = an_adjusted.sort_values(by='Amplitude (nA)', ascending=True).reset_index()
    an_adjusted = an_adjusted.drop(columns='index')
    confidentlist = []
    for ik in range(len(an_adjusted)):
        temppart = an_adjusted['Amplitude (nA)'].iloc[0:ik + 1]
        confidentlist.append((an_adjusted['Amplitude (nA)'][ik] - an_adjusted['Amplitude (nA)'][0]) / np.mean(temppart))
    an_adjusted['quality score'] = confidentlist
    if len(an_adjusted[an_adjusted['quality score'] >= 0.1]) >= 0:
        confi_threshold1 = an_adjusted[an_adjusted['quality score'] >= 0.1]['Amplitude (nA)'].min()
        confi_threshold = confi_threshold1

    else:
        confi_threshold = an_adjusted['Amplitude (nA)'].median()

    an_confident = an_adjusted[an_adjusted['Amplitude (nA)'] >= confi_threshold].reset_index()
    an_confident = an_confident.drop(columns='index')
    if direc_flag == 1:
        an_confident['Amplitude (nA)'] = -1*an_confident['Amplitude (nA)']
    an_confident = an_confident.sort_values(by='Peak Time (ms)', ascending=True).reset_index()
    an_confident = an_confident.drop(columns='index')

    for newindex_ in range(1, len(an_confident) + 1):
        confindex.append(newindex_)

    an_confident['Event index'] = confindex

    an_confident.to_csv(save_path + '/' + str(file_index) + '.csv', index=False)
    end_time = time.time()
    print('Processing finished within %ss, %s events detected.'%(round(end_time-start_time, 3),len(an_confident)))
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="AutoNanopore")
    parser.add_argument("--file_path", type=str, required=True, help="input the abf file path")
    parser.add_argument("--output_path", type=str, required=True, help="Output file path")
    parser.add_argument("--signal_direction", type=int, default=0, help="Direction of event detection, 0 for positive-going and 1 for negative-going")
    parser.add_argument("--theta", type=float, default=1.5, help="multiplication of amplitude IQR")
    parser.add_argument("--window_size", type=int, default=30, help="length of each segment (ms)")
    args = parser.parse_args()

    file_path = args.file_path
    output_path = args.output_path
    window_size = args.window_size
    theta = args.theta
    direct = args.signal_direction
    split_time = window_size / 1000
    direc_flag = direct 
    event_detection(file_path,args)