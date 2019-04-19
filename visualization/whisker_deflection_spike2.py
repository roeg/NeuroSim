import os.path
import numpy as np
import matplotlib.pyplot as plt
from neo.io import Spike2IO
from neo.core import AnalogSignal
from scipy.signal import butter, sosfilt
import quantities as pq

# file_groups = {'Test': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/200319_B1_control.txt'}
# fnames_test = (r'/Volumes/TOSHIBA/Muscimol/L6_control_expts/20032019/recording_pip/muscimol experiments/26062018_1211_rec_pip_before_musci_B2_trial1_Data.smr',)

# 200319
# file_groups = {'B2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/200319/200319_B2_control.txt',
#                'C2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/200319/200319_C2_control.txt',
#                'D3 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/200319/200319_D3_control.txt',
#                'B2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/200319/200319_B2_muscimol.txt',
#                'C2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/200319/200319_C2_muscimol.txt',
#                'D3 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/200319/200319_D3_muscimol.txt'}

# 220319
# file_groups = {'B2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/220319/220319_B2_control.txt',
#                'C2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/220319/220319_C2_control.txt',
#                'D2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/220319/220319_D2_control.txt',
#                'B2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/220319/220319_B2_muscimol.txt',
#                'C2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/220319/220319_C2_muscimol.txt',
#                'D2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/220319/220319_D2_muscimol.txt'}

# 030419
# file_groups = {'C1 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/030419/030419_C1_control.txt',
#                'D2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/030419/030419_D2_control.txt',
#                'E2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/030419/030419_E2_control.txt',
#                'C1 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/030419/030419_C1_muscimol.txt',
#                'D2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/030419/030419_D2_muscimol.txt',
#                'E2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/030419/030419_E2_muscimol.txt'}

# 020419
# file_groups = {'B2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/020419/020419_B2_control.txt',
               # 'C2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/020419/020419_C2_control.txt',
               # 'D2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/020419/020419_D2_control.txt',
               # 'B2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/020419/020419_B2_muscimol.txt',
               # 'C2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/020419/020419_C2_muscimol.txt',
               # 'D2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/020419/020419_D2_muscimol.txt'}

# 100419
# file_groups = {'C2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/100419/100419_C2_control.txt',
#                'D2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/100419/100419_D2_control.txt',
#                'E2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/100419/100419_E2_control.txt',
#                'C2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/100419/100419_C2_muscimol.txt',
#                'D2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/100419/100419_D2_muscimol.txt',
#                'E2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/100419/100419_E2_muscimol.txt'}

# 120419
file_groups = {'B2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/120419/120419_B2_control.txt',
               'C2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/120419/120419_C2_control.txt',
               'D2 control': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/120419/120419_D2_control.txt',
               'B2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/120419/120419_B2_muscimol.txt',
               'C2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/120419/120419_C2_muscimol.txt',
               'D2 muscimol': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_L6b/120419/120419_D2_muscimol.txt'}


def _load_group_fnames():
    file_group_names = {}
    for group in file_groups:
        fnames = []
        with open(file_groups[group], 'r') as f:
            for line in f:
                strip_line = line.strip()
                if not len(strip_line):
                    continue
                else:
                    fnames.append(strip_line)
        file_group_names[group] = fnames

    return file_group_names


def _load_whisker_deflection_signals_spike2(fname):
    # channel indices in Rajeev's recordings
    VM = 1
    LFP = 4
    WHISKER = 1
    r = Spike2IO(fname)
    blocks = r.read()
    analogsignals = blocks[0].segments[0].analogsignals
    events = blocks[0].segments[0].events
    return analogsignals[VM], analogsignals[LFP], events[WHISKER]


def _change_group_name_suffix(group, suffix):
    if '.' in suffix:
        offset = -4
    else:
        offset = -3
    fname = file_groups[group]
    pdf_name = fname[:offset] + suffix
    return pdf_name


def _save_psth(group, bins, psth):
    out_name = _change_group_name_suffix(group, '.csv')
    with open(out_name, 'w') as psth_file:
        header = 'Bin center (ms)\tAPs/stim/ms\n'
        psth_file.write(header)
        bins *= 1e3
        bin_center = 0.5*(bins[:-1] + bins[1:])
        for i in range(len(bin_center)):
            line = str(bin_center[i])
            line += '\t'
            line += str(psth[i])
            line += '\n'
            psth_file.write(line)


def _filter_trace(signal):
    filter_order = 3 # Daniel used 6th order with forward only filter
    fs = signal.sampling_rate.magnitude
    low = 300.0
    high = 5000.0
    sos = butter(filter_order, (low, high), btype='band', fs=fs, output='sos')
    print 'Filtering signal...'
    tmp_signal = signal.magnitude.flatten()
    filtered_signal_ = sosfilt(sos, tmp_signal)
    print '...done!'
    filtered_signal = AnalogSignal(filtered_signal_, units=signal.units, sampling_rate=signal.sampling_rate,
                                   t_start=signal.t_start)
    return filtered_signal


def _detect_spikes(signal, threshold):
    t_refract = 0.002*pq.sec # refractory period
    # threshold *= signal.units
    threshold_crossings = []
    tmp = signal.magnitude
    for i in range(1, len(tmp)):
        if tmp[i - 1] <= threshold and tmp[i] > threshold:
            threshold_crossings.append(signal.times[i])
    # throw out ISIs below refractory period
    ISIs = np.diff(threshold_crossings)
    # shift everything left by 1 index because ISI[i] is between spike[i] and spike[i+1]
    # but we want to throw out spike[i+1]
    threshold_crossings = np.roll(threshold_crossings, -1)
    threshold_crossings = threshold_crossings[np.where(ISIs > t_refract)]
    return np.roll(threshold_crossings, 1)


def view_recordings(fnames, detect_spikes=True):
    plt.ion()
    plt.show()

    # load and filter signal
    events = []
    spikes = []
    for i, name in enumerate(fnames):
        vm, lfp, whisker_deflections = _load_whisker_deflection_signals_spike2(name)
        vm_filtered = _filter_trace(vm)
        events.append(whisker_deflections)

        xlims = np.min(vm.times), np.max(vm.times)
        fig = plt.figure(i + 1)
        ax1 = plt.subplot(3, 1, 1)
        ax1.plot(vm.times, vm, 'k', linewidth=0.5)
        ax1.set_xlim(xlims)
        ax1.set_ylabel('Vm (mV)')
        ax2 = plt.subplot(3, 1, 2, sharex=ax1)
        ax2.plot(vm_filtered.times, vm_filtered, 'k', linewidth=0.5)
        helper_lines = [i*0.5 for i in range(1, int(np.max(vm_filtered.magnitude)/0.5))]
        for line_height in helper_lines:
            ax2.plot(xlims, [line_height, line_height], 'r--', linewidth=0.5)
        ax2.set_xlim(xlims)
        ax2.set_ylabel('Vm filt. (mV)')
        # ax3 = plt.subplot(5, 1, 4, sharex=ax1)
        # ax3.plot(lfp.times, lfp, 'b', linewidth=0.5)
        # ax3.set_xlim(xlims)
        # ax3.set_ylabel('LFP (mV)')
        # ax4 = plt.subplot(4, 1, 4, sharex=ax1)
        # ax4.eventplot(whisker_deflections, linewidths=0.5)
        # ax4.set_ylabel('Whisker')
        # ax4.set_xlabel('Time (s)')
        # ax4.set_xlim(xlims)

        plt.draw()
        plt.pause(0.1)

        if detect_spikes:
            threshold = None
            while threshold is None:
                threshold = raw_input('Spike threshold: ')
                try:
                    threshold = float(threshold)
                except:
                    threshold = None
                    print 'Please enter only a number (in mV)'

            spike_times = _detect_spikes(vm_filtered, threshold)
            spikes.append(spike_times)

            ax5 = plt.subplot(3, 1, 3, sharex=ax1)
            ax5.eventplot(spike_times, linewidths=0.5)
            ax5.set_ylabel('Spikes')
            ax5.set_xlim(xlims)
            plt.draw()
            plt.pause(0.1)

        raw_input('Continue...')
        plt.close(fig)

    plt.ioff()
    plt.show()

    return events, spikes


def event_aligned_spikes(events, spikes, group):
    t_before = -0.12 # CDK 2007
    t_after = 0.1

    all_aligned_times = []
    trial_aligned_times = []
    nr_trials = 0
    fig = plt.figure(1)
    ax1 = plt.subplot(2, 1, 1)
    ax1.set_title(group)
    for i in range(len(events)):
        file_events = events[i]
        for event in file_events:
            nr_trials += 1.0
            aligned_times_ = np.array(spikes[i]) - event.magnitude
            aligned_times = aligned_times_[np.where((aligned_times_ >= t_before)*(aligned_times_ <= t_after))]
            trial_aligned_times.append(aligned_times)
            all_aligned_times.extend(aligned_times)

    ax1.eventplot(trial_aligned_times, linewidths=1.0, colors='k')
    ax1.set_xlim([t_before, t_after])
    ax1.set_ylabel('Trial')

    bin_size = 0.001
    n_bins = int((t_after - t_before)/bin_size + 0.5) + 1
    bins = np.linspace(0, t_after - t_before, n_bins)
    bins += t_before
    psth, _ = np.histogram(all_aligned_times, bins=bins)
    psth = psth/nr_trials
    ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    ax2.step(bins[:-1], psth, where='post')
    ax2.set_xlim([t_before, t_after])
    ax2.set_ylabel('Spikes/stim/ms')
    ax2.set_xlabel('Time (s)')

    pdf_name = _change_group_name_suffix(group, '.pdf')
    plt.savefig(pdf_name)
    _save_psth(group, bins, psth)

    plt.show()


def event_aligned_lfp(fnames, group):
    # load all files and automatically extract event-aligned LFP
    # compute event-aligned mean and SE
    # display per file and total and save
    all_traces = []
    traces_per_file = {}
    lfp_window_duration = 0.1 # in s
    for i, name in enumerate(fnames):
        vm, lfp, whisker_deflections = _load_whisker_deflection_signals_spike2(name)
        whisker_deflection_samples = np.round(whisker_deflections.magnitude*lfp.sampling_rate.magnitude)
        lfp_window_samples = int(lfp_window_duration*lfp.sampling_rate.magnitude + 0.5)
        traces_per_file[name] = []
        for deflection in whisker_deflection_samples:
            deflection = int(deflection)
            snippet = lfp[deflection: deflection + lfp_window_samples]
            all_traces.append(snippet)
            traces_per_file[name].append(snippet)

    nr_subplots = len(traces_per_file) + 1
    plt.figure(1)
    ax1 = plt.subplot(nr_subplots, 1, 1)
    total_mean = np.mean(all_traces, axis=0)
    total_se = np.std(all_traces, axis=0)/np.sqrt(len(all_traces))
    time_axis = all_traces[0].times - np.min(all_traces[0].times)
    xlims = (np.min(time_axis), np.max(time_axis))
    ax1.plot(time_axis, total_mean, 'r')
    ax1.plot(time_axis, total_mean + 2*total_se, 'r--', linewidth=0.5)
    ax1.plot(time_axis, total_mean - 2*total_se, 'r--', linewidth=0.5)
    ax1.set_xlim(xlims)
    ax1.set_ylabel('LFP (mV)')
    ax1.set_title('%s mean LFP' % group)

    axes = []
    names = traces_per_file.keys()
    names.sort()
    for i, name in enumerate(names):
        file_traces = traces_per_file[name]
        file_mean = np.mean(file_traces, axis=0)
        file_se = np.std(file_traces, axis=0)/np.sqrt(len(file_traces))
        tmp_ax = plt.subplot(nr_subplots, 1, i + 2, sharex=ax1)
        tmp_ax.plot(time_axis, file_mean, 'k')
        tmp_ax.plot(time_axis, file_mean + 2*file_se, 'k--', linewidth=0.5)
        tmp_ax.plot(time_axis, file_mean - 2*file_se, 'k--', linewidth=0.5)
        tmp_ax.set_xlim(xlims)
        tmp_ax.set_ylabel('LFP (mV)')
        if i == len(traces_per_file) - 1:
            tmp_ax.set_xlabel('Time post-stim. (s)')
        title_str = os.path.split(name)[-1]
        tmp_ax.set_title(title_str)
        axes.append(tmp_ax)

    plt.savefig(_change_group_name_suffix(group, '_LFP.pdf'))
    plt.show()


def main():
    mode = raw_input('Spike analysis (S), LFP (L) or file viewer (V): ')
    if mode == 'S':
        file_group_names = _load_group_fnames()
        for group in file_group_names:
            fnames = file_group_names[group]
            events, spikes = view_recordings(fnames)
            event_aligned_spikes(events, spikes, group)
    elif mode == 'L':
        file_group_names = _load_group_fnames()
        for group in file_group_names:
            fnames = file_group_names[group]
            event_aligned_lfp(fnames, group)
    elif mode == 'V':
        fname = raw_input('Filename: ')
        view_recordings([fname], detect_spikes=False)


if __name__ == '__main__':
    main()