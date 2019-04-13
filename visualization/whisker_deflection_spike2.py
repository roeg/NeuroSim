import sys
import numpy as np
import matplotlib.pyplot as plt
from neo.io import Spike2IO
from neo.core import AnalogSignal
from scipy.signal import butter, sosfilt
import quantities as pq

file_groups = {'Test': '/Users/robert/Desktop/Neuron_Rev1/Muscimol_SupraGra/200319_B1_control.txt'}
fnames_test = (r'/Volumes/TOSHIBA/Muscimol/L6_control_expts/20032019/recording_pip/muscimol experiments/26062018_1211_rec_pip_before_musci_B2_trial1_Data.smr',)


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


def view_recordings(fnames):
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
        ax1 = plt.subplot(5, 1, 1)
        ax1.plot(vm.times, vm, 'k', linewidth=0.5)
        ax1.set_xlim(xlims)
        ax1.set_ylabel('Vm (mV)')
        ax2 = plt.subplot(5, 1, 2, sharex=ax1)
        ax2.plot(vm.times, vm_filtered, 'k', linewidth=0.5)
        ax2.set_xlim(xlims)
        ax2.set_ylabel('Vm filt. (mV)')
        ax3 = plt.subplot(5, 1, 4, sharex=ax1)
        ax3.plot(lfp.times, lfp, 'b', linewidth=0.5)
        ax3.set_xlim(xlims)
        ax3.set_ylabel('LFP (mV)')
        ax4 = plt.subplot(5, 1, 5, sharex=ax1)
        ax4.eventplot(whisker_deflections, linewidths=0.5)
        ax4.set_ylabel('Whisker')
        ax4.set_xlabel('Time (s)')
        ax4.set_xlim(xlims)

        plt.draw()
        plt.pause(0.1)

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

        ax5 = plt.subplot(5, 1, 3, sharex=ax1)
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


def event_aligned_spikes(events, spikes):
    t_before = -0.12 # CDK 2007
    t_after = 0.1

    all_aligned_times = []
    trial_aligned_times = []
    nr_trials = 0
    fig = plt.figure(1)
    ax1 = plt.subplot(2, 1, 1)
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
    n_bins = int((t_after - t_before)/bin_size + 0.5)
    bins = np.linspace(0, t_after - t_before, n_bins)
    bins += t_before
    psth, _ = np.histogram(all_aligned_times, bins=bins)
    psth = psth/nr_trials
    ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    ax2.step(bins[:-1], psth, where='post')
    ax2.set_xlim([t_before, t_after])
    ax2.set_ylabel('Spikes/stim/ms')
    ax2.set_xlabel('Time (s)')

    plt.show()

def main():
    file_group_names = _load_group_fnames()
    for group in file_group_names:
        fnames = file_group_names[group]
        events, spikes = view_recordings(fnames)
        event_aligned_spikes(events, spikes)

if __name__ == '__main__':
    main()