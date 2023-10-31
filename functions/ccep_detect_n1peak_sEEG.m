function [n1_peak_sample,n1_peak_amplitude,n1_peak_time] = ccep_detect_n1peak_sEEG(average_ccep,tt,params)

% Function for detecting the the N1 peaks of CCEPs

% INPUTS:
% - average_ccep
%   averaged signal for 1 stim/record pair
% - amplitude_thresh
%   Threshold factor that needs to be exceeded to detect a peak as
%   significant peak. This factor is multiplies the minimal pre simulation
%   standard deviation of 50 uV. amplitude_thresh of 3.4 recommended for
%   conservative algorithm (50 uV * 3.4 = 170 uV). 2.6 is recommended to
%   have a high sensitivity
% - n1_peak_range
%   End point (in ms) of the range in which the algorithm will search for the
%   peak of the N1. Algorithm will search between 10ms and the end point.
%   End point of 90 ms recommended for conservative algorithm (80ms has
%   similar performance.

% OUTPUTS:
% - n1_peak_sample
%   matrix[channels x stimulations] containing latency of the significant
%   n1 peaks that are detected.
% - n1_peak_amplitude
%   matrix[channels x stimulations] containing the amplitude of the
%   significant n1 peaks that are detected.

% Peakfinder function is used for peak detection, add to path

% This function is heavily based on the function and algorithm developed by
% Dorien van Blooijs during her master's thesis 'Improving the SPES protocol
% by automating ER and DR detection and evaluation of the spatial relation
% between ERs and DRs. Algorithm is validated during the thesis and used in
% publication (van Blooijs et al., 2018).

% Some slight differences:
% - While algorithm of van Blooijs detects both positive as negative peaks,
%   this algorithm only selects negative peaks
% - After validation the algorithm of van Blooijs used 125uV as amplitude threshold
%   and 100ms af n1 peak range, this algorithm performes best with 140 uV
%   and 80 ms. A reason for this possibly is the other ages of patients which are
%   used during validation. The first validation is done is only younger
%   (9yrs - 12 yrs.) patients, while later validation is also done in older
%   patients. Age seems to effect the characteristics of CCEPs and
%   therefore also the detection.
% - look at peak range....

% This version of the algorithm is validated in three patients (age 9, 21 and 50)
% Validation done by comparing the performance of the code and visual
% assesment. Algorithm optimized by setting different parameters and
% comparing their performances. (see parameters optimize function).

% For a conservative algorithm (high specificity of at least 95%), the following
% parameters are advised: (see validation matrices for performances with
% other parameters e.g. if you want a very sensitive algorithm)
%<<<<<<< HEAD
% - amplitude threshold of 140 uV (minSD * threshold = 50 uV * 2.8)
%   recommended
% - N1 peak range of (10 to) 70 ms is recommended
%=======
% - amplitude threshold of 170 uV (minSD * threshold = 50 uV * 3.4)
%   recommended
% - N1 peak range of (10 to) 90 ms is recommended
%>>>>>>> upstream/master

% FIXED PARAMETERS (that are validated by both van Blooijs & van der Aar):
% - sel = 20, which is how many samples around peak not considered as another peak
% - minSD = 50, minimum standard deviation (in uV) of prestimulus baseline

% original author: Dorien van Blooijs, UMC Utrecht, January 2018
% modified by: Jaap van der Aar, Dora Hermes, Dorien van Blooijs, Giulio Castegnaro, UMC Utrecht, 2019


amplitude_thresh    = params.amplitude_thresh;
n1_peak_range       = params.n1_peak_range;
peakSign            = params.peakSign; % 1 for positive -1 for negative
baseline_tt         = params.baseline_tt;%tt>-.5 & tt<-.020;
            
% take area before the stimulation of the new signal and calculate its SD
pre_stim_sd         = std(average_ccep(baseline_tt));

% if the pre_stim_sd is smaller that the minimally needed SD,
% which is validated as 50 uV, use this the minSD as pre_stim_sd
if pre_stim_sd < params.amplitudeThresh
    pre_stim_sd = params.amplitudeThresh;
end

% use peakfinder to find all positive and negative peaks and their
% amplitude.
% tt are the samples of the epoch based on the Fs and -2.5 to 2.5
% seconds sample of the total epoch
% As tt use first sample after timepoint 0
% till first sample after 0,5 seconds (rougly 1000 samples)
% sel = 20 , which is how many samples around a peak not considered as another peak

min_peak_tt         = 0.013;

% UPDATED: the range for finding peaks is set from 13ms:500ms
% post-stimulation.
[all_sampneg, all_amplneg] = ccep_peakfinder(average_ccep(find(tt>min_peak_tt,1):find(tt>0.5,1)),20,[],peakSign);

% If the first selected sample is a peak, this is not a real peak,
% so delete
all_amplneg(all_sampneg==1) = [];
all_sampneg(all_sampneg==1) = [];

% convert back timepoints based on tt, substract 1 because before
% the first sample after stimulation is taken
all_sampneg         = all_sampneg + find(tt>min_peak_tt,1) - 1;

% set the starting range in which the algorithm looks
% for a peak. Peak detection starts 13 ms after stimulation
n1_samples_start    = find(tt>min_peak_tt,1);

% find first sample that corresponds with the given n1
% peak range
n1_samples_end      = find(tt>n1_peak_range,1);

% for N1, first select the range in which the N1 could appear, and
% select the peaks found in this range
temp_n1_peaks_samp  = all_sampneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));
temp_n1_peaks_ampl  = all_amplneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));

% if peak(s) found, select biggest peak
if ~isempty(temp_n1_peaks_samp)
    max_n1_ampl         = find(abs(temp_n1_peaks_ampl) == max(abs(temp_n1_peaks_ampl)));
    n1_peak_sample      = temp_n1_peaks_samp(max_n1_ampl(1));
    n1_peak_amplitude   = temp_n1_peaks_ampl(max_n1_ampl(1));
    n1_peak_time        = tt(n1_peak_sample);
    % otherwise give the amplitude the value NaN
elseif isempty(temp_n1_peaks_samp)
    n1_peak_amplitude   = NaN;
    n1_peak_sample      = NaN;
    n1_peak_time        = NaN;
end

if peakSign==1 % if N1 < 0, it is deleted
    if temp_n1_peaks_ampl < 0
        n1_peak_sample  = NaN;
        n1_peak_amplitude = NaN;
        n1_peak_time    = NaN;
    end
elseif peakSign==-1 % if N1 > 0, it is deleted
    if temp_n1_peaks_ampl > 0
        n1_peak_sample  = NaN;
        n1_peak_amplitude = NaN;
        n1_peak_time    = NaN;
    end
end

% if the peak is not big enough to consider as a peak, assign NaN
if abs(n1_peak_amplitude) < amplitude_thresh * abs(pre_stim_sd)
    n1_peak_sample      = NaN;
    n1_peak_amplitude   = NaN;
    n1_peak_time        = NaN;
end

            

