function [n1_peak_sample,n1_peak_amplitude] = ccep_detect_n1peak_ECoG(average_ccep,good_channels,params)

% Function for detecting the the N1 peaks of CCEPs

% INPUTS:
% - average_ccep
%   array with averaged signals
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


amplitude_thresh = params.amplitude_thresh;
n1_peak_range = params.n1_peak_range;
epoch_prestim_length = params.epoch_prestim_length;
epoch_length = params.epoch_length;
srate = params.srate;

%% Script

% output in channels X stimulations X [latency amplitude]
n1_peak = NaN(size(average_ccep,1), size(average_ccep,2),2);

% for every averaged stimulation
for jj = 1:size(average_ccep,2)
    % for every channel
    for ii = 1:size(average_ccep,1)
        
        % if channel is a bad or not-recording electrode
        if ~ismember(ii,good_channels)
            n1_peak_sample = NaN;
            n1_peak_amplitude = NaN;
            
            % in other electrode
        else
            % create time struct
            tt = (1:epoch_length*srate) / ...
                srate - epoch_prestim_length;
            
            baseline_tt = tt>-2 & tt<-.1;
            
            if all(isnan(average_ccep(ii,jj,:))) % for example when stimulated
                n1_peak_sample = NaN;
                n1_peak_amplitude = NaN;
                
                % in other electrode
            else
                
                one_signal = squeeze(average_ccep(ii,jj,:));
                % testplot new signal: plot(tt,squeeze(one_signal))
                
                % take area before the stimulation of the new signal and calculate its SD
                pre_stim_sd = std(one_signal(baseline_tt));
                
                % if the pre_stim_sd is smaller that the minimally needed SD,
                % which is validated as 50 uV, use this the minSD as pre_stim_sd
                if pre_stim_sd < 50
                    pre_stim_sd = 50;
                end
                
                % use peakfinder to find all positive and negative peaks and their
                % amplitude.
                % tt are the samples of the epoch based on the Fs and -2.5 to 2.5
                % seconds sample of the total epoch
                % As tt use first sample after timepoint 0
                % till first sample after 0,5 seconds (rougly 1000 samples)
                % sel = 20 , which is how many samples around a peak not considered as another peak
                
                % UPDATED: the range is already set to 9ms:500ms
                % post-stimulation. In RESP0892, there is a clear switch
                % observed after stimulation. This switch is at 9.3ms, and
                % was until now detected as a N1 peak. With changing the
                % range here, this is not detected as N1 peak anymore.
                [all_sampneg, all_amplneg] = ccep_peakfinder(one_signal(find(tt>9/1000,1):find(tt>0.5,1)),20,[],-1);
                
                % If the first selected sample is a peak, this is not a real peak,
                % so delete
                all_amplneg(all_sampneg==1) = [];
                all_sampneg(all_sampneg==1) = [];
                
                % convert back timepoints based on tt, substract 1 because before
                % the first sample after stimulation is taken
                all_sampneg = all_sampneg + find(tt>9/1000,1) - 1;
                
                % set the starting range in which the algorithm looks
                % for a peak. At least 18 samples are necessary because
                % of the micromed amplifier does not record the
                % stimulated electrode before this. Peak detection
                % start 9 ms after stimulation, which is 19 samples
                % after stimulation
                n1_samples_start = find(tt>0.009,1);
                
                % find first sample that corresponds with the given n1
                % peak range
                n1_samples_end = find(tt>(n1_peak_range/1000),1);
                
                
                % for N1, first select the range in which the N1 could appear, and
                % select the peaks found in this range
                temp_n1_peaks_samp = all_sampneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));
                temp_n1_peaks_ampl = all_amplneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));
                
                % if peak(s) found, select biggest peak
                if ~isempty(temp_n1_peaks_samp)
                    max_n1_ampl = find(abs(temp_n1_peaks_ampl) == max(abs(temp_n1_peaks_ampl)));
                    n1_peak_sample = temp_n1_peaks_samp(max_n1_ampl(1));
                    n1_peak_amplitude = temp_n1_peaks_ampl(max_n1_ampl(1));
                    % otherwise give the amplitude the value NaN
                elseif isempty(temp_n1_peaks_samp)
                    n1_peak_amplitude = NaN;
                    n1_peak_sample = NaN;
                end
                
                % if N1 exceeds positive threshold, it is deleted
                if temp_n1_peaks_ampl > 0
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;
                end
                
                % when peak amplitude is saturated, it is deleted
                if abs(n1_peak_amplitude) > 3000
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;
                end
                
                % if the peak is not big enough to consider as a peak, assign NaN
                if abs(n1_peak_amplitude) < amplitude_thresh* abs(pre_stim_sd)
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;
                end
            end
            
        end
        % add properties to output frame
        n1_peak(ii,jj,1) = n1_peak_sample;
        n1_peak(ii,jj,2) = n1_peak_amplitude;
        
    end
end

% write n1_peak (sample and amplitude) to database
n1_peak_sample = n1_peak(:,:,1);
n1_peak_amplitude = n1_peak(:,:,2);

