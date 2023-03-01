% -*- MATLAB -*-
% ----------------
% first_line_match only works for files without extension
% see also: https://packagecontrol.io/packages/ApplySyntax
%
% Prominence Spiker
% (or BEAT IT ver.2)
%
% //__FeAR__//
%
% 21/02/2022
%
% An Octave/MATLAB script for the Multiparametric Description of Spiking Traces.
% Originally developed for the analysis of intracellular calcium oscillations in
% terms of basal levels, peak amplitude, and frequency values of spontaneously
% contracting IPS differentiated towards a muscular phenotype.
%
% NOTE: All along this script an important assumption upon the nature of the
% signals is made: within each experiment, all traces have a perfectly
% synchronous firing activity. Thus it follows that:
%   --> no synchronicity evaluation (e.g. by cross-correlation) has been
%       implemented, being it always equal to 1;
%   --> within each single experiment, the mean of the per-trace average
%       frequencies (i.e. 'freq_1') should be the same of the average frequency
%       of the mean trace (i.e. 'freq_2' and 'freq_3'). So a consensus strategy
%       can be adopted. When the three estimates are largely diverging,
%       asynchronicity may be suspected.
%
% QUICK START
% -----------
% 1. Start Octave
% 2. Move to the script-containing directory
%     @ DBIOS
%       cd 'D:\\UniTo Drive\\WORKS\\0010 - Ongoing\\Myo-IPS\\Custom Scripts'
%     @ HOME
%       cd 'E:\\UniTo Drive\\WORKS\\0010 - Ongoing\\Myo-IPS\\Custom Scripts'
% 3. Run the script beatit.m





% Preliminaries ----------------------------------------------------------------

% Do the cleaning
close all             % Close all figures
clear                 % Clear local variables
clear global          % Clear also global vars
diary off             % Close possible logs left open (because of crashed runs)

% Load packages
pkg load signal       % provides findpeaks() function
pkg load statistics   % provides boxplot() and normpdf() functions

% File separator for current platform
sep = filesep;

% Colors as globals
global magenta80  = [0.65   0.25    0.65];
global cyan80     = [0.10   0.75    0.75];
global darkBG     = [0.10   0.10    0.10];





% Settings ---------------------------------------------------------------------

delta = 750;          % Time window of interest (in seconds)
detrendFlag = true;   % Detrend data?
polyOrd = 4;          % Order of the polynomial to fit for detrending
promiFlag = true;     % Prominence-based peak filter?
diagnostics = false;  % Diagnostics plots (restricted to just 1 experiment)
exp_id = 4;           % Experiment to explore in case of diagnostics == true
fmt = 'png';          % Graphic output format
saveplot = true;      % Save plot to disk in fmt format?

% logical2str for logging
if (detrendFlag)
  detrendStr = ['order ', num2str(polyOrd)];
else
  detrendStr = 'OFF';
end
if (promiFlag) promiStr = 'ON'; else promiStr = 'OFF'; end
if (diagnostics) diagnosticStr = 'ON'; else diagnosticStr = 'OFF'; end





% Utility function definition --------------------------------------------------

function superPlot(t_vec, data, tit, xlab, ylab)
  
  global darkBG
  global cyan80
  c = size(data,2);
  
  figure, hold on

  if (c == 1)
    plot(t_vec, data, 'color', cyan80)
  else
    for i = 1:c
      p = plot(t_vec, data(:,i));

      % Color gradient ([~hue] .* [~brightness])
      set(p,'color',[(c-i)/(c-1), (i-1)/(c-1), 1] .* [0.65, 0.65, 0.65])
    end
  end
  set(gca,'color', darkBG) % Eighties cool!
  xlim([t_vec(1), t_vec(end)])
  title(tit), xlabel(xlab), ylabel(ylab), hold off
end


function quartiles = superBoxplot(data, tit, xlab, ylab)

  global magenta80
  global cyan80
  global darkBG
  c = size(data,2);

  % Use boxplot() to plot data distribution and to get the five-number summary
  % for each trace:
  %   Row 1:   Minimum
  %   Row 2:   1st quartile
  %   Row 3:   Median         --> ~ Basal value
  %   Row 4:   3rd quartile   --> ~ Noise threshold  
  %   Row 5:   Maximum        --> ~ Peak level
  %   Row 6:   Lower confidence limit for median
  %   Row 7:   Upper confidence limit for median
  %
  % 'data' is supposed to be a matrix with one column for each data set
  figure, hold on
  quartiles = boxplot(data);
  
  obj = get(gca,'children');
  % Objects meaning:
  %   1           = outliers
  %   2           = points beyond whiskers
  %   2+[1:c]     = median lines
  %   2+[1:c]+c   = upper whisker end ticks
  %   2+[1:c]+2*c = lower whisker end ticks
  %   2+[1:c]+3*c = upper whiskers
  %   2+[1:c]+4*c = lower whiskers
  %   2+(1:c)+5*c = boxes
  set(obj(1:end), 'color', cyan80, 'LineWidth', 1.0)
  set(obj(1), 'color', magenta80, 'LineWidth', 0.5)
  set(obj(2), 'color', magenta80, 'LineWidth', 0.5)
  set(obj(2+[1:c]), 'color', magenta80, 'LineWidth', 1.5)
  set(gca,'color', darkBG) % Eighties cool!

  % Superimpose mean and dispersion of median
  m = mean(quartiles(3,:));
  s = std(quartiles(3,:));
  plot([0 c+1], [m m], '-w')
  plot([0 c+1], m + [s s], '--w')
  plot([0 c+1], m - [s s], '--w')
  xlim([0 c+1]), title(tit), xlabel(xlab), ylabel(ylab), hold off

  % Change the stacking order to move white lines behind the boxes
  %obj = get(gca,'children'); % Now obj(1:3) are the new white lines
  %set(gca, 'children', [obj(4:end); obj(1:3)]) % This makes labels disappear...
end


% Get the name of a variable as a string!
function out = getVarName(var)
  out = inputname(1);
end


% Set precision (in terms of non-zero decimal places) of floating point numbers
% Useful for comparisons between floating point numbers!
function short = trunc(long, digit)
  short = round(long*10^digit)/(10^digit);
end


% Get current date and time in string format
function rightNow = getTime(form)

  switch form
    case 'long'   % To log
      string_form = 'dd-mmm-yyyy HH:MM:SS';
    case 'short'  % To file name
      string_form = 'yyyy-mm-dd_HH-MM-SS';
  end
  
  rightNow = datestr(fix(clock()), string_form);
end


% This function detect peaks with a prominence lower than 'minp' value and
% returns their index within 'locs' vector
%
% sig   [array]   -> 1D signal
% locs  [array]   -> peak vector (in sample unit)
% minp  [scalar]  -> minimum prominence
function upp = prominence(sig, locs, minp)

  pkn = size(locs,1); % Number of peaks

  % Find inter-peak minima (valleys)
  lp = [];
  lp(1) = sig(locs(1)) - min(sig(1:locs(1)));
  for k = 2:pkn
    % left prominence for the i-th peak
    lp(k) = sig(locs(k)) - min(sig(locs(k-1):locs(k)));
  end
  rp = [];
  for k = 1:pkn-1
    % right prominence for the i-th peak
    rp(k) = sig(locs(k)) - min(sig(locs(k):locs(k+1)));
  end
  rp(pkn) = sig(locs(pkn)) - min(sig(locs(pkn):end));

  % Interleave left and right prominences
  proseq = [];
  proseq([1:2:2*pkn]) = lp; % assign left prominences to odd indexes
  proseq([2:2:2*pkn]) = rp; % assign right prominences to even indexes

  logicPro = proseq < minp; % Underprominent peak cluster are 'islands' of ones

  csi = ceil(find(diff([0 logicPro 0]) == 1)/2); % Cluster start (peak indexes)
  cei = ceil((find(diff([0 logicPro 0]) == -1)-1)/2); % Cluster end (peak indexes)

  upp = []; % Under prominent peak indexes (to remove)
  if (size(csi,2) ~= size(cei,2))
    fprintf('\nSomething went wrong with prominence evaluation...');
    fprintf('\nLow prominence values have been detected somewhere,');
    fprintf('\nhowever, no peaks have been removed.\n');
  else
    for k = 1:size(csi,2)  
      clust = csi(k):cei(k);
      upp = [upp clust];

      if (length(clust) > 1)
        % Within underprominent peak clusters, retain only the highest peak
        [a b] = max(sig(locs(clust))); % [max index]
        upp(upp == clust(b)) = [];
      end
    end
  end
end


% Minimal demo-test for prominence() function
function prominence_Test()
  x = [1:0.1:40]';
  mu = [41 96 121 146 196 221 291 326]';
  coeff = [7 5 7 6 6 5.9 8 1];
  y = 0; % Sum of Gaussians
  for i = 1:length(coeff)
    y = y + coeff(i)*normpdf(x,x(mu(i)),0.9);
  end
  figure, plot(x,y), hold on
  [pks locs] = findpeaks(y);
  upp = prominence(y, locs, 1); display(upp)
  % Low prominent peak to be filtered out
  plot(x(locs(upp)), pks(upp), 'xr')
  xlabel('x'), ylabel('y'), hold off
end





% Start procedural scripting ---------------------------------------------------

% GUI Folder Selector
selpath = uigetdir(strcat(pwd,'\\..\\Data'), 'Select Group Folder');

% Create the output folder and move there
mkdir(selpath, 'Output');
startingPath = pwd;
cd([selpath, sep, 'Output']);

% List raw-data files
listing = dir(strcat(selpath,'\\*.txt')); % Structure array
nfiles = size(listing,1);

% Get basename
% Rawdata are supposed to live into folder named as the experimental group
[noUse1, group, noUse2] = fileparts(selpath);

% Log Command Window text to file
log_name = strcat('LogFile_', group, '_', getTime('short'), '.log');
diary(log_name)

% NOTE: strcat() removes trailing ASCII whitespace characters. To preserve them
%       character arrays also can be concatenated using square brackets.
fprintf(strcat('\n', ['Analysis started on: ', getTime('long')]))
fprintf(strcat('\n', ['Experimental Group:  ', group]))
fprintf(strcat('\n', ['Sample Size:         n=', num2str(nfiles)]))
fprintf(strcat('\n', ['Options:']))
fprintf(strcat('\n', [' |__ Delta Time  = ', num2str(delta)], ' s'))
fprintf(strcat('\n', [' |__ Detrend     = ', detrendStr]))
fprintf(strcat('\n', [' |__ Prominence  = ', promiStr]))
fprintf(strcat('\n', [' |__ Diagnostics = ', diagnosticStr]))
fprintf('\n--------------------------------------------------\n')

if (diagnostics)
  nfiles = 1;
end
% For-loop over files found in folder 'selpath'
for j = 1:nfiles

  if (diagnostics)
    j = exp_id;
  end

  filename = strcat(selpath, sep, listing(j).name);

  % Read tab-separated value files, skipping the first row (headings)
  data = dlmread(filename,'\t',1,0);

  % Get sampling time (in seconds)
  dt = mode(diff(data(:,1)));
  if (diagnostics)
    plot(diff(data(:,1)))
    xlim([1 size(data,1)-1]), ylim([0 3])
    title('Sampling time diagnostics plot')
    xlabel('time (s)'), ylabel('dt')
  end 

  % Get rid of the time column...
  data(:,1) = [];

  % Get data matrix size
  [r,c] = size(data);
  fprintf('\nExperiment %u/%u: %s', j, nfiles, listing(j).name)
  fprintf('\n |__ ROI number:    %u', c)
  fprintf('\n |__ Time points:   %u', r)
  fprintf('\n |__ Sampling time: %.4f s', dt)
  fprintf('\n |__ Total duration: %.2f s', (r-1)*dt)
  fprintf('\n\nData Preview:\n')
  disp(data(1:5, 1:5))

  % ...and build a perfectly-spaced time vector (starting from 0 s)
  t_vec = ([0:r-1]')*dt;

  % Plot raw data
  superPlot(t_vec, data, 'Raw Data', 'time (s)', 'ratio')
  if (saveplot)
    set(gcf, 'InvertHardCopy', 'off') % To print dark background
    figName = [num2str(j), '_1-Raws.', fmt];
    saveas(gcf, figName, fmt)
  end

  % Time window in cycles (last-step)
  lastep = round(delta/dt);

  % Trim data
  data = data(1:lastep,:);
  t_vec = t_vec(1:lastep);

  % Detrend data using a polyOrd-th-order polynomial fit, but keep DC component
  if (detrendFlag)
    data = detrend(data, polyOrd) + mean(data,1);
  end

  % Check possible negative values
  [r_neg, c_neg] = find(data < 0);
  data(:,c_neg) = [];
  if (~ isempty(c_neg))
    fprintf('\nWARNING: Data contains negative values!');
    fprintf('\n         Cannot run findpeaks() on the full dataset.');
    fprintf('\n         %u traces have been eliminated.', length(c_neg));
    fprintf('\n         Otherwise, try lowering polyOrd value...\n');
  end
  [r,c] = size(data);

  % Plot distribution of data points in time and get the five-number summary for
  % each trace:
  quartiles = superBoxplot(data, 'Spatio-Temporal Distribution of Ca^{2+} levels', ...
                           'ROI number', 'ratio');
  if (saveplot)
    set(gcf, 'InvertHardCopy', 'off') % To print dark background
    figName = [num2str(j), '_2-Boxes.', fmt];
    saveas(gcf, figName, fmt)
  end

  % FEATURE_1: Basal value #####################################################
    % mean level
      ave_basal = mean(quartiles(3,:));
    % homogeneity of basal levels (absolute and %-relative)
      homo_basal = std(quartiles(3,:));
      homo_basal_rel = 1e2*homo_basal/ave_basal;

  % Third quartile as the threshold level below which peaks are ignored
  % i.e. points below Q3 are just noise, while those above are the signal
  %thr = quartiles(4,:);
  thr = std(data,1) + mean(data,1); % 1 SD above the mean gives better results

  % Peak detection by findpeaks() function, trace by trace
  pks_detection = []; % This is important, to clean the array within the for loop!
  for i = 1:c
    [pks locs] = findpeaks(data(:,i), 'MinPeakHeight', thr(i));
    if (promiFlag)
      % Remove underprominent peaks
      upp = prominence(data(:,i), locs, (quartiles(5,i)-quartiles(3,i))/10);
      pks(upp) = [];
    end
    if (diagnostics)
      figure, hold on
      plot(t_vec, data(:,i))
      plot([t_vec(1) t_vec(end)], [thr(i) thr(i)], '-r')
      xlim([t_vec(1) t_vec(end)])
      title(['Trace ', num2str(i)])
      xlabel('time (s)'), ylabel('ratio')
    end
    pks_detection(i,:) = [median(pks), size(pks,1)];
  end

  % FEATURE_2a: Frequency (beating rate) #######################################
    % Under the assumption of 'always perfect synchronicity', here is the most
    % reliable estimate of the average beating rate (in mHz):
    event_mode = mode(pks_detection(:,2));
    reliability = 1e2*(sum(pks_detection(:,2) == event_mode)/c); % ~Consensus
    freq_1 = 1e3*(event_mode/delta);
    fprintf('\nAverage frequency (beating rate)')
    fprintf('\n |__ 1st estimate: %.2f mHz', freq_1)
    fprintf(' (%.1f %% reliability)', reliability)

  % Retrieve peak amplitudes (just from traces compliant with freq_1 estimate)
  ok_index = find(pks_detection(:,2) == event_mode);
  peakAmp = pks_detection(ok_index,1) - quartiles(3,ok_index)';

  % FEATURE_3: Peak amplitude ##################################################
    % average peak amplitude
      ave_peakAmp = mean(peakAmp);
    % homogeneity of peak amplitudes (absolute and %-relative)
      homo_peakAmp = std(peakAmp);
      homo_peakAmp_rel = 1e2*homo_peakAmp/ave_peakAmp;

  % Subtract basal (make median == 0)
  subTrace = data - quartiles(3,:);
  % Average trace
  aveTrace = median(subTrace,2);
  % NOTE: Octave's findpeak() only wants positive values! (re-add DC offset)
  offTrace = aveTrace + ave_basal;
  
  % offTrace third quartile
  %thrOff = quantile(offTrace, 0.75); % Third quartile
  thrOff = mean(offTrace) + std(offTrace); % 1 SD above the mean

  % Peak detection by findpeaks() function, on offTrace (average signal)
  [pks locs] = findpeaks(offTrace, 'MinPeakHeight', thrOff);
  if (promiFlag)
    % Remove underprominent peaks
    upp = prominence(offTrace, locs, ave_peakAmp/10);
    uppy = pks(upp);
    uppx = locs(upp);
    pks(upp) = [];
    locs(upp) = [];
  end

  % FEATURE_2b: Frequency (beating rate) #######################################
  % Measure frequency again...
    freq_2 = 1e3*size(pks,1)/delta;
    fprintf('\n |__ 2nd estimate: %.2f mHz', freq_2)
  
  % ..and plot peak detection on the average trace
  % with white crosses over underprominent peaks (already removed from computation) 
  superPlot(t_vec, offTrace, 'Average Trace Peak Detection', 'time (s)', 'ratio')
  hold on
  plot((locs-1)*dt, pks, 'o', 'color', magenta80, 'LineWidth', 1.5)
  plot((uppx-1)*dt, uppy, 'x', 'color', [1,1,1], 'LineWidth', 1.5)
  plot([t_vec(1), t_vec(end)], [thrOff, thrOff], '-w')
  hold off
  if (saveplot)
    set(gcf, 'InvertHardCopy', 'off') % To print dark background
    figName = [num2str(j), '_3-Peaks.', fmt];
    saveas(gcf, figName, fmt)
  end

  % FEATURE_2c: Frequency (beating rate) #######################################
    % Get mean frequency one more time!
    % NOTE:
    %   - When the trace rises above the threshold: diff(offTrace > thrOff) == 1
    %   - When the trace falls below the threshold: diff(offTrace > thrOff) == -1
    %   - Otherwise: diff(offTrace > thrOff) == 0
    % ...no risk of double-counting!
    up_crossing = sum(diff(offTrace > thrOff) == 1);
    freq_3 = 1e3*up_crossing/delta;
    fprintf('\n |__ 3rd estimate: %.2f mHz', freq_3)

  % Frequency Consensus Report
  fprintf('\n                   ----------')
  if (trunc(freq_1,2) == trunc(freq_2,2) && ...
      trunc(freq_1,2) == trunc(freq_3,2))
    fprintf('\n                   Full Consensus')
  elseif (trunc(freq_1,2) ~= trunc(freq_2,2) && ...
          trunc(freq_1,2) ~= trunc(freq_3,2) && ...
          trunc(freq_2,2) ~= trunc(freq_3,2))
    fprintf('\n                   WARNING: Total Discord!')
  else
    fprintf('\n                   WARNING: Minority Report!')
  end

  % Inter-Spike-Interval (ISI) Analysis
  ISI = diff(locs*dt);
  fprintf('\nISI')
  fprintf('\n |__ mean      = %.2f s', mean(ISI))
  fprintf('\n |__ mean^(-1) = %.2f mHz', 1e3/mean(ISI))
  % ...actually a forth way to get an estimate for frequency!

  % FEATURE_2d: Frequency (beating rate) #######################################
    % homogeneity of frequency (absolute and %-relative)
    homo_freq = std(ISI);
    homo_freq_rel = 1e-1*homo_freq*freq_1;
    % ...that is to say: 1e2*homo_freq/(1e-3*freq_1)^(-1)

  figure
  hist(ISI)
  colormap(summer());
  xlabel('ISI')
  ylabel('counts')
  if (saveplot)
    set(gcf, 'InvertHardCopy', 'off') % To print dark background
    figName = [num2str(j), '_4-ISIs.', fmt];
    saveas(gcf, figName, fmt)
  end

  % Store features in the output var
  out_matrix(j,:) = [ave_basal
                     homo_basal
                     homo_basal_rel
                     ave_peakAmp
                     homo_peakAmp
                     homo_peakAmp_rel
                     freq_1
                     homo_freq
                     homo_freq_rel]';
  
  fprintf('\n--------------------------------------------------\n')
end

fprintf('\n')
fprintf(['Analysis finished on: ', getTime('long')])
fprintf('\n--------------------------------------------------\n')

% Turn off the diary and save the log
diary off

% Currently unused
%headers = ['Basal_Average',
%           'Basal_Dispersion',
%           'Basal_RelativeSD',
%           'PeakAmplitude_Average',
%           'PeakAmplitude_Dispersion',
%           'PeakAmplitude_RelativeSD',
%           'Frequency_Average',
%           'Frequency_Dispersion',
%           'Frequency_RelativeSD']'; % Just to be a row...

out_name = strcat('Feature_Matrix_', group, '.txt');
dlmwrite(out_name, out_matrix, '\t')

% Back to the origin
cd(startingPath);
