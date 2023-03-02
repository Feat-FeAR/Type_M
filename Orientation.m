

% Requirements:
%
% 1. CircStat Toolbox for Circular Statistics needs to be present and stored in
%    'circ_path' folder (see below)
%
% 2. The Working Directory (WD) must be the folder containing both the script
%    (OrientScript.m) and the data (MasterMatrix_O.xlsx)
%
% 3. Data in MasterMatrix_O need to be re-saved as .xlsx to be properly loaded
%    with readmatrix()
%
% 4. MasterMatrix_O.xlsx must be organized as it follows:
%       - the first row is the header containing all the sample names
%       - the first column is the angle axis (from 89 to -90 by 1°-step)
%       - rows 182 and 183 (not considering the header row) contain the under-
%         and over-estimate of the sample size for each picture (i.e., the
%         number of cells in the field)
%



% Reset the Workspace
clc
clear all

% Detect the host machine
switch getenv('COMPUTERNAME')
   case 'SILVERLIFE-3'
      root = 'E:'; % @HOME
   case 'SKYLAKE2'
      root = 'D:'; % @DBIOS
end

% Check the Working Directory
% Get an error message and stop the script if WD is not as expected
datapath = '\UniTo Drive\WORKS\0010 - Ongoing\PMMA Orientation';
if ~strcmp(pwd, strcat(root, datapath))
	fprintf('\n')
	error('Wrong WD... aborted.')
	return
end

% Add the path of CircStat Toolbox to the 'MATLAB search path' list
circ_path = strcat(root,...
	'\\UniTo Drive\\Coding\\Z - Third-party Packages\\MATLAB - CircStat\\CircStat');
addpath(circ_path)

% Load Coherency-Weighted Orientation Histograms (MasterMatrix_O) as a table
filename = 'MasterMatrix_O';
cytor_table = readtable(strcat(root, datapath, '\\', filename, '.xlsx'),...
	'FileType', 'spreadsheet',...
	'ReadVariableNames', true,...
	'ReadRowNames', true);

% Set the Angular Domain (for binned data)
a = fliplr([-90:89]); % degs
axa = 2 * a; % axial data
arad = circ_ang2rad(axa); % rads
d = mode(-diff(arad)); % Bin width (for bias correction when using binned data)

% Samples and sample size
sampnames = cytor_table.Properties.VariableNames; % Column names (as cell array)
N = size(cytor_table, 2); % Number of images to be analyzed
n_doub = table2array(cytor_table(182:183,:)); % Estimates of (cell number)/image
%Choose one:
	%n = n_doub(1,:); % Under estimate
	%n = n_doub(2,:); % Over-estimate
	n = mean(n_doub,1); % Average of the previous



% Do the analysis for all the samples
for k = 1:N
	% Extract the k-th histogram from the master table and normalize it
	cyt = table2array(cytor_table(1:180,k))';
	cyt = (cyt/sum(cyt))*n(k);

	% Descriptive Statistics
	fprintf('\nAnalyzing image %d/%d: %s', k, N, sampnames{k})

	mu = circ_mean(arad, cyt);
	mudeg = circ_rad2ang(mu)/2;
	fprintf('\n\tMean direction = %.2f deg', mudeg);

	%me = circ_median(arad, cyt); % doesn't work... CircStat bug?
	%medeg = circ_rad2ang(me)/2;
	%fprintf('\n\tMedian direction = %.4f deg', medeg);

	[mw, mo_ind] = max(cyt);
	modeg = a(mo_ind);
	fprintf('\n\tMode direction = %.2f deg', modeg);

	r = circ_r(arad, cyt, d);
		rsum = r * mw;
		%zm = rsum * exp(i*mu); % Used with polar() function, but now unused...
	s = circ_var(arad, cyt, d);
	fprintf('\n\tResultant vector length = %.4f', r);
	fprintf('\n\tCircular variance = %.4f', s);

	% Inferential statistics - Hypothesis testing
	[pval, z] = circ_rtest(arad, cyt, d);
	pv = circ_vtest(arad, mu, cyt, d);
	fprintf('\n\tRayleigh test p-value = %.4f (n = %.1f)', pval, n(k));
	fprintf('\n\tV-test (against mean) p-value = %.4f (n = %.1f)', pv, n(k));
	fprintf('\n')

	% Create the output as an array
	out_array(k,:) = [n(k), mudeg, modeg, r, s, pval, pv];

	% Plot histograms
	if false % Linear Plots
		figure, hold on
		plot(a, cyt)
		xlim([-90, 89])
		yy = ylim;
		ylim([0, yy(2)])
		xlabel('Direction (deg)')
		plot([mudeg,mudeg], [0,yy(2)], '-r')
		%plot([medeg,medeg], [0,yy(2)], '-m')
		plot([modeg,modeg], [0,yy(2)], '-b')
		title(sampnames{k}, 'Interpreter', 'none')
		hold off
	end

	if false % 360° Circular Plots
		figure
		polarplot([pi arad], [cyt(end) cyt], 'k')
		hold on
		polarplot([0 mu], [0 rsum], 'r', 'linewidth', 1.5)
		title(sampnames{k}, 'Interpreter', 'none')
		hold off
	end

	if true % 180° Circular Plots
		figure
		polarplot([pi arad]/2, [cyt(end) cyt], 'k')
		hold on
		thetalim([-90,+90])
		polarplot([0 mu]/2, [0 rsum], 'r', 'linewidth', 1.5)
		title(sampnames{k}, 'Interpreter', 'none')

		% Also print p-values, in red if significant
		if pval < 0.05 col = '#A2142F';, else col = 'black';, end
		text((0.98)*pi, mw, strcat('$p_{Ray}=', num2str(pval), '$'),...
			'Interpreter', 'latex', 'FontSize', 12, 'Color', col)
		if pv < 0.05 col = '#A2142F';, else col = 'black';, end
		text((1.02)*pi, mw, strcat('$p_{V}=', num2str(pv), '$'),...
			'Interpreter', 'latex', 'FontSize', 12, 'Color', col)
		hold off
	end
end

% Convert output array into a table and save a copy
out_stats = {'n', 'mean', 'mode', 'R', 'S', 'Rayleigh_pval', 'Vtest_pval'};
out_table = array2table(out_array,...
	'RowNames', sampnames,...
	'VariableNames', out_stats);
% NOTE
% writetable() determines the file format based on the specified extension.
% Use .xls, .xlsm, or .xlsx for Excel spreadsheet files.
% The variable names of the table will become column headings in the first line
% of the file, by default.
writetable(out_table, 'out_table.xlsx', 'WriteRowNames', true)
