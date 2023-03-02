function MeanTrace = CalciumArea(filename,unit,mark,output)

%
%--------------------------------------------------------------------------------
% Script for the evaluation of the area under calcium signals 
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% CalciumArea(filename,unit,mark,output)
%
% INPUT       TYPE       EXAMPLE           MEANING
% -----       ----       -------           -------
% filename -> string  -> 'filename.csv' -> File Name
% unit     -> string  -> 's'            -> Time Unit: 's' or 'ms'
% mark     -> array   -> [300,310,510]  -> Marks* - Step Number
% output   -> boolean -> true           -> Print .eps Output Graphs
%
% OUTPUT      TYPE                         MEANING
% ------      ----                         -------
% -none-   -> plot                      -> 1 Plot - Sampling Quality Control
% -none-   -> plot                      -> 1 Plot - Trace Visualization
%
% *NOTE on marks: If the array features just one element, it is intended to be
%                 the time step of Agonist Administration. In this case the onset
%                 time-step will be automatically detected and the width of the
%                 Windows of Interest for the evaluation of the area is fixed
%                 by WindowOfInterest variable.
%                 If the array has two element the second one will overwrite
%                 the onset variable, thus forcing the starting point of the
%                 window of interest.
%                 If the array has three element also the width of the Window
%                 of Interest is user-defined.
%


% IMPORTANT
%
% This script was originally made to for the evaluation of the area under
% time-course traces of intracellular calcium concentration, whose analysis
% has been presented in 
% Distasi, Dionisi, et al. "The interaction of SiO2 nanoparticles with the
% neuronal cell membrane: activation of ionic channels and calcium influx",
% Nanomedicine (Lond.), 2019.
% See folder ...\UniTo Drive\WORKS\2019 - Article - SiO2 NPs III\Custom Scripts
% However, a somehow more advanced (more general and interactive) script with
% similar purposes was also implemented for the analysis of calcium traces from
% SiPM in the work
% Ruffinatti, Lomazzi et al. "Assessment of a Silicon-Photomultiplier-Based
% Platform for the Measurement of Intracellular Calcium Dynamics with Targeted
% Aequorin". ACS, 2020
% See folder: ...\WORKS\2020 - Article - SiPM-AEQ\Data\16 - Exp_8 - Dati
% See, in particular, `SiPM.m`, `SiPM2.m`, `SiPM3.m`, and `SiPM4.m` m-files.
%
% The day this script is to be used again, some integration between all these
% pieces of code is highly desirable!!


% Graphic parameters
s1 = 16; % X-Y TickLabel size
s2 = 19; % X-Y Label and text size
s3 = 24; % Title size

% Output print control - Default value = false = Print nothing
if (nargin < 4)
	output = false;
end

% Marker array control
if (nargin < 3 || length(mark) == 0)
	fprintf('\n\nWARNING: Array of marks cannot be empty!!\n');
	fprintf('\n\n');
	return
end

% Time unit control
if ~(strcmp(unit,'s') | strcmp(unit,'ms'))
	fprintf('\n\nWARNING: Invalid Time Unit\n');
	fprintf('\n\n');
	return
end

% Open data file
data = dlmread(filename);

% Extract time vector
timeVec = data(:,1);
timeVec = timeVec - timeVec(1); % Start from t=0s

% Sampling time mode
DtimeVec = diff(timeVec);
dt = mode(DtimeVec);
if (strcmp(unit,'ms'))
	fprintf('\n\nSampling Time (Mode) = %.2f ms',dt); % "%.2f" stands for "floating point number with only 2 decimal places"
else
	fprintf('\n\nSampling Time (Mode) = %.2f s',dt);
end

% Accepted time error for each sample
timeError = dt/2;

% Check if time steps are equal (evenly spaced)
if ~(isempty(find(abs(DtimeVec - dt) > timeError)))
	fprintf('\n\nWARNING: Not equal time steps!');
end

% From now on timeVec is measured in s
if (strcmp(unit,'ms'))
	timeVec = timeVec/1000;
end

% Print other general information
if (strcmp(unit,'ms'))
	fprintf('\n\nNyquist Frequency = %.2f Hz',1000/(2*dt));
else
	fprintf('\n\nNyquist Frequency = %.2f mHz',1000/(2*dt));
end
fprintf('\n\nTime Vector Length = %d',length(timeVec));
fprintf('\n\nLowest Frequency = %.2f mHz',1000/timeVec(end));
fprintf('\n\nOctaves Number = %d',floor(log2(length(timeVec))));
fprintf('\n');

% Data matrix
data(:,1) = [];

% Check user input
if (size(data,2) < 2)
	fprintf('\n\nWARNING: The number of traces is smaller than 2 or the column of time stamps is missing!\n');
	fprintf('\n\n');
	return; %Function breaks
end

% Store y-axis limits
massimo = max(max(data));
minimo = min(min(data));

% Trace selection: Dropped-to-0 traces
LinearIndex = find(data==0);
[i,j] = ind2sub(size(data),LinearIndex);

% Trace selection: Basal outliers evaluated over the first w time points 
w = 10;
BasalMean = mean(mean(data(1:w,:),1));
BMsigma = std(mean(data(1:w,:),1));
h = find(abs(mean(data(1:w,:),1) - BasalMean) > 3*BMsigma);

% Comprehensive index array of traces to be discarded
SelectOut = unique([j',h]);

%--------------------------------------------------------------------------------
% Mean Plot Computation
% Plot a single trace representing the mean value +/- SEM at each point in time
% after outliers and dropped-to-0 traces have been discarded
%--------------------------------------------------------------------------------

% Select Traces
SelectedData = data;
SelectedData(:,[SelectOut]) = [];

% Mean trace
y = mean(SelectedData,2);
MeanTrace = y;
sem = std(SelectedData,0,2)/sqrt(size(SelectedData,2));

% Linearly fit the basal signal
if (mark(1) > length(timeVec))
	fprintf('\n\nWARNING: Agonist Administration mark (1) exceeds time domain');
	fprintf('\n\n');
	return
end
coeff = polyfit(timeVec(1:mark(1)),y(1:mark(1)),1);
LinFit = polyval(coeff,timeVec);

% Find the onset of the response
onset = min(find(y(mark(1):end)-sem(mark(1):end) > LinFit(mark(1):end))) + mark(1);
if (isempty(onset) && length(mark)==1)
	fprintf('\n\nWARNING: The onset of the response is not detectable! ...it needs to be user-defined');
	fprintf('\n\n');
	return
end

% Array of marks
WindowOfInterest = 200; % Steps

switch (length(mark))
	
	case 1
		
		mark = [mark(1),onset,onset+WindowOfInterest];
		
	case 2
		
		mark = [mark(1),mark(2),mark(2)+WindowOfInterest];
		
	case 3
		
		% Do nothing: All marker are user-defined
		
	otherwise
		
		fprintf('\n\nWARNING: Too many time-marks in mark array!');
		fprintf('\n\n');
		return
		
end

if (mark(2) > length(timeVec))
	fprintf('\n\nWARNING: Response Onset mark (2) exceeds time domain');
	fprintf('\n\n');
	return
elseif (mark(3) > length(timeVec))
	fprintf('\n\nWARNING: Window of Interest End-mark (3) exceeds time domain');
	fprintf('\n\n');
	return
end
	
fprintf('\n\nOnset Step = %d',mark(2));
fprintf('\n\nWindow of Interest End-Step = %d',mark(3));
if (strcmp(unit,'ms'))
	fprintf('\n\nWindow of Interest Width = %.2f ms',(mark(3)-mark(2))*dt);
else
	fprintf('\n\nWindow of Interest Width = %.2f s',(mark(3)-mark(2))*dt);
end
fprintf('\n');

% Compute the area under the mean curve
AreaUnderMean = trapz(y(mark(2):mark(3))-LinFit(mark(2):mark(3)))*dt;
fprintf('\n\nArea under the Mean Curve = %.4f',AreaUnderMean);

% Plot Results
figure

subplot(2,1,2)
	hold on
	plot(timeVec,y,'-r','LineWidth',2);
	plot(timeVec,y+sem,'-b','LineWidth',1);
	plot(timeVec,y-sem,'-b','LineWidth',1);
	plot(timeVec,LinFit,'-c','LineWidth',1);
	
	xlim([timeVec(1),timeVec(end)])
	ylim([minimo,massimo])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	xlabel('Time (s)','FontSize',s2)
	ylabel('Ratio','FontSize',s2)
	title(['Selected Traces: ', num2str(size(SelectedData,2)), ' ROIs'],'FontSize',s3)
	%text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['legend'],'FontSize',s2,'Color','k')
	
	for l = 1:length(mark)
		plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	end
	
%--------------------------------------------------------------------------------
% Superimposition Plot Computation
% Plot all the original traces together with a Blue->Red color gradient
% Discarded traces are marked in green
%--------------------------------------------------------------------------------

subplot(2,1,1)
	
	counter = 1;
	for k = 1:size(data,2)
		
		y = data(:,k);
		
		hold on
		p = plot(timeVec,y,'-b','LineWidth',1);
		
		if (isempty(find(SelectOut == k)))
			
			% Blue->Red color gradient
			set(p,'Color',[(k-1)/(size(data,2)-1),0,(size(data,2)-k)/(size(data,2)-1)])
			
			% Linearly fit each basal signal
			coeff = polyfit(timeVec(1:mark(1)),y(1:mark(1)),1);
			LinFit = polyval(coeff,timeVec);
			
			% Compute the area under each curve
			AreaUnderCurve(counter) = trapz(y(mark(2):mark(3))-LinFit(mark(2):mark(3)))*dt;
			
			counter = counter+1;
			
		else
			% Green
			set(p,'Color',[0,1,0])
		end
	
	end
	
	xlim([timeVec(1),timeVec(end)])
	ylim([minimo,massimo])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	ylabel('Ratio','FontSize',s2)
	title(['All Traces: ', num2str(size(data,2)), ' ROIs'],'FontSize',s3)
	%text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['legend'],'FontSize',s2,'Color','k')
	
	for l = 1:length(mark)
		plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	end
	
	% Mean of the area under the curves
	%AreaUnderCurve' % Decomment this line to have the complete list of areas for statistical purposes
	MeanArea = mean(AreaUnderCurve);
	AreaSem = std(AreaUnderCurve)/sqrt(length(AreaUnderCurve));
	fprintf('\n\nMean Area under the Curves = %.4f +/- %.4f',MeanArea,AreaSem);
	
% Print output graph
if (output)
	print -depsc traceout.eps
end

%--------------------------------------------------------------------------------
% Sampling Quality Control
%--------------------------------------------------------------------------------

figure
plot(DtimeVec,'-b','LineWidth',1), hold on

xlim([1,length(timeVec)-1])
set(gca,'FontSize',s1,'XTick',unique([1,get(gca,'XTick'),length(timeVec)-1]));
xlabel('Time Step','FontSize',s2)
if (strcmp(unit,'ms'))
	ylabel('dt (ms)','FontSize',s2)
else
	ylabel('dt (s)','FontSize',s2)
end
title('Sampling Quality Control','FontSize',s3)
if (strcmp(unit,'ms'))
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['Sampling Time (Mode) = ',num2str(dt),' ms'],'FontSize',s2,'Color','k')
else
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['Sampling Time (Mode) = ',num2str(dt),' s'],'FontSize',s2,'Color','k')
end

for k = 1:length(mark) % If length(mark)==0 this does nothing
	plot([mark(k),mark(k)],[min(ylim),max(ylim)],'-r','LineWidth',1)
end


fprintf('\n\n');


%%------------------------------------------------------------------------------------------------------%%
%%---------------------------------------------END------------------------------------------------------%%
%%------------------------------------------------------------------------------------------------------%%