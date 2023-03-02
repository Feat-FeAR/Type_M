%---------------------------%
%---------------------------%
% PCA on Morpheus output    %
%---------------------------%
%---------------------------%

data = [MasterMatrix1;MasterMatrix2];
dataS = zscore(data);

a = size(MasterMatrix1,1);
b = size(MasterMatrix2,1);
k = size(dataS,2);

% Color to distinguish the groups
% Old MATLAB style
c1 = [0,0,1];
c2 = [1,0,0];
% New MATLAB style
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
% My style
c1 = [0 0.5 0.5];
c2 = [180/255 28/255 173/255];

c = [repmat(c1,a,1);repmat(c2,b,1)];

labels = {'Area','Perim','Major','Minor','Angle','Circ','Feret','FeretAngle','MinFeret','AR','Round','Solidity'};

fprintf('\nCovariance Matrix\n');
CovData = cov(data);
CovData
fprintf('\nCorrelation Matrix\n');
CorrData = corr(data);
CorrData

%-------------------------%
% Auto PCA - standardized %
%-------------------------%

% For matrix inputs, zscore() computes the z-scores of data using the mean and standard deviation along each column
[coeffS,scoreS,latentS,tsquaredS,explainedS,muS] = pca(dataS);

coeffS
latentS
scoreS(1:5,:)
explainedS
sum(explainedS(1:3))

% Inspect graphs
% Scree Plot
figure
bar(explainedS,'FaceColor',[0 0.4470 0.7410]);
xlim([0,k+1])
title('Scree Plot')
xlabel('Principal Component (PC)')
ylabel('Explained Variance (%)')

% Score Plot
maximum = max(max(scoreS));
minimum = min(min(scoreS));
OrthProj(scoreS(:,1:3),minimum,maximum,c)

figure
scatter3(scoreS(:,1),scoreS(:,2),scoreS(:,3),10,c,'filled')
title('3D - Score Plot of the First 3 PCs')
xlabel(['PC1 (',num2str(explainedS(1),3),'%)'])
ylabel(['PC2 (',num2str(explainedS(2),3),'%)'])
zlabel(['PC3 (',num2str(explainedS(3),3),'%)'])

% Biplot - Loading Plot
figure
biplot(coeffS(:,1:3),'Scores',scoreS(:,1:3),'varlabels',labels)

figure
format = { {'Marker','.','MarkerSize',10,'MarkerEdgeColor',c1}; {'Marker','.','MarkerSize',10,'MarkerEdgeColor',c2} };
biplotG(coeffS(:,1:2),scoreS(:,1:2),'Groups',c(:,1),'VarLabels',labels,'Format',format)


% Loading Analysis (la formula di combinazione delle variabili originali (descrittori) per la costruzione delle PC)
figure
plot([1:k],coeffS(:,1:3),'.-','MarkerSize', 20);
xlim([1,k])
set(gca,'xtick',[1:12],'xticklabel',labels,'XTickLabelRotation',45)
title('Loading Plot')
ylabel('Loading')
legend('PC1','PC2','PC3')

% Loading Analysis (la formula di combinazione delle PC per la ricostruzione delle variabili originali (descrittori))
figure
plot([1:k],coeffS(1:3,:),'.-','MarkerSize', 20); % Perché la trasposta di una matrice ortogonale è anche la sua inversa!
xlim([1,k])
set(gca,'xticklabel',{'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12'})
legend('Area','Perim','Major')

%-------------------------%
% Descriptor Spectra      %
%-------------------------%

% Original Standardized Spectra
figure
% Standardized Mean Spectra
subplot(2,1,1)
	%plot([1:k],dataS(1:10,:));
	m1 = mean(dataS(1:a,:),1);
	m2 = mean(dataS(1+a:a+b,:));
	sem1 = std(dataS(1:a,:))/sqrt(a);
	sem2 = std(dataS(1+a:a+b,:))/sqrt(b);
	hold on
	errorbar([1:k],m1,sem1,'Color',c1,'LineStyle','-','LineWidth',1,'Marker','.','MarkerSize',20);
	errorbar([1:k],m2,sem2,'Color',c2,'LineStyle','-','LineWidth',1,'Marker','.','MarkerSize',20);
	xlim([0.5,k+0.5])
	set(gca,'xtick',[1:12],'xticklabel',labels,'XTickLabelRotation',45)
	title('Standardized Mean Spectra')
	ylabel('Mean z-scores')
	legend('Group 1','Group 2') %%%
	hold off
	
% Cohen's d and marginal p-values
% Using data or dataS (standardized) for Choen's d computation and t-tests is the same
subplot(2,1,2)
	pooledSD = sqrt(((a-1)*(std(dataS(1:a,:)).^2)+(b-1)*(std(dataS(1+a:a+b,:)).^2))/(a+b-2));
	hold on
	cohen_d = abs(m2-m1)./pooledSD;
	plot([1:k],cohen_d,'Color','k','LineStyle','-','LineWidth',1,'Marker','.','MarkerSize',20);
	[h,p] = ttest2(data(1:a,:),data(1+a:a+b,:),'Vartype','unequal');
	text([1:k]+0.1,cohen_d-sign([diff(cohen_d),-1])*0.025,strsplit(num2str(p,3)));
	xlim([0.5,k+1])
	set(gca,'xtick',[1:12],'xticklabel',labels,'XTickLabelRotation',45)
	title('Mean Spectrum Distance')
	ylabel('Cohen''s {\itd} and {\itp}-values')
	hold off
	
% First 3 PCs Decomposition
figure
dataS3PCs = scoreS(:,1:3)*(coeffS(:,1:3))';
%plot([1:k],dataS3PCs(1:10,:));
plot([1:k],[mean(dataS3PCs(1:a,:),1);mean(dataS3PCs(1+a:a+b,:),1)],'.-','MarkerSize', 20);

% Residual - Last 4 PCs
figure
residualS = scoreS(:,9:12)*(coeffS(:,9:12))';
%plot([1:size(dataS,2)],residualS(1:10,:));
plot([1:k],[mean(residualS(1:a,:),1);mean(residualS(1+a:a+b,:),1)],'.-','MarkerSize', 20);

%-------------------------%
% Class Separation        %
%-------------------------%

% 2-class Linear separation by SVM
Mdl = fitclinear(scoreS(:,1:2), [zeros(a,1);ones(b,1)]);
figure, hold on
scatter(scoreS(:,1),scoreS(:,2),dotsize,c,'filled')
title('Linear Separation')
xlabel('Principal Component 1')
ylabel('Principal Component 2')
oldxlim = xlim;
oldylim = ylim;
plot(xlim,Mdl.Beta(1)*xlim+Mdl.Beta(2))
xlim(oldxlim)
ylim(oldylim)
hold off