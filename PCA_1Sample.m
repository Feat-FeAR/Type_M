%----------------------------------%
%----------------------------------%
% PCA for Dimensionality Reduction %
%----------------------------------%
%----------------------------------%

dataS = zscore(data);
[n, k] = size(dataS);

% Choose a color
% Old MATLAB style
c = [0,0,1];
c = [1,0,0];
% New MATLAB style
c = [0 0.4470 0.7410];
c = [0.8500 0.3250 0.0980];
% My style
c = [0 0.5 0.5];
c = [180/255 28/255 173/255];

labels = {'Sa', 'Sku', 'Smean', 'Sp', 'Sq', 'Ssk', 'Sv', 'Sz', 'Sdq', 'Sdr', 'Sal', 'Str'}

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
format = {{'Marker','.','MarkerSize',10,'MarkerEdgeColor',c}};
biplotG(coeffS(:,1:2),scoreS(:,1:2),'VarLabels',labels,'Format',format)


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
legend('Sa','Sku','Smean')



%-------------------------%
% Descriptor Spectra      %
%-------------------------%

% Import da PCA2SampleComparison.m, da riadattare...  %%%%%%%%%%%%
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

