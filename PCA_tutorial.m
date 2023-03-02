%-----------------------------%
%                             %
%  MATLAB script to test PCA  %
%                             %
%-----------------------------%

% This is a tutorial on PCA in the form of a MATLAB script.
%
% In this tutorial the following convention for matrices has been chosen: rows
% of 'data' matrix correspond to observations (n) and columns correspond to
% variables (m). This is the same convention used by MATLAB `pca()` function and
% the same you can find in my handwritten notes (see "Formalismo Trasposto")
% alongside the paper: Jonathon Shlens, 'A Tutorial on Principal Component
% Analysis', 2009.





%-------------------------------------------------------------------------------
%---- 1 sample -----------------------------------------------------------------

clear all

% Number of Variables
m = 3;
% Number of Observations
n = 500;
% Color for plots
c = 'blue';

% Seeds the random number generator for reproducible results
rng(101)
% Generate values from a multivariate normal distribution with specified mean
% vector and covariance matrix
mu = [3 2 4]; % The mean vector
% Define any square matrix (but it could also have been a 3x2 or 3x1 rectangular
% matrix) and multiply it by its transposed to make it symmetric and positive
% definite in order to be a plausible "covariance matrix"
S = [3 -4 2; -4 5 2.5; 2 2.5 1];
Sigma = S*S';
R = chol(Sigma); % Cholesky factorization
data = repmat(mu,n,1) + randn(n,m)*R;

% Inspect numerical data
fprintf('\nGenerative Covariance Matrix\n');
Sigma
fprintf('\nActual Covariance Matrix\n');
CovData = cov(data);
CovData
VAR = var(data,0,1);
VAR
Means = mean(data,1);
Means

% Inspect graphs (Orthogonal projection style)
maximum = max(max(data));
minimum = min(min(data));
OrthProj(data,minimum,maximum,c)
figure
scatter3(data(:,1),data(:,2),data(:,3))



%-------------------------------------------------------------------------------
%	Manual PCA 
%-------------------------------------------------------------------------------

ManData = data - repmat(mean(data,1),n,1); % Mean centering
ManCov = (1/(n-1))*ManData'*ManData % Compute covariance matrix
[EigenVectors,EigenValues] = eig(ManCov) % Compute eigenvectors and eigenvalues

% NOTE:
% EigenVectors is a matrix whose columns are the eigenvectors of the covariance
% matrix
round(ManCov*EigenVectors,10) == round(EigenVectors*EigenValues,10)
% EigenVectors is the change-of-basis matrix that diagonalizes the covariance
% matrix
round((EigenVectors')*ManCov*EigenVectors,10) == round(EigenValues,10)
% Columns of EigenVectors are the components of the new orthonormal basis with
% respect to the old one...
round(EigenVectors(:,1)'*EigenVectors,10) == [1 0 0]
% ...and then they are unitary and orthogonal (orthogonal matrix)
round((EigenVectors')*EigenVectors,10) == eye(3)
% Finally, eigenvectors are the principal components (LOADINGS) of our data

% Sort the variances in decreasing order
[junk, index] = sort(-diag(EigenValues));
EigenVectors = EigenVectors(:,index);

% Project the original dataset (compute the SCORES)
TransData = ManData*EigenVectors;

% Inspect new graphs
maximum = max(max(TransData));
minimum = min(min(TransData));
OrthProj(TransData,minimum,maximum,c)
figure
scatter3(TransData(:,1),TransData(:,2),TransData(:,3))



%-------------------------------------------------------------------------------
%	Auto PCA
%-------------------------------------------------------------------------------

% By default, `pca()` centers the data and uses the singular value decomposition
% (SVD) algorithm
[coeff,score,latent,tsquared,explained,mu] = pca(data);

% NOTE:
% coeff are the principal components of the data and are equal (except for sign)
% to the ranked EigenVectors
coeff
EigenVectors
% latent are the principal component variances and are equal to the ranked
% EigenValues
round(latent,10) == round(-sort(-diag(EigenValues)),10)
% scores are the representations of data in the principal component space and
% are equal (except for sign) to TransData
score(1:5,:)
TransData(1:5,:)
% explained is the percentage of the total variance explained by each principal
% component...
explained
% ...and is equal to:
(-sort(-diag(EigenValues))/sum(diag(EigenValues)))*100

% Inspect new graphs
OrthProj(score,minimum,maximum,c)
figure
scatter3(score(:,1),score(:,2),score(:,3))

% Columns of coeff matrix (EigenVectors) are the components of the new
% orthonormal basis respect to the old one...
% ...columns of inv(coeff)==coeff' are the components of the old orthonormal
% basis respect to the new one (place Data Cursor over blue dots!)
coeff'

% Scores and Loadings in the same space
figure
biplot(coeff(:,1:3),'Scores',score(:,1:3),'varlabels',{'v_1','v_2','v_3'})



%-------------------------------------------------------------------------------
%	Manual PCA - standardized
%-------------------------------------------------------------------------------

ManDataS = data - repmat(mean(data,1),n,1); % Mean centering
ManDataS = ManDataS./repmat(std(ManDataS,0,1),n,1); % Standardize the variance (z-score)
%ManDataS2 = zscore(data); % Alternatively

% Inspect graphs
maximum = max(max(ManDataS));
minimum = min(min(ManDataS));
OrthProj(ManDataS,minimum,maximum,c)
%figure
%scatter3(ManDataS(:,1),ManDataS(:,2),ManDataS(:,3))

ManCovS = (1/(n-1))*ManDataS'*ManDataS; % Compute covariance (correlation) matrix

% NOTE:
% all the following commands return the same value (correlation matrix)
ManCovS
cov(ManDataS)
cov(zscore(data))
corr(data)
corr(ManDataS)

[EigenVectorsS,EigenValuesS] = eig(ManCovS) % Compute eigenvectors and eigenvalues
% Also in this case EigenVectorsS is an orthogonal change-of-basis matrix (whose columns are unitary and orthogonal)
round((EigenVectorsS')*EigenVectorsS,10) == eye(3)

% Sort the variances in decreasing order
[junk, index] = sort(-diag(EigenValuesS));
EigenVectorsS = EigenVectorsS(:,index);

% Project the original dataset (compute the SCORES)
TransDataS = ManDataS*EigenVectorsS;

% Inspect new graphs
maximum = max(max(TransDataS));
minimum = min(min(TransDataS));
OrthProj(TransDataS,minimum,maximum,c)	
%figure
%scatter3(TransDataS(:,1),TransDataS(:,2),TransDataS(:,3))



%-------------------------------------------------------------------------------
%	Auto PCA - standardized
%-------------------------------------------------------------------------------

% For matrix inputs, zscore() computes the z-scores of data using the mean and
% standard deviation along each column
[coeffS,scoreS,latentS,tsquaredS,explainedS,muS] = pca(zscore(data));

% NOTE:
% coeff are the principal components of the data and are equal (except for sign)
% to the ranked EigenVectorsS
coeffS
EigenVectorsS
% latent are the principal component variances and are equal to the ranked
% EigenValuesS
round(latentS,10) == round(-sort(-diag(EigenValuesS)),10)
% scores are the representations of data in the principal component space and
% are equal (except for sign) to TransDataS
scoreS(1:5,:)
TransDataS(1:5,:)
% explained is the percentage of the total variance explained by each principal
% component...
explainedS
% ...and is equal to:
(-sort(-diag(EigenValuesS))/sum(diag(EigenValuesS)))*100

% Inspect new graphs
OrthProj(scoreS,minimum,maximum,c)
figure
scatter3(scoreS(:,1),scoreS(:,2),scoreS(:,3))

% Columns of coeffS matrix (EigenVectorsS) are the components of the new
% orthonormal basis respect to the old one...
% ...columns of inv(coeffS)==coeffS' are the components of the old orthonormal
% basis respect to the new one (place Data Cursor over blue dots!)
coeffS'

figure
biplot(coeffS(:,1:3),'Scores',scoreS(:,1:3),'varlabels',{'v_1','v_2','v_3'})





%-------------------------------------------------------------------------------
%---- 2 samples ----------------------------------------------------------------

clear all

% Number of Variables
m = 3;
% Number of Observations
n1 = 500;
n2 = 300;

% Seeds the random number generator for reproducible results
rng(101)
% Generate values from a multivariate normal distribution with specified mean
% vector and covariance matrix
mu1 = [3 2 4];
S1 = [3 -4 2; -4 5 2.5; 2 2.5 1]; % ...just any square matrix
Sigma1 = S1*S1'; % ...make it positive definite to be a "covariance matrix"
R1 = chol(Sigma1); % Cholesky factorization
data1 = repmat(mu1,n1,1) + randn(n1,m)*R1;

% Same covariance matrix
	mu2 = [20 -11 -30];
	Sigma2 = Sigma1;
% OR not
	mu2 = [20 -11 3];
	S2 = [5 6 -5; 7 -3.5 -2; -2.5 6 1]; % ...just another square matrix
	Sigma2 = S2*S2'; % ...make it positive definite to be a "covariance matrix"

R2 = chol(Sigma2); % Cholesky factorization
data2 = repmat(mu2,n2,1) + randn(n2,m)*R2;

% Inspect numerical data
fprintf('\nGenerative Covariance Matrix\n');
Sigma1
Sigma2
fprintf('\nActual Covariance Matrix\n');
CovData1 = cov(data1);
CovData1
CovData2 = cov(data2);
CovData2
Means1 = mean(data1,1);
Means1
Means2 = mean(data2,1);
Means2

% Joined data
data = [data1;data2];
 % Color to distinguish the two groups
c = [repmat([0,0,1],n1,1);repmat([1,0,0],n2,1)];

% Inspect graphs (Orthogonal projection style)
maximum = max(max(data));
minimum = min(min(data));
OrthProj(data,minimum,maximum,c)
figure
scatter3(data(:,1),data(:,2),data(:,3),5,c)



%-------------------------------------------------------------------------------
%	Auto PCA
%-------------------------------------------------------------------------------

% By default, pca() centers the data and uses the singular value decomposition
% (SVD) algorithm
[coeff,score,latent,tsquared,explained,mu] = pca(data);

coeff
latent
score(1:5,:)
explained

% Inspect new graphs
maximum = max(max(score));
minimum = min(min(score));
OrthProj(score,minimum,maximum,c)

% Scores and Loadings in the same space
figure
biplot(coeff(:,1:3),'Scores',score(:,1:3),'varlabels',{'v_1','v_2','v_3'})



%-------------------------------------------------------------------------------
%	Auto PCA - standardized
%-------------------------------------------------------------------------------

% For matrix inputs, zscore() computes the z-scores of data using the mean and
% standard deviation along each column
dataS = zscore(data);
% Inspect graphs
maximum = max(max(dataS));
minimum = min(min(dataS));
OrthProj(dataS,minimum,maximum,c)
figure
scatter3(dataS(:,1),dataS(:,2),dataS(:,3),5,c)

[coeffS,scoreS,latentS,tsquaredS,explainedS,muS] = pca(dataS);

coeffS
latentS
scoreS(1:5,:)
explainedS

% Inspect new graphs
maximum = max(max(scoreS));
minimum = min(min(scoreS));
OrthProj(scoreS,minimum,maximum,c)

figure
biplot(coeffS(:,1:3),'Scores',scoreS(:,1:3),'varlabels',{'v_1','v_2','v_3'})
