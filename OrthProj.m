% Orthographic projection

function OrthProj(x,minimum,maximum,c)

dotsize = 7;

if(size(x,2) ~= 3)
	fprintf('\nWARNING: Input data are not 3D or the matrix is not properly organized\n\n');
	return
end

figure

subplot(2,2,3)
	scatter(x(:,1),x(:,2),dotsize,c,'filled')
	xlim([minimum,maximum])
	ylim([minimum,maximum])
	set(gca,'Xdir','reverse','Ydir','reverse','XAxisLocation','top','YAxisLocation','right')
	ylabel('Component 2')
	%ylabel('Component 2 (31.8%)')
subplot(2,2,1)
	scatter(x(:,1),x(:,3),dotsize,c,'filled')
	xlim([minimum,maximum])
	ylim([minimum,maximum])
	set(gca,'Xdir','reverse','YAxisLocation','right')
	xlabel('Component 1')
	ylabel('Component 3')
	title('2D Score Plot of the First 3 PCs - Orthographic Projection')
	%xlabel('Component 1 (34.1%)')
	%ylabel('Component 3 (15.5%)')
subplot(2,2,2)
	scatter(x(:,2),x(:,3),dotsize,c,'filled')
	xlim([minimum,maximum])
	ylim([minimum,maximum])
	xlabel('Component 2')
	%xlabel('Component 2 (31.8%)')
	