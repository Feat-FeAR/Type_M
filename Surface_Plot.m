
data = [
28.72	1.64
22.60	8.51
43.74	2.92
26.49	4.73
60.13	5.46
70.58	4.68
];


ctrl = data([1 3 5],:);
OHP = data([2 4 6],:);

% OR transposed

data = data';
ctrl = data(:,[1 3 5]);
OHP = data(:,[2 4 6]);


figure

sc = surf(ctrl);
hold on
so = surf(OHP);

sc.EdgeColor = 'b';
so.EdgeColor = 'r';

sc.FaceColor = 'b';
so.FaceColor = 'r';

sc.FaceAlpha = 0.07;
so.FaceAlpha = 0.07;

xlabel('Time from Caps (min)')
ylabel('Time from OHP (h)')
zlabel('Firing Rate (60 \cdot Hz / ch)')

set(gca,'Xtick',[1 2])
xlim([0.9 2.1])

set(gca,'Ytick',[1 2 3])
ylim([0.9 3.1])

set(gca,'XTickLabel',[1 2]);
set(gca,'YTickLabel',[6 24 48]);

set(sc,'LineWidth',1)
set(so,'LineWidth',1)

