function timesignalgraph(X1, Y1)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 02-May-2018 09:02:10

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','XMinorTick','on',...
    'XTickLabel',{'0,0','0,1','0,2','0,3','0,4','0,5','0,6','0,7','0,8','0,9','1,0','1,1','1,2','1,3','1,4','1,5','1,6','1,7','1,8','1,9','2,0','2,1','2,2','2,3','2,4','2,5','2,6','2,7','2,8','2,9','3,0','3,1','3,2','3,3','3,4','3,5','3,6','3,7','3,8','3,9','4,0','4,1','4,2','4,3','4,4','4,5','4,6','4,7','4,8','4,9','5,0'},...
    'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5]);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 5]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-1 1]);
box(axes1,'on');
hold(axes1,'on');

% Create plot
plot(X1,Y1,'DisplayName','Pulso');

% Create xlabel
xlabel({'Tempo [s]'});

% Create ylabel
ylabel({'Amplitude [V]'});

% Create title
title({'Sinal no tempo'});

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.817038691589781 0.839535238072751 0.0543749994225798 0.0300046875618792]);

