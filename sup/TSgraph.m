function TSgraph(t1, S1)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 02-May-2018 09:17:38

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','XMinorTick','on',...
    'XTickLabel',{'0','0,2','0,4','0,6','0,8','1,0','1,2','1,4','1,6','1,8','2,0','2,2','2,4','2,6','2,8','3,0','3,2','3,4','3,6','3,8','4,0','4,2','4,4','4,6','4,8','5,0'},...
    'XTick',[0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5]);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 5]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-1 1]);
box(axes1,'on');
hold(axes1,'on');

% Create plot
plot(t1,S1,'DisplayName','Pulso');

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

