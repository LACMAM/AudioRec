%%%%
%
%   Esse arquivo irá gerar um atraso artifical entre canais de audio,
%   Os resultados desse atraso artificial serão utilizados para avaliar a
%   correlação cruzada de sinais de audio.
%
%
%
% https://www.mathworks.com/help/signal/ref/spectrogram.html?s_tid=doc_ta
%
%%%%

%Limpando a memória
clear; clc; close all;

%Leitura dos arquivos:
% [data0,faq0] = audioread('Marigold.mp3');
% [data1,faq1] = audioread('MATLABGrav#1.wav');
% [data2,faq2] = audioread('MATLABGrav#2.wav');
% [data3,faq3] = audioread('MATLABGrav#1.wav');
% [data4,faq4] = audioread('MATLABGrav#1.wav');
 [data1,faq1] = audioread('Pulso2.wav');
% [data1,faq1] = audioread('Chuva.flac');

%  data = decimate(


% Definição dos parâmetros da gravação
SO_npts = size(data1,1);
SO_chan = size(data1,2);
SO_dur = SO_npts/faq1;
dt = 1/faq1;

% Avaliando apenas um quinto do sinal (1 pulso)
data = data1(1:SO_npts/5,1);
Snpts = size(data,1);
Schan = size(data,2);
Sdur = Snpts/faq1;

% Criação de vetores de tempo para os gráficos
t = (0:Snpts-1)*dt;
% t_esp = (-Snpts+1:Snpts-1)*dt;

% Será feita a divisão do sinal em nd partes, e os atrasos a serem
% aplicados terão este tamanho.
ndiv = 3;                                   % Número de divisões que serão feitas;
Ldiv = floor(max(size(data))/ndiv);        % Número de pontos em cada divisão.
Tdiv = Sdur/ndiv;                           % Duração de cada divisão

% Gerando uma matriz com ndiv colunas e com o dobro de linhas do sinal original
% para alocação do sinal original e do sinal com os diversos atrasos
SM_npts = Snpts+Ldiv;               % Cálculo do número de pontos do sinal modificado
SM_data = zeros(SM_npts,ndiv);      % Geração da Matriz de sinais com atraso
pts=zeros(2,ndiv);                  % Vetor para alocação dos pontos iniciais
% SM_dur = SM_npts/faq1;
SM_d = Tdiv+Sdur;                  % Cálculo da duração do sinal modificação.

% Alocação do sinal na matriz com atraso
for div=1:ndiv
    pi=1+(div-1)*Ldiv;
    pf=pi+Snpts-1;
    pts(:,div)=[pi;pf];
    SM_data(pi:pf,div)=data(:,1);
end

% Gerando a correlação cruzada entre os sinais
[Sxy,lag] = xcorr(SM_data,'coeff');



t_esp = lag*dt;

% Definição dos parâmetros do gráfico e desenho.
mte = min(t_esp);
Mte = max(t_esp);
te_div = 15;
te_spc = (Mte-mte)/te_div;
te_tick = sort([mte:te_spc:Mte 0]);

% Sxy_tick = 
 %'XTickLabel',{}
grafSxy = figure();
exs = axes('Parent',grafSxy,'YGrid','on','XGrid','on','XMinorTick','on',...
           'XTick',te_tick);
xlabel({'Tempo [s]'});
ylabel({'Amplitude'});
title(['Correlações cruzadas de ',num2str(ndiv), ' sinais'])
hold on
plot(t_esp,Sxy);



%Definição do ta

% % % OUTRAS COISAS
% % SP_win = tukeywin(Ldiv,0.2);
% % spectrogram(data1(:,1),SP_win,0,Ldiv,faq1,'yaxis');