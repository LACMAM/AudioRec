%Tabula Rasa
clear;clc;close all;

% Abrindo arquivos
[data0,faq0] = audioread('G1C1.wav');
[data1,faq1] = audioread('G1C2.wav');
[data2,faq2] = audioread('G1C3.wav');

data = [data0 data1 data2];

% Obtendo os tamanhos limites dos vetores e parâmetros de gravação:
Ssiz = size(data,1);
Schn = size(data,2);
dt = 1/faq0; 

% Gerando os vetores temporais:
t = (0:Ssiz-1)/faq0;
t_esp = (-Ssiz+1:Ssiz-1)/faq0;

% Cálculo das correlações cruzadas
[S01] = xcorr(data1(),data0());
[S02] = xcorr(data2(),data0());

%Cálculo dos vetores de autocorrelação:
[S00,lag00] = xcorr(data0);
[S11,lag11] = xcorr(data1);
[S22,lag22] = xcorr(data2));

% Normalizando o vetor de correlação;
S00 = S00/faq0;
S11 = S11/faq1;
S22 = S22/faq2;