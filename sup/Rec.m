% Programa apra teste da gravação de audio com o Octave
% Limpando a memória
clear all
close all
clc

% Definição do sistema de aquisição de dados
DAqD = audiodevinfo();
DAqNames = {DAqD.input.Name};
DAqID = {DAqD.input.ID};

% Escolhendo o TASCAM US-800 como input device
DAq = 8;

% Definindo os parâmetros de gravação
T = 5.0;      % Duração da gravação
faq = 5000;  % Frequência de aquisição
canais = 2;   % Número de canais
nbits = 16;   % Número de bits

%Definições de plotagem
dt=1/faq;
t = 0:dt:T-dt;

%Definindo parâmetros da gravação e reprodução
recorder = audiorecorder(faq,nbits,canais);
%player = a

%Gravando
record(recorder,T);
sinal = getaudiodata(recorder);

%