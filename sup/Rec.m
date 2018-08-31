% Programa apra teste da grava��o de audio com o Octave
% Limpando a mem�ria
clear all
close all
clc

% Defini��o do sistema de aquisi��o de dados
DAqD = audiodevinfo();
DAqNames = {DAqD.input.Name};
DAqID = {DAqD.input.ID};

% Escolhendo o TASCAM US-800 como input device
DAq = 8;

% Definindo os par�metros de grava��o
T = 5.0;      % Dura��o da grava��o
faq = 5000;  % Frequ�ncia de aquisi��o
canais = 2;   % N�mero de canais
nbits = 16;   % N�mero de bits

%Defini��es de plotagem
dt=1/faq;
t = 0:dt:T-dt;

%Definindo par�metros da grava��o e reprodu��o
recorder = audiorecorder(faq,nbits,canais);
%player = a

%Gravando
record(recorder,T);
sinal = getaudiodata(recorder);

%