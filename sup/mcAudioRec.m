%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Adapta��o dos exemplo do MATLAB
%
%
%	https://www.mathworks.com/help/daq/multichannel-audio-input-and-output-1.html#bt_9lng-2
%	https://www.mathworks.com/help/matlab/ref/audiowrite.html
%   https://www.mathworks.com/help/daq/ref/range.html
%
%
%
%

clc
close all
clear all %#ok<CLSCR>



% Adapta��o exemplo de grava��o do MATLAB
d = daq.getDevices;
s = daq.createSession('directsound');
addAudioInputChannel(s,'Audio1',1:2);

% Obten��o da faq
faq = s.Rate;

% Defini��o da dura��o da grava��o
 s.DurationInSeconds = input('Dura��o: ');
% s.DurationInSeconds = 7;

% Inicio da Grava��o
disp('Inicio da grava��o')
[data,T]= startForeground(s);
disp('Fim da grava��o')

% Gerando gr�fico da Grava��o
plot(T,data);

%Adapta��o do exemplo de reprodu��o do MATLAB
dev = d(3);
noutchan = 2;
addAudioOutputChannel(s, dev.ID, 1:noutchan);
queueOutputData(s,data);
disp('Inicio da Reprodu��o')
rep=startForeground(s);
disp('Fim da Reprodu��o')


%Adapta��o do exemplo de cria��o de arquivo wave
grav = input('N�mero da Grava��o: ','s');
filename = ['MATLABGrav#',grav,'.wav'];
audiowrite(filename,data,s.Rate);