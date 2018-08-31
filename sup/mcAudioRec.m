%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Adaptação dos exemplo do MATLAB
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



% Adaptação exemplo de gravação do MATLAB
d = daq.getDevices;
s = daq.createSession('directsound');
addAudioInputChannel(s,'Audio1',1:2);

% Obtenção da faq
faq = s.Rate;

% Definição da duração da gravação
 s.DurationInSeconds = input('Duração: ');
% s.DurationInSeconds = 7;

% Inicio da Gravação
disp('Inicio da gravação')
[data,T]= startForeground(s);
disp('Fim da gravação')

% Gerando gráfico da Gravação
plot(T,data);

%Adaptação do exemplo de reprodução do MATLAB
dev = d(3);
noutchan = 2;
addAudioOutputChannel(s, dev.ID, 1:noutchan);
queueOutputData(s,data);
disp('Inicio da Reprodução')
rep=startForeground(s);
disp('Fim da Reprodução')


%Adaptação do exemplo de criação de arquivo wave
grav = input('Número da Gravação: ','s');
filename = ['MATLABGrav#',grav,'.wav'];
audiowrite(filename,data,s.Rate);