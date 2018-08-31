
% An�lise dos dados experimentais de "resposta impulsiva" das tubula��es

% Apagando a memoria do MATLAB
close all
clear all
clc
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       Lendo o Arquivo de medi��o      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definindo o caminho da medi��o e arquivo de medi��o
Pasta = 'C:\GSVL\0 - Mestrado\An�lise das Medi��es\TXTs'; % PC NOVO

% Defini��o do Nome da Medi��o

% Arq = 'Curva-90-Espuma';

% Arq = '(2014-08-14)MST-Aberta';
% Arq = '(2014-08-14)MST-Fechada';
% Arq = '(2014-08-14)MST-Anecoica';

% Arq = '(2014-06-26)MST-VGOFF';
% Arq = '(2014-06-26)MST-VGP3';
% Arq = '(2014-06-26)MST-VGP2';
% Arq = '(2014-06-26)MST-VGP1';
Arq = '(2014-06-26)MST-VGON';


% Arq = '(2014-08-01)MST-Silenciador';

% Arq = '(LTS) S2P25L500'; 
% Arq = '(LTS) S2P10L100'; 
% Arq = '(LTS) S8P10L500';
% Arq = '(LTS)S6P75L100';

% Arq = '(LTS)S4P10L100';
% Arq = '(LTS)S4P25L100';
% Arq = '(LTS)S4P50L100';
% Arq = '(LTS)S4P75L100';

% Arq = '(LTS)S4P10L250';
% Arq = '(LTS)S4P25L250';
% Arq = '(LTS)S4P50L250';
% Arq = '(LTS)S4P75L250';

% Arq = '(LTS)S4P10L500';
% Arq = '(LTS)S4P25L500';
% Arq = '(LTS)S4P50L500';
% Arq = '(LTS)S4P75L500';

% Fim da defini��o caminho e do arquivo de medi��o


% Importando os dados da medi��o
importfile([Pasta,'\',Arq,'.txt']);

% % Deletando poss�veis erros
% eNaN = isnan(data);
% data(eNaN) = 0;

% Obtendo os Par�metros de medi��o

    % Tamanho do vetor de medi��o
    [qpts,ncol] = size(data);
    nchan = ncol-1;
 
    %Obten��o do vetor temporal [s]
    t = data(1:qpts-1,1); 
    dt = mean(diff(t));       

    % Obten��o dos dados necess�rios para a FFT
    faq = 1/dt;                 % Frequ�ncia de Aquisi��o de dados

    % Vetor Espacial
    c = 343;            % Velocidade do Som no Fluido
    x = t*c/2;

% Definindo algumas propriedades dependentes do sistema
    
if (qpts < 6000)
    Lmax = 8;           % Definindo o comprimento  m�ximo de interesse 
    d = 50/1000;        % Di�metro Interno do Tubo
    s = 53/1000;        % Espa�amento entre microfones

else 
    Lmax = 120;
    d = 98.7/1000;        % Di�metro Interno do Tubo
    s = 40/1000;        % Espa�amento entre microfones
    x = x-min(x);
    [maxi,trig] = max(data(:,2:ncol));  %Para descontar o pr�-trigger
%     data(:,2:5) = data(:,2:ncol)/max(maxi);

    
end

dx = mean(diff(x));
nLmax = floor(Lmax/dx);

% Salvando as Medi��es
PMs = data(1:1:qpts-1,2:ncol);

% Limpando a mem�ria
clear data;

% Calculo das Frequ�ncias de Corte do Aparato Experimental
    fmp1 = (0.58*c)/d;          % Frequ�ncia de Corte do Tubo
    fmp2 = (0.45*c)/s;          % Frequ�ncia de Corte do Espa�amento dos microfones

% Definindo as frequ�ncias de interesse:
    fmax = min(fmp1,fmp2);
    fmin = 100;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                           %
            %       Processamento dos Dados medidos     %
            %                                           %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculo do N�vel de Press�o Sonora nos microfones
aux = PMs./(2*10^(-5));
NPS = 10*log10(aux.^2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Desenhando o Gr�fico da Medi��o no T  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 1 - Cria��o do Eixo X
xmin = round(min(x)*10)/10;
xmax = round(x(nLmax)*10)/10;
deltax = round((xmax-xmin)*10)/200;
vx = xmin:deltax:xmax;

tmin = (min(t)*10)/10;
tmax = (t(nLmax)*10)/10;
deltat = round((tmax-tmin)*10)/200;
tx = tmin:deltat:tmax;


% Desenhando o Gr�fico do Sinal no Tempo
figure('Position',[100 310 960 540])

if (qpts < 6000)
    
    %%%%%%%%%%%%%%
    %  FIGURA 1  %
    %%%%%%%%%%%%%%
    
    subplot(2,1,1)
%     plot(1000*t(1:nLmax,1),PMs(1:nLmax,1:3),'LineWidth',2);
%     set (gca,'XMinorGrid','on')
%     xlabel('Tempo [ms]','FontSize',14,'FontName','Verdana');

    plot(x(1:nLmax,1),PMs(1:nLmax,1:3),'LineWidth',1);
    set (gca, 'XTick',vx,'XMinorGrid','on')   
    xlabel('Dist�ncia Percorrida pelo Pulso [m]','FontSize',12,'FontName','Verdana');

    ylabel('Amplitude [Pa]','LineWidth',2,'FontSize',12,'FontName','Verdana');
    legend('Gerador de Sinais','Microfone 1','Microfone 2')
    axis tight
    grid on
    ylim([-1 1.5])

    %%%%%%%%%%%%%%
    %  FIGURA 2  %
    %%%%%%%%%%%%%%

    subplot(2,1,2)
    
    plot(x(1:nLmax,1),PMs(1:nLmax,4:5),'LineWidth',1);
    xlabel('Dist�ncia Percorrida pelo Pulso [m]','FontSize',12,'FontName','Verdana');
    set (gca, 'XTick',vx,'XMinorGrid','on')
    
%     plot(1000*t(1:nLmax,1),PMs(1:nLmax,4:5),'LineWidth',2);
%     xlabel('Tempo [ms]','FontSize',14,'FontName','Verdana');
%     set (gca,'XMinorGrid','on')

    ylabel('Amplitude [Pa]','LineWidth',2,'FontSize',12,'FontName','Verdana');
    legend('Microfone 3','Microfone 4')
    axis tight
    grid on
    ylim([-1 1.5])
    
else
    
    %%%%%%%%%%%%%%%
    % FIGURAS S#4 %
    %%%%%%%%%%%%%%%
    
    subplot(2,1,1)
    plot(x(1:nLmax+101,1),PMs(trig-100:nLmax+trig,1:2),'LineWidth',1);
    legend('Microfone 1','Microfone 2')
    xlabel('Dist�ncia Percorrida [m]','FontSize',12,'FontName','Verdana');
    ylabel('Amplitude [Pa]','LineWidth',2,'FontSize',12,'FontName','Verdana');
    axis tight
    set (gca, 'XTick',vx,'XMinorGrid','on')
    grid on
%     xlim([-0.5 Lmax])
%     ylim([-1.2 1.2])
    if (nchan > 2)
    subplot(2,1,2)
    plot(x(1:nLmax+101,1),PMs(trig-100:nLmax+trig,3:4),'LineWidth',1);
    legend('Microfone 3','Microfone 4')
    xlabel('Dist�ncia Percorrida [m]','FontSize',12,'FontName','Verdana');
    ylabel('Amplitude [Pa]','LineWidth',2,'FontSize',12,'FontName','Verdana');
    axis tight
    set (gca, 'XTick',vx,'XMinorGrid','on')
    grid on
%     xlim([-0.5 Lmax])
%     ylim([-1.2 1.2])
    end
end


% %Salvando a Imagem
% saveas(gcf,['Sinal - ',Arq], 'jpg')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Inicializa��o do Processamento    %
        %               da Medi��o              %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Definindo o tamanho da Janela de cada FFT
if (qpts < 6000)
LJan = 2;                                   % Tamanho da janela de an�lise em m
else LJan = 5;
end

TJan = LJan/c;                              % em segundos
npts = round(TJan/dt);                      % N�mero de Pontos dentro da Transformada de Fourier

% Configura��o da An�lise
zpf = 100;                                  % Porcentagem do sinal em zeropadding
olf = 10;                                   % Porcentagem do sinal sem overlapping (d� 10 cm de dEsp)
opts = floor((olf/100)*npts);               % Deslocamento do sinal (complemento do n�mero de Pontos de overlap)
zpts = floor((zpf/100)*npts);               % Pontos com zero (p/ zeropadding)
tpts = zpts+npts;                           % N�mero total de pontos

nFFTs = floor((length(PMs)-npts)/opts)+1;	% N�mero de blocos 

% Montagem do Vetor Espacial (Resolu��o espacial e limite de an�lise)
T = (opts/faq);                             % Resolu��o temporal da an�lise ]
dEsp = T*c;                                 % Resolu��o espacial da An�lise (FFTs calculadas de dEsp em dEsp)

binEsp = floor(2*Lmax/dEsp);                % Localiza o valor m�ximo a ser plotado no Eixo Espacial

% %Montagem do Vetor nos Gr�ficos
% 
% minEsp=0;                                   % Inicio do Gr�fico
% maxEsp = ceil(Lmax/dEsp);                   % Fim do Gr�fico
% deltaEsp = (maxEsp-minEsp)/40;                % Subdivis�ess
% vEsp = minEsp:deltaEsp:maxEsp;                      % Montagem do Vetor


% Detec��o da autom�tica posi��o do Pulso
[vMax,lMax] = max(abs(PMs));
bpulso = floor(lMax(:)/opts);

for knl=1:nchan
    if bpulso(knl) == 0
        bpulso(knl) = bpulso(knl)+1;
    end
end

    


CanalFixo = 1;                                  % Canal 1 � o gerador de Sinais

% Faixa "UTIL" dos vetores
Lu = floor(tpts/2);

% Gera��o do vetor de frequ�ncias
df = faq/tpts;                              % Calculo da resolu��o em frequ�ncias
% TvF = round(qpts/(2*nFFTs));              % Tamanho do meu vetor de frequ�ncias
vf = df*(0:npts/2-1)';                      % Gera��o do vetor de Frequ�ncia
% Defini��o de regi�o de interesse
binmin = floor(fmin/df)+1;
binmax = round(fmax/df)+1;

% Configura��o dos Graficos 3D
freq = df*(0:(Lu-1));
tempo = T*(0:(nFFTs-1));
espac = c*tempo/2;

% floor(min(espac)/10)*10;
% ceil(max(espac)/10)*10;

% Montagem da Janela
janela = zeros(tpts,1);
janela(1:npts) = tukeywin(npts,0.3);
% janela = hamming(tpts);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Processamento da Medi��o        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inicializa��o das Variaveis        
sft0 = zeros(tpts,nFFTs,nchan);                     %npts
ppm0 = zeros(tpts,nFFTs,nchan);                      %npts
ppm = zeros(tpts,nFFTs,nchan);                      %npts


% Separando os diversos canais em varios blocos.
for chan=1:nchan
    for j=1:nFFTs                                   % Preenchendo apenas com os "valores medidos" pulando o pr�-trigger
        pinicial = (j-1)*opts+1;                    % Ponto inicial do bloco
        pfinal = pinicial + npts;               	% Ponto final do bloco
        ppm0(1:npts,j,chan)=PMs(pinicial:pfinal-1,chan);
        ppm(:,j,chan) = ppm0(:,j,chan).*janela;
    end
end

sft0 = (sft0+fft(ppm))/npts;                       % Calculo das Transformadas de Fourrier.
                                            % Somando 1 para evitar divis�es por zero
% Avaliando a correla��o dos Blocos 

% % Determinando Aplitude e Fase das Transformadas de Fourier
AmpFFTs = zeros(size(sft0));                    % Inicializando a variavel
PhsFFTs = zeros(size(sft0));                    % Inicializando a variavel

for knl=1:nchan
    AmpFFTs(:,:,knl) = abs(sft0(:,:,knl));
    PhsFFTs(:,:,knl) = angle(sft0(:,:,knl));
end

AmpAux = AmpFFTs.*(2*(10^5));
AmpNPS = 10.*log10(AmpAux.^2);

% Truncando o Espectro para visualiza��o do Espectrograma
AmpNPS = AmpNPS(1:Lu,:,:);
PhsAngulo = PhsFFTs(1:Lu,:,:)*(180/pi);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             Espectrogramas            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 - Cria��o do Eixo X
Emin = round(min(espac)*10)/10;
Emax = round(espac(binEsp)*10)/10;
deltaE = round((Emax-Emin))/Lmax;
Ex = Emin:deltaE:Emax;

figure('Name',['Espectrograma dos Microfones da Medi��o ',Arq],'NumberTitle','off','PaperPositionMode','manual','Position',[100 100 1440 810])

if (nchan == 5)
    for knl=2:nchan    
    subplot(2,2,knl-1)
    surf(espac(1:binEsp+1),freq(binmin:binmax),AmpNPS(binmin:binmax,bpulso(knl):bpulso(knl)+binEsp,knl))
    view([45 45])
    title({['Espectrograma do Microfone ',num2str(knl-1)]},'FontSize',12,'FontName','Verdana')
    xlabel('Dist�ncia Percorrida pelo Pulso [m]','FontSize',12,'FontName','Verdana','Rotation',-20)
    ylabel('Frequ�ncia [Hz]','FontSize',12,'FontName','Verdana','Rotation',22)
    zlabel('NPS (dB ref 20\muPa)','FontSize',12,'FontName','Verdana')
    if qpts < 6000
    set (gca, 'XTick', Ex);
    end
    shading interp
    axis tight
    colorbar
    xlim([espac(1) Lmax])
    ylim([freq(binmin) freq(binmax)])
    zlim([30 110])
    grid on
    end
else
    for knl=1:nchan    
   subplot(2,2,knl)
    surf(espac(1:binEsp+1),freq(binmin:binmax),AmpNPS(binmin:binmax,bpulso(knl):bpulso(knl)+binEsp,knl))
    view([45 45])
    title({['Espectrograma do Microfone ',num2str(knl)]},'FontSize',12,'FontName','Verdana')
    xlabel('Dist�ncia Percorrida pelo Pulso [m]','FontSize',12,'FontName','Verdana','Rotation',-20)
    ylabel('Frequ�ncia [Hz]','FontSize',12,'FontName','Verdana','Rotation',22)
    zlabel('NPS (dB ref 20\muPa)','FontSize',12,'FontName','Verdana')
    if qpts < 6000
    set (gca, 'XTick', Ex);
    end
    shading interp
    axis tight
    colorbar
    xlim([espac(1) Lmax])    
    ylim([freq(binmin) freq(binmax)])
    grid on
    zlim([20 120])
    end
end          

%Salvando a Imagem
% saveas(gcf,['Espectrograma - ',Arq], 'jpg')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %               Reflet�ncia             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gerando os Espectros Cruzados
Cruzado = zeros(size(sft0));
Auto = zeros(tpts,nchan);
FRF = zeros(size(sft0));
Aux = zeros(size(sft0));
% Coher = zeros(size(sft0));

residue = 10^(-3);

for knl=1:nchan
    for jj=1:nFFTs
    Cruzado(:,jj,knl) = sft0(:,bpulso(knl),knl).*conj(sft0(:,jj,knl));
%     [Coher(:,jj,knl),fxy] = mscohere(ppm(:,bpulso(knl),knl),ppm(:,jj,knl));
    end
    Auto(:,knl) = Cruzado(:,bpulso(knl),knl);
    for jj=1:nFFTs
        Aux(:,jj,knl) = Cruzado(:,jj,knl)./(Auto(:,knl)+residue);
    end
end

ReRef = real(Aux);
ImRef = imag(Aux);
AmpRef = abs(Aux);
FasRef = unwrap(angle(Aux));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %        Gr�ficos de Reflet�ncia        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('Name',['Reflet�ncia dos Microfones ',Arq],'NumberTitle','off','Position',[100 100 1440 810],'PaperPositionMode','manual')

if (nchan == 5)
    for knl=2:nchan    
    subplot(2,2,knl-1)
%     figure('Name',['Reflet�ncia dos Canais ',num2str(knl)],'NumberTitle','off')
%     contour(espac(1:binEsp),freq(binmin:binmax),AmpRef(binmin:binmax,1:binEsp,knl),'LineWidth',2)
%     contour(espac(bpulso(knl):binEsp),freq(binmin:binmax),AmpRef(binmin:binmax,bpulso(knl):binEsp,knl))
    contour(espac(1:binEsp+1),freq(binmin:binmax),AmpRef(binmin:binmax,bpulso(knl):bpulso(knl)+binEsp,knl),'LineWidth',2,'LevelStep',0.05)
%     surfc(espac(bpulso(knl):binEsp),freq(binmin:binmax),AmpRef(binmin:binmax,bpulso(knl):binEsp,knl))
%     surfc(espac(:),freq(binmin:binmax),AmpRef(binmin:binmax,:,knl))
    title({['Reflet�ncia para Sinal do Microfone ',num2str(knl-1)]},'FontSize',12,'FontName','Verdana','LineWidth',2)
    xlabel('Dist�ncia Percorrida pelo Pulso [m]','FontSize',12,'FontName','Verdana')
    ylabel('Frequ�ncia [Hz]','FontSize',12,'FontName','Verdana')
    zlabel('Reflet�ncia','FontSize',12,'FontName','Verdana')
    % set (gca, 'XTick',vx,'XMinorGrid','on');
    shading interp
    axis tight
    colorbar
    xlim([espac(1) Lmax])
    ylim([freq(binmin) freq(binmax)])
    grid minor
%     zlim([0.01 1])  
    end
else
    for knl=1:nchan    
   subplot(2,2,knl)
%     surfc(espac(bpulso(knl):binEsp),freq(binmin:binmax),AmpRef(binmin:binmax,bpulso(knl):binEsp,knl))
%     contour(espac(1:binEsp),freq(binmin:binmax),AmpRef(binmin:binmax,1:binEsp,knl))
    contour(espac(1:binEsp+1),freq(:,:),AmpRef(1:size(freq,2),bpulso(knl):bpulso(knl)+binEsp,knl),'LineWidth',2,'LevelStep',0.05)
%     surf(espac(:,:),freq(binmin:binmax),AmpRef(binmin:binmax,:,knl))
    title({['Reflet�ncia para Sinal do Canal ',num2str(knl)]},'FontSize',12,'FontName','Verdana')
    xlabel('Dist�ncia Percorrida pelo Pulso [m]','FontSize',12,'FontName','Verdana')
    ylabel('Frequ�ncia [Hz]','FontSize',12,'FontName','Verdana')
    zlabel('Reflet�ncia','FontSize',12,'FontName','Verdana')
    % set (gca, 'XTick',vx,'XMinorGrid','on');
    shading interp
    axis tight
    colorbar
    xlim([espac(1) Lmax])
    ylim([freq(binmin) freq(binmax)])
    grid minor
%     zlim([0.01 1])  
    end
end

% %Salvando a Imagem
% set(gcf,'PaperPositionMode','manual');
% saveas(gcf,['Reflet�ncia - ',Arq], 'jpg')
% 
% 
% % Gravando o sinal dos microfones como sinal audivel
% wavwrite(PMs(:,1:nchan),faq,['Audivel - ',Arq,'.wav']) 
% save([Arq,' dEsp = ',num2str(LJan),' m - Dados.mat'])
