% reads an .s4p data file
% converts 4 port single ended s-param to 2 port mixed mode s-param
% (differential pairs)
% Plots forward transmission coeff in Freq domain
% Plots channel impulse response

clear all;
close all;
linestyle={'-b';'-r';'-g';'-m';'-k';'-c';'--b*';'--r*';'--g*'};
fsym = 10e9;
Tsym=1/fsym; %100e-12;	%%% Symbol Rate: e.g., Tsym = 1/fsym = 1/10 Gb/s
Ts=1e-12; %%% time step to sample, 1ps resolution measurement
nsym_short=500*100e-12/Tsym; 	%%% persistence of the impulse response


FileName = 'C:\Users\isabe\Documents\UBC_ELEC\ELEC403\project\ELEC403\peters_01_0605_B1_thru.s4p';

SingleEndedData = read(rfdata.data, FileName);
Freq = SingleEndedData.Freq;
% convert to 2-port mixed mode s-parameters (differential pair)
DifferentialSparams = s2sdd(SingleEndedData.S_Parameters);

% By default, S2SDD expects ports 1 & 3 to be inputs and ports 2 & 4 to be
% outputs. However if your data has ports 1 & 2 as inputs and ports 3 & 4
% as outputs, then use 2 as the second input argument to S2SDD to specify
% this alternate port arrangement. For example,
% DifferentialSparams = s2sdd(SingleEndedData.S_Parameters, 2);
rfdata.network('Data', DifferentialSparams, 'Freq', Freq);
thru_S11(1,:) = DifferentialSparams(1,1,:); thru_S21(1,:) = DifferentialSparams(2,1,:);
thru_S11_mag=20*log10(abs(thru_S11));
thru_S21_mag=20*log10(abs(thru_S21));


%figure; hold on;
%plot(Freq,thru_S21_mag,char(linestyle(1)));
%grid on;
%axis([0 15e9 -80 0]);
%ylabel('S21[dB]'); %

%%% Thru Channel %%%
H21=thru_S21;
imp=xfr_fn_to_imp(Freq',H21,Ts,Tsym); % Produces impulse response
imp_short=imp(1:floor(nsym_short*Tsym/Ts));
%figure; plot(imp,'b.-'); hold on; plot(imp_short,'r.-');
%ylabel('imp response'); legend('long','short');

% Sanity check the impulse response with fft
ch1=fft(imp);
ch1_freqs=(1/(Ts*size(ch1,2)))*(1:size(ch1,2))*1e-9; %(in GHz)

%figure;
%plot(ch1_freqs,20*log10(abs(ch1)),'-r',Freq'*1e-9,20*log10(abs(H21)),'-b');
%axis([0 15 -80 0]); grid on;
%ylabel('Freq response Check'); legend('rev','orig');

thru_imp_short=imp_short;
thru_imp=imp;

time=(1:size(imp,2))*1e-12;

% Impulse response plot
figure;
H=plot(time*1e9,imp*1e3,'b');
set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
%set(AX, 'XLim', [1.5 5]);
%set(AX, 'XTick', 1.5:Tsym*1e9:5);
set(AX, 'XLim', [0 5]);
set(AX, 'XTick', 0:1:5);
%set(AX, 'XTickLabel', {'0'; '5'; '10'; '15,-15'; '-10'; '-5'; '0'});
set(AX, 'YLim', [-0.5 8]); 
set(AX, 'YTick', 0:1:8);
%set(AX, 'YLim', [-0.55 0.55]); 
%set(AX, 'YTick', -0.5:0.1:0.5);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Impulse Response (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'B1 Backplane Channel Impulse Response','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
%L=legend('1','2','3','4',2);
%set(L, 'FontSize', [14]);
grid on;

% Forward Transmission Coefficient Frequency Plot
figure;
H=plot(ch1_freqs,20*log10(abs(ch1)),'-b',Freq'*1e-9,20*log10(abs(H21)),'-r');
set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
set(AX, 'XLim', [0 16]);
set(AX, 'XTick', 0:2:16);
%set(AX, 'XTickLabel', {'0'; '5'; '10'; '15,-15'; '-10'; '-5'; '0'});
set(AX, 'YLim', [-60 0]); 
set(AX, 'YTick', -60:10:0);
%set(AX, 'YLim', [-0.55 0.55]); 
%set(AX, 'YTick', -0.5:0.1:0.5);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Frequency (GHz)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'S_2_1 (dB)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'B1 Backplane Channel Frequency Response','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
L=legend('From Impulse Response','Measured');
set(L, 'FontSize', 14);
grid on;

ir=thru_imp_short; 
save ir_B1.mat ir;






