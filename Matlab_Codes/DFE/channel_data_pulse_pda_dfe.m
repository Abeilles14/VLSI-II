
% channel_data.m
% Sends data across a channel represented by an impulse response and produces an output wave

clear all;

bit_period=round(1e12/10e9);	% This is 10Gb/s with 1ps step time
opt_sample=1*(0e-2); % Used in plotting pulse response

% Load Channel Impulse Response
load ir_B1.mat;

tsample=1e-12;  % Impluse response has 1ps time step
sample_num=size(ir,2);  
ir_data=ir(1,:); 
scale_ir=1; % 1 for B12 (This is Vpp Differential)

sig_ir=ir_data*scale_ir;

% Plot Channel Impulse Response
%figure;
%plot(sig_ir);
%title('Channel Impulse Response');
time=(1:size(sig_ir,1))*1e-12;


% TX Equalization Function
eq_tap_number=3;    % Equalization Tap Number
precursor_number=0; % Number of Pre-Cursor Taps
num_bits=3;         % TX Equalization Resolution
% Calls TX Equalization Function
% taps output is un-quatized filter taps
% taps_quan output is filter taps quantized to equalizer resoltion
[taps,taps_quan] = tx_eq(sig_ir, bit_period, eq_tap_number, precursor_number, num_bits)


% For Pulse Response
 nt=110;
 m(1:20)=0; m(21)=1; m(22:110)=0;
% %m=-1*sign(m-0.5).^2+sign(m-0.5)+1;
% m_stream(1:20)=0; m_stream(21)=1; m_stream(22:27)=[0 0 0 1 0 1]; m_stream(28:100)=0;
% m_stream=-1*sign(m_stream-0.5).^2+sign(m_stream-0.5)+1;

% TX FIR Equalization Taps
% Set eq_taps equal to "taps" for infinite TX Eq resolution
% and "taps_quan" for finite TX Eq resolution
eq_taps = taps;  % or taps
m_fir = filter(eq_taps, 1, m);

%m_dr=reshape(repmat(m,bit_period,1),1,bit_period*size(m,2));
m_pulse = reshape(repmat(m,bit_period,1),1,bit_period*size(m,2));
m_dr = reshape(repmat(m_fir,bit_period,1),1,bit_period*size(m_fir,2));

% Note the 0.5 scale factor is to model a maximum logic swing of 0.5V per 
% symbol, or 1Vppd
data_channel=0.5*conv(sig_ir(1,:), m_dr(1:nt*bit_period));
data_channel_noeq=0.5*conv(sig_ir(1,:), m_pulse(1:nt*bit_period));


time_m=(1:size(m_pulse,2))*1e-12;
time_m_dr=(1:size(m_dr,2))*1e-12;

time_dc=(1:size(data_channel, 2))*1e-12;

%{
figure;
G=plot(time_m*1e9, m_pulse,'b');
hold on;
H=plot(time_m_dr*1e9, m_dr,'r');
hold off;
set(H, 'LineWidth', 1.0);
set(G, 'LineWidth', 1.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 5]);
set(AX, 'XTick', 0:0.5:5);

set(AX, 'YLim', [-0.5 1.5]); 
set(AX, 'YTick', 0:0.1:1);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (V)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Input Pulse With 3-Tap TX EQ','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
L=legend('No EQ', '3-Tap TX EQ');
set(L, 'FontSize', 14);
grid on;
%}

% 3 tap FIR response
figure;
H=plot(time_dc*1e9, 2*data_channel*1e3,'b');
set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 8]);
set(AX, 'XTick', 0:1:8);

set(AX, 'YLim', [-100 600]); 
set(AX, 'YTick', -100:100:600);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', '3-tap EQ Channel Pulse Response at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on;

save data_channel.mat data_channel; % Save Channel Output


% Uncomment below to plot pulse stuff
[max_data_ch, max_data_ch_idx] = max(abs(data_channel))
pulse_xaxis(1,:) = ((1:size(data_channel,2))-max_data_ch_idx)./bit_period;
sample_offset = opt_sample*bit_period;

% Take 10 pre-cursor, cursor, and 90 post-cursor samples
for i=1:101
    sample_points(i) = max_data_ch_idx + sample_offset + (i-11)*bit_period;
end
sample_values = data_channel(sample_points);
sample_points = (sample_points-max_data_ch_idx)./bit_period;

%{
figure;
plot(pulse_xaxis, 2*data_channel*1e3,'-b');
hold on;
stem(sample_points, 2*sample_values*1e3,'-or');
hold off;
set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [-4 14]);
set(AX, 'XTick', -4:2:14);

set(AX, 'YLim', [-100 600]); 
set(AX, 'YTick', -100:100:600);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', '3-tap EQ Channel Pulse Response at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on;
%}

% NO EQ
[max_data_ch_noeq, max_data_ch_idx_noeq] = max(abs(data_channel_noeq))
pulse_xaxis_noeq(1,:) = ((1:size(data_channel_noeq,2))-max_data_ch_idx_noeq)./bit_period;

% Take 10 pre-cursor, cursor, and 90 post-cursor samples
for i=1:101
    sample_points_noeq(i) = max_data_ch_idx_noeq + sample_offset + (i-11)*bit_period;
end
sample_values_noeq = data_channel_noeq(sample_points_noeq);
sample_points_noeq = (sample_points_noeq-max_data_ch_idx_noeq)./bit_period;


channel_delay = max_data_ch_idx-20*bit_period+1;

channel_delay_noeq = max_data_ch_idx_noeq-20*bit_period+1;

%%%%%%%%% DFE EQUALIZATION W/ 3-TAP %%%%%%%%
% Include DFE Equalization
dfe_tap_num = 3;
dfe_taps(1:dfe_tap_num) = sample_values(12:12+dfe_tap_num-1);
dfe_taps

dfe_taps_noeq(1:dfe_tap_num) = sample_values_noeq(12:12+dfe_tap_num-1);
dfe_taps_noeq

% Note this isn't a strtict DFE implementation - as I am not making a
% decision on the incoming data.  Rather, I am just using the known data
% that I transmitted, delay matching this with the channel data, and using 
% it to subtract the ISI after weighting with the tap values.  But, I think 
% it is good enough for these simulations.
m_dfe = filter(dfe_taps, 1, m);
m_dfe_dr = reshape(repmat(m_dfe,bit_period,1),1,bit_period*size(m_dfe,2));

m_dfe_noeq = filter(dfe_taps_noeq, 1, m);
m_dfe_dr_noeq = reshape(repmat(m_dfe_noeq,bit_period,1),1,bit_period*size(m_dfe_noeq,2));

time_dfe_dr=(1:size(m_dfe_dr, 2))*1e-12;

time_dfe_dr_noeq=(1:size(m_dfe_dr_noeq, 2))*1e-12;

figure;
%G=plot(time_m*1e9, m_pulse,'b');
%hold on;
%H=plot(time_m_dr*1e9, m_dr,'r');
K=plot(time_dfe_dr*1e9, m_dfe_dr, 'm');
%hold off;
%set(H, 'LineWidth', 1.0);
%set(G, 'LineWidth', 1.0);
set(K, 'LineWidth', 1.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 5]);
set(AX, 'XTick', 0:0.5:5);

set(AX, 'YLim', [-0.5 1.5]); 
set(AX, 'YTick', 0:0.1:1);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (V)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Input Pulse With 3-Tap TX EQ and DFE','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
%L=legend('No EQ', '3-Tap TX EQ');
%set(L, 'FontSize', 14);
grid on;


%data_channel = data_channel';
dfe_fb_offset = floor(bit_period/2); % Point at which the DFE taps are subtracted - can be anything from 0 to UI-1*time_step
data_channel_dfe = data_channel(channel_delay+dfe_fb_offset:channel_delay+dfe_fb_offset+size(m_dfe_dr,2)-1) - m_dfe_dr;

% Take samples for the pulse w/ DFE
clear pulse_xaxis;
[max_data_ch_dfe, max_data_ch_idx_dfe] = max(abs(data_channel_dfe));
pulse_xaxis(1,:)=((1:size(data_channel_dfe, 2))-max_data_ch_idx_dfe)./bit_period;
sample_offset = opt_sample*bit_period;

for k=1:101
    sample_points_dfe(k) = max_data_ch_idx_dfe + sample_offset + (k-11)*bit_period;
end
sample_values_dfe = data_channel_dfe(sample_points_dfe);
sample_points_dfe = (sample_points_dfe-max_data_ch_idx_dfe)./bit_period;

%%%%% NO EQv%%%%%%%%
%data_channel = data_channel';
data_channel_dfe_noeq = data_channel_noeq(channel_delay_noeq+dfe_fb_offset:channel_delay_noeq+dfe_fb_offset+size(m_dfe_dr_noeq,2)-1) - m_dfe_dr_noeq;

% Take samples for the pulse w/ DFE
clear pulse_xaxis_noeq;
[max_data_ch_dfe_noeq, max_data_ch_idx_dfe_noeq] = max(abs(data_channel_dfe_noeq));
pulse_xaxis_noeq(1,:)=((1:size(data_channel_dfe_noeq, 2))-max_data_ch_idx_dfe_noeq)./bit_period;

for k=1:101
    sample_points_dfe_noeq(k) = max_data_ch_idx_dfe_noeq + sample_offset + (k-11)*bit_period;
end
sample_values_dfe_noeq = data_channel_dfe_noeq(sample_points_dfe_noeq);
sample_points_dfe_noeq = (sample_points_dfe_noeq-max_data_ch_idx_dfe_noeq)./bit_period;

figure;
%J=plot(pulse_xaxis_noeq, 2*data_channel_dfe_noeq*1e3,'-b');
K=stem(sample_points_dfe_noeq, 2*sample_values_dfe_noeq*1e3,'-*b');
hold on;
%H=plot(pulse_xaxis, 2*data_channel_dfe*1e3,'-r');
%hold on;
I=stem(sample_points_dfe, 2*sample_values_dfe*1e3,'-*r'); 
hold off;
grid on;
%set(H, 'LineWidth', [2.0]);
set(I, 'LineWidth', [2.0]);
%set(J, 'LineWidth', [2.0]);
set(K, 'LineWidth', [2.0]);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', [16]);
set(AX, 'LineWidth', [2.0]);
set(AX, 'XLim', [-4 14]);
set(AX, 'XTick', -4:2:14);

set(AX, 'YLim', [-100 600]); 
set(AX, 'YTick', -100:100:600);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (UI)','FontName','utopia', 'FontSize', [20], 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', [20], 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', '3-Tap DFE Pulse Response','FontName','utopia', 'FontSize', [20], 'Color', [0 0 0]);
L=legend('No TX-FIR Filter', 'With 3-Tap TX Fitler');
set(L, 'FontSize', 14);


% Compute PDA
pda_eye_height=2*(sample_values(11)-(sum(abs(sample_values(1:10)))+sum(abs(sample_values(12:101)))))



