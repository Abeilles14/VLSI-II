% channel_data.m
% Sends data across a channel represented by an impulse response and produces an output wave

clear all;
%feature('javafigures',0);

bit_period=round(1e12/10e9);	% This is 10Gb/s with 1ps step time
opt_sample=1*(0e-2); % Used in plotting pulse response

% Load Channel Impulse Response
load ir_B1.mat;
%load C:\Users\Sam\Documents\research\tamu\electrical_io\channels\ir_server_clean.mat;

tsample=1e-12;  % Impluse response has 1ps time step
sample_num=size(ir,2);  
ir_data=ir(1, :); 
scale_ir=1; % 1 for B12 (This is Vpp Differential)

sig_ir=ir_data*scale_ir;

% Plot Channel Impulse Response
time=(1:size(sig_ir,1))*1e-12;

%{
figure;
H=plot(time*1e9, sig_ir*1e3,'b');
set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
set(AX, 'XLim', [0 8]);
set(AX, 'XTick', 0:0.5:8);

set(AX, 'YLim', [-100 800]); 
set(AX, 'YTick', -100:100:800);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Channel Pulse Response at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on;
%}

% TX Equalization Function
eq_tap_number=3;    % Equalization Tap Number
precursor_number=1; % Number of Pre-Cursor Taps
num_bits=3;         % TX Equalization Resolution
% Calls TX Equalization Function
% taps output is un-quatized filter taps
% taps_quan output is filter taps quantized to equalizer resoltion
[taps, taps_quan] = tx_eq(sig_ir, bit_period, eq_tap_number, precursor_number, num_bits)


% For Pulse Response
nt = 100;
input_pulse = zeros(1, nt);    % array of zeros for nt points
input_pulse(19) = 1; % pulse 1 UI long at 19

% TX FIR Equalization Taps
% Set eq_taps equal to "taps" for infinite TX Eq resolution
% and "taps_quan" for finite TX Eq resolution
eq_taps = [taps];  % or taps
m_fir = filter(eq_taps, 1, input_pulse);

%m_dr=reshape(repmat(m_fir,bit_period,1),1,bit_period*size(m_fir,2));

% reshape to bitrate
input_pulse_reshaped = reshape(repmat(input_pulse, bit_period, 1), 1, length(input_pulse) * bit_period);
pulse_response = conv(sig_ir(1,:), input_pulse_reshaped(1:nt*bit_period));
time_init = (1:size(pulse_response, 2))*1e-12;

input_pulse_fir = reshape(repmat(m_fir, bit_period, 1), 1, length(m_fir) * bit_period);
time_m_dr=(1:size(input_pulse_fir, 2))*1e-12;
time_tap = (1:size(input_pulse_fir, 2))*1e-12;

figure;
H=plot(time_m_dr*1e9, input_pulse_reshaped,'r');
hold on;
G=plot(time_tap*1e9 - 0.1, input_pulse_fir,'b');
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
set(AX, 'YTick', 0:0.5:1);
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

% Note the 0.5 scale factor is to model a maximum logic swing of 0.5V per 
% symbol, or 1Vppd
pulse_response_eq = conv(sig_ir(1,:), input_pulse_fir(1:nt*bit_period));
time_dc=(1:size(pulse_response_eq, 2))*1e-12;

%{
figure;
H=plot(time_dc*1e9, pulse_response_eq*1e3,'b');
set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 8]);
set(AX, 'XTick', 0:0.5:8);

set(AX, 'YLim', [-100 800]); 
set(AX, 'YTick', -100:100:800);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Channel Pulse Response at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on;
%}

save pulse_response_eq.mat pulse_response_eq; % Save Channel Output

% Uncomment below to plot pulse stuff
[max_data_ch, max_data_ch_idx]=max(abs(pulse_response_eq));
pulse_xaxis(1,:) = ((1:size(pulse_response_eq,1))-max_data_ch_idx)./bit_period;

% Take 10 pre-cursor, cursor, and 90 post-cursor samples
sample_offset = opt_sample*bit_period;
for i=1:101
    sample_points(i)=max_data_ch_idx+sample_offset+(i-11)*bit_period;
end

sample_values = pulse_response_eq(sample_points);
sample_points = (sample_points-max_data_ch_idx)./bit_period;
%figure; plot(2*data_channel,'-b'); hold on; stem(sample_points,2*sample_values,'-*r'); grid on;

figure;
H=plot(pulse_xaxis, pulse_response_eq*1e3,'-r');
hold on;
I=stem(sample_points, sample_values*1e3,'-*r');
hold off;

%set(H, 'LineWidth', [2.0]);
set(I, 'LineWidth', [2.0]);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', [16]);
set(AX, 'LineWidth', [2.0]);
set(AX, 'XLim', [-3 10]);
set(AX, 'XTick', -3:1:10);
set(AX, 'YLim', [-50 600]); 
set(AX, 'YTick', 0:100:600);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (UI)','FontName','utopia', 'FontSize', [20], 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', [20], 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Sampled 3-Tap Equalized Pulse Response','FontName','utopia', 'FontSize', [20], 'Color', [0 0 0]);
grid on;

% overlay pulse responses w and w/o EQ
% Load Channel Impulse Response


figure;
F = plot(time_init*1e9, pulse_response*1e3,'-r');
hold on;
G = plot(time_dc*1e9 - 0.1, pulse_response_eq*1e3,'-b');
hold off;

set(F, 'LineWidth', 2.0);
set(G, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
set(AX, 'XLim', [3 8]);
set(AX, 'XTick', 3:1:8);

set(AX, 'YLim', [-100 600]); 
set(AX, 'YTick', -100:100:600);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Equalized Pulse Response Comparison at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
L=legend('No EQ', '3-Tap TX EQ');
set(L, 'FontSize', 14);
grid on;


% Compute PDA
pda_eye_height=2*(sample_values(11)-(sum(abs(sample_values(1:10)))+sum(abs(sample_values(12:101)))))




