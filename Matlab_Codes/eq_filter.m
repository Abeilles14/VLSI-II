
% channel_data.m
% Sends data across a channel represented by an impulse response and produces an output wave

clear all;

% 1UI = 100ps
bit_length=round(1e12/10e9);	% This is 10Gb/s with 1ps step time
opt_sample=1*(0e-2); % Used in plotting pulse response

% Load Channel Impulse Response
load ir_B1.mat;

sample_num=size(ir,2);  
ir_data=ir(1,:); 
scale_ir=1; % 1 for B1 (This is Vpp Differential)

sig_ir=ir_data*scale_ir;

% For Pulse Response, create an input pulse from 0-1 for 1UI
nt = 100;
input_pulse = zeros(1, nt);    % array of zeros for nt points
input_pulse(19) = 1; % pulse 1 UI long at 19

% reshape to bitrate
input_pulse_reshaped = reshape(repmat(input_pulse, bit_length, 1), 1, length(input_pulse) * bit_length);

time1=(1:size(input_pulse_reshaped, 2))*1e-12;

% SINGLE BIT RESPONSE SBR
% pulse response = convolution of pulse with impluse
pulse_response = conv(sig_ir(1,:), input_pulse_reshaped(1:nt*bit_length));
time2=(1:size(pulse_response, 2))*1e-12;

% on pulse response, cursor seen at (4.130, 558.4164)
% find post and pre cursors, move 100ps steps from main cursor point
x = 1e9*time2;
y = 1e3*pulse_response;
max_cursor_y = max(y);
max_cursor_x = x(y == max_cursor_y);

% sample 4 pre-cursors and 10 post-cursors
step = 0.1;
x_interp_short = [max_cursor_x-4*step:step:max_cursor_x, max_cursor_x+step:step:max_cursor_x+10*step];
y_interp_short = interp1(x, y, x_interp_short);

% FLIPPED PULSE RESPONSE
flipped_pulse = fliplr(pulse_response);

y_flip = flipped_pulse;
max_y_flip = max(y_flip);
max_x_flip = x(y_flip == max_y_flip);

% take 60 post cursors, 40 pre-cursors
x_interp_long = [max_x_flip-60*step:step:max_x_flip, max_x_flip+step:step:max_x_flip+40*step];
y_interp_long = interp1(x, y_flip, x_interp_long);

% BIT PATTERN WORST CASE 0
% take flipped sampled pulse response, cursors > 0 = 1 else 0, cursor = 0
bits_worst0 = y_interp_long*1e3;
bits_worst0(bits_worst0 < 0) = 0;
bits_worst0(bits_worst0 > 0) = 1;
bits_worst0(x_interp_long == max_x_flip) = 0;   % set cursor to 0

% reshape to bitrate
bits_worst0_reshaped = reshape(repmat(bits_worst0, bit_length, 1), 1, length(bits_worst0) * bit_length);
time3=(1:size(bits_worst0_reshaped, 2))*1e-2;

% BIT PATTERN WORST CASE 1
% take flipped sampled pulse response, cursors > 0 = 0 else 1, cursor = 1
bits_worst1 = y_interp_long*1e3;
bits_worst1(bits_worst1 > 0) = 0;
bits_worst1(bits_worst1 < 0) = 1;
bits_worst1(x_interp_long == max_x_flip) = 1;   % set cursor to 1

% reshape to bitrate
bits_worst1_reshaped = reshape(repmat(bits_worst1, bit_length, 1), 1, length(bits_worst1) * bit_length);
time4=(1:size(bits_worst1_reshaped, 2))*1e-2;


bit0_response = conv(sig_ir(1,:), bits_worst0_reshaped(1:nt*bit_length));
time5=(1:size(bit0_response, 2))*1e-12;

bit1_response = conv(sig_ir(1,:), bits_worst1_reshaped(1:nt*bit_length));
time6=(1:size(bit1_response, 2))*1e-12;

% 3-TAP TX FIR EQUALIZATION FILTER
eq_taps = [1];      % need 3 taps
%eq_taps=[-0.046 0.78 -0.129 -0.0449]; % Random taps - not optimized for
%B12 channel
input_fir = filter(eq_taps, 1, input_pulse);

%m_dr=reshape(repmat(m,bit_period,1),1,bit_period*size(m,2));
pulse_fir = reshape(repmat(input_fir,bit_length,1), 1, bit_length*size(input_fir,2));

pulse_fir_response = 0.5*conv(sig_ir(:,1), input_fir(1:nt*bit_length));

figure;
plot(pulse_fir_response);
grid on;


%{
figure;
H=plot(time2*1e9, pulse_response*1e3,'b');
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
