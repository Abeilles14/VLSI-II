
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

%{
% Plot Channel Impulse Response
time=(1:size(sig_ir, 2))*1e-12;  % X axis
figure;
H=plot(time*1e9, ir*1e3,'b');

set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
set(AX, 'XLim', [1.5 5]);
set(AX, 'XTick', 1.5:(1/10e9)*1e9:5);

set(AX, 'YLim', [-0.5 8]); 
set(AX, 'YTick', 0:0.5:8);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Impulse Response (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Backplane Channel Impulse Response','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on;

% Generate Random Data
%nt=1e3;         %number of bits
%m=rand(1,nt+1);     %random numbers between 1 and zero, will be quantized later
%m=-1*sign(m-0.5).^2+sign(m-0.5)+1;

%}
% For Pulse Response, create an input pulse from 0-1 for 1UI
nt = 100;
input_pulse = zeros(1, nt);    % array of zeros for nt points
input_pulse(19) = 1; % pulse 1 UI long at 19

% reshape to bitrate
input_pulse_reshaped = reshape(repmat(input_pulse, bit_length, 1), 1, length(input_pulse) * bit_length);

time1=(1:size(input_pulse_reshaped, 2))*1e-12;
%{
figure;
H=plot(time1*1e9, input_pulse_reshaped,'b');
set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
set(AX, 'XLim', [0 5]);
set(AX, 'XTick', 0:0.5:5);

set(AX, 'YLim', [-0.5 1.5]); 
set(AX, 'YTick', -0.5:0.5:1.5);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (V)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Input Pulse 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on;

%}

% TX FIR Equalization Taps
%eq_taps=1;
%eq_taps=[-0.046 0.78 -0.129 -0.0449]; % Random taps - not optimized for
%B1 channel
%m_fir=filter(eq_taps,1,m);

%m_dr=reshape(repmat(m_fir,bit_length,1), 1, bit_length*size(m_fir,2));
%m_tx=m_dr;  % pulse

% convolution of impulse response and pulse
%data_channel=0.5*conv(sig_ir(:,1),m_dr(1:nt*bit_length));

% SINGLE BIT RESPONSE SBR
% pulse response = convolution of pulse with impluse
pulse_response = conv(sig_ir(1,:), input_pulse_reshaped(1:nt*bit_length));
time2=(1:size(pulse_response, 2))*1e-12;

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

%{
hold on;
plot(x_interp_short, y_interp_short, 'o', 'MarkerFaceColor', 'red');
hold off;
%}

% overlay pulse and pulse reponse on plot
%{
figure;
H=plot(time1*1e9, input_pulse_reshaped*1e3,'r');
hold on;
H=plot(time2*1e9, pulse_response*1e3,'b');
hold off;

set(H, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
set(AX, 'XLim', [0 8]);
set(AX, 'XTick', 0:0.5:8);

set(AX, 'YLim', [-100 1100]); 
set(AX, 'YTick', -100:100:1100);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ns)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Channel Pulse and Response at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on;

save pulse_response.mat pulse_response; % Save Channel Output
%}

% SAMPLED PULSE RESPONSE
%length(y_interp)

figure;
stem(-4:10, y_interp_short, "filled");

AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
set(AX, 'XLim', [-4 10]);
set(AX, 'XTick', -4:1:10);

set(AX, 'YLim', [-100 800]); 
set(AX, 'YTick', -100:100:800);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (UI)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Sampled Pulse Response 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on


% FLIPPED PULSE RESPONSE
flipped_pulse = fliplr(pulse_response);

y_flip = flipped_pulse;
max_y_flip = max(y_flip);
max_x_flip = x(y_flip == max_y_flip);

% take 60 post cursors, 40 pre-cursors
x_interp_long = [max_x_flip-60*step:step:max_x_flip, max_x_flip+step:step:max_x_flip+40*step];
y_interp_long = interp1(x, y_flip, x_interp_long);

%{
figure;
stem(0:100, y_interp_long*1e3, "filled");
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
set(AX, 'XLim', [0 100]);
set(AX, 'XTick', 0:10:100);

set(AX, 'YLim', [-100 800]); 
set(AX, 'YTick', -100:100:800);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (UI)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Sampled Flipped Pulse Response 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on

%}

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

%{
figure;
subplot(2,1,1);
H=plot(time3, bits_worst0_reshaped, 'b');
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 100]);
set(AX, 'XTick', 0:10:100);

set(AX, 'YLim', [-0.2 1.2]); 
set(AX, 'YTick', 0:1:1);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (UI)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Bit Logic','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Worst Case 0 Bit Pattern','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);

subplot(2,1,2);
H = plot(time4, bits_worst1_reshaped, 'b');
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 100]);
set(AX, 'XTick', 0:10:100);

set(AX, 'YLim', [-0.2 1.2]); 
set(AX, 'YTick', 0:1:1);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (UI)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Bit Logic','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Worst Case 1 Bit Pattern','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
%}


bit0_response = conv(sig_ir(1,:), bits_worst0_reshaped(1:nt*bit_length));
time5=(1:size(bit0_response, 2))*1e-12;

bit1_response = conv(sig_ir(1,:), bits_worst1_reshaped(1:nt*bit_length));
time6=(1:size(bit1_response, 2))*1e-12;

%{
figure;
plot(time5*1e9, bit0_response*1e3, 'b');
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 20]);
set(AX, 'XTick', 0:1:20);

set(AX, 'YLim', [-10 1100]); 
set(AX, 'YTick', 0:100:1100);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (UI)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Bit Logic','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Worst Case 0 Bit Response','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);

figure;
plot(time6*1e9, bit1_response*1e3, 'b');
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 20]);
set(AX, 'XTick', 0:1:20);

set(AX, 'YLim', [-10 1100]); 
set(AX, 'YTick', 0:100:1100);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (UI)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Bit Logic','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Worst Case 1 Bit Response','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);

%}


% EYE DIAGRAM AT 10 GB/S RANDOM AND WORST CASE

% RANDOM EYE DIAGRAM SAMPLES
rand_nt = 1e3;         % number of bits 1000
%rand_pulse = rand(1, nt+1);     %random numbers between 1 and zero, will be quantized later
rand_pulse = randi([0 1], 1, rand_nt);
input_rand_reshaped = reshape(repmat(rand_pulse, bit_length, 1), 1, length(rand_pulse) * bit_length);

rand_pulse_response = conv(sig_ir(1,:), input_rand_reshaped(1:rand_nt*bit_length));

j=1;
offset=144;

for (i=55:floor(size(rand_pulse_response, 2) / bit_length)-500)
    eye_data(:,j) = 2*rand_pulse_response(floor((bit_length*(i-1)))+offset: floor((bit_length*(i+1)))+offset);
    j=j+1;
end;

time_rand = 0:2*bit_length;     % 2UI 200ps

%{
figure;
H=plot(time_rand, eye_data*5*1e2, 'LineWidth', 1.0);

AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 12);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 200]);
set(AX, 'XTick', 0:20:200);
set(AX, 'YLim', [-100 1100]); 
set(AX, 'YTick', 0:200:1000);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ps)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Random Data Eye Diagram at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
grid on;

%}

% WORST CASE EYE DIAGRAM
j=1;
offset=144;

for (i=55:floor(size(bit0_response, 2) / bit_length)-500)
    eye_data0(:,j) = 2*bit0_response(floor((bit_length*(i-1)))+offset: floor((bit_length*(i+1)))+offset);
    j=j+1;
end;

j=1;
offset=144;

for (i=55:floor(size(bit1_response, 2) / bit_length)-500)
    eye_data1(:,j) = 2*bit1_response(floor((bit_length*(i-1)))+offset: floor((bit_length*(i+1)))+offset);
    j=j+1;
end;

time_bit = 0:2*bit_length;     % 2UI 200ps
height_x = [88 88];
height_y = [439.499 526.776];
width_x = [62 112];
width_y = [481.762 484.838];

%{
figure;
H=plot(time_bit, eye_data0*5*1e2, 'LineWidth', 1.0);
hold on;
plot(time_bit, eye_data1*5*1e2, 'LineWidth', 1.0);
h1 = plot(height_x, height_y, 'og-', 'LineWidth', 1.0);  % height
h2 = plot(width_x, width_y, 'or-', 'LineWidth', 1.0);  % width
hold off;

AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 12);
set(AX, 'LineWidth', 1.0);
set(AX, 'XLim', [0 200]);
set(AX, 'XTick', 0:20:200);
set(AX, 'YLim', [-100 1100]); 
set(AX, 'YTick', 0:200:1000);
set(AX, 'YColor', [0 0 0]);
HX = get(AX, 'xlabel');
set(HX, 'string', 'Time (ps)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
HY = get(AX, 'ylabel');
set(HY, 'string', 'Voltage (mV)','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
Htitle = get(AX, 'title');
set(Htitle, 'string', 'Worst Case Eye Diagram at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);

L=legend([h1(1), h2(1)], 'Eye Height = 87mV', 'Eye Width = 50ps');
set(L, 'FontSize', 14);
grid on;
%}