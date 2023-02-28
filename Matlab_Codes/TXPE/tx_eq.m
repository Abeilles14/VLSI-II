function [taps,taps_quan]=tx_eq(ir_in,bit_period,eq_tap_number,precursor_number,num_bits)




%%Include CTLE
%%Model Peaking Amplifier
%ctle.kacdc_fin=4;
%ctle.zero_fin=0.5e9;
%ctle.acbw_fin=10e9;
%ctle.gbw=25e9;
%ctle.gain_adjust_fac=1.0;
%ctle.format='pole/zero';
%ir_post_ctle=gen_ctle(ir_in,ctle,tsample);


% No RX CTLE Equalization
ir_process=ir_in;
% With RX CTLE Equalization
%ir_process=ir_post_ctle;

precursor_samples = 10;
pulse_response = conv(ir_process, ones(1,bit_period));
time_res=(1:size(pulse_response, 2));

[max_pulse_value, max_pulse_time] = max(pulse_response);
sample_times = [max_pulse_time-1*precursor_samples*bit_period:bit_period:max_pulse_time+30*bit_period];
sample_values = pulse_response(sample_times);

%{
figure;
F = plot(time_res*1e-3, pulse_response*1e3,'-b');
hold on;
stem(sample_times*1e-3, sample_values*1e3,'-or');
hold off;
set(F, 'LineWidth', 2.0);
AX=gca;
set(AX, 'FontName', 'utopia');
set(AX, 'FontSize', 14);
set(AX, 'LineWidth', 2.0);
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
set(Htitle, 'string', 'Channel Pulse Response at 10Gb/s','FontName','utopia', 'FontSize', 20, 'Color', [0 0 0]);
L=legend('No EQ', '3-Tap TX EQ');
set(L, 'FontSize', 14);
grid on;
%}
eq_taps=3;
precursor_number=1;


% Construct P matrix
P(:, 1) = sample_values';      % vertical vec
P(size(sample_values, 2) + 1:size(sample_values, 2) + eq_tap_number - 1, 1) = 0;
% add 2 zeros at end, col 1 of P matrix

for i=2:eq_tap_number-1  % 2:3-1 = 2:2
    P(1:i-1, i) = 0;     % set (1, 1) = 0
    size(sample_values, 2)
    P(i:size(sample_values, 2)+i-1, i);
    P(i:size(sample_values, 2)+i-1, i) = sample_values';
    P(size(sample_values, 2)+i:size(sample_values, 2)+eq_tap_number-1,1) = 0;
end;

P(1:eq_tap_number-1, eq_tap_number) = 0;
P(eq_tap_number:eq_tap_number+size(sample_values, 2)-1,eq_tap_number) = sample_values';

% Construct Y matrix
Ydes(1:precursor_samples-1, 1) = 0;
Ydes(precursor_samples:precursor_samples+precursor_number, 1) = 0;
Ydes(precursor_samples+precursor_number+1, 1) = 1;
Ydes(precursor_samples+precursor_number+2:size(P,1), 1) = 0;

W = (P'*P)^(-1)*P'*Ydes

Itot = 1;
taps = Itot*W/sum(abs(W));
taps_abs = abs(taps);
taps_sign = sign(taps)

partition = [1/(2*(2^num_bits-1)):1/(2^num_bits-1):1-1/(2*(2^num_bits-1))];
codebook = [0:1/(2^num_bits-1):1];
[index, abs_taps_quan] = quantiz(abs(taps), partition, codebook);
taps_quan = taps_sign'.*abs_taps_quan;
taps_quan(precursor_number+1) = sign(taps_quan(precursor_number+1))*(1-(sum(abs(taps_quan))-abs(taps_quan(precursor_number+1))));


