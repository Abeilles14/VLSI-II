function [ir_out]=gen_ctle(ir_in, ctle,time_step)
% function [ir_out]=gen_ctle(ir_in,A_dc,zeq,peq,gbw,time_step)

% if (strcmp(ctle.format, 'none(disable)') == 1)
%     ir_out = ir_in;
%     return;
% end

if ((strcmp(ctle.format, 'pole/zero') == 1) || ...
        (strcmp(ctle.format, 'none(disable)') == 1))

    kacdc = ctle.kacdc_fin;
    zeq = ctle.zero_fin;
    acbw = ctle.acbw_fin;
    gbw = ctle.gbw;
    gain_adj = ctle.gain_adjust_fac;


    A_dc=gbw/(acbw*kacdc);
    %account for gain adjustment (to limit range of CTLE output)
    A_dc=(A_dc*gain_adj);
    peq=zeq*kacdc;

    A_dc=max(A_dc,1e-3);
    zeq=max(zeq,1);
    peq=max(peq,1);

    gbw_rad=gbw*2*pi;
    zeq_rad=zeq*2*pi;
    peq_rad=peq*2*pi;

    ppar=(gbw_rad*zeq_rad)/(A_dc*peq_rad);

    a_tf=A_dc*peq_rad*ppar;
    b_tf=A_dc*peq_rad*ppar*zeq_rad;
    c_tf=zeq_rad;
    d_tf=(peq_rad+ppar)*zeq_rad;
    e_tf=peq_rad*ppar*zeq_rad;


    B_filt=[2*a_tf/time_step+b_tf;
        2*b_tf;
        b_tf-2*a_tf/time_step]';

    A_filt=[4*c_tf/time_step^2+2*d_tf/time_step+e_tf;
        2*e_tf-8*c_tf/time_step^2;
        4*c_tf/time_step^2-2*d_tf/time_step+e_tf]';

    B_filt=B_filt/A_filt(1);
    A_filt=A_filt/A_filt(1);

    ir_out=filter(B_filt,A_filt,ir_in);

elseif (strcmp(ctle.format, 'frequency table') == 1)
    fmax = ctle.fmax;
    %pick the appropriate CTLE response from the response set
    cidx = ctle.opt_freq_resp_idx;
    Hk = transpose(ctle.freq_resp_arr(:, cidx));

    Hk_ext = [Hk conj(Hk(end-1:-1:2))] ;

    hn_ext_ifft = real(ifft(Hk_ext, 'symmetric'));
    hn_ext_ifft(1) = (2 * hn_ext_ifft(1));
    %interpolate ifft
    N = length(ctle.f_arr);
    Ts = (1 / (2*fmax));
    t_arr = ((0 : (N-1)) * Ts); Tmax = t_arr(end); %10e-9; 
    treqd = (0 : time_step : Tmax);
    %interpolate step response instead of impulse response
    sr_ifft = cumsum(hn_ext_ifft(1:length(t_arr))); sr_ifft_t = t_arr;
    sr_interp = interp1(sr_ifft_t, sr_ifft, treqd, 'spline');
    ir_interp = diff(sr_interp);
    ir_ctle = ir_interp; clear ir_interp;
    %     ir_ctle = transpose((time_step / Ts) * interp1(t_arr, hn_ext_ifft(1:length(t_arr)), treqd, 'spline'));
    %scale by gain adjustment factor
    ir_ctle = (ir_ctle * ctle.gain_adjust_fac);
    ir_out = filter(ir_ctle, 1, ir_in);
    
else
    error('Error in CTLE format in gen_ctle!!!\n');
end
