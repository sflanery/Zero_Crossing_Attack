clc;clear;close all;

% setting parameters
%FCK: Three Resistors can be freely chosen, one has to be set such that
%net directional power flow is 0.
% 
%CHANGE TEMPS FOR ZERO CROSSING,USE TEMPERATURE FORMULAS IN ZC PAPER
k = 1.38e-23; %J/K
RAH = 46416; %Ohms
RAL = 278; %Ohms
RBL = 100;
%This is Bob's RH
RXH = (RAH * RBL) / RAL;
RXH = (RAH * RBL) / RAL
fB = 500; %Hz

%freely choose the RMS value of UAL
UAL_rms = 1;

TAL = UAL_rms^2/(4*k*RAL*fB); %UAL = UAL^2 = 1
TAH = (RAL/RAH) * TAL * (RBL *(RAH + RXH) + RAH * RXH + RAH^2) / (RAL^2 + RBL*(RAL+RXH) + RXH*RAL); %K
TBH = (RAL/RXH) * TAL * (RBL *(RAH + RXH) - RAH * RXH - RXH^2) / (RAL^2 + RBL*(RAL-RAH)- RAH*RAL);
TBL =(RAL/RBL) * TAL * (RBL *(RAH - RXH) - RAH * RXH + RBL^2) / (RAL^2 + RAL*(RXH-RAH)- RAH*RXH);

n = 1000; %sample size
iterations = 1000;

correct_guess = 'LH';
correct_guess_count = 0;

%Zero Cross Eve indexes and samples
Eve_index_LH = [];
Eve_sample_LH = [];
ms_zc_LH = [];

Eve_index_HL = [];
Eve_sample_HL = [];
ms_zc_HL = [];

ms_UwLH = [];
ms_IwLH = [];

ms_UwHL = [];
ms_IwHL = [];

for count = 1:1:iterations

    %noise voltages
    U_noise_AH = randn(1,n);
    U_noise_AH = sqrt(4*k*TAH*RAH*fB)*U_noise_AH;%/rms(U_noise_AH);
    
    U_noise_BH = randn(1,n);
    %Must change RBH, as this value is no longer arbitrary
    U_noise_BH = sqrt(4*k*TBH*RXH*fB)*U_noise_BH;%/rms(U_noise_BH);
    
    U_noise_AL = randn(1,n);
    U_noise_AL = sqrt(4*k*TAL*RAL*fB)*U_noise_AL;%/rms(U_noise_AL);
    
    U_noise_BL = randn(1,n);
    U_noise_BL = sqrt(4*k*TBL*RBL*fB)*U_noise_BL;%/rms(U_noise_BL);
    
    IwHL = (U_noise_AH-U_noise_BL)/(RAH+RBL);
    UwHL = IwHL.*RBL + U_noise_BL;
    
    IwLH = (U_noise_AL-U_noise_BH)/(RAL+RXH);
    UwLH = IwLH.*RXH + U_noise_BH;
    
    U2eff_LH = rms(UwLH)^2;
    I2eff_LH = rms(IwLH)^2;

    U2eff_HL = rms(UwHL)^2;
    I2eff_HL = rms(IwHL)^2;

    ms_UwLH = [ms_UwLH, U2eff_LH];
    ms_IwLH = [ms_IwLH, I2eff_LH];
    
    ms_UwHL = [ms_UwHL, U2eff_HL];
    ms_IwHL = [ms_IwHL, I2eff_HL];

    Eve_index_LH = [];
    Eve_sample_LH = [];

    Eve_index_HL = [];
    Eve_sample_HL = [];

    % Changed the lines of code, also looks sus
    for i = 1:1:n
        %Check if mean square voltage of BH and AL are the same and break (attempt single
        %bit exchange)
   
        if abs(UwLH(i))<= 1e-5
            Eve_index_LH = [Eve_index_LH,i];
            Eve_sample_LH = [Eve_sample_LH, IwLH(i)];
        end
        if abs(UwHL(i))<= 1e-5
            Eve_index_HL = [Eve_index_HL,i];
            Eve_sample_HL = [Eve_sample_HL, IwHL(i)];
        end
    end

    %RMS for HL and LH for each bit exchange
    ms_zc_LH = [ms_zc_LH, rms(Eve_sample_LH)^2];
    ms_zc_HL = [ms_zc_HL, rms(Eve_sample_HL)^2];
 
    %Power flow, should be approximately 0
    PHL = IwHL.*UwHL;
    PLH = IwLH.*UwLH;

    if ms_zc_HL(count) > ms_zc_LH(count)
        guess = 'HL';
    else
        guess = 'LH';
    end
    if (guess == correct_guess)
        correct_guess_count = correct_guess_count + 1;
    end
end
probability = correct_guess_count/iterations

figure;
%Line Plots
%Noise Voltage of AL V. BH
subplot(4,1,1);
plot(1:1:n,U_noise_AL,1:1:n,U_noise_BH);
title('Alice and Bob LH noise voltages');
xlim([0 100]);
%Noise Voltage of AH V BL (sample points determined by BL)
subplot(4,1,2)
plot(1:1:n,U_noise_AH,1:1:n,U_noise_BL);
title('Alice and Bob HL noise voltages');
xlim([0 100]);
%Wire current of AL + BH 
subplot(4,1,3);
plot(1:1:n,IwLH);
title('LH channel current');
xlim([0 100]);
%Wire current of AH + BL
subplot(4,1,4);
plot(1:1:n, IwHL);
title('HL channel current');
xlim([0,100]);

figure;
subplot(2,1,1);
plot(1:1:n,UwHL,Eve_index_HL,UwHL(Eve_index_HL),'*');
xlabel('UwHL');
subplot(2,1,2);
plot(1:1:n,IwHL,Eve_index_HL,IwHL(Eve_index_HL),'*');
ylabel('IwHL');

figure;
subplot(1,3,1);
histogram(ms_UwLH);
hold on;
histogram(ms_UwHL);
hold off;
legend('LH', 'HL')
xlabel('Uw^2');
subplot(1,3,2);
histogram(ms_IwLH);
hold on;
histogram(ms_IwHL);
hold off;
legend('LH', 'HL')
xlabel('Iw^2');
subplot(1,3,3);
histogram(ms_zc_LH);
hold on;
histogram(ms_zc_HL);
hold off;
legend('LH', 'HL')
xlabel('Uw_{zc}^2');
