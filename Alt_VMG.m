clc;clear;close all;
%CHANGE TEMPS FOR ZERO CROSSING,USE TEMPERATURE FORMULAS IN ZC PAPER

% setting parameters
k = 1.38e-23; %J/K
RAH = 46400; %Ohms
RAL = 278; %Ohms
RBH = 278;
RBL = 100;
fB = 500; %Hz

%freely choose the RMS value of UAL
%Are we still able to freely choose the RMS value of the VMG?
UAL_rms = 1;

%TAL held constant
TAL = UAL_rms^2/(4*k*RAL*fB); %UAL = UAL^2 = 1
TAH = (RAL/RAH) * TAL * (RBL *(RAH + RBH) + RAH * RBH + RAH^2) / (RAL^2 + RBL*(RAL+RBH) + RBH*RAL);
TBH = (RAL/RBH) * TAL * (RBL *(RAH + RBH) - RAH * RBH - RBH^2) / (RAL^2 + RBL*(RAL-RAH)- RAH*RAL);
TBL =(RAL/RBL) * TAL * (RBL *(RAH - RBH) - RAH * RBH + RBL^2) / (RAL^2 + RAL*(RBH-RAH)- RAH*RBH);

n = 10000; %sample size
iterations = 1000;

correct_guess = 'HL';
correct_guess_count = 0;

%initialize mean square voltages and currents
ms_UwLH = [];
ms_IwLH = [];
ms_zc_LH = [];

ms_UwHL = [];
ms_IwHL = [];
ms_zc_HL = [];

for count = 1:1:iterations
    %noise voltages
    U_noise_AH = randn(1,n);
    %U_noise_AH = sqrt(4*k*TAH*RAH*fB)*U_noise_AH/rms(U_noise_AH);
    U_noise_AH = sqrt(4*k*TAH*RAH*fB)*U_noise_AH;
    
    U_noise_BH = randn(1,n);
    %U_noise_BH = sqrt(4*k*TBH*RBH*fB)*U_noise_BH/rms(U_noise_BH);
    U_noise_BH = sqrt(4*k*TBH*RBH*fB)*U_noise_BH;
    
    U_noise_AL = randn(1,n);
    %U_noise_AL = sqrt(4*k*TAL*RAL*fB)*U_noise_AL/rms(U_noise_AL);
    U_noise_AL = sqrt(4*k*TAL*RAL*fB)*U_noise_AL;
    
    U_noise_BL = randn(1,n);
    %U_noise_BL = sqrt(4*k*TBL*RBL*fB)*U_noise_BL/rms(U_noise_BL);
    U_noise_BL = sqrt(4*k*TBL*RBL*fB)*U_noise_BL;
    
    IwHL = (U_noise_AH-U_noise_BL)/(RAH+RBL);
    UwHL = IwHL.*RBL + U_noise_BL;
    
    IwLH = (U_noise_AL-U_noise_BH)/(RAL+RBH);
    UwLH = IwLH.*RBH + U_noise_BH;  

    U2eff_LH = rms(UwLH)^2;
    I2eff_LH = rms(IwLH)^2;

    U2eff_HL = rms(UwHL)^2;
    I2eff_HL = rms(IwHL)^2;

    %Power flow for LH and HL
    PHL = IwHL.*UwHL;
    PLH = IwLH.*UwLH;

    ms_UwLH = [ms_UwLH, U2eff_LH];
    ms_IwLH = [ms_IwLH, I2eff_LH];
    
    ms_UwHL = [ms_UwHL, U2eff_HL];
    ms_IwHL = [ms_IwHL, I2eff_HL];

    Eve_index_LH = [];
    Eve_sample_LH = [];

    Eve_index_HL = [];
    Eve_sample_HL = [];

    %I did the switching here
    for i = 1:1:n
        if abs(UwLH(i))< 1e-5
            Eve_index_LH = [Eve_index_LH,i];
            %Eve_sample_LH = [Eve_sample_LH, UwLH(i)];
            Eve_sample_LH = [Eve_sample_LH, IwLH(i)];
            %Eve_sample_LH = [Eve_sample_LH, U_noise_BH(i)];
        end
        if abs(UwHL(i))< 1e-5
            Eve_index_HL = [Eve_index_HL,i];
            %Eve_sample_HL = [Eve_sample_HL, UwHL(i)];
            Eve_sample_HL = [Eve_sample_HL, IwHL(i)];
            %Eve_sample_HL = [Eve_sample_HL, U_noise_BL(i)];
        end
    end

    %RMS for HL and LH for each bit exchange
    ms_zc_LH = [ms_zc_LH, rms(Eve_sample_LH)^2];
    ms_zc_HL = [ms_zc_HL, rms(Eve_sample_HL)^2];

    %Powerflow for LH and HL, should not be 0
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
PHL = mean(PHL);
PLH = mean(PLH);
%Print out
PHL;
PLH;
U2eff_LH;
U2eff_HL;
I2eff_HL;
I2eff_LH;
CCC_LH = corrcoef(UwLH, IwLH);
CCC_HL = corrcoef(UwHL, IwHL);
Verify_CCC_LH = PLH/(rms(UwLH)*rms(IwLH));
Verify_CCC_HL = PHL/(rms(UwHL)*rms(IwHL));

probability = correct_guess_count/count

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