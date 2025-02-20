clc;clear;close all;

% GOAL: INSTEAD OF SAMPLING WHEN THE CURRENT IS ZER0 (VAnoise = VBnoise),
% SAMPLE when u_noise is 0
% setting parameters
k = 1.38e-23; %J/K
Teff = 1e15; %K
RAH = 10000; %Ohms
RAL = 1000; %Ohms
RBH = 10000;
RBL = 1000;
fB = 500; %Hz

n = 1000; %sample size, points that Eve could sample during a bit exchange
iterations = 1000; %bit exchanges

correct_guess = 'HL';
correct_guess_count = 0;

%initialize mean-square voltages and currents
ms_UwLH = [];
ms_IwLH = [];
ms_zc_LH = [];

ms_UwHL = [];
ms_IwHL = [];
ms_zc_HL = [];

for count = 1:1:iterations
    %noise voltages
    U_noise_AH = sqrt(4*k*Teff*RAH*fB)*randn(1,n);
    U_noise_AH = U_noise_AH/rms(U_noise_AH);
    
    U_noise_BH = sqrt(4*k*Teff*RBH*fB)*randn(1,n);
    U_noise_BH = U_noise_BH/rms(U_noise_BH);

    U_noise_AL = sqrt(4*k*Teff*RAL*fB)*randn(1,n);
    U_noise_AL = U_noise_AL/rms(U_noise_AL);

    U_noise_BL = sqrt(4*k*Teff*RBL*fB)*randn(1,n);
    U_noise_BL = U_noise_BL/rms(U_noise_BL);

    IwHL = (U_noise_AH-U_noise_BL)/(RAH+RBL);
    UwHL = IwHL.*RBL + U_noise_BL;

    IwLH = (U_noise_AL-U_noise_BH)/(RAL+RBH);
    UwLH = IwLH.*RBH + U_noise_BH;

    U2eff_LH = rms(UwLH)^2;
    I2eff_LH = rms(IwLH)^2;

    U2eff_HL = rms(UwHL)^2;
    I2eff_HL = rms(IwHL)^2;

    ms_UwLH = [ms_UwLH, U2eff_LH];
    ms_IwLH = [ms_IwLH, I2eff_LH];
    
    ms_UwHL = [ms_UwHL, U2eff_HL];
    ms_IwHL = [ms_IwHL, I2eff_HL];

    %Zero Cross Eve indexes, samples, and U^2 noise rms
    Eve_index_LH = [];
    Eve_sample_LH = [];
    
    Eve_index_HL = [];
    Eve_sample_HL = [];
   
    %Bit exchanges
    %for i = 1:1:n
        %if abs(IwLH(i)) < 1e-5  %Check if current is near 0
            %Eve_index_LH = [Eve_index_LH,i];
            %Eve_sample_LH = [Eve_sample_LH, UwLH(i)];
        %end
        %if abs(IwHL(i)) < 1e-5
            %Eve_index_HL = [Eve_index_HL,i];
            %Eve_sample_HL = [Eve_sample_HL, UwHL(i)];
        %end 
    %end

    %Bit exchanges for DUALITY
    for i = 1:1:n
        if abs(UwLH(i)) < 1e-1  %Check if voltage is near 0
            Eve_index_LH = [Eve_index_LH,i];
            Eve_sample_LH = [Eve_sample_LH, IwLH(i)];
        end
        if abs(UwHL(i)) < 1e-1
            Eve_index_HL = [Eve_index_HL,i];
            Eve_sample_HL = [Eve_sample_HL, IwHL(i)];
        end 
    end

    %Powerflow for LH and HL, should be approximately 0
    PHL = IwHL.*UwHL;
    PLH = IwLH.*UwLH;


    %RMS voltage for HL and LH for each bit exchange
    %FIXME: Generates the same values for a certain amount of iterations
    ms_zc_LH = [ms_zc_LH, rms(Eve_sample_LH)^2];
    ms_zc_HL = [ms_zc_HL, rms(Eve_sample_HL)^2];
 
    %Eve guesses resistor configuraiton 
    if ms_zc_HL(count) > ms_zc_LH(count)
        guess = 'HL';
    else
        guess = 'LH';
    end
    if (guess == correct_guess)
        correct_guess_count = correct_guess_count + 1;
    end
    %correct_guess_count
end
PHL = mean(PHL);
PLH = mean(PLH);

%Print out
probability = correct_guess_count / iterations
zc_LH = mean(Eve_sample_LH.^2);
zc_HL = rms(Eve_sample_HL)^2;

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
