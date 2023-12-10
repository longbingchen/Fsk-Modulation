function [snr_db] = snr_cal(Eb_N0,k,N)
% ----output: snr_db :信噪比，以db为单位
% ----input:  Eb_N0 :Eb/N0,bit能量与功率密度谱比值,db
%             k:调制符号携带bit数
%             N:Tsym/Tsamp,每个符号的采样点数
 

Es_N0_db = Eb_N0+10*log10(k);
Es_N0 = 10.^(Es_N0_db/10);
SNR = Es_N0/(0.5*N);  %real signal ,*0.5
% SNR = Es_N0/(N);
snr_db = 10*log10(SNR);
end

