function [snr_db] = snr_cal(Eb_N0,k,N)
% ----output: snr_db :����ȣ���dbΪ��λ
% ----input:  Eb_N0 :Eb/N0,bit�����빦���ܶ��ױ�ֵ,db
%             k:���Ʒ���Я��bit��
%             N:Tsym/Tsamp,ÿ�����ŵĲ�������
 

Es_N0_db = Eb_N0+10*log10(k);
Es_N0 = 10.^(Es_N0_db/10);
SNR = Es_N0/(0.5*N);  %real signal ,*0.5
% SNR = Es_N0/(N);
snr_db = 10*log10(SNR);
end

