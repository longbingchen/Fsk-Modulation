clc;
clear;
close all
%% parameter
fc = 2e6;%carrier frequency 10M
fs = 2e7; % Sampling frequency
T = 1/fs;  % 采样周期
N = 100;  % 每个符号采样点数
t = (0:N-1)*T; % T = 1/fs


%% 
M = 4;
k = log2(M);
s_m = 0:1:(M-1);
delta_f = 2*1/(2*N*T);  % △f_min = 1/(2*N*T)
fsk_set = eye(4);
%% 
bit_num = 20; % bit数
ber_test_data = fix(2*rand(bit_num,1)); %  bit data
ber_data = zeros(bit_num/k,1);
index = 1;
for i=1:k:bit_num
    for j = 0:1:k-1
        ber_data(index) =ber_data(index) + ber_test_data(i+j)*(2^(k-1-j));  % data 为M进制数
    end
    index = index + 1;
end
data = ber_data;

%%  4fsk信号空间 4维
f_1 = (2/N/T)*cos(2*pi*fc*t+ 2*pi*1*delta_f*t);
f_2 = (2/N/T)*cos(2*pi*fc*t+ 2*pi*2*delta_f*t);
f_3 = (2/N/T)* cos(2*pi*fc*t+ 2*pi*3*delta_f*t);
f_4 = (2/N/T)*cos(2*pi*fc*t+ 2*pi*4*delta_f*t);

s = zeros(length(t),length(data));

 for i = 1:1:length(data)	
    s(:,i) = cos(2*pi*fc*t + 2*pi*(data(i)+1)*delta_f*t);
 end

x = s(:);
%% channel
Eb_NO = 10;
snr_db = snr_cal(Eb_NO,k,N);
y = awgn(x,snr_db,'measured');    % awgn channel, 

figure(1);
subplot(2,1,2);
plot(1e6*(0:T:N*i*T-T),y); % after channel
title('信号时域图(经过信道后)');
xlabel('t(us)');
ylabel('y');

subplot(2,1,1);
plot(1e6*(0:T:N*i*T-T),x);
ylim([-1.5 1.5]);
title('信号时域图');
xlabel('t(us)');
ylabel('x');
%% power density spectrum
cor_x = xcorr(x); %计算序列的自相关函数
pds_x = fft(cor_x);
L = length(cor_x);
P2 = abs(pds_x/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

figure(4);
subplot(2,1,1);
plot(f/1e6,10*log10(P1));
xlim([1.6 3.8]);
ylim([-40 20]);
title('功率密度谱');
xlabel('f(MHZ)');
ylabel('幅度(dB)');

cor_y = xcorr(y); %计算序列的自相关函数
pds_y = fft(cor_y);
L = length(cor_y);
P2 = abs(pds_y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

figure(4);
subplot(2,1,2);
plot(f/1e6,10*log10(P1));
xlim([1.6 3.8]);
ylim([-20 20]);
title('功率密度谱(after channel Eb/N0 = 10dB)');
xlabel('f(MHZ)');
ylabel('幅度(dB)');
%% frequency domain
x_f = fft(x);
L = length(x);
P2 = abs(x_f/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

figure(2);
subplot(2,1,1);
plot(f/1e6,P1);
xlim([0 5]);
ylim([0 0.5]);
title('信号频域图');
xlabel('f(MHZ)');
ylabel('幅度');

%--------%

y_f = fft(y);
L = length(y);
P2 = abs(y_f/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

figure(2);
subplot(2,1,2);
plot(f/1e6,P1);
xlim([0 5]);
ylim([0 0.5]);
title('信号频域图(经过信道后)');
xlabel('f(MHZ)');
ylabel('幅度');

%% ber 
bit_error_table = bit_error_table_gen();
bit_error = 0;
Eb_N0_test = (-4:2:12);
L =100000;
ber = zeros(length(Eb_N0_test),L);
ber_result = zeros(length(Eb_N0_test),1);
j=1;
for m = 1:1:length(Eb_N0_test)
    for r = 1:1:L   %每个snr下 多次取平均值
    snr_db = snr_cal(Eb_N0_test(m),k,N);
    y = awgn(x,snr_db,'measured','dB');    % awgn channel, 
   
    bit_error = 0;
    r_k_1 = zeros(length(data),1);
    r_k_2 = zeros(length(data),1);
    r_k_3 = zeros(length(data),1);
    r_k_4 = zeros(length(data),1);
    Result = zeros(M,1);
    y_k = zeros(length(data),1); %判决输出
    for j = 1:1:length(data)
        r_k_1(j) = sum(y((j-1)*N+1:j*N).*f_1')*T;  % 解调输出
        r_k_2(j) = sum(y((j-1)*N+1:j*N).*f_2')*T;
        r_k_3(j) = sum(y((j-1)*N+1:j*N).*f_3')*T;  
        r_k_4(j) = sum(y((j-1)*N+1:j*N).*f_4')*T;  
            
        sm = [r_k_1(j),r_k_2(j),r_k_3(j),r_k_4(j)];
        res = fsk_set - sm; 
        for xx = 1:1:4       % 与每个点计算距离，取最小值
            Result(xx) = norm(res(xx,:));
        end
        [res_fsk,index] = min(Result);
        y_k(j) = s_m(index);
                bit_error = bit_error + bit_error_table(y_k(j)+1,ber_data(j) + 1);
    end
    ber(m,r)  = bit_error/bit_num;
    end
    ber_result(m) = sum(ber(m,:))/L;
end


%% 
figure(3);
semilogy(Eb_N0_test,ber_result,'red');
hold on;
theory_PM = M*qfunc(sqrt(k*10.^(Eb_N0_test/10))); %疏松上边界
theory_Pb = 2^(k-1)/(2^k -1)*theory_PM;
semilogy(Eb_N0_test,theory_Pb,'b--o');
grid on;
title('BER of 4-FSK');
xlabel('Eb/N0(dB)');
ylabel('误码率Pe');
legend('仿真结果','理论分析(疏松上边界)');









