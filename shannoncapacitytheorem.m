% % Shannon Capactity Theorem -  Pe < 10^-5

target_BER = 1e-5;  % Target BER value

% 1st: Get relation between energy efficiency and BER for different
% modulation types:

% A) M'ary PSK:

Mpsk_values = [4, 8, 16, 32, 64, 128]; % range of M values to plot
enEffdb = linspace(0, 30, 100); 

% I) Compute the bit error rate (BER) for each energy efficiency and M value

BER_psk = zeros(length(enEffdb), length(Mpsk_values));
for i = 1:length(enEffdb)
    for j = 1:length(Mpsk_values)
        M = Mpsk_values(j);
        BER_psk(i, j) = berawgn(enEffdb(i), 'psk', M, 'nondiff');
    end
end

% II) Plot the results
figure;
semilogy(enEffdb, BER_psk, 'LineWidth',1.5);
grid on;
xlabel('Energy Efficiency (dB)');
ylabel('Bit Error Rate (BER)');
title('Bit Error Rate vs. Energy Efficiency for M-ary PSK modulation');
legend('M = 4', 'M = 8', 'M = 16', 'M = 32', 'M = 64', 'M = 128');
ylim([10e-8,1]);
xlim([0,30]);

% III) store corresponding En/N0 values to BER = 10^-5

enEffdb_target_psk = zeros(1, length(Mpsk_values));  % Array to store enEffdb values

% Find the enEffdb value for each M value at the target BER
for j = 1:length(Mpsk_values)
    M = Mpsk_values(j);
    [~, idx] = min(abs(BER_psk(:, j) - target_BER));
    enEffdb_target_psk(j) = enEffdb(idx);
end


% B) M'ary FSK:

Mfsk_values = [ 4, 8, 16, 32]; % range of M values to plot
enEff=10.^(enEffdb/10);


% I) Compute the bit error rate (BER) for each energy efficiency and M value

BER_fsk = zeros(length(enEffdb), length(Mfsk_values));
for j = 1:length(Mfsk_values)
    M = Mfsk_values(j);
    BER_fsk(:,j)=((M-1)/2)*exp(-1*enEff*log2(M)/2);
end

% II) Plot the results
figure;
semilogy(enEffdb, BER_fsk, 'LineWidth',1.5);
grid on;
xlabel('Energy Efficiency (dB)');
ylabel('Bit Error Rate (BER)');
title('Bit Error Rate vs. Energy Efficiency for M-ary FSK modulation');
legend('M = 4', 'M = 8', 'M = 16', 'M = 32');
ylim([10e-8,1]);
xlim([0,20]);

% III) store corresponding En/N0 values to BER = 10^-5

enEffdb_target_fsk = zeros(1, length(Mfsk_values));  % Array to store enEffdb values

for j = 1:length(Mfsk_values)
    M = Mfsk_values(j);
    [~, idx] = min(abs(BER_fsk(:, j) - target_BER));
    enEffdb_target_fsk(j) = enEffdb(idx);
end

% C) QAM

Mqam_values=[4,8,16,32,64];

% I) Compute the bit error rate (BER) for each energy efficiency and M value

BER_qam = zeros(length(enEffdb), length(Mqam_values));
for j = 1:length(Mqam_values)
    M = Mqam_values(j);
    BER_qam(:, j) = 4*(1 - 1/sqrt(M)) * qfunc(sqrt(3*log2(M)*enEff/(M-1)));
end

% II) Plot the results
figure;
semilogy(enEffdb, BER_qam, 'LineWidth', 1.5);
grid on;
xlabel('Energy Efficiency (dB)');
ylabel('Bit Error Rate (BER)');
title('Bit Error Rate vs. Energy Efficiency for M-QAM modulation');
legend('M = 4', 'M = 16', 'M = 64', 'M = 256');
ylim([10e-8, 1]);
xlim([0, 20]);


% III) Find the enEffdb value for each M value at the target BER
enEffdb_target_QAM = zeros(1, length(Mqam_values));  % Array to store enEffdb values
for j = 1:length(Mqam_values)
    M = Mqam_values(j);
    [~, idx] = min(abs(BER_qam(:, j) - target_BER));
    enEffdb_target_QAM(j) = enEffdb(idx);
end


% 2nd: apply shannon Capacity theorem

figure;
k =-0.5:0.001:25; 
EbN0=(2.^k-1)./k;
EbN0_db = 10*log10(EbN0);
semilogy(EbN0_db,k, 'LineWidth', 1.5);
xlabel('E_b/N_o (dB)');
ylabel('Spectral Efficiency (\eta)');
title('Channel Capacity & Power efficiency limit')
hold on;
grid on; 
xlim([-2 30]);ylim([0.1 15]);
yL = get(gca,'YLim');
line([-1.59 -1.59],yL,'Color','r','LineStyle','--','LineWidth', 1.5);


% MPSK BW-Energy Tradeoff

for i = 1:length(Mpsk_values)
    bwEff = log2(Mpsk_values(i));
    EbN0_mpsk = enEffdb_target_psk(i);
    semilogy(EbN0_mpsk, bwEff, 'Marker', 'x', 'DisplayName', sprintf('%d-PSK',Mpsk_values(i)),'LineWidth', 1.5);
end


% QAM BW-Energy Tradeoff

for i = 1:length(Mqam_values)
    bwEff = log2(Mqam_values(i));
    EbN0_qam = enEffdb_target_QAM(i);
    semilogy(EbN0_qam, bwEff, 'Marker', '*', 'DisplayName', sprintf('%d-QAM',Mqam_values(i)),'LineWidth', 1.5);
end

% MPSK BW-Energy Tradeoff

for i = 1:length(Mfsk_values)
    bwEff = 2*log2(Mfsk_values(i))/Mfsk_values(i);
    EbN0_fsk = enEffdb_target_fsk(i);
    semilogy(EbN0_fsk, bwEff, 'Marker', 'o', 'DisplayName', sprintf('%d-FSK',Mfsk_values(i)),'LineWidth', 1.5);
end


% 
% % Display the results
% fprintf('Target BER: %.1e\n', target_BER);
% for j = 1:length(Mqam_values)
%     fprintf('M = %d, enEffdb = %.2f dB\n', Mqam_values(j), enEffdb_target_QAM(j));
% end
% 

