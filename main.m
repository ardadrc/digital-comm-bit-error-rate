% ELEC 365 - Fundamentals of Digital Communications %
% Spring 2024 %
% MATLAB Project %
% Arda DERİCİ - 2001020020251 %

% BER expression of this system over AWGN for a) P(1) = P(0) = 1/2

% SNR values in dB are determined
SNRdB= 0:1:15;
% Let Eb = 1 and using the equation SNR_dB = 10*log(Eb/N0), 16 numerical
% values of N0 are found
N0 = 1./(10.^(SNRdB/10));

% for each SNR value, number of bits is determined as 100 million
num_bits = 10^8;

% for each SNR value, 'num_bits' of bits are generated with randi()
% command. randi() is uniformly distributed for equal probability.
% 'bits' consist of bits 0 and 1.
bits = randi([0, 1], length(SNRdB), num_bits);

% in the analytical derivations, it is calculated as a1 = 1 and a2 = -1
a1 = 1;
a2 = -1;
% a matrix of size 16 x 100 million is generated 
ai = zeros(length(SNRdB), num_bits);
% if the bit is 1, then ai equals a1; if the bit is 0, then ai equals a2
ai(bits == 0) = a2;
ai(bits == 1) = a1;

% a normal distribution of the form N(0, σ=1) is generated as noise
% with randn() command
n0 = randn(length(SNRdB), num_bits);
% for each SNR value, a normal distribution of the form N(0, σ^2=N0)
% is generated. also, σ^2=N0 is found in the analytical derivations
% it is known that if Y = α.X, then μY=α.μX, σY^2=α^2.σX^2
% thus, the α value is found as the square root N0
% so that the variance of the n0 Gaussian noise can be
% at the N0 value found in the analytical derivations
n0 = sqrt(N0') .* n0;

% it is calculated z = ai + n0 for each SNR value
z = ai + n0;
% in the analytical derivations, it is calculated as γ = 0
gama = 0;

% a zeros matrix is genrated to hold
% the number of error bits for each SNR value
num_bit_err = zeros(length(SNRdB), 1);
% the decision part for error detection is made with a 'for' loop
for i = 1:length(SNRdB)
    % for each SNR value, the number of bits received correctly
    % at the receiver when bit 1 is transmitted is found 
    ctrl_for_bit1 = sum(ai(i, z(i, :) > gama) == 1);
    % for each SNR value, the number of bits received correctly
    % at the receiver when bit 0 is transmitted is found
    ctrl_for_bit0 = sum(ai(i,z(i,:) < gama) == -1);
    % the number of bit errors is calculated for each SNR value
    num_bit_err(i) = num_bits - (ctrl_for_bit1 + ctrl_for_bit0);
end

% it is the bit error rate (BER) expression in the analytical derivations
% it is used qfunc() command
theo_ber = qfunc(sqrt(1./N0));

% plotting process
figure;
% theoretical BER curve plotted
semilogy(SNRdB, theo_ber, 'LineWidth', 2);
hold on
% simulated BER curve plotted
semilogy(SNRdB, num_bit_err./num_bits, 'rh', 'LineWidth', 4);
legend('theoretical BER curve', 'simulated BER curve', ...
    'FontWeight', 'bold', 'FontSize', 14);
xlabel('SNR (dB)', 'FontWeight', 'bold', 'FontSize', 14);
xticks(SNRdB);
ylabel('Bit Error Probability (P_b)', 'FontWeight', 'bold', ...
    'FontSize', 14);
set(gca, 'FontSize', 14);
title('Theoretical and Simulated BER Curves for Equal Probability', ...
    'FontSize', 14, 'FontWeight','bold');
grid on;

%%

% BER expression of this system over AWGN for b) P(1) = 1/4, P(0) = 3/4

% SNR values in dB are determined
SNRdB= 0:1:15;
% Let Eb = 1 and using the equation SNR_dB = 10*log(Eb/N0), 16 numerical
% values of N0 are found
N0 = 1./(10.^(SNRdB/10));

% for each SNR value, number of bits is determined as 100 million
num_bits = 10^8;

% probabilities are determined for bits 1 and 0.
P_1 = 1/4;
P_0 = 3/4;
% random bits are generated for each SNR value with rand() command
bits = rand(length(SNRdB), num_bits);
% 0 and 1 bits are determined according to probabilities
% 'bits' consist of bits 0 and 1.
bits = bits < P_1;

% in the analytical derivations, it is calculated as a1 = 1 and a2 = -1
a1 = 1;
a2 = -1;
% a matrix of size 16 x 100 million is generated 
ai = zeros(length(SNRdB), num_bits);
% if the bit is 1, then ai equals a1; if the bit is 0, then ai equals a2
ai(bits == 0) = a2;
ai(bits == 1) = a1;

% a normal distribution of the form N(0, σ=1) is generated as noise
% with randn() command
n0 = randn(length(SNRdB), num_bits);
% for each SNR value, a normal distribution of the form N(0, σ^2=N0)
% is generated. also, σ^2=N0 is found in the analytical derivations
% it is known that if Y = α.X, then μY=α.μX, σY^2=α^2.σX^2
% thus, the α value is found as the square root N0
% so that the variance of the n0 Gaussian noise can be
% at the N0 value found in the analytical derivations
n0 = sqrt(N0') .* n0;

% it is calculated z = ai + n0 for each SNR value
z = ai + n0;
% in the analytical derivations, it is calculated as γ = 0.549*N0
gama = 0.549*N0;

% a zeros matrix is genrated to hold
% the number of error bits for each SNR value
num_bit_err = zeros(length(SNRdB), 1);
% the decision part for error detection is made with a 'for' loop
for i = 1:length(SNRdB)
    % for each SNR value, the number of bits received correctly
    % at the receiver when bit 1 is transmitted is found 
    ctrl_for_bit1 = sum(ai(i,z(i,:)>gama(i)) == 1);
    % for each SNR value, the number of bits received correctly
    % at the receiver when bit 0 is transmitted is found
    ctrl_for_bit0 = sum(ai(i,z(i,:)<gama(i)) == -1);
    % the number of bit errors is calculated for each SNR value
    num_bit_err(i) = num_bits - (ctrl_for_bit1 + ctrl_for_bit0);
end

% it is the bit error rate (BER) expression in the analytical derivations
% it is used qfunc() command
theo_ber = (1 - qfunc((0.549*N0 - 1) ./ sqrt(N0)))*1/4 + qfunc( ...
    (0.549*N0 + 1) ./ sqrt(N0))*3/4;

% plotting process
figure;
% theoretical BER curve plotted
semilogy(SNRdB, theo_ber, 'LineWidth', 2);
hold on
% simulated BER curve plotted
semilogy(SNRdB, num_bit_err./num_bits, 'rh', 'LineWidth', 4);
legend('theoretical BER curve', 'simulated BER curve', ...
    'FontWeight', 'bold', 'FontSize', 14);
xlabel('SNR (dB)', 'FontWeight', 'bold', 'FontSize', 14);
xticks(SNRdB);
ylabel('Bit Error Probability (P_b)', 'FontWeight', 'bold', ...
    'FontSize', 14);
set(gca, 'FontSize', 14);
title('Theoretical and Simulated BER Curves for P(1) = 1/4, P(0) = 3/4', ...
    'FontSize', 14, 'FontWeight','bold');
grid on;