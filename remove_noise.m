% delay = 100;
lengthMult = 8;
fc=5000;
f0 = 440;
M = 30; % number of taps for adaptive filter
u=.005;

[primary, fsV] = audioread("noise.m4a");
primary = primary';
primary = repmat(primary, 1, lengthMult); % extend

[b,a] = butter(5, 2*fc/fsV);
primary = filter(b,a,primary); % filter 
N = length(primary);

Atone = .05*cos(2*pi*f0/fsV*(0:N-1)); % add desired signal
primary = primary + Atone;


% 
% reference = [zeros(1,delay) primary]; % reference is a delayed version of primary
% primary = [primary zeros(1,delay)]; % pad the end of primary input to match length

[reference, fsN] = audioread("noise.m4a"); % 1x288704
reference = reference';
reference = repmat(reference, 1, lengthMult); % extend


% sound(primary, fsV);
% pause(7);


% primary = repmat(primary, 1, 10); % more time for filter to converge?
% reference = repmat(reference, 1, 10);




% R = zeros(M,M);
% for i = M:N
%     stk = reference(i:-1:i-M+1)';
%     R = R + (stk*stk');    
% end
% R=R/(N-M+1);






% spectrumNoise = fftshift(fft(reference));
% spectrumNoiseRange = -fsN/2 : fsN/length(spectrumNoise) : (fsN/2 - fsN/length(spectrumNoise));



%% LMS Algorithm

W = zeros(M,N);
y = zeros(1,N);

for i = M:N-1
    R = reference(i:-1:i-M+1)'; % [r[n] r[n-1] ... r[n-M+1]]'
    currW = W(:,i);
    
    y(i) = primary(i) - currW'*R;
    W(:,i+1) = currW + 2*u*y(i).*R;
%     W(:,i+1) = currW + (1/(norm(R).^2))*y(i).*R;
end



%% Plots

figure;
plot(W');


% figure;
% plot(dataNoise.^2);
% title("Noise");
% ylim([0 .3]);
% 
% figure;
% plot(dataVoice.^2);
% title("Voice recording");
% ylim([0 .3]);


% -2500 to 2500 hz is approximately where the majority of the power is
% concentrated
% figure;
% plot(spectrumNoiseRange, 20*log10(abs(spectrumNoise) / mean(abs(spectrumNoise))));
% title("Power Spectrum of Noise Recording");


% 
% figure; 
% spectrogram(y, 2000, 1990, 2000, fsN, 'yaxis');
% title("Output of Adaptive Filter");
% ylim([0 10]);

% 
% figure; 
% spectrogram(dataVoice, 20, 10, 20, fsV, 'yaxis');
% title("Voice recording");




%% Sounds
% 

% sound(primary, fsV);
% pause(7);

% sound(primary, fsV);
% sound(reference, fsV);



sound([primary(1:2*fsV) y], fsV);
% pause(10);



