function wvcfs=filterphase(Ve,Fs)
% filter and then transform to phase in each position

wvcfs=zeros(size(Ve,1),size(Ve,2),size(Ve,3));
dim=size(Ve); % dim=[N,N,L]
filtered_data=Ve;
general.fmin = 0.5; %Hz 0.5;%
general.fmax = 9;%5; %Hz
general.attenuation = 20;%10;%dB
general.order = 3; %4 did not work
sample_frequency=Fs; % sampling frequency 25 Hz
[b, a] = cheby2(general.order,general.attenuation,[general.fmin*2/sample_frequency,general.fmax*2/sample_frequency]);%Riedner 2007
for y = 1:dim(1)
    for x = 1:dim(2)    
        cache =  filtfilt(b,a,double(Ve(y,x,:)));
        mcache = mean(cache);
        filtered_data(y,x,:) = cache - mcache;
        wvcfs(y,x,:)=hilbert(filtered_data(y,x,:));
    end
end

