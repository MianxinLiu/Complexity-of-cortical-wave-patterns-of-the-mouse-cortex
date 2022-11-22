p=3;
trial=9
fmin=0.5;
fmax=12;
% cd( '\\158.182.15.58\test2\Junhao');
% % File name
fileName = ['result_of_p' num2str(p) '_trial' num2str(trial)];
load(fileName);
mask=ones(50,50);
ratioSequence=VolE(1:50,1:50,:);
T=1:size(ratioSequence,3);
for surroN=3:4
    %% Phase randomized surrogate
    %% fft
    %fft
    
    %% temporal filtering highpass
    general.fmin = 0.1; %Hz 0.2;%
    general.attenuation = 20;%10;%dB
    general.order = 3; %4 did not work
    %attention!!
    general.srate = 150; %sampling rate [Hz]
    [b, a] = cheby2(general.order,general.attenuation,general.fmin*2/general.srate,'high' );%Riedner 2007

    data=zeros(size(mask,1),size(mask,2),length(T));%(x,y,z)z represents total time duration
    for j=1:size(data,1)
        for i=1:size(data,2)
            cache=filtfilt(b,a,squeeze(double(ratioSequence(i,j,T))));
            mcache = mean(cache);
            data(i,j,:) = cache - mcache;
        end
    end

    data_s=data;
    for i=1:size(data,1)
        for j=1:size(data,2)
            if(mask(i,j)>0)
            Y = fft(squeeze(data(i,j,:)));   
            pp = randn(size(data,3)/2-1,1);
            r=abs(Y);
            rr=zeros(size(data,3),1);
            rr(2:size(data,3)/2)=pp;
            rr(size(data,3)/2+2:size(data,3))=-pp(end:-1:1);
            YX= r.*cos(angle(Y)+2*pi*rr)+r.*1i.*sin(angle(Y)+2*pi*rr);
            YX(Y==0)=0;%avoid Nan in atan(imag(Y)./real(Y), Y=0, no amplitude in this frequency
            X = ifft(YX,'symmetric');
            data_s(i,j,:)=X;
            end
        end
    end

    %% spatial filtering
    Mask_in=zeros(size(mask,1),size(mask,2));
    index=0;
    for r=1:size(mask,1)
        for c=1:size(mask,2)
            index=index+1;
            Mask_in(r,c)=index;                         
        end
    end
    mask0=mask;
    mask0(mask<1)=0;
    Mask_in0=Mask_in.*mask0;
    Mask_in0(Mask_in0==0)=[];
    mask_rand = randperm(length(Mask_in0));
    data_s_s=data_s;
    for i=1:length(Mask_in0)
        [r,c]=find(Mask_in==Mask_in0(i));
        [r1,c1]=find(Mask_in==Mask_in0(mask_rand(i)));    
        data_s_s(r,c,:)=data_s(r1,c1,:);
    end
    [vx,vy,vfs,wvcfs]=optical_flow_model_gp_surrogate(data_s_s,fmin,fmax);
    
    % run('mousePatternDetectionMain.m')

%     cd( '\\158.182.15.58\test2\Junhao\freq_0p5to12\surrogate');
    save( ['result_of_p' num2str(p) '_trial' num2str(trial) 'surrogate' num2str(surroN) '.mat' ],'data_s_s','mask','vx','vy','vfs','wvcfs');

end
