% test-code for VME
% authors: Mojtaba Nazari and S.Mahmoud Sakhaei
% mojtaba.nazari.21@gmail.com -- smsakhaei@nit.ac.ir
%
% M. Nazari, S. M. Sakhaei, Variational Mode Extraction: A New Efficient
% Method to Derive Respiratory Signals from ECG, IEEE Journal of Biomedical
% and Health Informatics, Vol. 22, No. 4, pp. 1059-1067, july 2018.
% http://dx.doi.org/10.1109/JBHI.2017.2734074

%% --------------- Start
close all
clear
clc

%--------------- General Preparation
alpha = 20000; % Alpha- compactness of mode constraint
T = 1000;
fs=1000;
t = (1:T)/T;
tau = 0; % tau- time-step of the dual ascent.
tol = 1e-7; % tol- tolerance of convergence criterion





%% -------------------------Example 1
ex1c1=1./(1.2+cos(2*pi*t));
ex1c2=1./(1.5+sin(2*pi*t));
ex1c3=cos(32*pi*t+.2*cos(64*pi*t));
f_ex1=(ex1c1+ex1c2.*ex1c3);%+.2*randn(size(ex1c1));
omega_init_ex1 =0;%initial guess for center frequency (Hz) of desired mode
% omega_init_ex1=0.002;
% omega_init_ex1=0.006;
% omega_init_ex1=0.009;


%% --------------------------------Example 2
ex2c1=2*cos(4*pi*t);
ex2c2=cos(30*pi*t).*(1+cos(2*pi*t))/2;
ex2c3=cos(80*pi*t).*(1+sin(2*pi*t))/2;
f_ex2=ex2c1+ex2c2+ex2c3;%+.2*randn(size(ex2c1));
omega_init_ex2=10;%initial guess for center frequency (Hz) of desired mode
% omega_init_ex2=0.013;
% omega_init_ex2=0.016;
% omega_init_ex2=0.019;


%% -------------------------------Example 3a
ex3ac2=2*cos(10*pi*t+10*pi*t.^2);
ex3ac3_1=cos(60*pi*t);
ex3ac3_1(round(length(t)/2)+1:end)=0;
ex3ac3_1=ex3ac3_1(1:500);
ex3ac3_2=cos(100*pi*(t)-10*pi);
ex3ac3_2(1:round(length(t)/2))=0;
ex3ac3_2=ex3ac3_2(501:1000);
ex3ac3=[ex3ac3_1 ex3ac3_2];
f_ex3a=ex3ac2+ex3ac3;%+.2*randn(size(ex3ac2));
omega_init_ex3a=6;%initial guess for center frequency (Hz) of desired mode
% omega_init_ex3a=0.009;
% omega_init_ex3a=0.012;
% omega_init_ex3a=0.015;


%% -----------------------------Example 3b
ex3bc2=2*cos(10*pi*t+10*pi*t.^2);
ex3bc3_1=cos(60*pi*t);
ex3bc3_1(round(length(t)/2)+1:end)=0;
ex3bc3_11=ex3bc3_1(1:500);
ex3bc3_2=cos(100*pi*(t)-10*pi);
ex3bc3_2(1:round(length(t)/2))=0;
ex3bc3_22=ex3bc3_2(501:1000);
ex3bc3=[ex3bc3_11 ex3bc3_22];
f_ex3b=ex3bc2+ex3bc3;%+.2*randn(size(ex3bc2));
omega_init_ex3b=26;%initial guess for center frequency (Hz) of desired mode
% omega_init_ex3b=0.029;
% omega_init_ex3b=0.032;
% omega_init_ex3b=0.035;



%% General Visualization

for i=1:5
    
    switch i
        
        case 1
            f=f_ex1;
            omega=omega_init_ex1;
            reference=ex1c1;
            name='Example 1';
            
        case 2
            f=f_ex2;
            omega=omega_init_ex2;
            reference=ex2c2;
            name='Example 2';
            
        case 3
            f=f_ex3a;
            omega=omega_init_ex3a;
            reference=ex3ac2;
            name='Example 3a';
            
        case 4
            f=f_ex3b;
            omega=omega_init_ex3b;
            reference=ex3bc3_1;
            name='Example 3b';
            
            
    end
    
    u=vme(f,alpha,omega,fs,tau,tol);
    
    figure
    plot(t,reference(end,:)./max(reference(end,:)),'r-.')
    hold on
    plot(t,u./max(u))
    title(name)
    xlabel('t(sec)')
    ylim([min(u./max(u)) max(u./max(u))])
    
end


%% -------------------------- Real Data Example (ECG)
% load('055m.mat');
% f_ECG=val(1,:);
% omega_init_ECG =0;
% f=f_ECG;
% omega=omega_init_ECG;
% 
% u=vme(f,alpha,omega,tau,tol);
% 
% 
% figure
% t=60*((1:length(f))/length(f));
% plot(t,val(end,:)./max(val(end,:)),'r-.')
% hold on
% plot(t,u./max(u))
% title('Example ECG')
% xlabel('t(sec)')

