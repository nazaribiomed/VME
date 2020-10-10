function [u_d,u_hatd,omega]=vme(signal,Alpha,omega_int,fs,tau,tol)
%% Variational Mode Extraction
% authors: Mojtaba Nazari and S.Mahmoud Sakhaei
% mojtaba.nazari.21@gmail.com -- smsakhaei@nit.ac.ir 
% Initial release 2016-8-10 (c) 2016
%
%
%
% Input and Parameters:
%
% signal  - the time domain input signal 
% Alpha   - compactness of mode constraint
% omega_int - initial guess of mode center-frequency (Hz)
% tau     - time-step of the dual ascent. set it to 0 in the presence of
%           high level of noise.
% tol     - tolerance of convergence criterion; typically around 1e-6
%
%
%
% Output:
% u       - the desired mode
% u_hatd   - spectra of the desired mode
% omega   - estimated mode center-frequency
%
%
%
%
% M. Nazari, S. M. Sakhaei, Variational Mode Extraction: A New Efficient
% Method to Derive Respiratory Signals from ECG, IEEE Journal of Biomedical
% and Health Informatics, Vol. 22, No. 4, pp. 1059-1067, july 2018.
%
% http://dx.doi.org/10.1109/JBHI.2017.2734074
%

%% ------------ Part 1: Start initializing
save_T = length(signal);

%--------------- Mirroring the signal to extend
T = save_T;
f_mir=zeros(1,T/2);
f_mir(1:T/2) = signal(T/2:-1:1);
f_mir(T/2+1:3*T/2) = signal;
f_mir(3*T/2+1:2*T) = signal(T:-1:T/2+1);
f = f_mir;

%------------- Time Domain (t -->> 0 to T)
T = length(f);
t = (1:T)/T;


udiff = tol+eps; %------ update step

%------------- Discretization of spectral domain 
omega_axis = t-0.5-1/T;


%-------- FFT of signal(and Hilbert transform concept=making it one-sided)
f_hat = fftshift((fft(f)));
f_hat_onesided = f_hat;
f_hat_onesided(1:T/2) = 0;

%------------ Max. number of iterations
N = 300;

%----------- Initializing omega_d
omega_d = zeros(N, 1);
omega_d(1) = omega_int/fs;


%----------- dual variables vector
lambda = zeros(N, length(omega_axis));

%---------- keeping changes of mode spectrum
u_hat_d = zeros(N, length(omega_axis));

%---------- keeping changes of residual spectrum
%F_r=zeros(N, length(omega_axis));

n = 1; %------------------ main loop counter


%% -------------- Part 2: Main loop for iterative updates

while ( udiff > tol &&  n < N ) 
    
%------------------ update ud 
    u_hat_d(n+1,:) = (f_hat_onesided+(u_hat_d(n,:).*...
        (Alpha.^2).*(omega_axis - omega_d(n)).^4)+lambda(n,:)/2)...
        ./((1+(Alpha.^2).*(omega_axis - omega_d(n)).^4).*...
        (1+2*Alpha^2.*(omega_axis - omega_d(n)).^4));
    
    
%------------------ update omega_d     
   omega_d(n+1) = (omega_axis(T/2+1:T)*(abs(u_hat_d(n+1, T/2+1:T)).^2)')...
        /sum(abs(u_hat_d(n+1,T/2+1:T)).^2);
    

% %------------------ update F_r
%  
%  F_r(n+1,:) = ((Alpha^2.*(omega_axis - omega_d(n+1)).^4).*...
%      (f_hat_onesided-(u_hat_d(n+1,:))))./(1+2*Alpha^2.*...
%       (omega_axis - omega_d(n+1)).^4);
    


%-----update lambda (dual ascent) ===> lambda = lambda + tau*(f(t)-(Ud+Fr))  

    lambda(n+1,:) = lambda(n,:) + (tau*(f_hat_onesided-(u_hat_d(n+1,:)+...
        ((Alpha^2.*(omega_axis - omega_d(n+1)).^4).*(f_hat_onesided-...
       (u_hat_d(n+1,:))))./(1+2*Alpha^2.*(omega_axis -omega_d(n+1)).^4))));  
    
   
    n = n+1;
    
    udiff = eps;
%------------------ 1st loop criterion    
    udiff = udiff + 1/T*(u_hat_d(n,:)-u_hat_d(n-1,:))*conj((u_hat_d(n,:)...
        -u_hat_d(n-1,:)))';
    
    udiff = abs(udiff);
    
end

%% ------------------ Part 3: Signal Reconstruction

N = min(N,n);
omega = omega_d(1:N,:);
u_hatd = zeros(T, 1);
u_hatd((T/2+1):T,:) = squeeze(u_hat_d(N,(T/2+1):T,:));
u_hatd((T/2+1):-1:2,:) = squeeze(conj(u_hat_d(N,(T/2+1):T,:)));
u_hatd(1,:) = conj(u_hatd(end,:));

u_d = zeros(1,length(t));
u_d(1,:)=real(ifft(ifftshift(u_hatd(:,1))));

%----------------- Remove mirror part
u_d = u_d(:,T/4+1:3*T/4);

u_hatd=fftshift(fft(u_d(1,:)))';
end
