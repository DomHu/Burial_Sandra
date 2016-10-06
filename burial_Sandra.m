function [y] = burial_Sandra(x)
%
% This function runs Sandra's burial.m routine using the input variables
% and returns the associated burial fluxes at diffent depths and times
%
% [y,Q_sim,STATES,FLUXES] = hymod_nse(param,rain,evap,flow)
% 
% Input:
% x = vector of model parameters (k1, f1,  KNH4, gammaNH4, gammaH2S)  - vector (1,5)
%
% Output:   - all scalar y(1) - y(16)
%     F(1)=por(I10)*w0*oc(I10)
%     FF10=f(I10)
%     F(2)=por(I100)*w0*oc(I100)
%     FF100=f(I100)
%     F(3)=por(I1000)*w0*oc(I1000)
%     FF1000=f(I1000)
%     F(4)=por(I10000)*w0*oc(I10000)
%     FF10000=f(I10000)
% 
%     dF(1)=por(11)*w0*oc(11)
%     dFF10=f(11)
%     dF(2)=por(101)*w0*oc(101)
%     dFF100=f(101)
%     dF(3)=por(1001)*w0*oc(1001)
%     dFF1000=f(1001)
%     dF(4)=por(10001)*w0*oc(10001)
%     dFF10000=f(10001)
%

% clear all
% close all

M = 7 ; % number of model parameters
x = x(:);
if ~isnumeric(x); error('input argument ''param'' must be numeric'); end
if length(x)~=M; error('input argument ''param'' must have %d components',M); end

plot_fig = false;
% % a=100; % [1e-1-1e5] yrs
% % nu=0.1;% [0.1-2]
% % w0=0.1;% [10 1e-4] cm/yr
% % Db0=25;% [0.1 30] cm2/yr
% % OC0=1; % [0.01-30] wt%
% % beta=1e-3; %[1e-5 1e-7] por(z)=por0*exp(-beta*z)
% % por0=0.75; % [0.5 0.9]

a=x(1);     % [1e-1-1e5] yrs
nu=x(2);    % [0.1-2]
w0=x(3);    % [10 1e-4] cm/yr
Db0=x(4);   % [0.1 30] cm2/yr
OC0=x(5);   % [0.01-30] wt%
beta=x(6);  %[1e-5 1e-7] por(z)=por0*exp(-beta*z)
por0=x(7);  % [0.5 0.9]

zmax=10000;
zbio=10; 



for z=1:zmax+1
    depth=z-1;
    por(z)=por0*exp(-beta*depth);
    if (depth<zbio)
        k=nu/a;
        oc(z)=OC0*exp((w0-sqrt(w0*w0+4.0*Db0*k))/(2.0*Db0)*depth);
        f(z)=oc(z)/OC0;
        age(z)=0;
    else
        age(z)=(depth+1.0/beta*por0*(exp(-beta*depth)-1))/(w0*(1-por0));
        oc(z)=oc(zbio)*(a/(a+age(z)))^nu;
        f(z)=oc(z)/OC0;
    end
end

I10=max(find(age<10));
I100=max(find(age<100));
I1000=max(find(age<1000));
I10000=max(find(age<10000));

y(1)=por(I10)*w0*oc(I10);
y(2)=f(I10);
y(3)=por(I100)*w0*oc(I100);
y(4)=f(I100);
y(5)=por(I1000)*w0*oc(I1000);
y(6)=f(I1000);
y(7)=por(I10000)*w0*oc(I10000);
y(8)=f(I10000);

y(9)=por(11)*w0*oc(11);
y(10)=f(11);
y(11)=por(101)*w0*oc(101);
y(12)=f(101);
y(13)=por(1001)*w0*oc(1001);
y(14)=f(1001);
y(15)=por(10001)*w0*oc(10001);
y(16)=f(10001);

% % %   Output
% % F(1)=por(I10)*w0*oc(I10);
% % FF10=f(I10);
% % F(2)=por(I100)*w0*oc(I100);
% % FF100=f(I100);
% % F(3)=por(I1000)*w0*oc(I1000);
% % FF1000=f(I1000);
% % F(4)=por(I10000)*w0*oc(I10000);
% % FF10000=f(I10000);
% % 
% % dF(1)=por(11)*w0*oc(11);
% % dFF10=f(11);
% % dF(2)=por(101)*w0*oc(101);
% % dFF100=f(101);
% % dF(3)=por(1001)*w0*oc(1001);
% % dFF1000=f(1001);
% % dF(4)=por(10001)*w0*oc(10001);
% % dFF10000=f(10001);


if(plot_fig)
    figure
    subplot(2,2,1)
    plot(oc,-(1:zmax+1))

    subplot(2,2,2)
    plot(f,-(1:zmax+1))

    subplot(2,2,3)
    plot(age,-(1:zmax+1))

    subplot(2,2,4)
    plot(age, f)

    figure
    subplot(2,2,1)
    bar(F)
    subplot(2,2,2)
    bar(dF)
end

end