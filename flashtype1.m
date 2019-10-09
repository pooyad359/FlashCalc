function [x,y,beta]=flashtype1(z,index,T,P)
%FLASHTYPE1 solve a flash calculation problem and determines the equilibrium
%   data for the components. Peng-Robinson EOS is used in calculation of
%   fugacity coefficients. This function use a library which contains the 
%   critical data CO2, CO, H2, H2O, C1-C30.
%   
%   [X,Y,PHI]=FLASHTYPE2(Z,INDEX,T,P)
%   Z: a vector which contains mole fraction of comonents in feed
%   INDEX: the index number of components
%   T: Temperature of flash calculation in Kelvin
%   P: Pressure of flash calculation in Bar
%   X: a vector which contains mole fraction of components in liquid phase
%   Y: a vector which contains mole fraction of components in vapor phase
%   BETA: amount of vapor phase to feed molar ratio
%
%   Index of Components in the library:
%
%   Index   Name
%   1       Carbon dioxide
%   2       Hydrogene
%   3       Water
%   4       Carbon monoxide
%   5       Methan
%   6       Ethan
%   7       Propane
%   8       i-butane
%   9       n-butane
%   10      i-pentane
%   11      n-pentane
%   12      n-hexane
%   13      n-heptane
%   14      n-octane
%   15      n-nonane
%   16      n-decane
%   17      nC11
%   18      nC12
%   19      nC13
%   20      nC14
%   21      nC15
%   22      nC16
%   23      nC17
%   24      nC18
%   25      nC19
%   26      nC20
%   27      nC21
%   28      nC22
%   29      nC23
%   30      nC24
%   31      nC25
%   32      nC26
%   33      nC27
%   34      nC28
%   35      nC29
%   36      nC30
%
%   EXAMPLE:
%   z=[.6;.4];
%   index=[8;11];    %I-BUTANE & N-PENTANE
%   T=300;  %Kelvine
%   P=2;    %Bar
%   [x,y,phi]=flashtype1(z,index,T,P)
%
%   Results:
%    x =
%        0.4456
%        0.5544
% 
%    y =
%        0.7854
%        0.2146
% 
%    phi =
%        0.4544
%
%   SEE ALSO FLASHTYPE2
if abs(sum(z)-1)>1e-8
    warning('Mole fractions normalaized. Sum of mole fractions must be equal to unity.')
    z=z/sum(z);
end
z=reshape(z,[length(z),1]);
index=reshape(index,[length(index) 1]);
load Critical_data Pcrit Tcrit omega
n=length(z);%calculate the number of components
Pc=Pcrit(index);
Tc=Tcrit(index);
w=omega(index);
%define a structure which contains equilibrium basic data
input.CriticalPressure=Pc;
input.CriticalTemperature=Tc;
input.AccentricFactor=w;
input.FeedMoleFraction=z;
input.Pressure=P;
input.Temperature=T;
try
    [x,y,tb]=flashtype2(z,index,P,0);
    if T<tb
        warning('Mixture is subcooled')
        x=z;
        y=zeros(n,1);
        beta=0;
        return
    end
    
catch ME
    warning(ME.identifier,'Unable to find Bubble point. The Mixture may be subcooled')
end
try
    [x,y,td]=flashtype2(z,index,P,1);
    if T>td
        warning('Mixture is superheated')
        y=z;
        x=zeros(n,1);
        beta=1;
        return
    end
catch ME
    warning(ME.identifier,'Unable to find Dew point. The Mixture may be superheated')
end
[x,y,beta]=idealflash(index,input);
k=kvalue(x,y,input);
if beta<0
    beta=.5;
end
while abs(flashtype1fun(beta,x,y,input))>1e-10
    beta=fzero(@(beta) flashtype1fun(beta,x,y,input),beta);
    x=z./(1-beta+k.*beta);
    y=x.*k;
    k=kvalue(x,y,input);
end
end

function f=flashtype1fun(beta,x,y,input)
z=input.FeedMoleFraction;
input.Vap2Feed=beta;
k=kvalue(x,y,input);
psi=z.*(k-1)./(1-beta+beta*k);
f=sum(psi);
end