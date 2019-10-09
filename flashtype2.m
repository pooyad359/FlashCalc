function [x,y,t]=flashtype2(z,index,P,beta)
%FLASHTYPE2 solve a flash calculation with unknown temperature problem and 
%   determines the equilibrium data for the components. Peng-Robinson EOS 
%   is used in calculation of fugacity coefficients. This function use a 
%   library which contains the critical data CO2, CO, H2, H2O, C1-C30.
%   
%   *Raghu Raman Algorithm mixed with MATLAB functions*
%
%   [X,Y,T]=FLASHTYPE2(Z,INDEX,P,BETA)
%   Z: a vector which contains mole fraction of comonents in feed
%   INDEX: the index number of components
%   P: Pressure of flash calculation in Bar
%   BETA: amount of vapor phase to feed molar ratio
%   X: a vector which contains mole fraction of components in liquid phase
%   Y: a vector which contains mole fraction of components in vapor phase
%   T: Temperature of flash calculation in Kelvin
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
%   Beta=.5
%   P=2;    %Bar
%   [x,y,phi]=flashtype2(z,index,P,Beta)
%
%   Results:
%     x =
% 
%         0.4281
%         0.5719
%
%     y =
% 
%         0.7719
%         0.2281 
% 
%     T =
% 
%       300.8614
%
%   SEE ALSO FLASHTYPE1
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
input.Vap2Feed=beta;
[x,y,t]=idealflash2(index,input,beta);
k=vaporpressure(index,t);
x=z./(1-beta+beta*k);
x=x/sum(x);
y=x.*k;
y=y/sum(y);
while abs(flashtype2fun(t,x,y,input))>1e-10
    t=fdzero(@(t) flashtype2fun(t,x,y,input),t);
    input.Temperature=t;
    k=kvalue(x,y,input);
    x=z./(1-beta+beta*k);
    y=k.*x;
end
end

function f=flashtype2fun(t,x,y,input)
input.Temperature=t;
k=kvalue(x,y,input);
z=input.FeedMoleFraction;
beta=input.Vap2Feed;
psi=z.*(k-1)./(1-beta+beta*k);
f=sum(psi);
end