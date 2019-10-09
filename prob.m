function prob
close all
z=[45.2,8.1,6.5,3.2,2.1,1.9,2.3,3.2,27.5]'/100;
Mw=[16.0425;30.069;44.0956;58.1222;58.1222;72.1488;72.1488;86.1754;142.282];
sum(z)
ind=[5:12,16]';
P=logspace(-1,2,20);
T=291.483;%K =65F
vf=zeros(size(P));% fraction of vapour phase
Vol_liq=zeros(size(P));% molar Volume of liquid phase Lit/mol
den_liq=zeros(size(P));% density of liquid phase kg/Lit
R=0.08314;% Bar.Lit/(mol.K)
for i=1:length(P)
    [x,y,b]=flashtype1(z,ind,T,P(i));
    zl=compressfac(ind,x,y,T,P(i));
    vol=zl.*R.*T./P(i);
    Vol_liq(i)=sum(x.*vol);
    den_liq(i)=sum(Mw.*x)./sum(vol.*x)./1000;
    vf(i)=b;
end
lf=1-vf;% fraction of liquid phase
gto=vf./lf;% gas to oil ratio
subplot(2,2,1)
plot(P,gto)
title('Gas to oil ratio vs Pressure')

subplot(2,2,2)
plot(P,Vol_liq)
title('molar volume of oil vs Pressure')

subplot(2,2,3)
plot(P,den_liq)
title('density of oil vs Pressure')

api=141.5./den_liq-131.5;
subplot(2,2,4)
plot(P,api)
title('API gravity of oil vs Pressure')