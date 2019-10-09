function [x,y,t]=idealflash2(index,input,beta)
%IDEALFLASH calculates the output of a ideal isothermal flash drum. in this
%   case for calculation of vapor pressure antoine equation has been used
%
%   [X,Y,T]=IDEALFLASH2(INDEX, INPUT,BETA)
%   INDEX is a vector which contains index of pure substances in the
%   mixture. the structure INPUt contains basic information for a flash
%   calculation such as Pressure, Temperature and mole fraction of
%   substances in the feed stream.
%
p=input.Pressure;
z=input.FeedMoleFraction;
n=length(index);
x=zeros(n,1);
y=zeros(n,1);
t0=boiling_point(index,p);
t0=sum(z.*t0);
option=optimset('Display','off');
t=fzero(@vle,t0,option);
if isnan(t)
    t0=boiling_point(index,p);
    tmax=max(t0);
    tmin=min(t0);
    t0=tmin+beta*(tmax-tmin);
    t=fzero(@vle,t0,option);
end
k=vaporpressure(index,t)/p;
for i=1:n
    x(i)=z(i)./(beta*k(i)+1-beta);
    y(i)=k(i)*x(i);
end

    function zz=vle(t)
        zz=0;
        k=vaporpressure(index,t)/p;
        for j=1:n
            zz=zz+z(j)*(k(j)-1)/(beta*k(j)+1-beta);
        end
    end
end