function [x,y,phi]=idealflash(index,input)
%IDEALFLASH calculates the output of a ideal isothermal flash drum. in this
%   case for calculation of vapor pressure antoine equation has been used
%
%   [X,Y,PHI]=IDEALFLASH(INDEX, INPUT)
%   INDEX is a vector which contains index of pure substances in the
%   mixture. the structure INPUt contains basic information for a flash
%   calculation such as Pressure, Temperature and mole fraction of
%   substances in the feed stream.

p=input.Pressure;
z=input.FeedMoleFraction;
t=input.Temperature;
n=length(index);
x=zeros(n,1);
y=zeros(n,1);
k=vaporpressure(index,t)/p;
option=optimset('Display','off');
phi=fsolve(@vle,.5,option);
for i=1:n
    x(i)=z(i)./(phi*k(i)+1-phi);
    y(i)=k(i).*x(i);
end

    function zz=vle(phi)
        zz=0;
        for j=1:n
            zz=zz+z(j)*(k(j)-1)/(phi*k(j)+1-phi);
        end
    end
end
