function tb=boiling_point(ind,P)
% VAPORPRESSURE calculates the vapor pressure of a pure component using
% modified Antoine vapor pressure model
%   VAPORPRESSURE(IND, TEMP)
%   open the antoine coefficient data file and calculate vapor pressure at
%   temperature equals TEMP (Kelvin) for component with index number IND.
%   Vapor pressure unit is Bar
%
%   CAUTION: for every pure component antoine equation have a validation
%   temperature interval. if the input temperature is out of the interval
%   the result could be inappropriate.

tb=zeros(size(ind));
P=P*100;%convert Bar to kPa
load Antoine ant_coef
for i=1:length(ind)
    a=ant_coef(ind(i),:);
    t0=a(2)/(log(P)-a(1))-a(3);
    t1=t0;
    if boil(a,t0)>0
        while boil(a,t1)>0
            t1=t1*1.1;
        end
    else
        while boil(a,t1)<0
            t1=t1/1.1;
        end
    end
    tb(i)=fzero(@(temp) boil(a,temp),[t0,t1]);
end

function f=boil(a,temp)
    f=log(P)-(a(1)+a(2)/(temp+a(3))+a(4)*log(temp)+a(5)*temp^a(6));
end

end