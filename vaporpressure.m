function vp=vaporpressure(ind,temp)
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

vp=zeros(size(ind));
load Antoine ant_coef
for i=1:length(ind)
    a=ant_coef(ind(i),:);
    vp(i)=a(1)+a(2)/(temp+a(3))+a(4)*log(temp)+a(5)*temp^a(6);
    vp(i)=exp(vp(i));
end
%convert kPa to Bar
vp=vp/100;