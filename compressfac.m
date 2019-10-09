function Zl=compressfac(ind,xmol,ymol,t,P)
%COMPRESSFAC finds the compressability factor using Peng-Robinson Equation 
%    of State.
n=length(xmol);
index=reshape(ind,[length(ind) 1]);
load Critical_data Pcrit Tcrit omega
Pc=Pcrit(index);
Tc=Tcrit(index);
w=omega(index);

% reduced temperature and pressure
Tre=t./Tc;
Pre=P./Pc;
% m, a, Ai, Bi, Aij, A, B for the PR EOS
m=0.37464 + 1.54226.*w-0.26992.*w.^2; 
a=(1+m.*(1-Tre.^0.5)).^2;
Ap=0.45724.*a.*Pre./Tre.^2;
Bp=0.07780.*Pre./Tre;
Ab=zeros(n,n);
for i=1:n
    for j=1:n
        Ab(i,j)=(Ap(i)*Ap(j))^0.5;
    end
end
Av=0;
for i=1:n
    for j=1:n
        Av=Av+ymol(i)*ymol(j)*Ab(i,j);
    end
end
Bv=0;
for i=1:n
    Bv=Bv+ymol(i)*Bp(i);
end
Bl=0;
for i=1:n
    Bl=Bl+xmol(i)*Bp(i);
end
Al=0;
for i=1:n
    for j=1:n
        Al=Al+xmol(i)*xmol(j)*Ab(i,j);
    end
end
Alsum=zeros(n,1);
for i=1:n
    for j=1:n
        Alsum(i)=Alsum(i)+xmol(j)*Ab(i,j);
    end
end
Avsum=zeros(n,1);
for i=1:n
    for j=1:n
        Avsum(i)=Avsum(i)+ymol(j)*Ab(i,j);
    end
end
% liquid and gas phase compressibility factors
Zv=max(roots([1 -1+Bv Av-3*Bv^2-2*Bv -Av*Bv+Bv^2+Bv^3])); 
Zl=min(roots([1 -1+Bl Al-3*Bl^2-2*Bl -Al*Bl+Bl^2+Bl^3]));