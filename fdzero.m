function x=fdzero(funh,x0)
x=x0;
n=0;
while abs(funh(x))>1e-12
    n=n+1;
    if n>500
        error('Iteration number exceeds its maximum')
    end
    if x==0
        x1=-.01;
        x2=.01;
    else
        x1=x*.999;
        x2=x*1.001;
    end
    f1=funh(x1);
    f2=funh(x2);
    x=x-funh(x)*(x2-x1)/(f2-f1);
end