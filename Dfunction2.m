function D=Dfunction2(D0,t,b,a,tx)
    D=D0+[(b(1)/(abs(a)*sqrt(pi)))*exp(-((t-tx(1))/a)^2),0;0,(b(2)/(abs(a)*sqrt(pi)))*exp(-((t-tx(2))/a)^2)];
end
