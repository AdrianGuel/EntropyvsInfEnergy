function D=Dfunction4(D0,t,b,a,tx)
    D=D0+[(b(1)/(abs(a)*sqrt(pi)))*exp(-((t-tx(1))/a)^2),0,0,0;0,(b(2)/(abs(a)*sqrt(pi)))*exp(-((t-tx(2))/a)^2),0,0;
        0,0,(b(3)/(abs(a)*sqrt(pi)))*exp(-((t-tx(1))/a)^2),0;0,0,0,(b(4)/(abs(a)*sqrt(pi)))*exp(-((t-tx(2))/a)^2)];
end