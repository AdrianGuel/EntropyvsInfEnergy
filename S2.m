function r=S2(S)
    x=eig(S);
    n=length(S);
    r=0;
        for j=1:n
            for i=1:j-1
                r=r+x(i)*x(j);
            end
        end
end