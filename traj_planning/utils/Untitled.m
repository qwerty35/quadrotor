    r=3;
    
    C = zeros(Npoly+1);
    for i = 1:r
        C(i,:) = timevector(0,Npoly,r-i);
        C(i+r,:) = timevector(ts(end), Npoly, r-i);
    end