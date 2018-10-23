t = 0.2;
T = 1;
       
[QptT,Apt0] = QpAp(t,T,t,0,3);
[Qp0t,AptT] = QpAp(0,t,t,T,3);
[Qp0T,Ap0t] = QpAp(0,T,0,t,3);

Qnzt = Qp0t(4:end,4:end);
QnztT = QptT(4:end,4:end);

Qnzt * [10*t^-3 6*t^-2 3*t^-1; -15*t^-4 -8*t^-3 -3*t^-2; 6*t^-5 3*t^-4 t^-3;] %* [T^3-t^3 3*(T^4-t^4) 6*(T^5-t^5); -3*(T^2-t^2) -8*(T^3-t^3) -15*(T^4-t^4); 3*(T-t) 6*(T^2-t^2) 10*(T^3-t^3);]

Ntest = 100;
result = zeros(Ntest);
for i = 1:Ntest
    for j = 1:Ntest
        sum = 0;
        for k = 0:i-1
            sum = sum + ((-1)^k)/(factorial(k)*factorial(j+i-1-k));
        end
        result(i,j) = sum - ((-1)^(i-1))/(factorial(i-1)*factorial(j-1)*(i+j-1));
    end
end