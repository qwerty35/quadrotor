function tv = timevector(t,n,d)
% return d derivated nth order polynomial time vector
%
% input:  t    time
%         n    order of polynomial
%         d    number of derivatives

% output: tv   time vector
% if d = 0, [t^n t^n-1 ... t^2 t 1]
% if d = 1, [n*t^(n+1) (n-1)*t^(n-2) ...2t 1 0]


tv = ones(1,n+1);

for k = 1:d
    tv = polyder(tv);
end
tv(n+1) = (d == 0); % rescale to 1x(n+1)

for i = 1:n+1
    if n-(i-1)-d < 0
        tv(i) = 0;
    else
        tv(i) = tv(i) * t^(n-(i-1)-d);
    end
end

end

