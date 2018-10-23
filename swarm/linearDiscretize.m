function [A,B] = linearDiscretize(s, u, params)
% s      - 16 x 1, state vector = [x, y, z, xd, yd, zd, qw, qx, qy, qz, p, q, r, F, M1, M2]
% u      - 3 x 1

    h = 0.00001;
    s_dim = size(s,1);
    u_dim = 3;
    
    A = zeros(s_dim);
    B = zeros(s_dim,u_dim);
%     c = zeros(s_dim,1);
    
    sr = s;
    sl = s;
    for i = 1:s_dim
        sr(i) = s(i) + h;
        sl(i) = s(i) - h;
        A(:,i) = (quadEOM_CO_readonly(sr,u,params)-quadEOM_CO_readonly(sl,u,params))/(2*h);
        sr = s;
        sl = s;
    end
    
    ur = u;
    ul = u;
    for i = 1:u_dim
        ur(i) = u(i) + h;
        ul(i) = u(i) - h;
        B(:,i) = (quadEOM_CO_readonly(s,ur,params)-quadEOM_CO_readonly(s,ul,params))/(2*h);
        ur = u;
        ul = u;
    end
end
%     temp = zeros(s_dim+u_dim+1);
%     temp(1:s_dim,1:s_dim) = A;
%     temp(1:s_dim, s_dim+1:s_dim+u_dim) = B;
%     temp = expm(t*temp);
%     
%     F = temp(1:s_dim, 1:s_dim);
%     G = temp(1:s_dim, s_dim+1:s_dim+u_dim);
% end
