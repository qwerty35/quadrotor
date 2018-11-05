function collide = collision_inter_check(p, map, margin)
collide = 0;
if(size(p,1) <= 1)
    return;
end
% p(:,3) = p(:,3)/3; % scale z-axis by 3 to make it ellipsoid
% dis = pdist(p);
% if min(dis) < 2*margin
%     collide = 1;
% end
for i = 1:size(p,1)
    for j = i+1 : size(p,1)
        if norm(p(i,:)-p(j,:)) < 2*margin
            collide = 1;
        end
    end
end