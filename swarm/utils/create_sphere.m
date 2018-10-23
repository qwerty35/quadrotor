function points = create_sphere(NUM_POINTS, XYRADIUS, ZRADIUS)
    dlong = pi * (3-sqrt(5));
    dz = 2/NUM_POINTS;
    longt = 0;
    z = 1 - dz/2;
    points = zeros(3, NUM_POINTS);

    for i= 1:NUM_POINTS
        points(:,i) = [2*XYRADIUS*cos(longt)*sqrt(1-z*z); 2*XYRADIUS*sin(longt)*sqrt(1-z*z); 2*ZRADIUS*z];
        z = z - dz;
        longt = longt + dlong;
    end
end






