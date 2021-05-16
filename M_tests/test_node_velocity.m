close all
a = 10*rand-5;
b = 2*rand-1;
alp = 2.0;
xint = roots([(1+a^2),2*a*b,b^2-alp^2]);
yint = a*xint+b;
velx1 = alp./((1+a^2)*xint + a*b)

velx2 = sign(yint).*alp./(a*sqrt(alp^2 - xint.^2) + sign(yint).*xint)

fprintf('max abs diff = %g\n',max(abs(velx1-velx2)))
xpts = linspace(-alp,alp,100);
lineypts = a*xpts+b;
circUpPts = sqrt(alp^2 - xpts.^2);
circDwPts = -circUpPts;
figure
plot(xpts,lineypts,'k')
hold on
plot(xpts,circUpPts,'r',xpts,circDwPts,'r')
plot(xint,yint,'bo')
axis equal