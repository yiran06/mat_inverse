% plot thk model
function plt_model_rbh(thk,dns1,vp1,vs1,linetype)
tmp = cumsum(thk);
z1 = [0;tmp];
z2 = [tmp;tmp(end)+20];

NL = length(thk);

z  = zeros(1 + 2*NL + 1,1);
dns= z; vp=z; vs=z;

z(1) = z1(1);
z(2:2:end-1) = z2(1:end-1);
z(3:2:end-1) = z2(1:end-1);
z(end) = z2(end);

dns(1) = dns1(1);
dns(2:2:end-1) = dns1(1:end-1);
dns(3:2:end-1) = dns1(1:end-1);
dns(end)= dns1(end);

vp(1) = vp1(1);
vp(2:2:end-1) = vp1(1:end-1);
vp(3:2:end-1) = vp1(1:end-1);
vp(end)= vp1(end);

vs(1) = vs1(1);
vs(2:2:end-1) = vs1(1:end-1);
vs(3:2:end-1) = vs1(1:end-1);
vs(end)= vs1(end);

plot(vs,z,['r',linetype]);
hold on;
plot(vp,z,['b',linetype]);
plot(dns,z,['g',linetype]);
set(gca,'YDir','reverse');
hold off;
title('Input model');
legend('vs','vp','rho');

end