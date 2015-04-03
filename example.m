% 1-10 sec Love wave disperion from CVM model
function example
close all;
addpath('/home/yma/Codes/Surf.Codes/mat_disperse');

itheo = 0;
disp_theo  = 'disp_theo.txt';
model_theo = 'CVM_basin.mdl';

%%%
% creat 1~10 s dispersion curve from CVM basin model
if itheo    
    % create synthetic dispersion curve
    period = 1:0.25:10;
    period = period(:);
    nfreq= length(period);
    
    z = 0:1:40; % top of the layers
    
    % create the model
    [thk,dns,vp,vs] = create_model(z);
    
    % plot and write the model
    figure(1)
    plt_model_rbh(thk,dns,vp,vs,'-');
    write_model_rbh(model_theo,thk,vp,vs,dns);
    dlmwrite('period.txt',period,'delimiter',' ');
    
    % calculate the theoretical dispersion curve
    % Love for now!
    [vr,~,~,z,zdvrvs,zdvrrho]= mat_disperse(thk,dns,vp,vs,1./period,'L');
    [T_rbh,vr_rbh] = read_rbh;  % read benchmark from "rbh"
    
    % plot and write the disperion
    figure(2);
    plot(T_rbh,vr_rbh,'ko-');
    hold on;
    plot(period, vr,'ro-');
    xlabel('T'); ylabel('C');
    legend('rbh','this code');
    title('Dispersion curve');
    hold off;
    
    dlmwrite(disp_theo,[period(:) vr(:)],'delimiter',' ','precision',5);
    
    % plot the partial derivatives
    period1 = 1:2:8;
    [~,id]  = ismember(period1,period);
    z1 = z(:,id);
    zdvrvs  = squeeze(zdvrvs);
    zdvrrho = squeeze(zdvrrho);
    zdvrvs1 = zdvrvs(id,:);
    zdvrrho1= zdvrrho(id,:);
    figure(3)
    for i=1:length(period1)
        subplot(121)
        plot(zdvrvs1(i,:),z1(:,i))
        set(gca,'YDir','reverse');
        ylim([0 10]);
        hold all;
        subplot(122)
        plot(zdvrrho1(i,:),z1(:,i))
        set(gca,'YDir','reverse');
        ylim([0 10]);
        hold all;
    end
    subplot(121)
    legend('1','3','5','7');
    title('dc/d\beta');
    subplot(122)
    legend('1','3','5','7');
    title('dc/d\rho');
end
%%%

if ~itheo
    data   = load(disp_theo);
    period = data(:,1);
    freq   = 1./period;
    nfreq  = length(period);
    vr = data(:,2);
end

% add noise to the dispersion curve
sigma = 0.02 + zeros(nfreq,1);
err   = sigma.*randn(nfreq,1);
vr_exp= vr + err;


figure(2)
hold on;
plot(period,vr_exp,'b.-');
hold off;

dlmwrite(['disp_exp.txt_',num2str(sigma(1))],[period vr(:) vr_exp(:)],...
    'delimiter',' ','precision',5);

%%%
% start the inversion!
% create initial model
% parameter for the inversion
maxiter = 10;  
mu = 10;
tol_vs = 0.01;

% initial model, 15 km
z1 = 0:1:15;
[thk1,dns1,vp1,vs1] = create_model(z1);
% make a constant shift
dns1 = 120/100 * dns1;
vp1  = 120/100 * vp1;
vs1  = 120/100 * vs1;

[niter,vr_iter,vp_iter,vs_iter,dns_iter] = mat_inverse('L',freq,vr_exp,sigma,thk1,vp1,vs1,dns1,maxiter, mu, tol_vs);

fprintf('Total number of iterations: %d \n',niter);

figure(2)
hold on;
plot( period, vr_iter(:,niter),' ro-');

% plot the model
[thk0,dns0,vp0,vs0] = read_model_rbh(model_theo);

figure(4)
plt_model_rbh(thk0,dns0,vp0,vs0,'-');
hold on;
plt_model_rbh(thk1,dns_iter(:,niter),vp_iter(:,niter),vs_iter(:,niter),'o-');
hold off;

end

function [T,C] = read_rbh
% read love for now
file = 'rbh/SLEGN.ASC';
fid = fopen(file);
fgetl(fid);
m = textscan(fid,'%d %d %f %f %f %f %f %f');
fclose(fid);
T = m{3};
C = m{5};
end


function [thk,rho,vp,vs] = create_model(z)
% z is top of each layer, except the last value
% create the mid points
z1   = z(1:end-1);
z2   = z(2:end);
zmid = (z1 + z2)./2; 

vp_file = '/export/nobackup/yma/LASSIE/CVMh_ref/vp_basin.grd';
vs_file = '/export/nobackup/yma/LASSIE/CVMh_ref/vs_basin.grd';
rho_file= '/export/nobackup/yma/LASSIE/CVMh_ref/rho_basin.grd';

vp_data = load(vp_file);
vs_data = load(vs_file);
rho_data= load(rho_file);

vp = interp1(vp_data(:,1),vp_data(:,2),zmid);
vs = interp1(vs_data(:,1),vs_data(:,2),zmid);
rho= interp1(rho_data(:,1),rho_data(:,2),zmid);

thk= z2-z1;
thk= thk(:);
vp = [vp(:);vp(end)];
vs = [vs(:);vs(end)];
rho= [rho(:);rho(end)];  
end


