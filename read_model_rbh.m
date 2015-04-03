function [thk,rho,vp,vs] = read_model_rbh(file)
fid = fopen(file);
for i=1:12
    fgetl(fid);
end
m = textscan(fid,repmat('%f',1,10));
thk = m{1};
 vp = m{2};
 vs = m{3};
rho = m{4};
thk(end) = [];