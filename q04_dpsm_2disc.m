clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter
mm=1e-3; kHz=1e3; j=1i;
c=340; rho=1.293;
f=20*kHz; w=2*pi*f; k=w/c;

a0=50*mm; %[m] disc radius
v0=1; %[m/s] disc velocity
height=2*c/f; %上円盤高さ

%observing plane
x=[-2*a0:a0/80:2*a0];
z=[0:height/120:height];
[X,Z] = meshgrid(x,z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disk division
rs = a0/10;
[lw_disc, lw_norms]=def_circ3(rs, a0, [0 0 0], [0 0 1]); %下面円盤
[up_disc, up_norms]=def_circ3(rs, a0, [0 0 height], [0 0 -1]); %上面円盤

pos_tgt=[lw_disc; up_disc];
norms  =[lw_norms;up_norms];
pos_src = pos_tgt - rs*norms ; %未知音源は放射面より -rsだけ下げた位置に配置する
figure(1); 
plot3(pos_tgt(:,1), pos_tgt(:,2), pos_tgt(:,3), 'o', pos_src(:,1), pos_src(:,2), pos_src(:,3), 'x');
xlabel('x'); ylabel('y'); zlabel('z');

[lw_num, ~] = size(lw_disc);
[up_num, ~] = size(up_disc);
[num_src, ~] = size(pos_tgt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine amplitude
dps_M = zeros(num_src, num_src);
dps_V = [v0*ones(lw_num,1); zeros(up_num,1)];

for si=1:num_src %source loop
    for ti=1:num_src %target loop
        
        %src -> tgt への位置ベクトル
        pos_st = pos_tgt(ti,:) - pos_src(si,:);
        R = norm(pos_st);
        G = exp(-j*k*R)/R; %音圧グリーン関数
        
        beta = j*k + 1/R;
        rdotn = dot(norms(ti,:), pos_st/R); %位置ベクトルと界面法線の内積
        dps_M(ti,si) = rdotn/(j*w*rho) * beta * G; %法線粒子速度のグリーン関数
    end
end
dps_A = dps_M\dps_V;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integral
p=zeros(size(X));
for si=1:num_src
    r_x = X-pos_src(si,1);
    r_y = 0-pos_src(si,2);
    r_z = Z-pos_src(si,3);    
    R = sqrt(r_x.^2 + r_y.^2 + r_z.^2);
    G = exp(-j*k*R)./R;
    p=p+dps_A(si)*G;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
figure(2); surf(x,z,abs(p));
shading interp;
view(2); axis tight;
xlabel('x'); ylabel('z');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

