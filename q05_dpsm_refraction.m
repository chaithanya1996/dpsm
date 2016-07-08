clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter
mm=1e-3; kHz=1e3; j=1i;
f=20*kHz; w=2*pi*f;

c1=340; rho1=1.293;  k1=w/c1; lamb1=c1/f;
c2=1500; rho2=998;  k2=w/c2; lamb2=c2/f;

a0=30*mm; %[m] 放射源
v0=1; %[m/s] disc velocity
rs = lamb1/10;

len=16*lamb1;
height= (33/8)*lamb1; %放射源高さ
tilt =  10;     %[deg] 放射源傾き

%observing plane
x=[-len/2+len/4:rs:len/2+len/4];
y=[-2*height:rs/2:height];
[X,Y] = meshgrid(x,y);
is_m1 = (Y >= 0);
is_m2 = (Y < 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%source placement
ang = tilt/180*pi;
snorm = [sin(ang) -cos(ang)  0];
scent = -height * [sin(ang) -cos(ang)  0];
[pos_s, norms_s] = def_line3(rs, 2*a0, scent, snorm); %0,0,0に放射する線音源
[pos_a, norms_a] = def_line3(rs, len, [len/4 0 0], [0 1 0]);  %空気側放射界面
[pos_w, norms_w] = def_line3(rs, len, [len/4 0 0], [0 -1 0]); %水側放射界面

pos_tgt=[pos_s; pos_a; pos_w];
norms  =[norms_s;norms_a;norms_w];
pos_src = pos_tgt - rs*norms ; %未知音源は放射面より -rsだけ下げた位置に配置する
figure(1);
plot(pos_tgt(:,1), pos_tgt(:,2), 'o', pos_src(:,1), pos_src(:,2), 'x'); axis equal;
xlabel('x'); ylabel('y');

[s_num, ~] = size(pos_s);
[a_num, ~] = size(pos_a);
[w_num, ~] = size(pos_w);
[num_src, ~] = size(pos_tgt);

kind = [1*ones(s_num,1); 2*ones(a_num,1); 3*ones(w_num,1)];  %伝播媒質種別
k    = [k1*ones(s_num+a_num,1); k2*ones(w_num,1)]; %媒質で波数が異なる
jwr   = j*w*[rho1*ones(s_num+a_num,1); rho2*ones(w_num,1)]; %媒質で密度が異なる

M_coef = [1 1 0; 1 1 -1; 0 0 0];
G_coef = [0 0 0; 0 0 0; 1 1 -1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine amplitude
dps_M = zeros(num_src, num_src);
dps_V = [v0*ones(s_num,1); zeros(a_num+w_num,1)];

for si=1:num_src %source loop
    for ti=1:num_src %target loop
        
        %src -> tgt への位置ベクトル
        pos_st = pos_tgt(ti,:) - pos_src(si,:);
        R = norm(pos_st);
        G = besselh(0, 2, k(si)*R);
        
        rdotn = dot(norms(ti,:), pos_st/R); %位置ベクトルと界面法線の内積
        M = rdotn/jwr(si) * k(si) * besselh(1, 2, k(si)*R); %法線粒子速度のグリーン関数
        
        tk = kind(ti); sk = kind(si);
        dps_M(ti,si) = M_coef(tk,sk)*M + G_coef(tk,sk)* G;
    end
end
figure(2);
surf(log(abs(dps_M))); view(2); axis ij; 
axis tight;
dps_A = dps_M\dps_V;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integral
p=zeros(size(X));
for si=1:num_src
    r_x = X-pos_src(si,1);
    r_y = Y-pos_src(si,2);
    r_z = 0-pos_src(si,3);
    R = sqrt(r_x.^2 + r_y.^2 + r_z.^2);
    G = besselh(0, 2, k(si)*R);
    
    rdotn = r_y./R;
    M = rdotn/jwr(si) * k(si) .* besselh(1, 2, k(si)*R); 
   
    sk = kind(si);
    prop = is_m1 .* ((sk==1)|(sk==2)) + is_m2 .* (sk==3);
    p=p+dps_A(si)*(prop.*G);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
figure(3); 
surf(x,y,real(p)); view(2);
shading interp;
axis equal;
xlabel('x'); ylabel('y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
