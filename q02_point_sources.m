clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter
mm=1e-3; kHz=1e3; j=1i;
c=340; rho=1.293;
f=20*kHz; w=2*pi*f; k=w/c;

a0=100*mm; %[m] disc radius
v0=1; %[m/s] disc velocity

%observing plane
x=[-300:5:300]*mm;
z=[5:5:1000]*mm;
[Z,X] = meshgrid(z,x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disk division
rs = a0/40;
[pos,norms]=def_circ3(2*rs, a0, [0 0 0], [0 0 1]);
figure(1); plot3(pos(:,1), pos(:,2), pos(:,3), 'o');
[num_src, ~] = size(pos);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integral
p=zeros(size(X));
for si=1:num_src
    dS = pi*rs^2;
    r_x = X-pos(si,1);
    r_y = 0-pos(si,2);
    r_z = Z-pos(si,3);    
    R = sqrt(r_x.^2 + r_y.^2 + r_z.^2);

    A=j*w*rho/(2*pi)*v0*dS;
    p=p+A*exp(-j*k*R)./R;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
figure(2); surf(z,x,abs(p));
shading interp;
view(2); axis tight;
xlabel('z'); ylabel('x');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

