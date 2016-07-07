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
da=a0/20;
a=[da:da:a0];
dtheta=pi/36;
theta=[0:dtheta:2*pi-dtheta];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integral
p=zeros(size(X));
for ai=1:length(a)
    for ti=1:length(theta)
        dS = a(ai)*da*dtheta;
        a_x = a(ai)*cos(theta(ti));
        a_y = a(ai)*sin(theta(ti));
        R = sqrt((X-a_x).^2 + (a_y).^2 + Z.^2);
        A=j*w*rho/(2*pi)*v0*dS;
        p=p+A*exp(-j*k*R)./R;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
surf(z,x,abs(p));
shading interp;
view(2); axis tight;
xlabel('z'); ylabel('x');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
