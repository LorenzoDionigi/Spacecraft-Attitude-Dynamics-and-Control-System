clear;
%e=0.05;
h_p=800;
h_a=800;
i=0*pi/180;
theta0=0;
P_E=670;


I_x=0.093296;
I_y=0.161261916;
I_z=0.160253916;
J=[I_x 0 0; 0 I_y 0; 0 0 I_z];
mu=astroConstants(13);
h=400;
r=astroConstants(23);
Rp=h_p+r;
Ra=h_a+r;
e=(Ra-Rp)/(Ra+Rp);
a=(Rp+Ra)/2;
n=sqrt(mu/a^3);
alpha1=0.9;
alpha2=0.1;
L=0.3/2;
d=0.05/2;
kState_obs=0.9;
alphaState_obs=-0.9;
gamma=30*pi/180;
w=[1 1 1 1]';
w0_x=-0.008328;
w0_y=0.004962;
w0_z=-0.0004756;
T=(2*pi)*sqrt(a^3/mu);
om=0;
k1SL=8;
k2SL=0.09;
k1=12;
k2=0.09;
kp=0.09;
k=5;
kd=12;

I=I_x+I_y+I_z;
kpx=kp*(I_x/I);
kpy=kp*(I_y/I);
kpz=kp*(I_z/I);
kdx=kd*(I_x/I);
kdy=kd*(I_y/I);
kdz=kd*(I_z/I);
RW=0.002;
RWmin=0;
Thruster_Power=0.03;
Tmin=0.03/3;
MIB=5*1e-3;
MMS=0.037;




mu_s=astroConstants(4);
a_s=149600000;
e_s=0.0167086;
i_s=1.578690*180/pi;
T_s=(2*pi)*sqrt(a_s^3/mu_s);
n_s=sqrt(mu_s/a_s^3);
eps=23.45*pi/180;
err=10*pi/180;
ax=rand(1,1)*2*pi;
ay=rand(1,1)*2*pi;
az=rand(1,1)*2*pi;
%A0=[0 -az ay; az 0 -ax; -ay ax 0  ];
%A0=A*(3/2)-(A*A'*A)/2;
%A0=[0 1 0; -1 0 0; 0 0 1];
%A0=eye(3);
Az=[cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
Ax=[1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
Ay=[cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
%A0=Ax*Ay*Az;
A0=[-0.2505 -0.7848 0.5668; 0.7978 0.1644 0.5801; -0.5485 0.5975 0.585];
%A0=rand(3);
%A0=[-0.1085 0.9922 0.06211; -0.8902 -0.1248 0.4382; 0.4426 -0.007732 0.8967];
Angle0=[pi/2 0 0];
q4=0.5*(1+A0(1,1)+A0(2,2)+A0(3,3))^(1/2);
q1=(1/(4*q4))*(A0(2,3)-A0(3,2));
q2=(1/(4*q4))*(A0(3,1)-A0(1,3));
q3=(1/(4*q4))*(A0(1,2)-A0(2,1));
q0=[q1 q2 q3 q4];
t=2*T;
B=[cos(n*t) sin(n*t) 0; -sin(n*t) cos(n*t) 0; 0 0 1];
c=299792458;
F_e=1358;
F_earth=670;
P=F_e/c;
P_E=F_earth/c;
ps_sc=0.5;
pd_sc=0.1;
ps_sp=0.1;
pd_sp=0.1;
c_m = 8.215*(1e-3)*[1;0;0];

L2 = 0.2; %Body Length [m]
L3 = 0.1; %Body width [m]
L1 = 0.3; %Body height [m]
A1 = L1*L2;
A2 = A1;
A3 = L1*L3;
A4 = A3;
A5 = L2*L3;
A6 = A5;
A7=A1;
A8=A1;
A9=A1;
A10=A1;

r1 = [0;0;-L3/2] + c_m;
r2 = [0; 0 ;L3/2] + c_m;
r3 = [0; L2/2;0] + c_m;
r4 = [0; -L2/2;0] + c_m;
r5 =  [-L1/2; 0 ;0] + c_m;
r6 = [L1/2; 0;0]+ c_m;
r7=[-L1/2; 0 ; L3/2 + L1/2] + c_m;
r8=r7;
r9=[-L1/2; 0 ; -L3/2 - L1/2];
r10=r9;


length_SP = 0.3;
height_SP = 0.2;
A_vect = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10]';


g10=-29615;
g11=-1728;
h11=5186;
an=(11.5*pi)/180; %inclinazione asse terreste


%MAGNETIC FIELD
n_B =13; %magnetic field order
% tabulate coefficient
coeff_str = importdata('MAGNETIC_COEFFICIENTS.txt');
coeff = coeff_str.data;
gh = coeff_str.textdata;
gh = gh(5:end,1);

gh_vect = coeff(:,27);
h_vect = [];
g_vect = [];

for i = 1:length(gh_vect)
    if  isequal(gh{i},'g')
        g_vect = [g_vect; gh_vect(i)];
        if  isequal(gh{i+1}, gh{i})
            h_vect = [h_vect; 0];
        end
    else
        h_vect = [h_vect; gh_vect(i)];
    end
end
h_vect = h_vect'.*1e-9;
g_vect =  g_vect'.*1e-9;

% P_MAT_0 = ones(n_b+1,n_b+1);
% dP_MAT_0 = zeros(n_b+1,n_b+1);
% p_vect00 = 1;
% K_nm = [];
% 
% for j = 1 : n_b
%     for k = 0 : j
%         if j == 1 
%             K_nm = [K_nm; 0];
%         else
%             K_nm = [K_nm; (((j-1)^2) - k^2)/((2*j-1)*(2*j-3))];
%         end
%     end
% end
H0test=sqrt(g10^2+g11^2+h11^2);
thetam=acos(g10/H0test);
phim=atan(h11/g11);

m=[0.01 0.05 0.01];
w_E=(15.04*pi/180)/3600; 