clc;
clear all;
close all;
%% %%%%%%%%%% Geometric parameters %%%%%%%%%%%%
thickness_count =0;


cases = [10 5;14 7;18 9;10 7.5;14 10.5;18 13.5;10 10;14 14;18 18];
do_b = 15; %bucket outer diameter
thickness = 0.03; %bucket thickness
l_b =15; %length of the bucket skirt
Z_lid = 0.001; % thickness of the lid of bucket
n = 32; %number of strips in a single circumference
n_l = 12; % number of rings in the bottom
n_q =16;
eccentricity = 31.5;
th1=0;
Deadload = 10000; %kN
thickness_count = thickness_count + 1; 
l = 0.0;
%% %%%%%%%%%% Soil parameters %%%%%%%%%%%%
% input_gamma = 19;%kN/m3 %input soil
gamma =10.27; %19-9.8;  %input_gamma-9.8; % input_gamma;% specific soil weigh %kN/m3 
S_gamma = 0.6;
e_max =0;
e_min =0;
D_r = 0.6 ; %input soil
m_c = 161.42*D_r^2+199.8*D_r+36.877; % 80-150 soft sand % 150-250 medium sand % 250-400 stiff sand  (rewrite based on dr) 
A_c = 0.3474*D_r^2+0.4222*D_r+0.328; 
fi_crit = 32.5;% input soil
m =3; % for finding the FI_preak %m is either 0 or 3
delta = 0.75*32.5;% 0.75*33; % for the field test 0.75*33     % centrifuge test 0.65~0.7 *33 % kim = 21 wang =20.8
%%
ro_b = do_b/2; %outer radious of the bucket
di_b =do_b-2*thickness; %bucket inner diameter
ri_b = di_b/2; %outer radious of the bucket
r_pile = di_b/2 + thickness/2; %from the center of buckt to the center of the thickness
r_b = abs(di_b-do_b)/4; %radius of the pile
penetration =l_b; %fully penetrated
circle_split = linspace(0,360,n+1);
R_split = linspace(0,r_pile,n_q);
%% Code initiation
QQ = zeros(1,10000);
kkk =0;
abc=0;
cc1=0;
% full penetration action
    abc=abc+1;
 cc1=cc1+1;
   if abs(l_b - penetration)<=0.01
       FP =1;
   else
        FP=0;
   end
count =0;

%%  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ dealload effect with no tilt angle $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  %%
dz=0;
SUMz=0;
while true
    if SUMz>Deadload
    dz=dz-0.0000125;
    else
     dz=dz+0.0000125;
    end
for dd=0
   q=[0 -dd];
phi = q(1);     % Rotation around x-axis
theta = q(2);   % Rotation around y-axis
psi = 0;     % Rotation around z-axis
j=0;
for k=0:l_b/n_l:l_b-(l_b/n_l)
for i=2:n+1 
    j=j+1;
    first_s_end(:,j) = [l*cos(th1)+ro_b*cosd(circle_split(i)) l*sin(th1)+ro_b*sind(circle_split(i)) -k-(l_b/n_l)/2]'; %  bucket points on the skirt  w. r. t. the body frame
    % notice: the stars are in the middle of layers.
end
end
% w. r. t.  fixed frame
A = zeros(size(first_s_end));
A(3,:) = -dz;
first_s_end_fixed = A + first_s_end;


j=0;
for i=2:n+1 
    j=j+1;
    Pile_first_s_end(:,j) = [l*cos(th1)+r_pile*cosd(circle_split(i)) l*sin(th1)+r_pile*sind(circle_split(i)) -l_b]'; % tip of the bucket points  w. r. t the body frame
end
B = zeros(size(Pile_first_s_end));
B(3,:) = -dz;
Pile_first_s_end_fixed = Pile_first_s_end + B;
j=0;
jj = 0;
co =0;
R_sec = R_split(2)-R_split(1);
for jj =2:size(R_split,2)
j=0;
% w.r.t body frame
for i=2:n+1 
    j=j+1;
r_Pile_first_s_end(:,j+co*n) = [l*cos(th1)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)) l*sin(th1)+(R_split(jj)-R_sec/2)*sind(circle_split(i)) -Z_lid]';
end
co = co+1;
end
j=0;
C = zeros(size(r_Pile_first_s_end));
C(3,:) = -dz;
r_Pile_first_s_end_fixed = r_Pile_first_s_end+C;
z= ones(n_l,n);
kk=1;
for p=1:n_l
    z(p,:) = z(p,:)*l_b/n_l*kk;
    kk=kk+1;
end
nn = size(z, 1);
mm = size(z, 2);
z_values = zeros(1, nn * mm); % z_value is the depth of each node located in the center of elements % bucket lid excluded

for i = 1:nn * mm
    row = floor((i - 1) / mm) + 1;
    col = mod(i - 1, mm) + 1;
    z_values(i) = z(row, col);
end
z_values = z_values-(l_b/(2*n_l)); %point the center of layers

%% activating the springs
active_z_values = zeros(size(z_values));
for i=1:n_l
if penetration <= i*l_b/n_l && penetration > (i-1)*l_b/n_l
        active_z_values(length(active_z_values)-i*n+1:end) = z_values(1:i*n);
end
end
FI_peak = fi_finder(gamma,active_z_values,D_r,fi_crit,m);
fi = FI_peak;
fi_end = fi(end);
displacement = (abs(first_s_end_fixed(3, :)) - abs(first_s_end(3, :))) * 1e3;  
%% (t-z curve) %%%%
[~,~,P_i_out] = coeffinder_deadload(first_s_end,first_s_end_fixed,active_z_values,l_b,n_l,fi,n,do_b,gamma); %KN
[tz_in,tz_out]=tzcurve(P_i_out,fi,gamma,do_b,di_b,active_z_values,displacement,l_b,n_l,n,fi_crit,delta,D_r,e_max,e_min);
KKz1_in = tz_in;
KKz1_out = tz_out;
KKz1 =tz_out+tz_in;
%% (q-z curve) %%%%
[qb1_1] = qzcurve2(fi_end,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb2_1] = qzcurve4(fi(1),Z_lid,r_Pile_first_s_end, r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1);
    
end
    abs(Deadload - SUMz)
 if abs(Deadload - SUMz)<Deadload/100
    break
 end

end 


%%  $$$$$$$$$$$$$$$$$$$$$$$$$   bucket rotation angle change step by step $$$$$$$$$$$$$$$$$$$$$$$$$ %%
Q = [0,0,0]';     %initial acenter of rotation
for dd = [0.25 0.5 0.75 1.0 1.5 2.0 3.0 3.5]*pi/180
q=[0 -dd];
phi = q(1);      % Rotation around x-axis
theta = q(2);   % Rotation around y-axis
psi = 0;           % Rotation around z-axis
%%%%%%%%%%%Suction arrangement angle and geometrica vectors%%%%%%%%%%%
tower_tip = [0,0,eccentricity]';
lid_center = [0,0,0]';

tower_tip_fixed = transformation(tower_tip,Q,[phi,theta,psi]);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(1)^2);
SUMz =0; 
count = count + 1;
ddd(count) = dd*180/pi;    
%%%%%%%%%%%
z= ones(n_l,n);
kk=1; 
for p=1:n_l
    z(p,:) = z(p,:)*l_b/n_l*kk;
    kk=kk+1;
end
nn = size(z, 1);
mm = size(z, 2);
z_values = zeros(1, nn * mm);
for i = 1:nn * mm
    row = floor((i - 1) / mm) + 1;
    col = mod(i - 1, mm) + 1;
    z_values(i) = z(row, col);
end
z_values = z_values-(l_b/(2*n_l));
active_z_values = zeros(size(z_values));
for i=1:n_l
if penetration <= i*l_b/n_l && penetration > (i-1)*l_b/n_l
        active_z_values(length(active_z_values)-i*n+1:end) = z_values(1:i*n);
end
end

ee=0;
ee_tar = eccentricity;
nnn=0;
cc=0;
mmm=0;
QQQ = zeros(1,10000);
ccc =0;
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  Outer While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %% 
while  abs((ee)-ee_tar)>0.001*ee_tar %while condition to find the z-axis of the center of rotation
    tower_tip_fixed = transformation(tower_tip,Q,[phi,theta,psi]);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
ccc = ccc+1;
 if ((ee)-ee_tar)>0.0
 mmm=mmm-0.005;
 else
 mmm=mmm+0.005;
 end
if cc==1
 Q = [0, 0, -(0.8*l_b)+mmm]';
else
 Q = [(0+nnn)*cosd(angle_indicator), (0+nnn)*sind(angle_indicator), -(0.8*l_b)+mmm]';
end
if rem(ccc,5)==0
fprintf('\rdegree: %0.3f  | ((ee)-ee_tar): %0.2f |(SUMz-Deadload): %0.3f | Qx: %0.3f| Qy: %0.2f | Qz: %0.3f', dd*180/pi, CON2,CON1,Q(1),Q(2),Q(3));  
end
QQQ(ccc+2)=Q(3);
if abs(abs(QQQ(ccc+2))-abs(QQQ(ccc)))==0
     break
end

%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  inner While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %% 
 while abs(SUMz-Deadload)>(Deadload/1500) 
CON1 = abs(SUMz-Deadload); %track the condition
tower_tip_fixed = transformation(tower_tip,Q,[phi,theta,psi]);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
cc=cc+1;
if SUMz-Deadload>0
nnn=nnn+0.005;
else
nnn=nnn-0.005;
end
Q = [(0+nnn)*cosd(angle_indicator), (0+nnn)*sind(angle_indicator), Q(3)]';% iterative increasing or decreasing the x-axis of the rotation center
fprintf('\rdegree: %0.3f  | ((ee)-ee_tar): %0.3f |(Deadload - SUMz): %0.3f | Qx: %0.3f| Qy: %0.2f | Qz: %0.3f', dd*180/pi,((ee)-ee_tar) ,CON1,Q(1),Q(2),Q(3));

QQ(cc+2)=Q(1);
if abs(abs(QQ(cc+2))-abs(QQ(cc)))==0
     break
 end

 [first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Pile_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,phi,theta,psi,Z_lid,l_b,l,ro_b,ri_b,r_pile);
 %% py and tz forces inner while
[KKx1,KKy1,~] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed_in(3, :)),l_b,n_l,fi,n,do_b,gamma,abs(Q(3)),A_c,m_c); %KN

displacement = (abs(first_s_end_fixed(3, :)) - abs(first_s_end(3, :))) * 1e3; %1e3 is to make it to mm
displacement_in = (abs(first_s_end_fixed_in(3, :)) - abs(first_s_end_in(3, :))) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi,n,do_b,gamma,Q(3),A_c,m_c); %KN
[~,tz_out]=tzcurve(P_i_out,fi,gamma,do_b,di_b,abs(first_s_end_fixed(3, :)),displacement,l_b,n_l,n,fi_crit,delta,D_r,e_max,e_min);
[tz_in,~]=tzcurve(P_i_out,fi,gamma,do_b,di_b,abs(first_s_end_fixed_in(3, :)),displacement_in,l_b,n_l,n,fi_crit,delta,D_r,e_max,e_min);

KKz1_in = tz_in;
KKz1_out = tz_out;
KKz1 =KKz1_in+KKz1_out;
%% qz forces inner while
[qb1_1] = qzcurve2(fi_end,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb2_1] = qzcurve4(fi_end,Z_lid,r_Pile_first_s_end, r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1);

 end
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  End of inner While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %% 
j=0;
fi=FI_peak;
 [first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Piles_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,phi,theta,psi,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[KKx1,KKy1,~] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi,n,do_b,gamma,abs(Q(3)),A_c,m_c); %KN
% KKX1(count,:) = KKx1;
% KKY1(count,:) = KKy1;
SUMx = sum(KKx1);
SUMy = sum(KKy1);
%% py and tz forces outer while
displacement = (abs(first_s_end_fixed(3, :)) - abs(first_s_end(3, :))) * 1e3; %1e3 is to make it to mm
displacement_in = (abs(first_s_end_fixed_in(3, :)) - abs(first_s_end_in(3, :))) * 1e3; %1e3 is to make it to mm
[~,~,P_i_out] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi,n,do_b,gamma,Q(3),A_c,m_c); %KN
[~,tz_out]=tzcurve(P_i_out,fi,gamma,do_b,di_b,abs(first_s_end_fixed(3, :)),displacement,l_b,n_l,n,fi_crit,delta,D_r,e_max,e_min);
[tz_in,~]=tzcurve(P_i_out,fi,gamma,do_b,di_b,abs(first_s_end_fixed_in(3, :)),displacement_in,l_b,n_l,n,fi_crit,delta,D_r,e_max,e_min);
KKz1_in = tz_in;
KKz1_out = tz_out;
KKz1 =tz_out+tz_in;
% KKZ1(count,:) = tz_out+tz_in;
%% qz forces outer while
[qb1_1] = qzcurve2(fi_end,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb2_1] = qzcurve4(fi_end,Z_lid,r_Pile_first_s_end, r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1);
%% py moments outer while
forces_py = [KKx1;KKy1;0*KKz1]';
positions = first_s_end';
moments_py = cross(positions, forces_py);
total_moment_py = sum(moments_py);
total_moment_py = sqrt(total_moment_py(1)^2+total_moment_py(2)^2+total_moment_py(3)^2);
for i=1:n_l
    pc(i)= n*i;
end

%% tz moments outer while
forces_tz_out = [0*KKx1;0*KKy1;tz_out]';
forces_tz_in = [0*KKx1;0*KKy1;tz_in]';

positions_out = first_s_end';
positions_in = first_s_end_in';

moments_tz_out = cross(positions_out, forces_tz_out);
moments_tz_in = cross(positions_in, forces_tz_in);

total_moment_tz = sum(moments_tz_in)+sum(moments_tz_out);
total_moment_tz = sqrt(total_moment_tz(1)^2+total_moment_tz(2)^2+total_moment_tz(3)^2);


%% qz lid moments outer while

positions_pile = Pile_first_s_end;

XX = zeros(size(qb1_1));
YY = zeros(size(qb1_1));
force_vectors = [XX;YY;qb1_1]';
moments_qz_1 = cross(positions_pile', force_vectors);
total_moment_qz_1 = sum(moments_qz_1);
total_moment_qz_1 = sqrt(total_moment_qz_1(1)^2+total_moment_qz_1(2)^2+total_moment_qz_1(3)^2);

%% qz tip moments outer while
r_positions_pile = r_Pile_first_s_end;

XXX = zeros(size(qb2_1));
YYY = zeros(size(qb2_1));
force_vectors_2 = [XXX;YYY;qb2_1]';
if FP==1
moments_qz_2 = cross(r_positions_pile', force_vectors_2);
else
moments_qz_2 = cross(r_positions_pile', force_vectors_2);
end
total_moment_qz_2 = sum(moments_qz_2);
total_moment_qz_2 = sqrt(total_moment_qz_2(1)^2+total_moment_qz_2(2)^2+total_moment_qz_2(3)^2);
%% outer While condition calculation
MOMENT = (total_moment_qz_2+total_moment_qz_1+total_moment_tz+total_moment_py);
F_Hx = sum(KKx1);
F_Hy = sum(KKy1);
F_H  = sqrt (F_Hx^2+F_Hy^2);
ee = abs((MOMENT)/F_H);

CON2=abs((ee)-ee_tar); 
end 
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  End of outer While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %% 
tower_tip_fixed = transformation(tower_tip,Q,[phi,theta,psi]);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(1)^2);
j=0;
[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Piles_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,phi,theta,psi,Z_lid,l_b,l,ro_b,ri_b,r_pile);

%% qz  moments 
[KKx1,KKy1,~] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi,n,do_b,gamma,abs(Q(3)),A_c,m_c); %KN
% KKX1(count,:) = KKx1;
% KKY1(count,:) = KKy1;
SUMx = sum(KKx1);
SUMy = sum(KKy1);

displacement = (abs(first_s_end_fixed(3, :)) - abs(first_s_end(3, :))) * 1e3; %1e3 is to make it to mm
displacement_in = (abs(first_s_end_fixed_in(3, :)) - abs(first_s_end_in(3, :))) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out] = coeffinder(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi,n,do_b,gamma,Q(3),A_c,m_c); %KN
[~,tz_out]=tzcurve(P_i_out,fi,gamma,do_b,di_b,abs(first_s_end_fixed(3, :)),displacement,l_b,n_l,n,fi_crit,delta,D_r,e_max,e_min);
[tz_in,~]=tzcurve(P_i_out,fi,gamma,do_b,di_b,abs(first_s_end_fixed_in(3, :)),displacement_in,l_b,n_l,n,fi_crit,delta,D_r,e_max,e_min);
KKz1_in = tz_in;
KKz1_out = tz_out;
KKz1 =tz_out+tz_in;
KKZ1(count,:) = tz_out+tz_in;
%% qz  moments 
[qb1_1] = qzcurve2(fi_end,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb2_1] = qzcurve4(fi_end,Z_lid,r_Pile_first_s_end, r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1);
%% py  moments 
forces_py = [KKx1;KKy1;0*KKz1]';
positions = first_s_end';
moments_py = cross(positions, forces_py);
total_moment_py = sum(moments_py);
total_moment_py = sqrt(total_moment_py(1)^2+total_moment_py(2)^2+total_moment_py(3)^2);

for i=1:n_l
    pc(i)= n*i;
end

%% tz  moments 
forces_tz_out = [0*KKx1;0*KKy1;tz_out]';
forces_tz_in = [0*KKx1;0*KKy1;tz_in]';
positions_out = first_s_end';
positions_in = first_s_end_in';
moments_tz_out = cross(positions_out, forces_tz_out);
moments_tz_in = cross(positions_in, forces_tz_in);
total_moment_tz = sum(moments_tz_in)+sum(moments_tz_out);
total_moment_tz = sqrt(total_moment_tz(1)^2+total_moment_tz(2)^2+total_moment_tz(3)^2);


%% qz tip moments 
positions_pile = Pile_first_s_end;
XX = zeros(size(qb1_1));
YY = zeros(size(qb1_1));
force_vectors = [XX;YY;qb1_1]';
moments_qz_1 = cross(positions_pile', force_vectors);
total_moment_qz_1 = sum(moments_qz_1);
total_moment_qz_1 = sqrt(total_moment_qz_1(1)^2+total_moment_qz_1(2)^2+total_moment_qz_1(3)^2);

%% qz lid moments 
r_positions_pile = r_Pile_first_s_end;
XXX = zeros(size(qb2_1));
YYY = zeros(size(qb2_1));
force_vectors_2 = [XXX;YYY;qb2_1]';
if FP==1
moments_qz_2 = cross(r_positions_pile', force_vectors_2);
else
moments_qz_2 = cross(r_positions_pile', force_vectors_2);
end
total_moment_qz_2 = sum(moments_qz_2);
total_moment_qz_2 = sqrt(total_moment_qz_2(1)^2+total_moment_qz_2(2)^2+total_moment_qz_2(3)^2);
%% total moment calculation
MOMENT(:,count) = (total_moment_qz_2+total_moment_qz_1+total_moment_tz+total_moment_py); % moment arount center of lid
F_Hx = sum(KKx1);
F_Hy = sum(KKy1);
F_HH(:,count)  = sqrt (F_Hx^2+F_Hy^2)/1000;
M_1(:,count)  = F_HH(thickness_count,count)  *eccentricity;
M_2(:,count)  = MOMENT(thickness_count,count)/1000;%+Deadload*Q(1);
lid_center_fixed(:,count) = transformation(lid_center,Q,[phi,theta,psi]);
del(:,count) = lid_center_fixed(1,count);
figure(46)
plot(del(thickness_count,count) ,F_HH(thickness_count,count),'*b','LineWidth',5)
hold on
grid on
colors2 = [
    0.111, 0.111, 0.111;    % Dark Grey
    0.8500, 0.3250, 0.0980; % Orange
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330; % Light Blue
    0.9290, 0.6940, 0.1250; % Yellow
    0.6350, 0.0780, 0.1840  % Dark Red
];
figure(44)


figure(4)
plot(dd*180/pi,M_1(count),'*r','LineWidth',5)
hold on
plot(dd*180/pi,M_2(count),'*b','LineWidth',5)

grid on
hold on
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',22)
xlabel('Rotation (degree)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Moment Load (MNm)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
figure(107)
title('CENTER OF ROTATION')
scatter3(Q(1),Q(2),Q(3),'o','LineWidth',8)
xlabel('X', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Y', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
Zlabel('Z', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')

hold on
% 
py_force(:,count) = sum(reshape(KKx1,[n,n_l])',2);
tz_forcee = reshape(KKz1_out,[n,n_l])'
tz_force_left(:,count)  = sum(tz_forcee(:,[4:9]),2);
tz_force_right(:,count) = sum(tz_forcee(:,[1:2]),2)+sum(tz_forcee(:,[10:12]),2);

inner_right(:,count) = sum(KKz1_in(pc)+KKz1_in(pc-11)+KKz1_in(pc-10)+KKz1_in(pc-1)+KKz1_in(pc-2))/gamma/do_b^3;
inner_left(:,count) = sum(KKz1_in(pc-8)+KKz1_in(pc-7)+KKz1_in(pc-6)+KKz1_in(pc-5)+KKz1_in(pc-4))/gamma/do_b^3;
outer_left(:,count) = sum(KKz1_out(pc-8)+KKz1_out(pc-7)+KKz1_out(pc-6)+KKz1_out(pc-5)+KKz1_out(pc-4))/gamma/do_b^3;
outer_right(:,count) = sum(KKz1_out(pc)+KKz1_out(pc-11)+KKz1_out(pc-10)+KKz1_out(pc-1)+KKz1_out(pc-2))/gamma/do_b^3;
% figure(77) 
% plot(dd*180/pi,inner_right(count), '>b', 'LineWidth', 6); %inner right
% hold on
% plot(dd*180/pi,inner_left(thickness_count,count), '>r', 'LineWidth', 6); %inner left
% hold on
% %% 
% plot(dd*180/pi,outer_left(thickness_count,count), 'or', 'LineWidth', 6); %outer left
% hold on
% %%
% plot(dd*180/pi,outer_right(thickness_count,count), 'ob', 'LineWidth', 6); %outer right
% hold on
% set(gca, 'XScale', 'log');  % Set the x-axis to logarithmic scale
% hold on;
% grid on;
% xlabel('Angle (degrees)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% ylabel('Skirt wall friction (F_t_z/\gamma^\prime/D^3) ','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% % set(gca,'TickLabelInterpreter','latex');
% % set(gca,'fontweight','bold','fontsize',22)
% % set(gca, 'YAxisLocation', 'left');
% legend({'Inner right side friction','Inner left side friction','Outer left side friction','Outer right side friction'},'fontsize',16,'interpreter','latex','location','northeast')

%  CENTE_OF_ROTATION(3*thickness_count-2:3*thickness_count,count) = [Q(1)/do_b,Q(2)/do_b,Q(3)/l_b]';
% 
% Share_py(thickness_count,count) = total_moment_py/MOMENT(thickness_count,count)*100% figure(3)
% Share_tz(thickness_count,count) = total_moment_tz/MOMENT(thickness_count,count)*100% figure(3)
% Share_qz_tip(thickness_count,count) = (total_moment_qz_1)/MOMENT(thickness_count,count)*100% figure(3)
% Share_qz_lid(thickness_count,count) = (total_moment_qz_2)/MOMENT(thickness_count,count)*100% figure(3)

%  plot(dd*180/pi,total_moment_tz/total_moment(count)*100,'-*b','LineWidth',5)
%   plot(dd*180/pi,total_moment_py/total_moment(count)*100,'-*r','LineWidth',5)
%     plot(dd*180/pi,2*(total_moment_qz_1+total_moment_qz_2)/total_moment(count)*100,'-*k','LineWidth',5)
% %       plot(dd*180/pi,total_moments_qz_2(2)/total_moment(count)*100,'-*m','LineWidth',5)   
 hold on


end

