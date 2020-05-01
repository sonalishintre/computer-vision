%% Sonali Shintre
%%ID 7620
%Assignment 3 
%Computer vision (CSC I6716)
%% question i
% Generate data of a "virtual” 3D cube 

cnt = 1;
% plane : Xw = 1
for i=0.2:0.2:0.8
 for j=0.2:0.2:0.8
   Pw(cnt,:) = [1 i j];
   cnt = cnt + 1;
 end
end

% plane : Yw = 1
for i=0.2:0.2:0.8
 for j=0.2:0.2:0.8
   Pw(cnt,:) = [i 1 j];
   cnt = cnt + 1;
 end
end

% plane : Xw = 0
for i=0.2:0.2:0.8
 for j=0.2:0.2:0.8
   Pw(cnt,:) = [0 i j];
   cnt = cnt + 1;
 end
end

% plane : Yw = 0
for i=0.2:0.2:0.8
 for j=0.2:0.2:0.8
   Pw(cnt,:) = [i 0 j];
   cnt = cnt + 1;
 end
end

% plane : Zw = 1
for i=0.2:0.2:0.8
 for j=0.2:0.2:0.8
   Pw(cnt,:) = [i j 1];
   cnt = cnt + 1;
 end
end

% plane : Zw = 0
for i=0.2:0.2:0.8
 for j=0.2:0.2:0.8,
   Pw(cnt,:) = [i j 0];
   cnt = cnt + 1;
 end
end

N = cnt;
% so here Pw is a 1x1x1 m3 cube, 96 points are generated and each surface have 16 points.

%% Question ii
%Design a "virtual” camera with known intrinsic parameters 
% set the intrinsic parameters
f = 0.016;
Ox = 256;
Oy = 256;

Sx = 0.0088/512.0;
Sy = 0.0066/512.0;

Fx = f/Sx;
Fy = f/Sy;

% 4 meters away
T = [0 0 4]';  

gamma = 0.0*pi/180.0; 
Rr = [ [cos(gamma) -sin(gamma) 0];
       [sin(gamma) cos(gamma)  0];
       [  0          0         1]; ];

beta = 0.0*pi/180.0;

Rb = [ [cos(beta) 0 -sin(beta)];
       [0         1       0];
       [sin(beta) 0  cos(beta)]; ];

% tilt angle of 30 degree
alpha = 30.0*pi/180.0;  
Ra = [ [1      0                0];
       [0   cos(alpha)  -sin(alpha)];
       [0   sin(alpha)   cos(alpha)]; ];

R = Ra*Rb*Rr;

% Generate Image coordinates
% surface Xw = 1
cnt = 1;
for cnt = 1:1:16,
   Pc(cnt,:) = (R*Pw(cnt,:)' + T)';
   p(cnt,:)  = [(Ox - Fx*Pc(cnt,1)/Pc(cnt,3)) (Oy - Fy*Pc(cnt,2)/Pc(cnt,3))];
end
plot(p(:,1), p(:,2), 'r+');
axis([0 512 0 512]);
grid;
hold

% surface Yw = 1
for cnt = 17:1:32
   Pc(cnt,:) = (R*Pw(cnt,:)' + T)';
   p(cnt,:)  = [(Ox - Fx*Pc(cnt,1)/Pc(cnt,3)) (Oy - Fy*Pc(cnt,2)/Pc(cnt,3))];
end
plot(p(17:32,1), p(17:32,2), 'g+');
axis([0 512 0 512]);
%%plot3(Pc(:,1), Pc(:,2), Pc(:,3), '+');
grid;

% surface Xw = 0
for cnt = 33:1:48
   Pc(cnt,:) = (R*Pw(cnt,:)' + T)';
   p(cnt,:)  = [(Ox - Fx*Pc(cnt,1)/Pc(cnt,3)) (Oy - Fy*Pc(cnt,2)/Pc(cnt,3))];
end
plot(p(33:48,1), p(33:48,2), 'r+');
axis([0 512 0 512]);
grid;


% surface Yw = 1
for cnt = 49:1:64
   Pc(cnt,:) = (R*Pw(cnt,:)' + T)';
   p(cnt,:)  = [(Ox - Fx*Pc(cnt,1)/Pc(cnt,3)) (Oy - Fy*Pc(cnt,2)/Pc(cnt,3))];
end
plot(p(49:64,1), p(49:64,2), 'g+');
axis([0 512 0 512]);
%%plot3(Pc(:,1), Pc(:,2), Pc(:,3), '+');
grid;


% surface Zw = 1
for cnt = 65:1:80
   Pc(cnt,:) = (R*Pw(cnt,:)' + T)';
   p(cnt,:)  = [(Ox - Fx*Pc(cnt,1)/Pc(cnt,3)) (Oy - Fy*Pc(cnt,2)/Pc(cnt,3))];
end
plot(p(65:80,1), p(65:80,2), 'b+');
axis([0 512 0 512]);
grid;


% surface Zw = 0
for cnt = 81:1:96
   Pc(cnt,:) = (R*Pw(cnt,:)' + T)';
   p(cnt,:)  = [(Ox - Fx*Pc(cnt,1)/Pc(cnt,3)) (Oy - Fy*Pc(cnt,2)/Pc(cnt,3))];
end
plot(p(81:96,1), p(81:96,2), 'b+');
%%plot3(Pc(:,1), Pc(:,2), Pc(:,3), '+');
grid;
hold


%% Question 3
% Direction calibration method 

% first, we know the image center Ox = 256; Oy = 256;
p2(:,1)=p(:,1)-Ox;
p2(:,2)=p(:,2)-Oy;


for cnt = 1:1:96
    A(cnt,:) = [p2(cnt,1)*Pw(cnt,1) p2(cnt,1)*Pw(cnt,2) p2(cnt,1)*Pw(cnt,3) p2(cnt,1) -p2(cnt,2)*Pw(cnt,1) -p2(cnt,2)*Pw(cnt,2) -p2(cnt,2)*Pw(cnt,3) -p2(cnt,2)];
end

[U,D,V] = svd(A);
v=(V(:,8))';
r_ab = sqrt(v(1)^2+v(2)^2+v(3)^2);
ar=sqrt(v(5)^2+v(6)^2+v(7)^2)/r_ab;

%assume the sign s is positive
v2 = v/r_ab;
R2 = ([v2(1) v2(2) v2(3)])';
Ty = v2(4);
R1 =([v2(5)/ar v2(6)/ar v2(7)/ar])';
Tx1 = v2(7)/ar;
R3 = cross(R1,R2);
% test sign of s use the first point.
X = R1(1)*Pw(1,1) + R1(2)*Pw(1,2) + R1(3)*Pw(1,3) + Tx1;
% X = -1 and x = 1 so s is negative.
R_c = ([-R1 -R2 R3])';

% calculate Tz, Fx and Fy
for cnt = 1:1:96
    A2(cnt,:) = [p2(cnt,1) R_c(1,:)*(Pw(cnt,:))'+Tx1 -p2(cnt,1)*R_c(3,:)*(Pw(cnt,:))'];
    A3(cnt,:) = [p2(cnt,1) R_c(1,:)*(Pw(cnt,:))'+Tx1];
    b1(cnt) = -p2(cnt,1)*R_c(3,:)*(Pw(cnt,:))';
end
b2 = b1';

result = inv(A3'*A3)*A3'*b2;
Tz = result(1);
fx = result(2);
fy = fx/ar;

% so the parameters are all found.

%% question 3(iii)
% add some noises.add 0.1 mm random error to 3D points and 0.5 pixel random error to 2D points.
Pw = Pw + 0.0001*randn(96,3);
p_n= p + 0.5*randn(96,2);

p2(:,1)=p_n(:,1)-Ox;
p2(:,2)=p_n(:,2)-Oy;


for cnt = 1:1:96
    A(cnt,:) = [p2(cnt,1)*Pw(cnt,1) p2(cnt,1)*Pw(cnt,2) p2(cnt,1)*Pw(cnt,3) p2(cnt,1) -p2(cnt,2)*Pw(cnt,1) -p2(cnt,2)*Pw(cnt,2) -p2(cnt,2)*Pw(cnt,3) -p2(cnt,2)];
end

[U,D,V] = svd(A);
v=(V(:,8))';
r_ab = sqrt(v(1)^2+v(2)^2+v(3)^2);
ar_n=sqrt(v(5)^2+v(6)^2+v(7)^2)/r_ab;

%assume the sign s is positive
v2 = v/r_ab;
R2 = ([v2(1) v2(2) v2(3)])';
Ty_n = v2(4);
R1 =([v2(5)/ar_n v2(6)/ar_n v2(7)/ar_n])';
Tx = v2(7)/ar_n;
R3 = cross(R1,R2);
% test sign of s use the first point.
X = R1(1)*Pw(1,1) + R1(2)*Pw(1,2) + R1(3)*Pw(1,3) + Tx;
% X = -1 and x = 1 so s is negative.
R_c_n = ([-R1 -R2 R3])';

% calculate Tz, Fx and Fy

for cnt = 1:1:96
    A2(cnt,:) = [p2(cnt,1) R_c_n(1,:)*(Pw(cnt,:))'+Tx -p2(cnt,1)*R_c_n(3,:)*(Pw(cnt,:))'];
    A3(cnt,:) = [p2(cnt,1) R_c_n(1,:)*(Pw(cnt,:))'+Tx];
    b1(cnt) = -p2(cnt,1)*R_c_n(3,:)*(Pw(cnt,:))';
end
b2 = b1';


result = inv(A3'*A3)*A3'*b2;
Tz_n = result(1);
fx_n = result(2);
fy_n = fx/ar_n;

























