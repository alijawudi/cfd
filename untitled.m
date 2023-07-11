clear all
close all
clc
%% velocity relaxation at coeficients
clc; clear; close all;
tic
L = 1; % Length of Cavity
nH = 63; % Horizontal Number of Grids
nV = nH; % Vertical Number of 
dx = L/nH;
dy = L/nV;
Re = 100;
%% Under Relaxation Factors
alpha_p = 0.3;
alpha_u = 0.7;
alpha_v = alpha_u;
R_Continuty = 10^-6;
Beta1 = 1.85;
%% Nodes locations in y Direction for u
y = zeros(1,nV+2);
y(2:end-1) = L/nV/2:L/nV:L-L/nV/2;
y(end) = L;
u0 = zeros(nV+2,nH+1);
u0(1,:) = 1; 
v0 = zeros(nV+1,nH+2);
p0 = zeros(nV+2,nH+2);
D_x = 1/Re/dx*dy; 
D_y = 1/Re/dy*dx;
% u = [0 0 0;0 2 0;0 1 0];
% v = [0 0 0 0; 0 2 1 0;0 0 0 0];
% p = [1 2 3 4];
% d_U = [1;2];
% d_V = [2 2];
k = 0; %% Number of iteration
% Initial Guess
un = u0;
vn = v0;
pn = p0;
pp = p0;
bp = p0;
du = zeros(nV+2,nH+1);
dv = zeros(nV+1,nH+2);
R_TDMA = 10^-6; %% TDMA Convergance Criteria
Residual_Continuty = 1;
R0 = 0;
Ap_l =zeros(nH+2,6);
Ap = zeros(nH+2,6,nV+2);
kr = 0;
params1 = [dx dy D_x D_y];
niteration = 40000;
iteration = zeros(1 , niteration-6);
residualS = zeros(1 , niteration-6);
% Simple Loop
while true
k = k+1;
if k>6
iteration(k-6) = k;
residualS(k) = Residual_Continuty;
end
%% Step 1: Solve Discritized Momentum Equations
%% Solve Discritized Momentum Equation in X Direction
R_TDMA = 10^-3; %% TDMA Convergance Criteria
% if Residual_Continuty <0.16 && Residual_Continuty >0.15
% alpha_u = 0.6;
% alpha_p = 0.3;
% aloha_v =0.6;
% end
u = un;
v = vn;
p = pn;
% j_u = 0;
%% Solve Discritized Momentum Equation in X Direction
for j = 2:nH
 for i = 2:nV+1
 % F: Connvective Flux Per unit mass
 Fe = dy*(u(i,j)+u(i,j+1))/2;
 Fw = dy*(u(i,j)+u(i,j-1))/2;
 Fn = dx*(v(i-1,j)+v(i-1,j+1))/2;
 Fs = dx*(v(i,j)+v(i,j+1))/2;
 % D: Diffusive Conductance
 if i==2
 Dn = 2*D_y;
 Ds = D_y;
 elseif i==nV+1
Ds = 2*D_y;
 Dn = D_y;
 else
 Dn = D_y;
 Ds = D_y;
 end
 De = D_x;
 Dw = D_x;
 % Coefficients of Discritized Equations
 ae = De + max(-Fe,0);
 aw = Dw + max(Fw,0);
 an = Dn + max(-Fn,0);
 as = Ds + max(Fs,0);
 ap = ae+aw+an+as+ ((Fe-Fw) + (Fn-Fs));
 du(i,j) = alpha_u*dx/ap;
 % caculate of u Velocity
 un(i,j) = alpha_u*(1/ap*(ae*u(i,j+1) + aw*u(i,j-1) + an*u(i-1,j) + as*u(i+1,j)+(p(i,j)-p(i,j+1))*dy))+(1-alpha_u)*u(i,j);
 end
end
%% Solve Discritized Momentum Equation in Y Direction
for j = 2:nH+1
 for i = 2:nV
 
 % F: Connvective Flux Per unit mass
 Fe = dy*(u(i,j)+u(i+1,j))/2;
 Fw = dy*(u(i,j-1)+u(i+1,j-1))/2;
 Fn = dx*(v(i,j)+v(i-1,j))/2;
 Fs = dx*(v(i,j)+v(i+1,j))/2;
 
 % D: Diffusive Conductance
 Dn = D_y;
 Ds = D_y;
 if j==2
 De = D_x;
 Dw = 2*D_x;
 elseif j== nH+1
 De = 2*D_x;
 Dw = D_x;
 else
 Dw = D_x;
 De = D_x;
 end
% Coefficients of Discritized Equations
 ae = De + max(-Fe,0);
 aw = Dw + max(Fw,0);
 an = Dn + max(-Fn,0);
 as = Ds + max(Fs,0);
 ap = ae+aw+an+as+ ((Fe-Fw) + (Fn-Fs));
 dv(i,j) = alpha_v*dy/ap;
 % caculate of v Velocity
 vn(i,j) = (alpha_v)*(1/ap*(ae*v(i,j+1) + aw*v(i,j-1) +an*v(i-1,j) + as*v(i+1,j)+ (p(i+1,j)-p(i,j))*dx))+(1-alpha_v)*v(i,j);
 end
end
%% Pressure Correction Equations
pp = p0;
% Pprime Discritized Equations Coefficient
for i=2:nV+1
 for j=2:nH+1
 ae = dy*du(i,j);
 aw = dy*du(i,j-1);
 an = dx*dv(i-1,j);
 as = dx*dv(i,j);
 ap = ae+aw+an+as;
%  ap = ap-as*(Beta1-1);
 bprime = dy*(un(i,j-1) - un(i,j)) + dx*(vn(i,j)-vn(i-1,j));
 Ap_l(j,:) = [ap ae aw an as bprime];
 bp(i,j) = bprime;
 end
 Ap(:,:,i) = Ap_l;
end
% Solving pprime with TDMA Line by Line Methode
pp0 = pp;
for i=2:nV+1
 A = Ap(:,:,i);
 A = A(2:end-1,:);
 A1 = A(:,1:3);
 An = A(:,4:5);
 su = A(:,6);
 
 su = su+An(:,1).*(pp(i-1,2:end-1))'+An(:,2).*(pp0(i+1,2:end-1)-(Beta1-1)*pp0(i,2:end-1))';
 A = [A1 su];
 c = A;
 c(:,1) = A(:,3);
c(:,2) = A(:,1);
 c(:,3) = A(:,2);
% TDMA %%%%%%%%%%%%%%%%%%%%
 A = zeros(1,nH+1);
 Cp = zeros(1,nH+1); 
 for j=2:nH+1
 alphaj = c(j-1,3);
 Dj = c(j-1,2);
 Cj = c(j-1,4);
 betaj = c(j-1,1);
 A(j) = alphaj/(Dj-betaj*A(j-1));
 Cp(j) = (betaj*Cp(j-1)+Cj)/(Dj-betaj*A(j-1));
 end
 A(1) = [];
 Cp(1) = [];
 phi = zeros(1,nH+1);
 for j=nH:-1:1
 phi(j) = A(j)*phi(j+1)+Cp(j);
 end
 phi(end) = [];
 pp(i,2:end-1) = phi;
 
 
end
pn(2:end-1,2:end-1) = p(2:end-1,2:end-1) + alpha_p*pp(2:end-1,2:end-1);
%% u Velocity Update
for j=2:nH
 for i=2:nV+1
 un(i,j) = un(i,j)+du(i,j)*(pp(i,j)-pp(i,j+1));
% un(i,j) = alpha_u*un(i,j)+(1-alpha_u)*u(i,j);
 end
end
%% v Velocity Update
for j=2:nH+1
 for i=2:nV
 vn(i,j) = vn(i,j)+dv(i,j)*(pp(i+1,j)-pp(i,j));
% vn(i,j) = alpha_v*vn(i,j)+(1-alpha_v)*v(i,j);
 end
end
%% Residuals
Residual_P = max(max(abs(pn-p)));
Residual_U = max(max(abs(un-u)))/max(max(u));
Reseidual_V = max(max(abs(vn-v)));
if k<6R0 = R0+sum(sum(abs(bp)));
%  if k==5
%  R0 = R0/5;
%  end
elseif k>5
 Residual_Continuty = sum(sum(abs(bp)))/R0; 
end
disp(['R_Continuty : ',num2str(Residual_Continuty)]);
if Residual_Continuty >1000
 break
end
%% Convergence Condition
if Residual_Continuty<R_Continuty
 break
end
%% Dispalying Residuals
% disp('iteration : ',num2str(k),'Continuty : ',num2str(Residual_Continuty),...
% 'u: ',num2str(Residual_U),'v: ',num2str(Reseidual_V),...
% 'p: ',num2str(Residual_P))
end
ux = u(:,64);
ux_mid = flipud(ux);
plot(ux_mid,y)
n = length(iteration)
semilogy(iteration,residualS(1:n))
x = 0:L/nV:1;
[X_u,Y_u] = meshgrid(x,y);
u_contour = flipud(u);
figure(1)
contourf(X_u,Y_u,u_contour,20)
title('u velocity contour')
[X_v,Y_v] = meshgrid(y,x);
v_contour = flipud(v);
figure(2)
contourf(X_v,Y_v,v_contour,20)
title('v velocity contour')
[X_p,Y_p] = meshgrid(y,y);
p_contour = flipud(p);
figure(3)
contourf(X_p,Y_p,p_contour,25)
title('pressure contour')
[X_s,Y_s] = meshgrid(x,y);
figure(4)
streamline(X_s,Y_s,u,v)
title('stream line contour')
xlabel('x');
ylabel('y');
colorbar
streamslice(X_u(1:end-1,:),Y_u(1:end-1,:),u_contour(1:end-1,:),v_contour(:,1:end-1),5)
ali1=[1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0
]
ali2=[1.00000 0.48223 0.46120 0.45992 0.46036 0.33556 0.20087 0.08183 -0.03039 -0.07404 -0.22855 -0.33050 -0.40435  -0.43643 -0.42901 0 ]
for i = 1:nH
    for j = 1:nV
        u_final(i,j) = 0.5*(un(i,j) + un(i+1,j));
        v_final(i,j) = 0.5*(vn(i,j) + vn(i,j+1));
    end
end
yCenter = dx/2:dx:L-dx/2;
figure
subplot(2,2,1)
plot( u_final(:,(nH+1)/2), 1-yCenter, 'r','LineWidth', 2 )
subplot(2,2,2)
plot(iteration,Residual_Continuty , 'ko','MarkerFaceColor', 'k')
xCenter=dx:dx:L-dx/2;
subplot(2,2,3)
plot(xCenter, v_final((nV+1)/2,2:nV), 'r','LineWidth', 2 )
yCenter = dx/2:dx:L-dx/2;
subplot(2,2,4)
plot( u_final(:,(nH+1)/2), 1-yCenter, 'r','LineWidth', 2 )