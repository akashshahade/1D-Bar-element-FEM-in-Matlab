%-------------------------------------------------------------------%
% Software by -  Akash Shahade
% Title - Program for FE Analysis of 1D bar element With axial Load.
%-------------------------------------------------------------------%
clear;
clc;
fprintf(' \nFEM ANALYSIS OF 1-D BAR WITH UNIFORM CROSS-SECTION \n\n');

%------------------------------------------------------%
% PRE PROCESSING
%------------------------------------------------------%
l = 1000;                     % Length of bar (mm)
area=1000;                    % Area of bar (mm^2)
E=2.1*10^5;                   % Young's Modulus (Mpa)
u = 0.3;                      % Poisson's Ratio
dof=1;                        % Degree of Freedom per node
nnode=3;                      % Number of Nodes
z=nnode*dof;                  % Size of Stiffness Matrix
k=zeros(z,z);                 % Zero K matrix
f=zeros(z,1);                 % Zero Force matrix
U=zeros(nnode,1);             % zero Displacement matrix
nel=2;                        % number of elements
le=l/nel;                     % Length of each element
 
fixed_node = 1;

kel=((area*E)/le)*[1 -1;-1 1]; % Elemental stiffness matrix


%------------------------------------------------------%
%Assembly of Global Matrix
%------------------------------------------------------%

for i =1:nel
 range = (i:i+1);
 k(range,range) = k(range,range) + kel;
end
fprintf('________________________________________________\n\n');
fprintf('Global Stiffness Matrix:\n');
disp(k);

%------------------------------------------------------%
% Applied load
%------------------------------------------------------%

p=100*10^3;
f(3,1)=p;
fprintf('________________________________________________\n\n');
fprintf('Global Load Vector:\n');
disp(f);
fprintf('________________________________________________\n\n');

%------------------------------------------------------%
% Applying Boundary Conditions
%------------------------------------------------------%

kg1=k(:,(z-fixed_node):end);
kg=kg1((z-fixed_node):end,:);
fg=f((z-fixed_node):end,:);

%------------------------------------------------------%
% SOLVE
%------------------------------------------------------%
U=kg\fg ;

%------------------------------------------------------%
% POST PROCESSING
%------------------------------------------------------%

Ut=[0;U];
e1=([-1 1]*([Ut(1) Ut(2)]).')/le;
S1=E*e1;
e2=([-1 1]*([Ut(2) Ut(3)]).')/le;
S2=E*e2;
fprintf('NODAL DISPLACEMENTS:\n');
Node=[1;2;3];
Displacement = [Ut(1,1);Ut(2,1);Ut(3,1)];
T=table(Node,Displacement);
disp(T);
fprintf('________________________________________________\n\n');
fprintf('\n STRESS IN ELEMENTS:\n')
Element=[1;2];
Stress=[S1;S2];
M=table(Element,Stress);
disp(M);
fprintf('________________________________________________\n');
fprintf('END OF PROGRAM.\n')
