%length = mm , F = N, time = s
% E = 0.2 MPa, radius = 24.2/2, a = EA, c = 0, f = parabolic x(L-x) loading, 
% L = 300 mm 

% This is the code to obtain FEA solution to problem

%-d/dx(a du/dx) + cu - f = 0, 0<x<L
%with BCs
%u(0) = 0; u(1) = 0;


%Input parameters 
a = 1011; c = -1; % Material parameters
L = 1; %Geometric parameters
N = 3; p = 2; he = L/N;  %Finite element method parameters
%f = -x^2; % Distributed force function

n = (p-1)*N+1; %n --> Total number of global nodes in the problem
QG = zeros(n,1); fG = zeros(n,1); kG = zeros(n,n);

EssentialBCNodes = [1]; %Nodes on which essential BCs are prescribed
EssentialBCs = [0]; %Values of essential BCs

NaturalBCNodes = [n]; %Nodes on which non-trivial natural BCs are prescribed
NaturalBCs = [-50000]; %Values of non-trivial natural BCs

%Assembly of element matrices

for e=1:N 
  xe = he*(e-1)+he/2; ae=a; ce=c; fe = 0; %Determing the Center of the element e to determine fe on that elementt
  [feM, keM] = EM_321(ae, ce, fe, he, p);
  for i = 1:p
    I = (e-1)*(p-1)+i;
    fG(I) = fG(I) + feM(i);
    for j = 1:p
      J = (e-1)*(p-1)+j;
      kG(I,J) = kG(I,J) + keM(i,j);
    end
  end
end

UG = zeros(n,1); s = length(EssentialBCNodes);
for i = 1:s    %Assigning the BC value in the solution matrix initialized 
  UG(EssentialBCNodes(i)) = EssentialBCs(i);
end

fC = fG - kG*UG;
fC(EssentialBCNodes) = []; %Removing the rows in fC to create condesed form 

UCNodes = 1:n; UCNodes(EssentialBCNodes)=[]; 

for i=1:length(NaturalBCNodes)
  QG(NaturalBCNodes(i)) = NaturalBCs(i);
end
QC = QG;
QC(EssentialBCNodes)=[];


kC = kG;
kC(EssentialBCNodes, :)=[];
kC(:,EssentialBCNodes)=[];

UC = kC\(fC+QC);
UG(UCNodes) = UC;

if (p==2)
    xG = 0:he:L;
else
    xG = 0:he/2:L;
end
 
NoPts = 100;    %no of points in a element to determine the solution with in the element
xList = zeros((NoPts-1)*N+1,1); 
uList = zeros((NoPts-1)*N+1,1);
count = 1;

for e = 1:N
    for i = 1:(NoPts-1)
        
        xe = (e-1)*he;
        xep1 = xe+he;

        xbar = (i-1)/(NoPts-1);
        xList(count) = xe + he*xbar;

        if (p==2)
            uList(count) = UG(e)*(1-xbar) + UG(e+1)*(xbar);
        end

        if (p==3)
            uList(count) = UG(2*e-1)*(1-xbar)*(1-2*xbar) + ...
                UG(2*e)*(xbar)*(1-xbar)*4 + ...
                UG(2*e+1)*xbar*(1-2*xbar)*(-1);
        end
        count = count + 1;


    end
end
xList(count) = L;
uList(count) = UG(n);


plot(xList,uList,'k', xG, UG, 'xk', 'linewidth',2);
xlabel('x'); ylabel('u');
set(gca,'fontsize', 18);

hold on;