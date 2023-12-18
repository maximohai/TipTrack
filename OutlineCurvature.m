%OutlineCurvature.m

%Rico Rojas, 2/21/21

%This program calcuates the meridional curvature, normal vectors, and
%normal angle of the outline of a tip-growing cell.  Curvature is
%calculated K=dPhi/ds, where Phimid is the angle from the x axis and s is the 
%arclength.  End points are interpolated from neighboring points.  
%This is an update of curvature6.m.

%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%X: row or column vector of x coordinates of the outline
%Y: row or column vector of y coordinates of the outline
%smooth: number of neighboring points over which to average, should be an odd number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%K: Meridional curvature of outline.  The sign convention is such that K is
%   positive if the outline is curving clockwise as the coordinate index
%   increases.  Column vector.
%N: nx X 2 (where nx is the number of coordinates) matrix with the X and Y
%   components of unit normal vectors to the outlines.
%Phi: Angles between the normals and the x axis.  Column vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [K,Phi,N] = OutlineCurvature(X,Y,smooth)

nx=length(X);%Number of coordinates

%Reshape X and Y into column vectors
X=reshape(X,[nx,1]);
Y=reshape(Y,[nx,1]);

%Pre-allocate matrices
Phi=zeros(nx,1);%Angles between normals and x-axis
KRaw=zeros(nx,1);%Meridional curvature before smoothing
N=zeros(nx,2);%Meridional curvature after smoothing

%Calculate Normal Vectors
dX=diff(X);
dY=diff(Y);
dS=sqrt(dX.^2+dY.^2);%Arc-distances between coordinates

Nx=dY./dS;%Normal vectors at midpoints between coordinates
Ny=-dX./dS;

Phimid=acos(Ny)-pi/2;

%Calculate curvature
dphi=-diff(Phimid);
KRaw(2:nx-1)=dphi./((dS(1:nx-2)+dS(2:nx-1))./2);
KRaw(1)=KRaw(2);
KRaw(nx)=KRaw(nx-1);

%Calculate Phi at original coordinates
Phi(2:nx-1)=Phimid(1:nx-2)+dphi/2;
Phi(1)=Phi(2);
Phi(nx)=Phi(nx-1);

%Calculate normals at original coordinates
N(:,1)=cos(Phi);
N(:,2)=-sin(Phi);

%Smooth curvature
K=movingaverage(KRaw',smooth)';


