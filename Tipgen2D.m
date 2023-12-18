%Tipgen2D
%Generates X,Y tip morphology from profile of curvature data.
%Curvature profile is first interpolated to yield equally spaced points.

function [X,Y,Sint,Kint,theta] = Tipgen2D(S,K)

[S,IS,~]=unique(S); 
%returns S without repetitions and in sorted order
%IS is the indices of the sorted S vector with no repetions
K=K(IS);
%Index curvature to match dimensions of sorted S vector

ls=length(S);
hls=round(ls/2); %median index - but not pole index necessarily

ms=min(abs([S(1) S(end)])); %absolute of minimum S, which is the most negative arclength
sp=ms/(hls-1); %scale minimum S by median(?)

%Define the arc length limits
Slim1=-floor(-S(1)/sp)*sp;%new higher S minimum, i.e Slim1 > -ms, but lower absolute value
Slim2=floor(S(end)/sp)*sp;%new lower S maximum 

Sint=[Slim1:sp:Slim2];%new arc length vector that is evenly spaced at sp intervals

Kint=interp1(S,K,Sint');% estimating new curvatures

%figure,plot(S,K)
%figure,plot(Sint,Kint),pause

%interpolated curvature values using linear interpolation  taking S as point, 
%K as values corresponding to given S and Sint' as coordinates to me map
%new estimated K onto

Kmid=(Kint(1:end-1)+Kint(2:end))/2; %Averaging between Kint values with the indices next to each other
%reduces Kint by one, but basically Kmid== Kint
%Why Kmid then??

%sp = 0.75; %empirical estimate that seems to give right shapes
dthet=Kmid*sp; %sp;%change in theta - biggest at pole, why multiply with sp?? What is sp?
theta=zeros(size(Kint)); %size of kint = size(kmid)+1

pole=-Slim1/sp+1;%correct pole index(at least for image 1)

%Angle theta - looks right to me with asigmoidal curve of increasing
%absolute angles away from pole
theta(1:pole-1)=-flipud(cumsum(flipud(dthet(1:pole-1))));
theta(pole+1:end)=cumsum(dthet(pole:end));
%figure,plot(theta)

xincs=sp*sin(theta); %sigmoidal function of Sint, peak at sint= 0 aka pole
yincs=sp*cos(theta); %Parabola function of Sint
%what is this sp???

X=zeros(size(Kint));
Y=zeros(size(Kint));
%figure,plot(yincs),pause
X(pole+1:end)=-cumsum(xincs(pole+1:end));
X(1:pole)=flipud(cumsum(flipud(xincs(1:pole))));%cumulative sum produces a parabola
Y(pole+1:end)=cumsum(yincs(pole+1:end));
Y(1:pole)=-flipud(cumsum(flipud(yincs(1:pole))));%cumulative sum produces a sigmoid

%lx = length(X); %added parameter