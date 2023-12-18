%TipTrack
%
%7/29/21
%Rico Rojas
%
%Tracks the outline of tip growing cells from phase contrast image stacks.  
%
%INSTRUCTIONS FOR USE:
%Rotate and crop image stacks such that cell is growing from left to right
%and is rougly centered vertically. Save phase image stack in a directory
%by itself.
%
%OUTPUT:
%X: X coordinates of cell outlines.  
%Y: Y coorinates of cell outlines
%K: Curvature of cell outlines
%Kav: Average curvature of cell outlines
%v: Elongation rate of cell during movie
%s: Arcdistance from pole

%Calls on the following m-files:
%norm16bit.m
%OutlineCurvature.m
%movingaverage.m
%fig2pretty.m

clear
close all
addpath('/Users/maximohairwe/Desktop/MATLAB/Achlya/TipVideos')
%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basename='03172022_LatB80_1_31';%Name of image stack to be tracked

dirname=['/Users/maximohairwe/Desktop/MATLAB/Achlya/TipVideos/' basename];%Directory in which image stack is stored.

ppix=0.25;%Percentage of pixels to saturate at min and max values (default=0.25).
sigma=5;%Standard deviation of Gaussian filter with which to smooth image in edge detector (Default=5)
min_grad=0.75;%Minimum gradient threshold for cell outline (Default=0.75.  Must be between 0 and 1)
minL=350;%Minimum length threshold for cell outline (Default=300)
cs=0.001;%Parameter for smoothing spline (Default=0.001)
s1=10;%Number of pixels over which to normal vectors to cell outline (Default=10)
R1=16;%Radii of structural elements used in image operations(Default=16)
R2=2;%(Default=2)
R3=5;%(Default=5)
R4=7;%%(Default=7)
NPath=51;%Number of paths to test to test for polar path (Default=51)
s2=11;%Smoothing paramter used for calculating polar path (Default=11)
s3=30;%Number of pixels over which to calculate curvature (Default=10)

lscale= 0.1548;%um/pixel (60X oil objective)
%lscale =0.09;%100 x
tscale=30;%Frame rate in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
curdir=cd;
cd(dirname);%Point to directory in which file is saved
directory=dir('*.tif');%List filenames of images
T=length(directory);%Determine the number of frames
%T= T-2; %for tip 8
%T= T-10; %for tip 3

%T=29;%LatB_3_h
%T=10;%LatB_2
T =31;%LatB_1 
%T=29;%LatB_3_h

%Pre-allocate matrices
imagename=directory(1).name;
I=imread(imagename);
All_Outlines=zeros(size(I));

[IM,IN]=size(I);
IM_half=round(IM/2);
[IX,IY]=meshgrid(1:IN,1:IM);
IY=flipud(IY);


%Detect cell outlines
f = waitbar(0,'Detecting cell outlines'); %progress bar

for t=1:T
    
f = waitbar(t/T,f,'Detecting cell outlines');
%try
    %Load images
    imagename=directory(t).name;
    I=imread(imagename);
    Ic=imcomplement(I);
    [ImM,ImN]=size(I);
    
    %Normalize images
    I=norm16bit(I,ppix);
    
    %Calculate background
    [imcounts,bins]=imhist(I);
    [imcounts,idx]=sort(imcounts);
    bins=bins(idx);
    bg=bins(end);
    I2=imadjust(I,[bg/65535 1],[]);
    
    %Calculate image gradient
    %This is the weakest part of the code.  The reason I calculated Gdir is to
    %use it as a criteria for selecting the right outline from the edge
    %detector below, but I don't think the code will be robust to things like
    %cell drift, or if the cells are off center.
    
    [Gmag,Gdir]=imgradient(I2);
    out_cri=zeros(IM,IN);
    
    out_cri(Gdir>100 | Gdir < -100)=1;
    out_cri(Gdir>80 & Gdir <100 & IY<IM_half)=1;
    out_cri(Gdir<-80 & Gdir >-100 & IY>IM_half)=1;
    
    %Find edges
    [ed,thresh2]=edge(I,'canny',[],sigma);
    
    %Bridge small gaps
    ed=bwmorph(ed,'bridge'); 
    
    %Find cell outline based on length of outline and gradient direction
    cc=bwconncomp(ed,8);
    stats=regionprops(ed,out_cri,'Area','PixelList','MeanIntensity');
    Pixel_Cell={stats.PixelList};
    idx=find([stats.Area]>minL & [stats.MeanIntensity]>min_grad);
    ed2=ismember(labelmatrix(cc),idx);
    ed3=ed2;
    
    %Extend outlines to x=0 if they do not already
    end_points=bwmorph(ed2,'endpoints');
    [epy,epx]=find(end_points);
    epx_dist=epx(1:2);
    epy_dist=epy(1:2);
    cell_lim_1=max(epx(1:2));
    
    if cell_lim_1>2
        ed2(epy_dist(1),2:epx_dist(1))=1;
        ed2(epy_dist(2),2:epx_dist(2))=1;
    end
    
    %Watershed cell
    peaks=bwdist(ed2);
    
    SE1=strel('disk',R1,4);
    peaks_erode=imerode(peaks,SE1);
    peaks_recon=imreconstruct(peaks_erode,peaks);
    peaks_dilate=imdilate(peaks_recon,SE1);
    peaks_recon_2=imreconstruct(imcomplement(peaks_dilate),imcomplement(peaks_recon));
    peaks_recon_2=imcomplement(peaks_recon_2);
    peaks_recon_com=imcomplement(peaks_recon_2);
    
    boundary=zeros(IM,IN);
    boundary(1,:)=1;
    boundary(end,:)=1;
    boundary(:,end)=1;
    
    reg_max=imregionalmax(peaks_recon_2);
    valleys=imimposemin(peaks_recon_com,boundary | reg_max);
    
    ws=watershed(valleys,8);
    
    SE2=strel('disk',R2,4);
    cell_ws=zeros(IM,IN);
    cell_ws(ws>=2)=1;
    cell_ws=imclose(cell_ws,SE2);
    
    ws(cell_ws==1)=2;
    cell_body=ws==2;
    
    figure,imshowpair(I,ed2)
    %This is also for if the tracked cell did not extend to the left side of
    %the image
    SE3=ones(R3,1);
    cell_lim=1;
    if cell_lim_1>2 
        cell_body_dil=imdilate(cell_body,SE3);
        cropped_ed=cell_body_dil & ed3;
        end_points=bwmorph(cropped_ed,'endpoints');
        [epy,epx]=find(end_points);
        cell_lim=max(epx(1:2))+R4-1;
    end
    
    %Find boundary of cell
    SE4=strel('disk',R4,4);%Erode cell so that outline corresonds to actual cell boundary from images
    cell_body=imerode(cell_body,SE4);
    cell_boun_cell=bwboundaries(cell_body,4,'noholes');
    cell_boun=cell_boun_cell{1};
    
    Xb=cell_boun(:,2);
    Yb=cell_boun(:,1);
    
    Yb(Xb<3)=0;%Erase outline at left (distal) edge of cell
    Xb(Xb<3)=0;
    Xb=Xb(Xb>0);
    Yb=Yb(Yb>0);
    
    cell_out=bwmorph(cell_body,'remove');
    cell_out(:,1:cell_lim)=0;
    All_Outlines(cell_out)=1;
    
    %Shift outlines so that they start at the upper left point
    indx=find(Xb,1);
    Xb=circshift(Xb',1-indx)';
    Yb=circshift(Yb',1-indx)';
    Xb=flipud(Xb);
    Yb=flipud(Yb);
    Xb=Xb';
    Yb=Yb';
    
    %Calculate arcdistance
    dX=diff(Xb);
    dY=diff(Yb);
    ds=sqrt(dX.^2+dY.^2);
    s=[0 cumsum(ds)];
    
    %Smooth outline
    Xb=csaps(s,Xb,cs,s)';
    Yb=csaps(s,Yb,cs,s)';
    
    %Re-spline to make arc-spacing constant (1 pixel)
    snew=[0:s(end)];
    
    Xb=spline(s,Xb,snew);
    Yb=spline(s,Yb,snew);
    Nout(t)=length(Xb);

    %Re-Calculate arcdistance after smoothing Xb/Yb
    dX=diff(Xb);%Xi+1 - Xi
    dY=diff(Yb);
    ds=sqrt(dX.^2+dY.^2); %approximates arclength if very small
    s =[0 cumsum(ds)]; %smoothened arc length

    
    %Save outlines
    X(Nout(t),t)=0;
    Y(Nout(t),t)=0;
    S(Nout(t),t) =0;

    X(1:Nout(t),t)=Xb;
    Y(1:Nout(t),t)=Yb;
    S(1:Nout(t),t) =s;
    
    %Calculate radius
    radius(t) = max(Yb) - min(Yb);
    
    %Calculate curvature
    Kout(Nout(t),t)=0;
    Kout(1:Nout(t),t)=OutlineCurvature(Xb,Yb,s3);

    %Find preliminary pole
    %Horizontal pole
    [tip_x(t),tip_ind(t)]=max(Xb); %Assumption: pole is the most horizontally extended point
    tip_y(t) = Yb(tip_ind(t));

    %Curvature pole
    [tip_k(t),tip_indk(t)] = max(Kout(1:Nout(t),t));

    %Center Arclength around pole
    %pole(t) = tip_indk(t);%Use curvature pole
    pole(t) = tip_ind(t);%Horizontal pole

    S(1:Nout(t),t) = S(1:Nout(t),t)- S(pole(t),t);

end 

%T = T-skip;
%Remove all zero columns
Kout(:,all(Kout == 0)) = [];
%S(:,all(S == 0)) = [];
X(:,all(X == 0)) = [];
Y(:,all(Y == 0)) = [];
Nout = nonzeros(Nout);
T = length(Nout);

%Calculate preliminary growth rate
X_disp=diff(tip_x);
X_GR=X_disp*lscale/tscale;
time=[0:T-2]*tscale;

%Find normal paths
for a=-T:-2
    try
        t=-a;
        f=waitbar((T-t)/T,f,'Finding polar axis');
        
        X1=X(1:Nout(t-1),t-1);%Cell outline and time point t-1
        Y1=Y(1:Nout(t-1),t-1);
        
        X2=X(1:Nout(t),t);%Cell outline and time point t
        Y2=Y(1:Nout(t),t);
        
        X2mid=(X2(s1+1:end)+X2(1:end-s1))/2;
        Y2mid=(Y2(s1+1:end)+Y2(1:end-s1))/2;
        
        %Find unit tangent vectors to outline at time point t
        dX2=X2(s1+1:end)-X2(1:end-s1);
        dY2=Y2(s1+1:end)-Y2(1:end-s1);
        Tangents=[dX2';dY2']';
        Lengths=sqrt(Tangents(:,1).^2+Tangents(:,2).^2);
        Lengths=[Lengths,Lengths];
        Tangents=Tangents./Lengths;
        
        %Calculate unit normal vectors to outline at time point t
        rotmat=[0 1;-1 0];
        Norms=Tangents*rotmat;
        
        %Calculate slope and y-intercept of lines normal to outline at time point t
        Slopes_2=Norms(:,2)./Norms(:,1);
        Y0_2=Y2mid-Slopes_2.*X2mid;
        
        %Find lines orthogonal these lines that pass through each point on cell
        %outline in time point t-1 (??)
        OrthSlopes=-Norms(:,1)./Norms(:,2);
        
        Y0_1_mat=ones(length(OrthSlopes),1)*Y1'-OrthSlopes*X1';
        Y0_2_mat=Y0_2*ones(1,length(X1));
        
        Slopes_mat=Slopes_2*ones(1,length(X1));
        OrthSlopes_mat=OrthSlopes*ones(1,length(X1));
        
        %Find points where the lines intersect
        X_intersect=-(Y0_1_mat-Y0_2_mat)./(OrthSlopes_mat-Slopes_mat);
        Y_intersect=Slopes_mat.*X_intersect+Y0_2_mat;
        
        %Use these intersection points to calculate the intersection between lines
        %normal to the outline at time point t and cell outline at time point t-1
        Y1_mat=ones(length(Y2mid),1)*Y1';
        X1_mat=ones(length(X2mid),1)*X1';
        
        X2_mat=X2mid*ones(1,length(X1));
        Y2_mat=Y2mid*ones(1,length(Y1));
        
        dists=sqrt((Y_intersect-Y1_mat).^2+(X_intersect-X1_mat).^2);
        dists2=sqrt((Y2_mat-Y1_mat).^2+(X2_mat-X1_mat).^2);
        
        %Eliminate intersection points that are too far away. 
        dist_thresh=2*X_disp(t-1);
        dists(dists2>dist_thresh)=Inf;
        [~,min_ind]=min(dists');
        
        Pix_buf=ceil(s1/2);
        
        %Save connectivity
        Connectivity(Nout(t),t)=0;
        Connectivity(Pix_buf+1:Nout(t)-Pix_buf,t)=min_ind';
    catch
        continue;
    end

end
close(f)

%Find polar path
Path_i=zeros(NPath,T);
Path_i(:,T)=[tip_ind(T)-floor(NPath/2):tip_ind(T)+floor(NPath/2)]';

%Fill in gaps in paths
for a=-T+1:-1
t=-a;

RawPath=Path_i(:,t+1);
Path_i(:,t)=Connectivity(RawPath,t+1);
Path_i(Path_i(:,t)==1,t)=Path_i(Path_i(:,t)==1,t+1);
Path_i(Path_i(:,t)==0,t)=Path_i(Path_i(:,t)==0,t+1);
end

%Redundant code(?) - only difference is it finds pole_K- which is for the
%pole of the 51 paths
XYsize=size(X);
for i=1:NPath
    Path_Test=Path_i(i,:);
    Kshift=zeros(2*XYsize(1),T);
    Kshift(1:XYsize(1),:)=Kout;
    for t=1:T
        Kshift(:,t)=circshift(Kshift(:,t),XYsize(1)-Path_Test(t));
    end
    Kshift(Kshift==0)=NaN;
    Pole_K(i)=mean(Kshift(XYsize(1),:));%mean curvature at the tip across all frames for a given path(?)
    Kshiftav=mean(Kshift',"omitnan")';
end

%Calculate polar path as the normal path that maximizes the path integral
%of meridional curvature
Pole_K=movingaverage(Pole_K,s2);
[~,pole]=max(Pole_K);

%Calculate average curvature profile across time
Polar_Path=Path_i(pole,:);
X_shift=zeros(2*XYsize(1),T);
Y_shift=zeros(2*XYsize(1),T);
Kshift=zeros(2*XYsize(1),T);

X_shift(1:XYsize(1),:)=X;
Y_shift(1:XYsize(1),:)=Y;
Kshift(1:XYsize(1),:)=Kout;
for t=1:T
    X_shift(:,t)=circshift(X_shift(:,t),XYsize(1)-Polar_Path(t));
    Y_shift(:,t)=circshift(Y_shift(:,t),XYsize(1)-Polar_Path(t));
    Kshift(:,t)=circshift(Kshift(:,t),XYsize(1)-Polar_Path(t));
end

Kshift(Kshift==0)=NaN;
Kav=mean(Kshift',"omitnan")';
Kav(isnan(Kav))=0;

First_ind=find(Kav,1,'first');
Last_ind=find(Kav,1,'last');

X_shift=X_shift(First_ind:Last_ind,:);
Y_shift=Y_shift(First_ind:Last_ind,:);
K=Kshift(First_ind:Last_ind,:);

Kav=Kav(First_ind:Last_ind);
%Pole=XYsize(1)-First_ind+1;


%Create time and arclength vectors
Pole = tip_ind(end); %pole at most horizontal point
%[Pole_val,Pole] = max(Kav); %- pole at max curvature
%Pole = 1154;%tip7
%Pole = 1156;%tip8
%Pole = 1152; %tip10
%Pole = 1161;%tip11
%Pole = 1097;%tip14
%Pole = 1115; %tip4
%Pole = 1307;%tip5
%Pole =1300; %tip8
%Pole= 1287;%tip10
%Pole=1231;%tip11
%Pole=1094;%tip12
%s = s - s(Pole);
time=[0:T-1];

%Calculate elongation rate
Xpole=X_shift(Pole,:);
Ypole=Y_shift(Pole,:);
v=sqrt(diff(Xpole).^2+diff(Ypole).^2);

%Dimsionalize variables
s=s*lscale;
S = S*lscale;
time=time*tscale/60;
Kout=Kout/lscale;
Kav=Kav/lscale;
v= v*lscale/tscale*60;

%Plot results
figure,hold on
for t=1:T
    plot(X(1:Nout(t),t),-Y(1:Nout(t),t),'-k','Linewidth',0.5)
end
plot(X_shift(Pole,:),-Y_shift(Pole,:),'-r')
axis equal
fig2pretty

figure,plot(S(:,end),Kav)
xlim([-50 50])
xlabel('Arclength (\mu m)')
ylabel('Meridional Curvature (1/\mum)')
fig2pretty

figure,plot(time(1:end-1),v)
xlabel('Time (m)')
ylabel('Elongation Rate (\mum/m)')
fig2pretty


%Generate average outline
[Xav,Yav,Sav,Kav_int,theta_av] = Tipgen2D(S(:,end),Kav);

figure,
plot(Xav,-Yav)
axis off, axis equal
fig2pretty

cd('/Users/maximohairwe/Desktop/MATLAB/Achlya/TipVideos/' )
save('03172022_LatB_1.mat')

toc
