%Fit Tip Shapes to morphospace

clear
close all


load("06142022_Control7.mat")
%Pre-allocate matrices
kj_coords = NaN(2,T);
kj_eucdist = NaN(2,T);

%Load simulated profiles
simShapes = importdata('SimulatedShapes.mat');
simulated_size = size(simShapes);

points =75;
cutoff = 0.9;


%%%%Fit average outline%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Symmetrize average outline
[Ksym_av,Ssym_av] = TipSymmetrize(Kav_int,Sav');

[Xsym_av,Ysym_av,Ssym_av,Ksym_int,theta_av] = Tipgen2D(Ssym_av,Ksym_av);

%Define apical tip 
lower_cut = find(theta_av > -cutoff * pi/2,1);
upper_cut = find(theta_av < cutoff * pi/2,1,"last");

Xsym_cut = Xsym_av(lower_cut:upper_cut);
Ysym_cut = Ysym_av(lower_cut:upper_cut);
Ssym_cut = Ssym_av(lower_cut:upper_cut); 

%Normalize for Radius      
radius = max(Ysym_cut) - min(Ysym_cut);
Ysym_cut = Ysym_cut./radius;
Xsym_cut = Xsym_cut./radius;  
Ssym_cut = Ssym_cut./radius; 
    
%Interpolate outlines to 75 points
exptXY = NaN(points*2,1);
SXY = NaN(points,1);

[exptXY,SXY] = interpolate_raw(Xsym_cut,Ysym_cut,Ssym_cut,points);

%Fit profile to morphospace
for j = 1:simulated_size(2)
        for k = 1:simulated_size(1)
           
           simXY = simShapes{k,j}';
            
           %Method 1: euclidean distance between points in two vectors
           norm_method1(k,j) = norm(exptXY - simXY); %l2 norm
            
           %Method 2: Minimum distance between point and line
           %distance between simulated point and line between two experimental point
           min_d = distPoint2Line(exptXY', simXY);
           norm_method2(k,j) = norm(min_d,2);
            %method 2 usually much better fit                    
        end    
end

%store coordinates of best fit for average profile
min_norm1 = min(min(norm_method1));
[k_eucdist_av,j_eucdist_av] = find(norm_method1 == min_norm1); %[row, col]    

min_norm2  = min(min(norm_method2));
[k_coords_av,j_coords_av] = find(norm_method2 == min_norm2); %[row, col]

figure,
hold on
plot(Xsym_cut,-Ysym_cut)
%scatter(exptXY(1:points),-exptXY(points+1:end))
%Method 2 simulated outlines
k = k_coords_av;j=j_coords_av;
plot(simShapes{k,j}(1:points), -simShapes{k,j}(points+1:end))
%Euclidiean fits
k = k_eucdist_av;j=j_eucdist_av;
plot(simShapes{k,j}(1:points), -simShapes{k,j}(points+1:end))

legend('expt','','simulated-method2','eucdist')

fig2pretty
axis equal
title("Average outline fit")


%%%%%%Fit each outline individually %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:T
     %Symmetrize the individual outlines
    [Ksymm,Ssymm,lout] = TipSymmetrize(Kout(1:Nout(t),t),S(1:Nout(t),t));
    Ksym(1:lout,t) = Ksymm;
    Ssym(1:lout,t) = Ssymm; 
    [lx(t),Xint(1:lx(t),t),Yint(1:lx(t),t),Sint(1:lx(t),t),Kint(1:lx(t),t),theta(1:lx(t),t)] = Tipgen2D_annotated(Ssymm,Ksymm);

    %Define the apical tip
    lower_cut = find(theta(1:lx(t),t) > -cutoff * pi/2,1);
    upper_cut = find(theta(1:lx(t),t) < cutoff * pi/2,1,"last");
    
    Xcut = Xint(lower_cut:upper_cut,t);
    Ycut = Yint(lower_cut:upper_cut,t);
    Scut = Sint(lower_cut:upper_cut,t); 

    
    
    %Normalize for Radius      
    radius = max(Ycut) - min(Ycut);
    Ycut = Ycut./radius;
    Xcut = Xcut./radius;  
    Scut = Scut./radius; 
        
    %Interpolate outlines to 75 points
    exptXY = NaN(points*2,1);
    SXY = NaN(points,1);
    
    [exptXY,SXY] = interpolate_raw(Xcut,Ycut,Scut,points);
    
  
    %Fit profile to morphospace
    for j = 1:simulated_size(2)
            for k = 1:simulated_size(1)
               
               simXY = simShapes{k,j}';
                
               %Method 1: euclidean distance between points in two vectors
               norm_method1(k,j) = norm(exptXY - simXY); %l2 norm
                
               %Method 2: Minimum distance between point and line
               %distance between simulated point and line between two experimental point
               min_d = distPoint2Line(exptXY', simXY);
               norm_method2(k,j) = norm(min_d,2);
                %method 2 usually much better fit                    
            end    
    end

    %Save fit coordinates and profiles
    %minimum error fit
    min_norm1 = min(min(norm_method1));
    [kj_eucdist(1,t),kj_eucdist(2,t)] = find(norm_method1 == min_norm1); %[row, col]    

    min_norm2  = min(min(norm_method2));
    [kj_coords(1,t),kj_coords(2,t)] = find(norm_method2 == min_norm2); %[row, col]

    figure,
    hold on
    %plot(Xint(1:lx(t),t),-Yint(1:lx(t),t))
    plot(Xcut,-Ycut)
    %Method 2 simulated outlines
    k = kj_coords(1,t);j=kj_coords(2,t);
    plot(simShapes{k,j}(1:points), -simShapes{k,j}(points+1:end))
    %Euclidiean fits
    k = kj_eucdist(1,t);j=kj_eucdist(2,t);
    plot(simShapes{k,j}(1:points), -simShapes{k,j}(points+1:end))
  
    legend('expt-cut','simulated-method2','eucdist')
    %ylim([-50 50])
    fig2pretty
    axis equal

%     
%     for t = 1:T
%        figure,
%        plot(Xint(1:lx(t),t),-Yint(1:lx(t),t))
%        fig2pretty
%        axis equal
%     end
%    



end

%Save fitted profiles
% cd('/Users/maximohairwe/Desktop/MATLAB/Achlya/TipVideos')
% save("03162022_LatB_2_fit")

%Map onto Morphospace
Parameters = importdata("ParameterMatrices.mat");
fields = fieldnames(Parameters);
Arat_mat = Parameters.(fields{1});
lrat_mat = Parameters.(fields{2});

labels = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','1','2','3','4','5','6','7','8','9'];

%Euclidean fits
figure,
hold on
for n =1:T
    
    k1 = kj_eucdist(1,n); j1=kj_eucdist(2,n);
    scatter(lrat_mat(k1,j1),Arat_mat(k1,j1),'+')
    text(lrat_mat(k1,j1),Arat_mat(k1,j1),labels(n))

end
scatter(lrat_mat(k_eucdist_av,j_eucdist_av),Arat_mat(k_eucdist_av,j_eucdist_av),'o')

xlabel("lrat")
ylabel("Arat")

title("LatB_1 euclidean fits")
xlim([1 10])
ylim([0 0.52])
xline(5)
yline(0.05)
fig2pretty

%Method2 fits
figure,
hold on
for n =1:T

    k2 = kj_coords(1,n);j2=kj_coords(2,n);
    scatter(lrat_mat(k2,j2),Arat_mat(k2,j2),'o')
    text(lrat_mat(k2,j2),Arat_mat(k2,j2),labels(n))
end

scatter(lrat_mat(k_coords_av,j_coords_av),Arat_mat(k_coords_av,j_coords_av),'+')

xlabel("lrat")
ylabel("Arat")

title("LatB_1 Method 2 fits")
xlim([1 10])
ylim([0 0.52])
xline(5)
yline(0.05)
fig2pretty


%%%Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [expt_outline_interp,S_interp] = interpolate_raw(X_expt,Y_expt,S_expt,num)
%Reinterpolate raw outlines for fitting with simulated outlines
%X_expt = X coordinates form expt outlines
%Y_expt = Y coordinates from expt outlines
%S_expt = Arc length of expt outlines
%num = number of points to reinterpolate outlines to (default=75)

    pole = find(S_expt==0, 1, 'first');
    %last = find(isnan(S_expt),1,'first') - 1; --if input is matrix columns
    last = length(S_expt);
    
    %arc length of each symmetrical half
    Shalf = abs(S_expt(1) - S_expt(pole));
    sp = 1/floor(num/2);                      %spacing between points for each symmetrical half
    ds = sp*Shalf;                            %sp*mag is interpolated ds
    new_pole = ceil(num/2);                   %redefine pole position
   
    %interpolate arclength to num points
    %S_interp = NaN(num,1);
    S_interp(1:new_pole) = [S_expt(1):ds:S_expt(pole)]; 
    S_interp(new_pole:num) = [S_expt(pole):ds:S_expt(last)]; %constant arclength between points
     
    
    %Interpolate coordinates to num points
    
    expt_outline_interp(1:num) = interp1(S_expt(1:last),X_expt(1:last),S_interp);
    
    expt_outline_interp(num+1:2*num) = interp1(S_expt(1:last),Y_expt(1:last),S_interp);
     
 end