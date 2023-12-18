%Fit experimental outlines to simulated outlines by finding shorested
%distance between a point and nearest line segment
function[min_dist] = distPoint2Line(outline,simpoints)
%INPUT
%outline: 1 X 2*num_points array - experimental outline(treated as line)
%simpoints: 1 X 2 * num_points array - simulated outlines(treated as points)
%OUTPUT
%mindist: 1 X (num_points-1) array - minimum(perpendicular distance between 1 simpoints
%        and line between two outline points, pole of simpoints is excluded
    num_points = length(outline)/2;
    
    xdiffs = diff(outline(1:num_points));
    ydiffs = diff(outline(num_points+1:num_points*2));
    
    %distance between consecutive points on outline
    euc_dist = (xdiffs.^2 + ydiffs.^2).^0.5;

    before_pole = floor(num_points/2);
    after_pole = before_pole +2;
    
    Xexpminussim = cat(1,(outline(1:before_pole) - simpoints(1:before_pole)),(outline(after_pole:num_points) - simpoints(after_pole:num_points)));
    Yexpminussim = cat(1,(outline(num_points+1:before_pole+num_points) - simpoints(num_points+1:before_pole+num_points)), (outline(after_pole+num_points:num_points*2) - simpoints(after_pole+num_points:num_points*2)));

    %distance between simulated point and line between two experimental points
    min_dist = abs(((xdiffs.*Yexpminussim) - (Xexpminussim .* ydiffs)))./euc_dist;
end