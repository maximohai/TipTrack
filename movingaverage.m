function y=movingaverage(x,s);
%

[mx,lx]=size(x);
y=zeros(mx,lx);

if s>lx
    s=round(lx/2);
end

for i=1:lx
        s2=floor(s/2);
    if i-s2<=0
        y(:,i)=mean(x(:,1:i+s2)',"omitnan")';
    elseif i+s2>lx
        y(:,i)=mean(x(:,i-s2:end)',"omitnan")';
    else
        y(:,i)=mean(x(:,i-s2:i+s2)',"omitnan")';
    end
end

%Update: nanmeans were repalced by mean(,"omitnan")