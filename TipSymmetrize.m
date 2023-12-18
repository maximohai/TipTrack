%Symmetrize the curvature and arc length

function[Ksym,Ssym,len] = TipSymmetrize(K,S)
%assumption: pole is the point of highest curvature
[~,Pole_index] = max(K);

%Symmetrize the curvature profile
last = find(K,1,'last'); % last nonzero element
Kbefore_pole = K(1:(Pole_index-1));
Kafter_pole = K((Pole_index+1):last);
Sbefore_pole = S(1:(Pole_index-1));
Safter_pole = S((Pole_index+1):last);

%Half before pole bigger than half after pole
if length(Kbefore_pole) > length(Kafter_pole)
   %Shorten half before pole
   Kbefore_pole = K((Pole_index - length(Kafter_pole)):(Pole_index -1));
   Sbefore_pole = S((Pole_index - length(Kafter_pole)):(Pole_index -1));

   newPoleindex = length(Kbefore_pole)+1;
   %symmetrized pole index
   Ksym(newPoleindex) = K(Pole_index);
   %Ssym(newPoleindex) = S(Pole_index);%zero
   Ssym(newPoleindex) = 0;%zero

   %first half of symmetrized curvature before pole
   Ksym(1:newPoleindex-1) = (Kbefore_pole + flipud(Kafter_pole))/2;
   Ssym(1:newPoleindex-1) = (Sbefore_pole - flipud(Safter_pole))/2;
   %second half of symmetrized curvature after pole
   Ksym(newPoleindex+1: newPoleindex + length(Kafter_pole)) = (flipud(Kbefore_pole) + Kafter_pole)/2;
   Ssym(newPoleindex+1: newPoleindex + length(Safter_pole)) = (-flipud(Sbefore_pole) + Safter_pole)/2;

%Half before pole smaller than half after pole   
   elseif length(Kbefore_pole) <= length(Kafter_pole)
   %Shorten half after pole
   Kafter_pole = K((Pole_index+1): (Pole_index + length(Kbefore_pole)));
   Safter_pole = S((Pole_index+1): (Pole_index + length(Sbefore_pole)));

   newPoleindex = Pole_index;
   %symmetrized pole
   Ksym(newPoleindex) = K(Pole_index);
   %Ssym(newPoleindex) = S(Pole_index);%zero
   Ssym(newPoleindex) = 0;%zero
   
   %Symmetrized half before pole
   Ksym(1:newPoleindex-1) = (Kbefore_pole + flipud(Kafter_pole))/2;
   Ssym(1:newPoleindex-1) = (Sbefore_pole -flipud(Safter_pole))/2;
   %Symmetrized half after pole
   Ksym(newPoleindex+1: newPoleindex + length(Kafter_pole)) = (flipud(Kbefore_pole) + Kafter_pole)/2;
   Ssym(newPoleindex+1: newPoleindex + length(Safter_pole)) = (-flipud(Sbefore_pole) + Safter_pole)/2;

end

len = length(Ksym);
end