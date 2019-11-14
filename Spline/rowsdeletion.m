
s=who
for i=1:length(s)
    
     dd= eval(s{i} )
     
     
     
%      j= 3
%     dd(j)=[];
%      
     j= 8
    dd(j)=[];
   
%      j= 6
%       dd(j)=[];
%        
% %       j= 8
%     dd(j)=[];
% %    
%      j= 7
%       dd(j)=[];
%     
   
      
      
    
 eval([  s{i} '= dd '])
    
    
        
    
end

clear s dd i j
