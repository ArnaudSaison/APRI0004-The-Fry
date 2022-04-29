CG_empty = [2 : 0.01 : 4];
empty_weight = 1121.77 ; 
k=0;
for i=1:length(CG_empty)
[min_CG, max_CG] = center_gravity(empty_weight, CG_empty(i));
[ Kn_max , hn_min_CG, xn_min_CG] = stability(min_CG);
[ Kn_min , hn_max_CG, xn_max_CG] = stability(max_CG);
    
    if 0<Kn_min && Kn_min>0.05 && 0<Kn_max && Kn_max<0.2
        
        if k==0 
        CG_empty_min = CG_empty(i) ;
        k=k+1 ;
        end 
        CG_empty_max = CG_empty(i) ;
    end
end
fprintf('CG of the empty plane must be between %d\n and %d\n', CG_empty_min, CG_empty_max) 

