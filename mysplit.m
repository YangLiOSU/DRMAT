function [o]=mysplit(s,e, Scur)
global N p  Nxterm
global Sall Kall List
global len

if ((e-s)+1) < 2*len+1
    o      =-100;
    return;
end

IC_best=1e100; % initiate the information criterion
s1_best=1;
s2_best=1;
o_best=1;
Sall_cur_best=1;

for i=(s+len-1):(e-len)
    s1=mytsreg(s, i);
    s2=mytsreg(i+1, e);
    Sall_cur = Sall-Scur+s1+s2; 
    Kall_cur = Kall+Nxterm;
    IC_cur  = N*log(det(Sall_cur/N))+(Kall_cur*p+(p*(p+1)./2))*log(N); % BIC
    if (IC_cur < IC_best)
        IC_best=IC_cur;
        o_best=i; 
        s1_best=s1;
        s2_best=s2;
        Sall_cur_best=Sall_cur;
    end
end

IC_all  = N*log(det(Sall/N))+(Kall*p+(p*(p+1)./2))*log(N);

if (IC_best< IC_all)
   Sall = Sall_cur_best;
   Kall = Kall_cur;
   List = [List, o_best];
   [Oleft ]=mysplit(s,o_best,  s1_best);    
   [Oright]=mysplit(o_best+1,e, s2_best);
    o    = o_best; 
    return;
else
    o      =-100;
    return;
end
