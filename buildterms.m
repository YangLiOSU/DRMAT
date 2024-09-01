function Xall=buildterms(X, List)

if isempty(List)
    Xall=X;
    return;
end

List = sort(List);
N=size(X,1);
K=size(X,2);
Xall=[];
for i=1:(numel(List)+1)
    Xtmp=X;
    if i==1
        s=1;
    else
        s=List(i-1);
    end
    
    if i==(numel(List)+1)
        e=N;
    else
        e=List(i)-1;
    end
     
   Xtmp(1:(s-1),:)=0;
   Xtmp((e+1):N,:)=0;
   Xall=[Xall, Xtmp];
end

