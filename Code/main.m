%%
clc
clear
close all

%%
root_dir = '/Users/liyang/Documents/Study/DRMAT/';

data_dir = sprintf('%sData/', root_dir);
addpath(data_dir)
load('Brazil.mat')
load('Ref.mat')

code_dir = sprintf('%sCode/',root_dir);
addpath(code_dir)

%%
M1 = 0;                  % detected disturbances disagree withreference
M2 = 0;                  % total number of detected disturbances
M3 = 0;                  % reference disturbances disagree with detection
M4 = size(Ref_year,1);   % total number of reference disturbances

band_fin = band;

%%
tic
for indx = 1:size(id,1)
    indx
    band1  = band_fin(indx,:,:);
    band1 = reshape (band1,[size(band,2),size(band,3)]);
    ndvi1 = ndvi(indx,:)';
    id1 = id(indx);
    time1  = t;
    T=time1;
    Y=band1;
    ymd1 = ymd;

    %%%%%%%%%%%%% get rid of NaN
    band1_nan = sum(band1,2);
    nan_indx = isnan(band1_nan);
    idx = find(nan_indx==0);  % non NaN

    band1 = band1(idx,:);
    ndvi1 = ndvi1(idx,:);
    time1 = time1(idx,:);
    T     = T(idx,:);
    Y     = Y(idx,:);
    ymd1  = ymd1(idx,:);

    %%%%%%%%%%%%%
    global Yg Xg Tg N p  Nxterm
    global  Sall Kall List;

    Yg = Y;
    Tg = T;

    N=length(T);
    p=size(Y,2);

    Xg=[Tg-Tg+1, Tg, sin(2*pi*Tg),cos(2*pi*Tg), ...
        sin(4*pi*Tg),cos(4*pi*Tg),...
        sin(6*pi*Tg),cos(6*pi*Tg)];
    %     Xg=zscore(Xg);

    global XS XT  XX YY
    global len
    len=6;

    XT = Xg(:,1:2);
    XS = Xg(:,3:end);
    %%
    XX=[XT XS];
    YY=Yg;
    Scur = mytsreg(1,N);
    Sall = Scur;
    Kall = size(XX,2);
    List = [];
    Nxterm=size(XX,2);
    o=mysplit(1,N, Scur);


    ListS=List;
    ListT=List;
    shouldContinue=1;

    i=0;
    while (shouldContinue)
        i=i+1;
        ListSOld=ListS;
        ListTOld=ListT;

        Xt=buildterms(XT, ListT);
        Xs=buildterms(XS, ListS);

        if  mod(i,2)==0
            len=40;
            YY=dtrend(Yg,Xt,Xs);
            XX=XS;
        else
            len=20;
            YY=dseason(Yg,Xt,Xs);
            XX=XT;
        end

        Scur = mytsreg(1,N);
        Sall = Scur;
        Kall = size(XX,2);
        List = [];
        Nxterm = size(XX,2);
        o = mysplit(1,N, Scur);
        if  mod(i,2)== 0
            ListS = List;
        else
            ListT = List;
        end

        if i<=1
            continue
        end

        if length(ListT)~=length(ListTOld) || length(ListS)~=length(ListSOld)
            shouldContinue=1;
        else
            if sum(sort(ListT)-sort(ListTOld))==0 && sum(sort(ListS)-sort(ListSOld))==0
                shouldContinue=0;
            else
                shouldContinue=1;
            end
        end

    end


    %% validation
    % reference data
    ref_idx = find(Ref_year(:,1)==id1);
    Ref_year_cur = Ref_year(ref_idx,2);

    % num to year
    ListT = ymd1(ListT);
    ListS = ymd1(ListS);

    [ListS ia] = setdiff(ListS,ListT);  % remove the breakpoints in S that have already been recorded in T
    if(isempty(ListS))
        ListS = [];
    end

    ListAll = [ListT,ListS];
    ListAll = sort(ListAll);

    %%
    clear i j
    % commission error
    for i = 1:size(ListAll,2)
        if ismember(ListAll(i),Ref_year_cur)==1
            M2 = M2+1;
        else
            M1 = M1+1;
            M2 = M2+1;
        end
    end
    % omission error
    for j = 1:size(Ref_year_cur,1)
        if ismember(Ref_year_cur(j),ListAll)==0
            M3 = M3+1;
        end
    end
    % indx
end
% commission error = detected disturbances disagree withreference(M1)/total number of detected disturbances(M2)
comission = M1./M2;
% omission error = reference disturbances disagree with detection(M3)/total number of reference disturbances(M4)
omission = M3./M4;
F1 = (((1-comission).*(1-omission))./(2-comission-omission)).*2;
toc
