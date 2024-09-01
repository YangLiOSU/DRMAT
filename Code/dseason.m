function Ye = dseason(Y, Xt, Xs)

    X=[Xs Xt];
    if  (nargin<4)
        ridge=0.001;
    end

    XtX=X'*X;
    XtY=X'*Y;

    diagRidge    = diag(XtX);
    diagRidge(:) = ridge;

    XtX = XtX + diag(diagRidge);
    B   = linsolve(XtX,XtY);

    Ye=Y-Xs*B(1:size(Xs,2),:);
    
end