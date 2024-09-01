function Ye = dtrend(Y, Xt, Xs)

    X=[Xt Xs];
    if  (nargin<4)
        ridge=0.001;
    end

    XtX=X'*X;
    XtY=X'*Y;

    diagRidge    = diag(XtX);
    diagRidge(:) = ridge;

    XtX = XtX + diag(diagRidge);
    B   = linsolve(XtX,XtY);

    Ye = Y-Xt*B(1:size(Xt,2),:);
    
end