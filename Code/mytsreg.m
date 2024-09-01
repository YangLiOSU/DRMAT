function SIG = mytsreg(s, e, ridge)

global YY XX

X=XX(s:e,:);
Y=YY(s:e,:);

if  (nargin<3)
    ridge=0.001;
end

XtX=X'*X;
XtY=X'*Y;

diagRidge    = diag(XtX);
diagRidge(:) = ridge;

XtX = XtX + diag(diagRidge);
B   = linsolve(XtX,XtY);

E=X*B-Y;
SIG=E'*E;
end