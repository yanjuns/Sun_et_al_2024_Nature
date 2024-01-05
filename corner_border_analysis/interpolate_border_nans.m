function M=interpolate_border_nans(M0)
M=M0;
M(1:end,1)=clean_nans(M0(1:end,1));
M(1:end,end)=clean_nans(M0(1:end,end));
M(1,1:end)=clean_nans(M0(1,1:end));
M(end,1:end)=clean_nans(M0(end,1:end));
function Z=clean_nans(Z0)
X=1:length(Z0);
X0=X;
Zcopy=Z0;
aux=isnan(Z0);
X0(aux)=[];
Z0(aux)=[];
if(length(X0)>2)
    Z=interp1(X0,Z0,X,'spline','extrap');
else
    Z=Zcopy;
end
