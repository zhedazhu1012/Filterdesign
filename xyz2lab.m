function Lab=xyz2lab(XYZ,XYZn)

% XYZ2LAB: calculates L,a,b from XYZ
% according to CIE Publication 15.2
%
% Example: LAB=XYZ2Lab(XYZ,XYZn)
% where XYZn is the reference white.
%
% If no argument is supplied for the reference white,
% D50 perfect diffuser is assumed. 
%
%   Colour Engineering Toolbox
%   author:    � Phil Green
%   version:   1.1
%   date:  	   16-02-2002
%   book:      http://www.wileyeurope.com/WileyCDA/WileyTitle/productCd-0471486884.html
%   web:       http://www.digitalcolour.org

if nargin<2
   rwhite=d(50);
else
   rwhite=XYZn;
end

Xn=rwhite(1);Yn=rwhite(2);Zn=rwhite(3); 
X=XYZ(:,1);Y=XYZ(:,2);Z=XYZ(:,3);

% normalize to media white
Yrel=Y/Yn;Xrel=X/Xn;Zrel=Z/Zn;
fY=Yrel.^(1/3);fy=fY;
fX=Xrel.^(1/3);
fZ=Zrel.^(1/3);

% Check for T/Tn <= 0.008856
p=find(Yrel<=0.008856);
q=find(Xrel<=0.008856);
r=find(Zrel<=0.008856);

% calculate f(T) for T/Tn <= 0.008856
fY(p)=7.787*Yrel(p)+16/116;
fX(q)=7.787*Xrel(q)+16/116;
fZ(r)=7.787*Zrel(r)+16/116;

% calculate L,a,b
L=116*fy-16;
L(p)=903.3*Yrel(p);
a=(fX-fY)*500;
b=(fY-fZ)*200;

Lab=[L,a,b];
