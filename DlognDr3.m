% this function takes a bunch of particles, rp, and a grid, r, and casts
% the particles onto the grid to make a density distribution, then
% calculates (dn/dr)/n=d(log(n))/dr.  

% It uses a local least-squares fit at each grid point

function [dlogndr, dens1, dens2,vr1,vr2]=DlognDr2(Np1,weight1,rp1,vp_1,Np2,weight2,rp2,vp_2,r)
global logn dens1 

% takes in two particle arrays, converts them to radial density functions
% n(r), then adds them to get the total density and compute d/dr(log(n))
% note: the radial grid is cell-centered with ghost points at both ends

dr = r(2)-r(1) ; % assume a uniform radial grid, cell-center
Vol = 4/3*pi*((r+.5*dr).^3-(r-.5*dr).^3);  % cell volumes

% estimate the widths of the particle distributions
rrms1 = sqrt( mean(rp1.^2) ) ./ sqrt(3) ;
rrms2 = sqrt( mean(rp2.^2) ) ./ sqrt(3) ;

% cast density to the 2 nearest neighbor grid points using linear
% interpolation.

Nr=length(r);

% #####################  start particle type 1 #######################
dens1=0*r;vr1=0*r;
for n=1:Np1
    
    % find the cell a particle is in
    x=rp1(n)/dr;
    i = ceil(x+0.5)+1; % the index of the cell-centered grid point that is just
                 % above where the particle is
    i=min(i,Nr+1);
    f = (x-i+1.5); % fraction of the way the particle is from the lower edge
                   % of the cell to the upper edge.
    f = (x-i+2.5) ; %sdb fix
                 
%    fprintf(' rp %g i %g r(i) %g  r(i-1) %g f %g \n',rp1(n)/dr,i,r(i)/dr,r(i-1)/dr,f)
%    pause
                 
    % quadratic weighting, 2, linear weighting ,1 
    q=1;

     f=1 - (r(i)^q-rp1(n)^q)/(r(i)^q-(r(i)-dr)^q);
    % cast density to the grid point below i and to i. Do something tricky
    % take care of the case that i=1, so i-1 is zero, which is not a grid
    % point
    
    dens1(i) = dens1(i)+f;
    vr1(i) = vr1(i)+f*vp_1(n);
    ibelow=i-1;
    %  ibelow=max(i-1,1);
    if ibelow>1
      dens1(ibelow)= dens1(ibelow) + 1-f;
      vr1(ibelow)= vr1(ibelow) + (1-f)*vp_1(n);
    end
    
% check for casting trouble
%     if f<0 | f>1
%         fprintf(' rp1 %g f %g x %g i %g \n',rp1(n),f,x,i)
%         pause
%     end
%     

    
end


% convert these integer counting densities to densities by dividing by
% the cell volumes

% average the velocities in their cells
vr1=vr1./(dens1+1e-12);  % dens1 is just a particle count at this stage,
                         % divide by it to average
% load the ghost point below r=0
vr1(1)=-vr1(2);
dens1=dens1./Vol*weight1; % convert the dens1 particle count into a density


% the first few points are noisy - smooth them

ibad=2; % ibad is the last point to smooth. Use the next delta points to make a fit.
[~,dum1] = min( (r - rrms1).^2 ) ;
% dum1 = min( dum1, length(r) ) ;
dum1 = ceil(dum1/3) +1;
delta=max(3,dum1); % adjust these if you don't like the way the density looks near r=0


x=r(ibad+1:ibad+delta).^2/dr^2; % density is a function of r^2
y=dens1(ibad+1:ibad+delta);
p=polyfit(x,y,1); % fit
for i=2:ibad+1
    dens1(i)=polyval(p,r(i).^2/dr^2);
end
% load the ghost point below r=0
dens1(1)=dens1(2);
% now make the integral of the density be equal to the number of particles


Npdens = 4*pi*sum(dens1(2:end).*r(2:end).^2)*dr;
if Npdens > 0
    scale = weight1*Np1/Npdens;
    dens1=dens1*scale;
end

% extrapolate the velocity beyond it's loaded range
[vmax,imax]=max(vr1);
i2=max(floor(.8*imax),1);
i1=max(floor(.6*imax),1);
p=polyfit(r(i1:i2),vr1(i1:i2),1);
vr1(imax:end)=polyval(p,r(imax:end));

% #####################  end particle type 1 #######################

% #####################  start particle type 2 #######################
dens2=0*r;vr2=0*r;
for n=1:Np2
    
    % find the cell a particle is in
    x = rp2(n)/dr;
    i = ceil(x+0.5)+1; % the index of the cell-centered grid point that is just
                 % above where the particle is
    i = min(i,Nr);
    f = (x-i+1.5); % fraction of the way the particle is from the lower edge
                 % of the cell to the upper edge.
                 
%    fprintf(' rp %g i %g r(i) %g  r(i-1) %g f %g \n',rp(n)/dr,i,r(i)/dr,r(i-1)/dr,f)
%    pause
                 
    % quadratic weighting, 2, linear weighting ,1 
    q = 1;

     f=1 - (r(i)^q-rp2(n)^q)/(r(i)^q-(r(i)-dr)^q);
    % cast density to the grid point below i and to i. Do something tricky
    % take care of the case that i=1, so i-1 is zero, which is not a grid
    % point
    
    dens2(i) = dens2(i)+f;
    vr2(i) = vr2(i)+f*vp_2(n);
    ibelow=i-1;
    %  ibelow=max(i-1,1);
    if ibelow>1
      dens2(ibelow)= dens2(ibelow) + 1-f; 
      vr2(ibelow)= vr2(ibelow) + (1-f)*vp_2(n);
    end
    

    
end


% convert these integer counting densities to densities by dividing by
% the cell volumes

vr2=vr2./(dens2+1e-12); % dens is a particle count at this stage
% load the ghost point
vr2(1)=-vr2(2);

dens2=dens2./Vol*weight2; % convert particle count to density

% the first few points are noisy - smooth them

ibad=2; % ibad is the last point to smooth. Use the next delta points to make a fit.
[~,dum2] = min( (r - rrms2).^2 ) ;
dum2 = ceil(dum2/3) +1;
delta=max(3,dum2); % adjust these if you don't like the way the density looks near r=0


x=r(ibad+1:ibad+delta).^2/dr^2; % density is a function of r^2
y=dens2(ibad+1:ibad+delta);

p=polyfit(x,y,1); % fit
for i=2:ibad+1
    dens2(i)=polyval(p,r(i).^2/dr^2);
end
% load the ghost point below r=0
dens2(1)=dens2(2);
% now make the integral of the density be equal to the number of particles

Npdens = 4*pi*sum(dens2(2:end).*r(2:end).^2)*dr;
if Npdens > 0
    scale = weight2*Np2/Npdens;
    dens2=dens2*scale;
end

% extrapolate the velocity beyond it's loaded range

[vmax,imax]=max(vr2);
i2=max(floor(.8*imax),1);
i1=max(floor(.6*imax),1);

% added by SDB 2/6/2018
i2 = max(i2, i1+1);

p=polyfit(r(i1:i2),vr2(i1:i2),1);
vr2(imax:end)=polyval(p,r(imax:end));


% #####################  end particle type 2 #######################

%smooth the vr1 and vr2 at the origin
[~,dumv] = min( (r - rrms1).^2 ) ;
dumv = ceil(dumv/3) +1;
delta=max(3,dumv); % adjust these if you don't like the way the density looks near r=0

qv = mean( vr1(dumv+2:dumv+6) ./ r(dumv+2:dumv+6) ) ;
vr1(1:dumv+1) = r(1:dumv+1) .* qv ;
qv2 = mean( vr2(dumv+2:dumv+6) ./ r(dumv+2:dumv+6) ) ;
vr2(1:dumv+1) = r(1:dumv+1) .* qv2 ;

% now compute the total density

dens = dens1 + dens2;

% if the density is zero it's going to cause trouble in dlogn/dr
% find all of the zero values and put them in the mask array

mask=dens==0;

% find the first point at which dens=0 (mask=0)



% now compute d(log(n))/dr

logn=log(dens+mask); % this will make log(n) be correct when mask=0
                     % and be zero when mask=1
                     
% find the first mask=1 point
Nr=length(r);
jfirst=0;
j=min(rrms1,rrms2)/2/dr; % don't start at the first grid point because late in time there might
                         % not be any particles near the origin
j=ceil(j);
j=max(j,2);
while jfirst==0
    j=j+1;
    
    if  j<Nr & mask(j)==1 % don't let this index run out the end of the grid
       jfirst=j;
    elseif j==Nr
        jfirst=Nr;
    end
      
end


% define N to be a few grid points earlier than jfirst
N=jfirst-6;

%    
% r(N)
% plot(r,logn)
% disp('1')
% pause
%         

% it should be a function of r^2; try a polynomial first

y=logn(1:N);
x=r(1:N).^2/r(end)^2;

% loop over the points and fit a straight line in x at each point
dlogndr=0*x;
stride=3;
for j=stride:N-1 % don't do the first point and the last point
    j1=j-stride;j2=j+stride;
    j1=max(j1,2);
    j2=min(j2,N);
    q=x(j1:j2)-x(j);
    p=polyfit(q,y(j1:j2),1);
    dlogndr(j) = p(1); % actually dlogn/dx, x=r^2, fixed later
end


%  fprintf('dum1 %g \t dum2 %g \t %g \t %g \n',dum1, dum2, length(r), length(dlogndr));



% smooth

% close all
% plot(dlogndr,'r-.')
% hold on

p=mean(dlogndr(stride:dum2));
dlogndr(1:dum2)=p;
dlogndr(N)=2*dlogndr(N-1)-dlogndr(N-2);



% dum2
% dlogndr(1:dum2)
% 
% plot(dlogndr,'b-.')
% disp('$$')
% pause(5)


for is=1:2
   dlogndr(stride:N-1)=.25*dlogndr(stride-1:N-2)+.5*dlogndr(stride:N-1)+.25*dlogndr(stride+1:N);
   p=mean(dlogndr(stride:dum2));
   dlogndr(1:dum2)=p;
   dlogndr(N)=2*dlogndr(N-1)-dlogndr(N-2);
end

% dum2
% figure(5);
% plot(r,vr1,'b-',r,vr2,'r-')
% pause(eps)

% chain rule, converting from d(log(nj))/dx to d(log(n))/dr
dlogndr(2:N) = 2*r(2:N)/r(end)^2.*dlogndr(2:N);
dlogndr(1) = -dlogndr(2);

% this only loads dlogndr up to where logn starts to be infinite/Nan
% extrapolate linearly the rest of the way, or at least clean up
% the end if there are particles in every cell

j1 = N-1-3*dum2;
j2 = N-1 ;
jmax = length(r);
%fprintf('N %g j1 %g j2 %g dum2 %g \n',N,j1,j2,dum2)


p1=polyfit(r(j1:j2), dlogndr(j1:j2), 1);
dlogndr(j1:jmax) = polyval(p1,r(j1:jmax));

% plot(dlogndr)

% smooth the tail

for is=1:2
   dlogndr(j1:jmax-1)=.25*dlogndr(j1-1:jmax-2)+.5*dlogndr(j1:jmax-1)+.25*dlogndr(j1+1:jmax);
   dlogndr(jmax)=2*dlogndr(jmax-1)-dlogndr(jmax-2);
end



return

    
     