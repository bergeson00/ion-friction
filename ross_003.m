% ross_003.m -- code from Ross on 2/2/2018, edited by SDB to include
%   - self-collisions
%   - options for saving the histogram data
%   - on-the-fly vz data
%   - lots of different logLambda values

clear;
close all;

loglambdaArray = [3 1 0]; % --> the "real" log lambda is half of this value
% loglambdaArray = [3] ;

save_flag = 0 ;

NNmax = length(loglambdaArray) ;
for NN = 1 : NNmax

% global logn dens1

rng default

myConstants

% Plasma parameters for species 1
Te1_0 = 96 ;            % initial electron temperature in Kelvin 96 K == 100 cm-1
mi1   = mp * 40 ;       % the ion mass in kg
r01_0 = 2.99e-4 ;       % initial rms width of the plasma in meters
t01   = 0e-6 ;          % the time during which the Ca plasma expands before the Yb ionization
tau1  = sqrt( mi1 * r01_0^2 / kb / Te1_0 ) ;% characteristic plasma expansion time
Te1 = Te1_0 / ( 1 + t01^2 / tau1^2 );       % Te after initial expansion
r01 = r01_0 * sqrt( 1 + t01^2 / tau1^2 ) ;  % rms size after initial expansion
dvdr1 = t01 / tau1^2 / ( 1 + t01^2 / tau1^2 ) ;      % the coefficient of plasma expansion
Ti1 = 2 ; % initial ion temperature, species 1

% Plasma parameters for species 2
Te2   = 96 ;
mi2   = mp * 174 ;      % the ion mass in kg
r02   = 3.71e-4 ;       % initial rms width of the plasma in meters
% r02 = r01 ;
tau2  = sqrt( mi2 * r02^2 / kb / Te2 ) ; % characteristic plasma expansion time
dvdr2 = sqrt( kb * Te2 / mi1 ) / r02 ;   % the coefficient of plasma expansion
Ti2 = 2 ; % initial ion temperature, species 1

% Computed properties : mu, b0
reduced_mass = mi1 * mi2 / ( mi1 + mi2 ) ;
r0fac = sqrt(pi/2)*e^2/(4*pi*e0*.5*reduced_mass);

% Computational parameters
N1real = 2e5 ;
N2real = 3.82e5 ;
% N2real = N1real ;
loglambda = loglambdaArray(NN) ; 
% the number of simulation particles used in the calculation
n_particles_1 = 2e6 ;   
n_particles_2 = 2e6 ;
% calculate the simulation particle weights
weight1 = N1real/n_particles_1;
weight2 = N2real/n_particles_2;

% set the radial grid to be cell-centered with ghost points at both ends
% so that interp1 will be happy using our particles
num_of_cells  = 100 ;   % number of radial grid points
rmax = 15*r01 ;         % maximum size of the grid
dr   = rmax / num_of_cells ;
r    = -dr/2:dr:rmax+dr/2 ; % we need a ghost point beyond rmax for particles

vr1=0*r;vr2=vr1;  % define the radial fluid velocity arrays

tmax = 5e-6;            % maximum time

dtlimit = dr / sqrt(kb*Te1_0/min(mi1,mi2)) ;
fprintf('courant time step limit = %g \n',dtlimit);

% dt = input('  enter the time step: ');
dt = 0.1e-7 ;
num_of_steps = floor( tmax / dt );

t    = (dt:dt:tmax)' ;
rms_size_1  = zeros([num_of_steps,1]);
rms_size_2  = zeros([num_of_steps,1]);
rms_size_3  = zeros([num_of_steps,1]);


% load particles, randn uses the probability distribution
% exp(-x^2/2)/sqrt(2*pi)
x = randn([n_particles_1,1]).*r01 ;
y = randn([n_particles_1,1]).*r01 ;
z = randn([n_particles_1,1]).*r01 ;
rp_1 = sqrt(x.^2 + y.^2 + z.^2) ;
% Because the Ca plasma has expanded, we need to calculate the velocity at
% the expansion time. ( or not )
% aws = (3 / 4 / pi / 1.8e16 )^(0.3333) ;
% Tdih = 2/3 * e^2 / 4 / pi / e0 / aws / kb / 2.3 ;
% vdih1 = sqrt(kb * Tdih / mp / 40 ) ;
% vp_1 = dvdr1 .* rp_1 ...
%    + vdih1 .* sqrt(randn(size(rp_1)).^2 + randn(size(rp_1)).^2 + randn(size(rp_1)).^2) ;
vp_1 = dvdr1 .* rp_1 ;

x = randn([n_particles_2,1]).*r02 ;
y = randn([n_particles_2,1]).*r02 ;
z = randn([n_particles_2,1]).*r02 ;
rp_2 = sqrt(x.^2 + y.^2 + z.^2) ;
vp_2 = zeros([n_particles_2,1]);
% vdih2 = sqrt(kb * Tdih / mp / 174 ) ;
% vp_2 = dvdr2 .* rp_2 ...
%    + vdih2 .* sqrt(randn(size(rp_2)).^2 + randn(size(rp_2)).^2 + randn(size(rp_2)).^2) ;

% Now we can move the particles. We are calculating the field and then
% letting the field move the particles. We are going to try to add in a
% little friction between species. The friction should be proportional to
% the density-weighted velocity difference between species.
options   = optimset();
num = 1e12;

Te = ( N1real * Te1 + N2real * Te2 ) / ( N1real + N2real ) ;
Te_array = zeros([num_of_steps,1]);
rms_vp1  = zeros([num_of_steps,1]);
rms_vp2  = zeros([num_of_steps,1]);

dlogndr=0*r;
a1 = 0;
a2 = 0;

% find the initial energy


h2 = figure(2);
set(h2,'position',[2000 10 1700 1200])

h3 = figure(3) ;
set(h3, 'position', [50 100 600 400]);

h4 = figure(4) ;
set(h4, 'position', [950 100 900 900]) ;

qnum = 0 ;

for n = 1 : num_of_steps

    rp_1 = min( rp_1 + vp_1*dt, rmax ) ; % let the particles pile up at rmax
    rms_size_1(n) = sqrt( mean( rp_1.^2 ) ) / sqrt(3) ;
    
    rp_2 = min( rp_2 + vp_2*dt, rmax ) ; % let the particles pile up at rmax
    rms_size_2(n) = sqrt( mean( rp_2.^2 ) ) / sqrt(3) ;
    
   [dlogndr,dens1,dens2,vr1,vr2]=DlognDr3(n_particles_1,weight1,rp_1,vp_1,n_particles_2,weight2,rp_2,vp_2,r);
   

    % use energy conservation to estimate Te_now
    
    Te_now = ( ...
        N1real * ( Te1_0 - mi1 * mean( vp_1.^2 ) / ( 3 * kb ) ) + ...
        N2real * ( Te2   - mi2 * mean( vp_2.^2 ) / ( 3 * kb ) ) ...
        ) / ( N1real + N2real ) ;
    
    if Te_now < 0
        Te_now = Te_array(n-1);
    end
    
    Te_array(n) = Te_now ;
%     Te_now
        
    %---------------------------
    % To calculate the friction between species use the density of each
    % species and the fluid velocity of the species with which the particle
    % is colliding
    
    % Use F12 = m1*n2*pi*r0^2/2*abs(v1-v2)*(v1-v2)
    
    % Build a1=F12/m1*dt, a friction force array for each particle
    
    
    vother     = vr2;
    vrel       = vp_1 - interp1( r, vother, rp_1 ) ;
    vthrel     = sqrt( vrel.^2 + Ti1 * kb / mi1 + Ti2 * kb / mi2 ) ; % add the thermal sideways component, roughly, to the relative speed
    vrelself   = vp_1 - interp1( r, vr1, rp_1 ) ;
    vthrelself = sqrt( vrelself.^2 + Ti1 * kb / mi1 ) ;
    a1 = -loglambda * mi2 / (mi1+mi2) * interp1( r, dens2, rp_1 ).*r0fac^2./abs(vthrel).^3.*vrel*dt...
         -loglambda * 0.5 * interp1( r, dens1, rp_1 ).*r0fac^2./abs(vthrelself).^3.*vrelself*dt;

    vother     = vr1;
    vrel       = vp_2 - interp1( r, vother, rp_2) + 1e-6 ;
    vthrel     = sqrt( vrel.^2 + Ti1 * kb / mi1 + Ti2 * kb / mi2 ) ;
    vrelself   = vp_2 - interp1( r, vr2, rp_2 ) + 1e-6 ;
    vthrelself = sqrt( vrelself.^2 + Ti2 * kb / mi2 ) ;
    a2 = -loglambda * mi1 / (mi1+mi2) * interp1( r, dens1, rp_2 ).*r0fac^2./abs(vthrel).^3.*vrel*dt...
         -loglambda * 0.5 * interp1( r, dens2, rp_2 ).*r0fac^2./abs(vthrelself).^3.*vrelself*dt;
% 
%     a1 = a1 * 0 ;
%     a2 = a2 * 0 ;
    % do the particle velocity update
    
    Kbefore = sum(.5*mi1*vp_1.^2+.5*mi2*vp_2.^2);
   
    vp_1 = vp_1 + a1;
    vp_2 = vp_2 + a2;
    
    Kafter = sum(.5*mi1*vp_1.^2+.5*mi2*vp_2.^2);
    
    % now check on how badly we messed up energy conservation for the ions
    
    % dE = (Kafter-Kbefore)*2/(Kafter+Kbefore)
    
    
    
    
    
    
    
    
    %-------------------------------------
    
    
    a1 = -kb*Te_now/mi1*dlogndr * dt;
    a2 = -kb*Te_now/mi2*dlogndr * dt;
    
    % let the particles change their velocities because of the ambipolar
    % field
    vp_1 = vp_1 + interp1(r, a1,rp_1,'linear');
    vp_2 = vp_2 + interp1(r, a2,rp_2,'linear');
    
    rms_vp1(n) = sqrt( 0.3333 * mean( vp_1.^2 ) ) ;
    rms_vp2(n) = sqrt( 0.3333 * mean( vp_2.^2 ) ) ;
    
    if ( mod(n,20) == 0 || n == 1 )
        figure(1)
        subplot(211)
        yexact1=weight1*length(rp_1)/(2*pi*rms_size_1(n)^2)^1.5*exp(-.5*r.^2./rms_size_1(n)^2 );
        yexact2=weight2*length(rp_2)/(2*pi*rms_size_2(n)^2)^1.5*exp(-.5*r.^2./rms_size_2(n)^2 );
        plot( r, dens1,'k', r,dens2, 'b',r,yexact1,'r--')
        title(num2str( n / num_of_steps ) );
        subplot(212)
        plot(r,dlogndr);
        
        figure(2)
        vhist = linspace(0,500,501);
        [y,x]=histcounts(vp_1,vhist);
        xplot = (x(2:end) + x(1:end-1))./2 ;
        plot(xplot,y)
        [a,b]=max(y);
%         hold on;
%         plot( xplot(2*b) .* [1 1], [0 1].*a,'k');
%         hold off;
%         xin = xplot(1:2*b) ;
%         yin = y(1:2*b) ;
        xin = xplot(1:end) ;
        yin = y(1:end) ;
        options = optimset() ;
        w = [ x(b) a/x(b)^2 ] ;
        W = fminsearch(@expRfitOne,w,options,xin,yin) ;
        yfit = expRone(W,xplot);
        hold on;
        plot(xplot,yfit);
        hold off;
        qnum = qnum+1;
        time_array(qnum) = t(n) ;
        w1array(qnum) = abs( W(1) ) ;
        w2array(qnum) = abs( W(2) ) ;
        sum(dens1.*r.^2.* dr)
        
        figure(3)
        plot(rp_1,vp_1,'.', rp_2,vp_2,'.')
        set(gca, 'xlim',[0 15e-4], 'ylim', [0 200] );
        
        figure(4)
        delta_vz = 2.5 ;
        max_vz = ceil( max(vp_1)/delta_vz + 50) * delta_vz ;
        vz_bins = -max_vz : delta_vz : max_vz ;
        vz = vz_bins .* 0 ;
        num_bins = length( vz_bins ) ;
        num_angles = 200 ;
        for nn = 1 : n_particles_1
            vz_dum = ( 1 - 2 .* rand( [ num_angles , 1 ] ) ) .*  vp_1(nn) ;
            for mm = 1 : num_angles 
                bb = round( vz_dum(mm) / delta_vz + (num_bins+1)/2 ) ;
                vz(bb) = vz(bb) + 1 ;
            end
        end
        lorentz = 1 ./ (1 + (vz_bins ./ 4.4).^2 ) ;
        vz_final = conv(lorentz,vz,'same') ;
%         vz_final = vz ;
        plot(vz_bins, vz_final) ;
        vz_final_rms(qnum) = sqrt( sum( vz_bins.^2 .* vz_final ) / sum( vz_final ) ) ;
        w = [ vz_final(qnum) max(vz_final)  ] ;
        options = optimset();
        W = fminsearch(@expFunctionFit,w, options, vz_bins, vz_final) ;
        hold on;
        plot(vz_bins, expFunction(W,vz_bins));
        hold off;
        gausfit(qnum) = abs( W(1) ) ;
        pause(0.5) ;
        
        fname0 = ['vhist_20180328_0_', num2str(loglambda),'_',num2str(time_array(qnum)),'.mat'];
        if save_flag == 1
            save(fname0,'xin','yin','rp_1','vp_1','rp_2','vp_2','vz_bins','vz_final','vz');
        end
        pause(eps);
    end
    
%     fprintf('%g \t %g \t %g \t %g \t %g \t %g \n \n',Te_now, mean(a1), mean(a2), mean(vp_1), mean(vp_2), Eperp_time(n))
    
end

figure('position',[50 50 1200 600]) ;

t1 = t + t01 ;
subplot(121)
title('species 1');
plot( t, rms_size_1, 'b-', t, r01_0 .* sqrt(1 + t1.^2 ./ tau1^2), 'r--', ...
    t, rms_size_2, 'k-', t, r02   .* sqrt(1 +  t.^2 ./ tau2^2), 'g--', 'linewidth', 2 )
xlabel('time (s)');
ylabel('r_{rms} (m)');
legend('ca sim','ca mod','yb sim','yb mod','location','northwest');

subplot(122)
plot(t, rms_vp1, 'b-', t, t1 .* r01_0 ./ tau1^2 ./ sqrt(1 + t1.^2 ./ tau1^2), 'r--', ...
     t, rms_vp2, 'k-', t, t  .* r02   ./ tau2^2 ./ sqrt(1 +  t.^2 ./ tau2^2), 'g--', 'linewidth', 2 )
% hold on;
% plot(time_array, gausfit,'co','markerfacecolor','c');
% hold off;
xlabel('time (s)');
ylabel('vrms (m/s)');
legend('ca sim','ca mod','yb sim','yb mod','location','northwest');


figure;
subplot(211); 
plot(time_array,w1array, time_array, w2array)

subplot(212);
plot(time_array,w1array./w2array, time_array,w2array./w1array)
set(gca, 'ylim',[0 1]);

end
