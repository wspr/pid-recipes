function [k0,L,T,a,tau] = plant_fotd(Gmodel,time,nt)

if nargin < 3
  nt = [0.05 0.63];
end

[y,t] = step(Gmodel,time);

%% Characterising FOTD part 1

k0 = dcgain(Gmodel);
p1 = find(y>nt(1)*k0,1,'first');
p2 = find(y>nt(2)*k0,1,'first');

y1 = y(p1);
y2 = y(p2);
t1 = t(p1);
t2 = t(p2);

fprintf('Points on the time trace are:\n')
fprintf('    %2.0f%%: (%1.2f, %1.2f)\n',nt(1)*100,t1,y1)
fprintf('   %2.0f%%: (%1.2f, %1.2f)\n',nt(2)*100,t2,y2)

%% Characterising FOTD part 2
%
% Doing it programmatically of course.
%
% We don't actually use the tangent equation, but I wrote it out originally
% to derive the other points around it.

% tangent_eq = @(x) y(p1)+(x-t(p1))*( y(p2)-y(p1) ) / ( t(p2) - t(p1) );
x0 = t(p1) - y(p1)/( y(p2)-y(p1) )*( t(p2) - t(p1) );
xk = t(p1) + (k0 - y(p1))/( y(p2)-y(p1) )*( t(p2) - t(p1) );

plot(t,y,'linewidth',3)

%plot(t([p1,p2]),y([p1,p2]),'r.-','markersize',40)

plot([0,x0,xk,t(end)], [0, 0,k0,k0], 'k-')
plot([x0,x0,xk], [0, k0,k0], 'k-')
plot([0,t2,t2], [y2,y2,k0], 'k-')
plot([0, x0, t2, xk], k0*[1,1,1,1], 'k-','markersize',20)

%text(x0/2,k0,'~$L$~','interpreter','latex','edgecolor','black','backgroundcolor','white','horizontalalignment','center')
%text((x0+t2)/2,k0,'~$T$~','interpreter','latex','edgecolor','black','backgroundcolor','white','horizontalalignment','center')

axis tight
xlabel('Time, s')
ylabel('Amplitude')

L = x0;
T = t2-x0;
a = k0*L/T;
tau = L/(L+T);
plot([L, L+T],[0, k0],'-k')

box on
