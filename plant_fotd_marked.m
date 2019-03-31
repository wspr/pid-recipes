function [k0,L,T,tau] = plant_fotd(Gmodel,time,nt)

if nargin < 3
  nt = [0.05 0.63];
end

[y,t] = step(Gmodel,time);

%% Characterising FOTD part 1

k0 = dcgain(Gmodel);

fprintf('Points on the time trace are:\n')

N = numel(nt);
for ii = 1:N
  pp(ii) = find(y>nt(ii)*k0,1,'first');
  yy(ii) = y(pp(ii));
  tt(ii) = t(pp(ii));
  fprintf('    %2.0f%%: (%1.2f, %1.2f)\n',nt(ii)*100,tt(ii),yy(ii))
end

%% Characterising FOTD part 2

tangent_eq = @(x) y(pp(1))+(x-t(pp(1)))*( y(pp(2))-y(pp(1)) ) / ( t(pp(2)) - t(pp(1)) );
ya = tangent_eq(0);
t0 = t(pp(1)) - y(pp(1))/( y(pp(2))-y(pp(1)) )*( t(pp(2)) - t(pp(1)) );
tk = t(pp(1)) + (k0 - y(pp(1)))/( y(pp(2))-y(pp(1)) )*( t(pp(2)) - t(pp(1)) );

L = t0;
T  = tt(2)-t0;
T2 = tt(3)-t0;
a  = k0*L/T;
a2 = k0*L/T2;
tau = L/(L+T);

KK1 = -1/ya
KK2 = 1/a
(KK2-KK1)/KK1
(KK1-KK2)/KK2

%%

plot(t,y,'linewidth',3,'color',[0 0.8 0.2])


%plot(t([p1,p2]),y([p1,p2]),'r.-','markersize',40)

plot([0,t(end)], [0, 0], 'k-')   % origin line
plot([0,t(end)], [k0, k0], 'k-') % ss line

plot([0, L+T],[-a, k0],'-k') % diag "real" a
plot([0, L+T2],[-a2, k0],'-k') % diag "wrong" a
plot([0,tk], [ya, k0], 'k-') % diag "fake" a

%plot([0,t0,t0,tk], [ya,0, k0,k0], 'k-')
%plot([0, t0, tt(2), tk], k0*[1,1,1,1], 'k-','markersize',20)

%text(x0/2,k0,'~$L$~','interpreter','latex','edgecolor','black','backgroundcolor','white','horizontalalignment','center')
%text((x0+t2)/2,k0,'~$T$~','interpreter','latex','edgecolor','black','backgroundcolor','white','horizontalalignment','center')

axis tight
xlabel('Time, s')
ylabel('Amplitude')


box on

mark_point([tk,k0],'A')
mark_point([L,0],'D')

for ii = N:-1:2
  plot([0,tt(ii),tt(ii)], [yy(ii),yy(ii),k0], 'k-')
  mark_point([tt(ii),yy(ii)],'E'-ii)
end

mark_point([0,ya],'E',0.1,-0.02)
mark_point([0,-a2],'F',0.1,-0.07)
mark_point([0,-a],'G',0.1,-0.13)


end

function mark_point(p,s,xx,yy)

if nargin < 4
  xx = +0.07;
  yy = -0.07;
end

gg = 0.5*[1 1 1];

plot(p(1)+[0 xx],p(2)+[0 yy],'color',gg)
plot(p(1),p(2),'.','markersize',25,'color','black')

text(p(1)+xx,p(2)+yy,sprintf('%s: (%2.2f, %2.2f)',s,p(1),p(2)),...
  'userdata',sprintf('matlabfrag:\\smallbox{%s: $(%2.2f, %2.2f)$}',s,p(1),p(2)),...
...  'edgecolor',gg,'backgroundcolor','white',...
  'horizontalalignment','left','verticalalignment','middle')

end

