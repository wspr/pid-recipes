%% PID Example Tuning Laws

close all
clear all
s = tf('s');

k = 0.8;
L = 0.1;
T = 0.3;
a = k*L/T;
N = 20;
GG = k/(T*s+1)*exp(-s*L);

figure(99); clf; hold on
[kk,LL,TT] = plant_fotd(GG,linspace(0,2,1000),[0.01 0.63]);
ylim([0 1])

disp('Actual vs identified:')
disp([k kk; L LL; T TT])

md = {'ZN','CHR_DR_0','CHR_DR_20','CHR_SP_0','CHR_SP_20','CC','WJC','KT1','KT2'};

H = pid_recipe(md,kk,LL,TT,N);
CL = cell(size(H));

for ii = 1:numel(md)
  HH = H{ii};
  if iscell(HH), HH = HH{2}; end
  CL{ii} = feedback(GG*HH,+1);
end

figure(11); clf; hold on
step(GG)
for ii = 1:numel(md)
  step(CL{ii});
end
legend(['OL',md])

colourplot