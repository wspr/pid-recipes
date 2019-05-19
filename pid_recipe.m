function [H,g] = pid_tune(method,k,L,T,N)

if ~iscell(method)
  method = {method};
end

if nargin < 5
  N = 20;
end

M = numel(method);
H = cell(1,M);

for ii = 1:M
  [H{ii}, g.(method{ii})] = method_pid(method{ii},k,L,T,N);
end

if M == 1
  H = H{1};
end

end

%% %%%%%%%%%%%%

function [H,g] = method_pid(M,k,L,T,N)

a   = k*L/T;
tau = L/(L+T);
R   = a/L;

switch M
  case 'ZN'
    Kp = 1.2/a;
    Ti = 2*L;
    Td = L/2;
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td);
    H = pidstd(Kp,Ti,Td,N);
    
  case 'CHR_DR_0'
    Kp = 0.95/a;
    Ti = 2.4*L;
    Td = 0.42*L;
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td);
    H = pidstd(Kp,Ti,Td,N);
    
  case 'CHR_DR_20'
    Kp = 1.2/a;
    Ti = 2*T;
    Td = 0.42*L;
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td);
    H = pidstd(Kp,Ti,Td,N);
    
  case 'CHR_SP_0'
    Kp = 0.6/a;
    Ti = T;
    Td = 0.5*L;
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td);
    H = pidstd(Kp,Ti,Td,N);
    
  case 'CHR_SP_20'
    Kp = 0.95/a;
    Ti = 1.4*T;
    Td = 0.47*L;
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td);
    H = pidstd(Kp,Ti,Td,N);
    
  case 'CC'
    Kp = 1.35/a*(1+0.18*tau/(1-tau));
    Ti = (2.5-2*tau)/(1-0.39*tau)*L;
    Td = 0.37*(1-tau)/(1-0.81*tau)*L;
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td);
    H = pidstd(Kp,Ti,Td,N);
    
  case 'WJC'
    Kp = (0.7303+0.5307*T/L)*(T+0.5*L)/k/(L+T);
    Ti = T+0.5*L;
    Td = L*T/(2*T+L);
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td);
    H = pidstd(Kp,Ti,Td,N);
    
  case 'FCL' % "Simplified IMC-PID tuning rules"
    
    Kp = 1/(2*a);
    if T/L > 3
      Ti = 5*L;
    else
      Ti = T;
    end
    Td = 0.5*L;
    
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td);
    H = pidstd(Kp,Ti,Td,N);
    
  case 'KT1'
    % Ms = 1.4
    aKp = 3.80*exp(-8.40*tau+7.3*tau^2);
    TiL = 5.20*exp(-2.50*tau-1.4*tau^2);
    TdL = 0.89*exp(-0.37*tau-4.1*tau^2);
    
    Kp = aKp/a;
    Ti = TiL*L;
    Td = TdL*L;
    
    b  = 0.22*exp(0.65*tau+0.051*tau^2);
    
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td,'b',b);
    H = pidstd2(Kp,Ti,Td,N,b,0);

  case 'KT2'
    % Ms = 1.4
    aKp = 8.4*exp(-9.6*tau+9.8*tau^2);
    TiL = 3.2*exp(-1.5*tau-0.93*tau^2);
    TdL = 0.86*exp(-1.9*tau-0.44*tau^2);
    
    Kp = aKp/a;
    Ti = TiL*L;
    Td = TdL*L;
    
    b  = 0.22*exp(0.65*tau+0.051*tau^2);
    
    g = struct('Kp',Kp,'Ti',Ti,'Td',Td,'b',b);
    H = pidstd2(Kp,Ti,Td,N,b,0);

  otherwise
    error('Method unknown')

end

end

