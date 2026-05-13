clear

%% simulating bacteria-phage dynamics
% phage population can be free or bound to bacterial hosts

r = 0.9/6; % Atf growth rate
K = 1e8; % Atf local carrying capacity
m = 0.02; % bacterial death rate
lf = 3e-7; % leaching factor per cell density
xp = 3e-4; % production rate
xc = 1e4; % agent consumption rate during leaching
cMIC = 5; % MIC for metal toxicity from NMC

b0 = 2e6; % initial bacterial density
c0 = 1; % initial NMC cathode material

dt = 0.001; % simulation time step
tf = 96; % total simulation time
trng = 0:dt:tf;
Nt = length(trng);
tIntro = [2 4 8 12 24 48];
NI = length(tIntro);

b = zeros(NI,Nt); % bacteria population
x = zeros(NI,Nt); % conc. of compound responsible for leaching
c = zeros(NI,Nt); % NMC cathode material
cf = zeros(NI,1); % NMC after 48 hours of leaching

b(:,1) = b0; % initial Atf density

for cI = 1:NI
    ct = 0;
    for t = trng(1:Nt-1)
        ct = ct+1;

        b(cI,ct+1) = b(cI,ct) + dt*(r*(t/(t+4))*(1-c(cI,ct)/cMIC)*(1 - 1/K*b(cI,ct)) - m)*b(cI,ct);
        c(cI,ct+1) = max(c(cI,ct) - dt*lf*x(cI,ct)*c(cI,ct),0);
        x(cI,ct+1) = max(x(cI,ct) + dt*xp*b(cI,ct) + xc*(c(cI,ct+1)-c(cI,ct)),0);
        if abs(t-tIntro(cI))<=dt
            c(cI,ct+1) = c0; % NMC at introduction
        end
        if abs(t-tIntro(cI)-48)<=dt
            cf(cI) = c(cI,ct+1); % NMC after 48 hours of leaching
        end
    end
end

newcolors = [0, 0, 0; 1, 0, 0; 0, 0, 1; 0, 0.8, 0; 1, 0.2, 1; 0.8, 0.4, 0.8];
figure
plot(trng, b)
xlim([0 96])
ylim([0 1e8])
colororder(newcolors)
legend('2 hr','4 hr','8 hr','12 hr','24 hr','48 hr')
xlabel('Time (hours)')
ylabel('Atf density (cells/ml)')

figure
plot(trng, c)
xlim([0 96])
ylim([0 1.1])
colororder(newcolors)
legend('2 hr','4 hr','8 hr','12 hr','24 hr','48 hr')
xlabel('Time (hours)')
ylabel('Fraction of initial NMC')

figure
plot(trng, x)
xlim([0 96])
% ylim([0 1.1])
colororder(newcolors)
legend('2 hr','4 hr','8 hr','12 hr','24 hr','48 hr')
xlabel('Time (hours)')
ylabel('Concentration of leaching agent')

figure
bp = bar(1:NI, 100*(1-cf),'FaceColor','flat');
bp.CData = newcolors;
set(gca,'XTickLabel',tIntro)
ylim([0 110])
xlabel('Time of NMC introduction (hours)')
ylabel('Leaching efficiency (%)')
