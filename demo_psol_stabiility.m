%% Demo for linear stability of periodically forced problem
clear
ddebiftool_path([pwd(),'/ddebiftool_25_03_07']);
zeta = 0.01; % zeta
delta = 1500; % delta or: 150
epsilon = 150; % epsilon
tau = 2*pi; % Time delay tau
b = 0.5;
T = 2*pi;
%% time dependent matrices nx x nx x nt
A=@(t)reshape([...
    0*t;...
    -delta-epsilon.*cos(2*pi./T.*t);...
    1+0*t;...
    -2*zeta+0*t;...
    ...
    0*t;...
    b+0*t;...
    0*t;...
    0*t],...
    2,2,2,numel(t));
delays=[0,tau];
degree=200;
intervals=1;
tic;
stability=dde_period_floquet(A,delays,T,...
    'degree',degree,'intervals',intervals,...
    'geteigenfuncs',true,'max_number_of_eigenvalues',5,'eigmatrix','sparse')
tspent=toc
%% check residual of linear DDE
tfine=linspace(0,1,10000);
err=dde_floquet_residual(stability.mu,stability.eigenfuncs,tfine,A,delays)
%% Plot eigenfunctions
lw={'linewidth',2};
figure(1);clf;ax=gca;hold(ax,'on');
profile=dde_coll_eva(stability.eigenfuncs(1),tfine);
plot(ax,tfine*T,real(profile),'-',lw{:});
set(ax,'ColorOrderIndex',1)
plot(ax,tfine*T,imag(profile),'--',lw{:});
set(gca,lw{:},'Fontsize',16)
ltx={'interpreter','latex'};
xlabel('time',ltx{:});
ylabel('Re $z(t)$, Im $z(t)$',ltx{:});
%%
exportgraphics(figure(1),'DelayedMathieuFlqouet1.pdf')