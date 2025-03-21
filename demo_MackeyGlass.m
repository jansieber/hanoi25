%% Minimalistic DDE-Biftool demo Mackey-Glass Equation 
%
% <http://www.scholarpedia.org/article/Mackey-Glass_equation>
%
% The Mackey-Glass equation is given by
% 
% $$x'(t)=\beta \frac{x(t-\tau)}{1+x(t-\tau)^n}-\gamma x(t)$$
% 
% Parameters are (in this order) |beta|, |n|, |tau| (|gamma| is not part of
% parameter vector).
% See manual for extensive comments.
%% load DDE-Biftool into path
clear
ddebiftool_path([pwd(),'/ddebiftool_25_03_07']);
format compact
clr=lines();
%% Define problem: rhs f and its directional derivative df
[ib,in,itau,ig]=deal(1,2,3,4); % define parameter indices
f= @(x,b,n,g)b.*x(2,:)./(1+x(2,:).^n)-g.*x(1,:); % r.h.s., all vectorized for speed-up
df1=@(x,n,dx,dn)1./(x.^n+1).^2.*(dx.*(1+x.^n-n.*x.^n)-dn.*x.*x.^n.*log(x));
df=@(x,b,n,g,dx,db,dn,dg)db.*x(2,:)./(1+x(2,:).^n)+b.*df1(x(2,:),n,dx(2,:),dn)-dg.*x(1,:)-g.*dx(1,:);
%% Convert r.h.s. to DDE-Biftool f(x,p) format,
% note that the parameters are row vectors for historic reasons
funcs=set_funcs('sys_rhs',{1,[ib,in,ig],f},'sys_tau',@()itau,...
    'sys_dirderi',{df},'x_vectorized',true,'p_vectorized',true); % set up problem
%% Set initial guess for state and parameter
% set boundaries, max stepsize, and perform initial correction.
[b,n,tau,g]=deal(2, 10, 0,1);   % initial parameters
x0=(b/g-1)^(1/n);               % initial non-trivial equilibrium
bounds={'max_bound',[itau,2;ib,5],'max_step',[0,0.3]};
[eqbr,suc]=SetupStst(funcs,'x',x0,'parameter',[b,n,tau,g],...
    'step',0.1,'contpar',itau,bounds{:})
%% Compute, find stability and bifurcations of non-trivial equilibria 
disp('Nontrivial equilibria');
figure(1);clf;ax1=gca;xlabel(ax1,'tau');ylabel(ax1,'x (eq)');set(ax1,'FontSize',16);
eqbr=br_contn(funcs,eqbr,10,'ax',ax1);
[eqbr_wbif,nunst_eqs,eqbifs,ihopf]=MonitorChange(funcs,eqbr);
fprintf('Hopf bifurcation near point %d, tau=%g\n',ihopf(1),eqbr_wbif.point(ihopf(1)).parameter(itau));
%% Illustrate eigenvalues of example point
pt=eqbr_wbif.point(ihopf);
[A,tau]=dde_stst_linearize_rhs(funcs,pt);
pt.stability=dde_stst_eig_cheb(A,tau,'max_number_of_eigenvalues',50,'min_number_of_eigenvalues',50);
figure(2);clf;ax2a=gca;
plot(ax2a,real(pt.stability.l0),imag(pt.stability.l0),'o','LineWidth',2);
grid(ax2a,'on');xlabel(ax2a,'Re');ylabel(ax2a,'Im');
set(ax2a,'FontSize',16,'LineWidth',1,'xlim',[-10,1]);
xline(ax2a,0,'LineWidth',1)
%% Plot initial one-parameter bifuration diagram
% Note that stability boundaries are not accurate after br_contn and
% br_stable as these functions do not refine the branch. Use MonitorChange or
% LocateSpecialPoints for refinement.
figure(5);clf;ax5=gca;hold(ax5,'on');
Plot2dBranch(eqbr_wbif,'ax',ax5);
plot(eqbr_wbif.point(ihopf(1)).parameter(itau),eqbr_wbif.point(ihopf(1)).x,'ko',...
    'DisplayName','Hopf','MarkerFaceColor','k')
legend('Location','southeast')
xlabel(ax5,'tau');ylabel(ax5,'x (eq)');set(ax5,'FontSize',16);
%% Continue Hopf bifurcation in two parameters b and tau
% starting from point where stability change was detected
% (eqbr.point(ihopf)).
[hopfbr,suc]=SetupHopf(funcs,eqbr_wbif,ihopf(1),...
    'contpar',[ib,itau],'dir',itau,'step',1e-1)
figure(2);clf;ax2=gca; xlabel('b');ylabel('tau');set(ax2,'FontSize',16);
hopfbr=br_contn(funcs,hopfbr,30,'ax',ax2);
hopfbr=br_rvers(hopfbr);
hopfbr=br_contn(funcs,hopfbr,30,'ax',ax2);
[hopfbr,nunst_h]=br_stabl(funcs,hopfbr);
%% Re-plot Hopf bifurcation in 2 parameter plane
figure(3);clf;ax3=gca;
Plot2dBranch(hopfbr,'ax',ax3);
xlabel(ax3,'b');ylabel(ax3,'tau');set(ax3,'FontSize',16);
%% Branch off toward periodic orbits
[per_orb,suc]=SetupPsol(funcs,eqbr_wbif,ihopf,'intervals',20,'degree',4,'max_step',[itau,0.5],...
    'plot_measure',{@(p)p.parameter(itau),@(p)max(p.profile)});
%% Plot initial periodic orbit to illustrate discretization
figure(4);clf;ax4=gca;
plot(ax4,per_orb.point(2).mesh,per_orb.point(2).profile,'o-',...
    per_orb.point(1).mesh,per_orb.point(1).profile,'o-');
xlabel(ax4,'t/period');ylabel(ax4,'x(t)');set(ax4,'FontSize',16);
%% Continuation with "live plot" of periodic orbits
% Floquet multipliers are determined in call to br_stabl. The number of
% unstable Floquet multipliers, nunst_per, indicates that stability is lost
% for increasing tau.
figure(1);ylabel(ax1,'x(eq), max(x(t))')
per_orb=br_contn(funcs,per_orb,60,'ax',ax1);
[per_orb_wbif,nunst_per,perbifs,ipd]=MonitorChange(funcs,per_orb);
%% plot all profiles, highlighting first period doubling
figure(4);clf;ax4=gca;hold(ax4,'on');
for i=1:length(per_orb_wbif.point)
    pt=per_orb_wbif.point(i);
    plot(ax4,pt.mesh,pt.profile,'r-');
end
pt=per_orb_wbif.point(ipd);
plot(ax4,pt.mesh,pt.profile,'k-','linewidth',3);
xlabel(ax4,'t/period');ylabel(ax4,'x(t)');set(ax4,'FontSize',16);
%% Illustrate Floquet multipliers of example periodic orbit
ptper=per_orb_wbif.point(ipd);
ptper.stability=dde_psol_eig(funcs,ptper,per_orb_wbif.method.stability,...
    'geteigenfuncs',true,'eigmatrix','sparse');
figure(6);clf;tiledlayout(1,2);nexttile; ax6a=gca;hold(ax6a,'on');
plot(ax6a,real(ptper.stability.mu),imag(ptper.stability.mu),'o','LineWidth',2);
s=linspace(0,2*pi,100);
plot(ax6a,cos(s),sin(s),'k-','LineWidth',2);
grid(ax6a,'on');xlabel(ax6a,'Re');ylabel(ax6a,'Im');
set(ax6a,'FontSize',16,'LineWidth',1,'xlim',[-1.1,1.1],'DataAspectRatio',[1,1,1],...
    'PlotBoxAspectRatio',[1,1,1]);
nexttile;ax6b=gca;hold(ax6b,'on');grid(ax6b,'on');
ef=ptper.stability.eigenfuncs(1:2);
plot(ax6b,ef(1).mesh,[ef(1).profile;ef(2).profile],'linewidth',2)
xlabel(ax6b,'t/period');ylabel(ax6b,'$\delta_x(t)$ (eigenfunction)','Interpreter','latex');set(ax4,'FontSize',16);
set(ax6b,'FontSize',16,'LineWidth',1);
%% Check that first stability change is a period doubling
per_orb_wbif.point(ipd).stability.mu(1:3)
%% Add periodic orbits to one-parameter bifuration diagram
figure(5);
Plot2dBranch(per_orb_wbif,'ax',ax5);
plot(ax5,per_orb_wbif.point(ipd).parameter(itau),max(per_orb_wbif.point(ipd).profile),'o','color',clr(1,:),...
    'DisplayName','PD1','MarkerFaceColor',clr(1,:));
%% Find period doubling bifurcations in two parameters
[pdfuncs,pdbr1,suc]=SetupPeriodDoubling(funcs,per_orb_wbif,ipd,...
    'contpar',[ib,itau],'dir',ib,'step',1e-1,bounds{:},'plot_measure',[])
%% Continuation
figure(2);
pdbr1=br_contn(pdfuncs,pdbr1,30,'ax',ax2);
pdbr1=br_rvers(pdbr1);
pdbr1=br_contn(pdfuncs,pdbr1,30,'ax',ax2);
[pdbr1,nunst_pd]=br_stabl(pdfuncs,pdbr1);
%% Branch off at period doubling 
% (Solutions at far end get inaccurate.)
[per2,suc]=DoublePsol(funcs,per_orb_wbif,ipd);
figure(1);
per2=br_contn(funcs,per2,60,'ax',ax1);
[per_orb2_wbif,nunst_p2,perbifs2,ipd2]=MonitorChange(funcs,per2);
%% Plot final one-parameter diagram
% Note that stability boundaries are not accurate as we have not refined
% the branch there. Use MonitorChange or LocateSpecialPoints for
% refinement.)
figure(5);
Plot2dBranch(per_orb2_wbif,'ax',ax5);
plot(per_orb2_wbif.point(ipd2(1)).parameter(itau),max(per_orb2_wbif.point(ipd2(1)).profile),'o','color',clr(4,:),...
    'DisplayName','PD2','MarkerFaceColor',clr(4,:));
legend(ax5,'Location','southeast');
%% Find secondary period doubling and track in two parameters
[pd2funcs,pdbr2,suc]=SetupPeriodDoubling(funcs,per_orb2_wbif,ipd2,...
    'contpar',[ib,itau],'dir',ib,'step',1e-1,'plot_measure',[],bounds{:});
figure(2);
pdbr2=br_contn(pd2funcs,pdbr2,30,'ax',ax2);
pdbr2=br_rvers(pdbr2);
pdbr2=br_contn(pd2funcs,pdbr2,30,'ax',ax2);
[pdbr2,nunst_pd2]=br_stabl(pd2funcs,pdbr2);
%% Two-parameter bifurcation plot cleaned-up
% pass on pdfuncs for period doublings to help identify type of bifurcation
% (without this hint, pdbr1 and pdbr2 look just like periodic orbits)
figure(3);hold(ax3,'on');
Plot2dBranch(pdbr1,'funcs',pd2funcs,'DisplayName','PD1');
Plot2dBranch(pdbr2,'funcs',pd2funcs,'DisplayName','PD2');
%% plot all profiles of period-2 orbits, highlighting period doubling
figure(4);clf;ax4=gca;hold(ax4,'on');
for i=1:length(per_orb2_wbif.point)
    pt=per_orb2_wbif.point(i);
    plot3(ax4,pt.mesh,pt.parameter(itau)+0*pt.mesh,pt.profile,'r-');
end
pt=per_orb2_wbif.point(ipd2);
plot3(ax4,pt.mesh,pt.parameter(itau)+0*pt.mesh,pt.profile,'k-','linewidth',3);
xlabel(ax4,'t/period');ylabel(ax4,'tau');zlabel(ax4,'x(t)');set(ax4,'FontSize',16);
grid(ax4,'on');view(ax4,[5,30])