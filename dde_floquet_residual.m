function err=dde_floquet_residual(mu,ef,tfine,Afun,delays)
[nev,nd,nt]=deal(length(mu),length(delays),length(tfine));
period=ef(1).period;
tdel=(tfine(ones(nd,1),:)-repmat(delays(:)/period,1,nt))';
tdel=tdel(:).';
At=permute(Afun(tfine*period),[1,2,4,3]);
err=NaN(1,nev);
for i=1:nev
    err(i)=loc_floquet_residual(mu(i),ef(i),At,tdel,tfine,period);
end
end
function err=loc_floquet_residual(mu,ef,At,tdel,tfine,period)
[nd,nt,nx]=deal(size(At,4),length(tfine),size(ef.profile,1));
mupow=repmat(mu.^floor(tdel),nx,1);
efx=dde_coll_eva(ef,mod(tdel,1));
efdxdt=dde_coll_eva(ef,tfine,'diff',1)/period;
exfmu=efx.*mupow;
Aexfmu=pagemtimes(At,reshape(exfmu,nx,1,nt,nd));
sumAexfmu=reshape(sum(Aexfmu,4),nx,nt);
res=[sumAexfmu-efdxdt,ef.profile(:,1)*mu-ef.profile(:,end)];
err=norm(res(:),'inf');
end
