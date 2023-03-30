% La fonction d'intégration pour rk23
% * Le paramètre f est une fonction, par exemple
% f=@ode_ligne
% * Le paramètre t est ici un vecteur contenant
% les instants auxquels on souhaite calculer la solution
% * Le paramètre y0 est l'état initial
% * Le tableau y en résultat comprend N lignes
% où N est le nombre d'instants sur lesquels
% la solution est estimée. Chaque ligne de ce
% tableau est l'état à l'instant correspondant
function [t,y]=myode23(f,t,y0)
N=length(t);
y=zeros(N,length(y0));
k1=zeros(1,length(y0));k2=k1;
y(1,:)=y0(:);
for n=2:N
    yprec=y(n-1,:);tprec=t(n-1);
    h=t(n)-tprec;
    k1(:)=f(tprec,yprec);
    k2(:)=f(tprec+h,yprec+h*k1);
    y(n,:)=yprec+(h/2)*(k1+k2);
end