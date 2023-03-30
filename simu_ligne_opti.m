% Programme à compléter pour la simulation numérique
% d'une ligne de trasmission.
function simu_ligne_opti(V,Zc,N,D,dt,methode)

% Les paramètres sont ceux donnés par défaut
% si on tape dans la fenêtre de commande Matlab :
% >> simu_ligne_a_modifier
    % Il est possible aussi de passer d'autres paramètres
    % en arguments du programme :
    % >> simu_ligne_a_modifier(2.4e8,120,2000,500,1e-9,@myode23)
if nargin==0
    V=2.3e8;Zc=100;
    N=2000;D=500;dt=10^(-9);
    methode=@myode45;
end
    
% Fermeture de toutes les fenêtres graphiques  
% préalablement créées
close all

% Paramètres de fonctionnement.
C=1/(V*Zc);L=Zc/V;
T=4*D/V;
Re=sqrt(L/C);Rs=sqrt(L/C);dx=D/N;
disp(['V*dt/dx = ' num2str(V*dt/dx)])
Ldx=L*dx;Cdx=C*dx;

% La variable "tliste" sera un vecteur ligne contenant
% 0, dt,2*dt, ..., n*dt jusqu'à ce que n*dt dépasse T
tliste=0:dt:T;

% Vecteur d'état initial
y0=zeros(2*N,1);

% Lancement de l'approximation numérique
[~,y]=methode(@(t,y) ode_ligne(t,y,Ldx,Cdx,N,Re,Rs),tliste,y0);

% Visualisation des résultats
% (On ne représente ici que Nv=1000 éléments régulièrement répartis
% de la solution pour accélérer l'affichage)
Nv=2000;Nt=length(tliste);
for n=1:Nv
    plot((1:N)*dx,y(round(1+(Nt-1)*(n-1)/(Nv-1)),1:N))
    axis([dx N*dx -2 2])
    drawnow
end

% La fonction décrivant l'équation dynamique
% * t est l'instant 
% * y est l'état à instant t
% * yp est la dérivée du vecteur d'état à ce même instant
function yp=ode_ligne(t,y,Ldx,Cdx,N,Re,Rs)
% Les paramètres qui suivent t et y sont définis dans le programme principal

% La partie suivante est à compléter
yp=zeros(2*N,1);
%-------------- partie à compléter --------------
% Équivalent à remplir la matrice A et B mais avec moins d'opération
% élémentaires nécessaires.
for i=1:2*N
    if i < N
        yp(i) = (-y(N+i+1)+y(N+i))/Cdx;
    elseif i == N
        yp(i) = (-(y(N)/Rs)+y(2*N))/Cdx;
    elseif i == N+1
        yp(i) = (-y(1)+entree(t)-Re*y(N+1))/Ldx;
    elseif i > N+1
        % index out of range result in 0
        yp(i) = (y(i-N-1)-y(i-N))/Ldx;
    end
end
%------------------------------------------------
end
% La fonction décrivant la source de tension
% à l'instant t
function src=entree(t)
    Te=200e-9;
    if t<Te
        src=(1-cos(2*pi*t/Te))/2;
    else
        src=0;
    end
end
end