%%% 3.3.1.2

model = createpde;

% les valeurs des rayons r1 et r2 peuvent être exprimées en unité réduite
r1 = 1;
r2 = 0.2;

% Bords de la figure
C1 = [1;0;0;r1];
C2 = [1;0;0;r2];

% Définition des paramètres de decsg
geom = [C1,C2];
sf = 'C1-C2';
ns = char('C1', 'C2')';

% Utilisation de decsg et geometryFromEdges
g = decsg(geom, sf, ns);
geometryFromEdges(model,g);

% tracé de la description géométrique obtenue
figure,pdegplot(model,'EdgeLabels','on'),drawnow

%%% 3.3.1.3

% paramètres de specifyCoefficient
m = 0; d = 0; c = 1; a = 0; f = 1;
specifyCoefficients(model,'m',m,'d',d,'c',c,'a',a,'f',f);

%%% 3.3.1.4

% les valeurs des potentiel V1 et V2 sont exprimés en volts
V1=10;V2=20;

% appel à la fonction applyBoundaryCondition, avec les potentiels imposés
applyBoundaryCondition(model,'dirichlet','Edge',1:4,'u',V1);
applyBoundaryCondition(model,'dirichlet','Edge',5:8,'u',V2);

%%% 3.3.1.5

% différents tests
Hmax=0.1;
generateMesh(model,'Hmax',Hmax,'GeometricOrder','linear');
%generateMesh(model,'Hmax',Hmax,'GeometricOrder','quadratic');
% le tableau [1 2 3 4] contient les numéros de segment représentant
% le contour du disque de rayon r1
%generateMesh(model,'Hmax',Hmax,'Hedge',{[1 2 3 4],Hmax/5});

figure
pdemesh(model)

%%% 3.3.2.1

results = solvepde(model);

%%% 3.3.2.2

u_n = results.NodalSolution;
disp(length(u_n))
figure,pdeplot(model,'XYData',u_n,'ZData',u_n)
colormap(jet)
drawnow

% Comparaison
x=results.Mesh.Nodes(1,:)';
y=results.Mesh.Nodes(2,:)';
r=sqrt(x.^2+y.^2);

V0 = (V1-V2)/log(r1/r2);
r0 = r2*10^(-V2/V0);
u_a=V0*(log(r)-(log(r2)-V2/V0));

erreur=abs(u_n-u_a);
disp(max(erreur));
figure,pdeplot(model,'XYData',erreur,'ZData',erreur)
colormap(jet)
drawnow

%%% 3.4.2

Hmax=logspace(-0.5,-1.5,100);
tab_erreur=zeros(length(Hmax),1);
tab_numel=zeros(length(Hmax),1);
for n=1:length(Hmax)
 generateMesh(model,'Hmax',Hmax(n),'Hedge',{[1 2 3 4],Hmax(n)/5});
 results = solvepde(model);
 u_n = results.NodalSolution;
 disp(length(u_n))
 tab_numel(n)=length(u_n);
 x=results.Mesh.Nodes(1,:)';
 y=results.Mesh.Nodes(2,:)';
 r=sqrt(x.^2+y.^2);
 u_a=V0*(log(r)-log(r0));
 erreur=abs(u_n-u_a);
 tab_erreur(n)=max(erreur);
end
figure,semilogx(Hmax,20*log10(tab_erreur),Hmax,40*log10(Hmax)-20)
grid on
figure,semilogx(Hmax,20*log10(tab_numel),Hmax,-40*log10(Hmax)+26)
grid on
