clear all;
close all;

E_mem = 3.35e9;                % Young modulus of membrane
nu_mem = 0.33;                % Poisson ratio of membrane
rho_mem = 900;                % mass density of membrane
G_mem = E_mem/(2*(1+nu_mem)); % Shear modulus of membrane

E_susp = 1.4e7;               % Young modulus of suspension
nu_susp = 0.4999;             % Poisson ratio of suspension
rho_susp = 1100;              % mass density of suspension
G_susp = E_mem/(2*(1+nu_mem));   % Shear modulus of suspension

Loudspeaker = struct( 'Diaphragm_radius',        0.045     ,...
                      'Diaphragm_height',        0.03      ,...
                      'DustCap_radius',          0.0225    ,...
                      'DustCap_sphere_radius',   0.03      ,...
                      'VoiceCoil_radius',        0.015     ,...
                      'SuspensionWidth',         10E-3     ,...
                      'MembraneThickness',       0.7E-3    );

res = 40;
[ NiHu_mesh ] = create_loudspeaker( Loudspeaker, res );
plot_mesh(NiHu_mesh);
%%
FEnode = [NiHu_mesh.Nodes(:,1), zeros(size(NiHu_mesh.Nodes,1),3), NiHu_mesh.Nodes(:,2:4)];
FEel0 = [
    inf abs('quad4')
    NiHu_mesh.Elements(:,[5:8 3:4]);
    ];
[node_suspension, element_suspension] = mesh_select(NiHu_mesh, sprintf('r > (%g+0.000)', Loudspeaker.Diaphragm_radius), 'ID', 'any');
FEel0(element_suspension+1,5:6) = 2*FEel0(element_suspension+1,5:6);

FEelt = [];
femesh('AddSel');
openFEM_mesh = femesh('Model');
openFEM_mesh.pl = [
    %   MatID type                       E      nu        rho    G    eta
    1     fe_mat('m_elastic','SI',1)    E_mem   nu_mem  rho_mem  G_mem   0;     % membrane and dust cap
    2     fe_mat('m_elastic','SI',1)    E_susp nu_susp  rho_susp G_susp  0;     % suspension
    ];
openFEM_mesh.il = [
    %   ProID type                     f d 0 w
    1     fe_mat('p_shell','SI',1)     0 1 0 Loudspeaker.MembraneThickness % membrane  property
    2     fe_mat('p_shell','SI',1)     0 1 0 Loudspeaker.MembraneThickness % suspension property
    ];


% openFEM model, FE mass and stiffness matrices
DOF = feutil('GetDOF', openFEM_mesh);
[m, k, mdof] = fe_mk(openFEM_mesh, 'Options', [0 2]);

%% Nodes of boundary conditions and excitation
node_per = mesh_select(NiHu_mesh,sprintf( 'r >= %g-0.0015',...
    Loudspeaker.Diaphragm_radius+Loudspeaker.SuspensionWidth+1E-3), 'ID', 'any');


node_exc1 = mesh_select(NiHu_mesh, sprintf('r < (%g+0.002)', Loudspeaker.VoiceCoil_radius), 'ID', 'any');
node_exc2 = find(NiHu_mesh.Nodes(:,4) < min(NiHu_mesh.Nodes(:,4))+0.002 );
node_exc = intersect(node_exc1,node_exc2);
%node_exc = mesh_select(NiHu_mesh, sprintf('r < (%g+0.002)', Loudspeaker.Diaphragm_radius-1E-2), 'ID', 'any');
nodes = union(node_per,node_exc);
plot_mesh(NiHu_mesh,'node',nodes,ones(length(nodes),1))
shading interp
light
lighting phong
drawnow
caxis([-1,1])
%% Implement boundary conditions
dof_per_zeroed = node_per + 3e-2;
DOFind = fe_c(DOF, dof_per_zeroed, 'ind', 2);

m2 = m(DOFind,DOFind);
k2 = k(DOFind,DOFind);
freeDOF = DOF(DOFind);

% Implement excitation with load vector
DOFexcitation = fe_c(freeDOF, node_exc+3e-2, 'ind', 1);
F = zeros(length(freeDOF),1);
F(DOFexcitation) = 1;
%% Solve system of linear equations
f = 10;
omega = 2*pi*f;
Left_side = k2 - omega^2*m2;
v(DOFind) = Left_side\F;
v = v(:);
zind = round(mod(DOF,1)*100) == 3;

plot_mesh(NiHu_mesh, 'node', floor(DOF(zind)), (v(zind)));
colorbar
AVD = sum( abs( v(zind) ) );
%%
[NiHu_mesh2, v_z2] = post_process_mesh(NiHu_mesh, v, DOF, zind, Loudspeaker);
%%
animate_solution( NiHu_mesh, omega, 1*v, 15, DOF, zind);
%% Map to rect
N = 1024;
total = Loudspeaker.Diaphragm_radius + Loudspeaker.SuspensionWidth;
x = linspace(-total,total,N);
y = linspace(-total,total,N);
[X, Y] = meshgrid(x,y);
Fn = scatteredInterpolant(NiHu_mesh2.Nodes(:,2),NiHu_mesh2.Nodes(:,3),v_z2);
r0 = sqrt( X.^2 + Y.^2 );

Vq = Fn(X,Y);
Vq(isnan(Vq)) = 0;
hx = sum(Vq)';
surf(X,Y,Vq);
xlabel('x [m]');ylabel('y [m]');
shading interp;axis equal;axis tight;
%%

dx = x(2)-x(1);
temp = [zeros(4096*30,1);hx;zeros(4096*30,1)];
Nx = length(temp);
kx = 2*pi*( (-Nx/2:Nx/2-1)/(Nx*dx) );
k0 = omega/343.1;
sinfi = kx/k0;
Vkx = fftshift(fft(fftshift(temp)));

sinfi2 = sinfi((abs(sinfi)<=1));
Vfi = Vkx((abs(sinfi)<=1));
figure
polar(asin(sinfi2),abs(Vfi'));
%xlim([-500,500])