clear all;
close all;

E_sp = 1.5e9;                % Young modulus of membrane
nu_sp = 0.33;                % Poisson ratio of membrane
rho_sp = 900;                % mass density of membrane
G_sp = E_sp/(2*(1+nu_sp));   % Shear modulus of membrane

Loudspeaker = struct( 'Diaphragm_radius',        0.045     ,...
    'Diaphragm_height',        0.03      ,...
    'DustCap_radius',          0.025    ,...
    'DustCap_sphere_radius',   0.05      ,...
    'VoiceCoil_radius',        0.015     ,...
    'SuspensionWidth',         10E-3     ,...
    'MembraneThickness',       0.7E-3    );

res = 60;
[ NiHu_mesh ] = create_loudspeaker( Loudspeaker, res);
%[ NiHu_mesh ] = create_circle( Loudspeaker.Diaphragm_radius + Loudspeaker.SuspensionWidth, res);
%
plot_mesh(NiHu_mesh);
drawnow
%%
FEnode = [NiHu_mesh.Nodes(:,1), zeros(size(NiHu_mesh.Nodes,1),3), NiHu_mesh.Nodes(:,2:4)];
FEel0 = [
    inf abs('quad4')
    NiHu_mesh.Elements(:,[5:8 3:4]);
    ];
FEelt = [];
femesh('AddSel');
openFEM_mesh = femesh('Model');
openFEM_mesh.pl = [
    %   MatID type                       E    nu    rho    G    eta
    1     fe_mat('m_elastic','SI',1)    E_sp nu_sp  rho_sp G_sp  0
    ];
openFEM_mesh.il = [
    %   ProID type                     f d 0 w
    1     fe_mat('p_shell','SI',1) 0 1 0 Loudspeaker.MembraneThickness % plate property
    ];

% openFEM model, FE mass and stiffness matrices
DOF = feutil('GetDOF', openFEM_mesh);
[m, k, mdof] = fe_mk(openFEM_mesh, 'Options', [0 2]);

%% Nodes of boundary conditions and excitation
node_per = mesh_select(NiHu_mesh,sprintf( 'r >= %g-0.0015',...
    Loudspeaker.Diaphragm_radius+Loudspeaker.SuspensionWidth+1E-3), 'ID', 'any');


node_exc1 = mesh_select(NiHu_mesh, sprintf('r < (%g+0.02)', Loudspeaker.VoiceCoil_radius), 'ID', 'any');
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
%% Calculate transfer function
f = (0:20:5000)';
omega = 2*pi*f;
transfer = zeros(size(f));
for n = 1:length(f)
    fprintf('Processing frequency bin no. %d out of %d \n', n , length(f));
    Left_side = k2 - omega(n)^2*m2;
    v(DOFind) = Left_side\F;
    v = v(:);
    zind = round(mod(DOF,1)*100) == 3;
    transfer(n) = sum( abs( v(zind) ) );
end

plot( f, 20*log10( transfer ) );