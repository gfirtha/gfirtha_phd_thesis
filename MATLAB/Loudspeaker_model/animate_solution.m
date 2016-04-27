function [] = animate_solution( NiHu_mesh,omega, v, T, DOF, zind)%, aviobj)
orig = NiHu_mesh.Nodes(:,4);
node_ind = floor(DOF(zind));
v_z = v(zind);
for i = 1:10*T
cla
NiHu_mesh.Nodes(:,4) = orig + v_z*sin(2*pi*i/T);
%plot_mesh(NiHu_mesh, v_z);
plot_mesh(NiHu_mesh, 'node', node_ind, v_z);
zlim([-0.085 0.085]);
shading interp
light
lighting phong
drawnow
%print('-dpng', sprintf('animation/%03d.png',i),'-r200')
%F(i) = getframe(gcf);
end

end

