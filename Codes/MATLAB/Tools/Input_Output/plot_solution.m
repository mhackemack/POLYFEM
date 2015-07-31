function plot_solution(mesh, DoF, FE, x, basis_name)

if nargin == 4
    switch(mesh.Dimension)
        case(1)
            plot_1D_solution(mesh, DoF, FE, x);
        case(2)
            plot_2D_solution(mesh, DoF, FE, x);
        case(3)
            plot_3D_solution(mesh, DoF, FE, x);
    end
elseif nargin == 5
    switch(mesh.Dimension)
        case(1)
            plot_1D_solution(mesh, DoF, FE, x, basis_name);
        case(2)
            plot_2D_solution(mesh, DoF, FE, x, basis_name);
        case(3)
            plot_3D_solution(mesh, DoF, FE, x, basis_name);
    end
end
    