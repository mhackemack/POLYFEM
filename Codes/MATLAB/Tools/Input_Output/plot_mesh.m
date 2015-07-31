function plot_mesh(varargin)

if varargin{1}.Dimension == 1
    plot_1D_mesh(varargin{:});
elseif varargin{1}.Dimension == 2
    plot_2D_mesh(varargin{:});
elseif varargin{1}.Dimension == 3
    plot_3D_mesh(varargin{:});
end