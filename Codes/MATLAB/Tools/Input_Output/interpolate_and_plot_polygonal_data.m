function interpolate_and_plot_polygonal_data(v, xin, yin, data, flag)
if nargin < 5, flag = 'surface'; end
[mx,nx] = size(xin);
[my,ny] = size(yin);
[md,nd] = size(data);

if nx > 1 && ny > 1 && nd > 1
    x=xin;
    y=yin;
    z = remove_values(v,x,y,data);
elseif nx == 1 && nd ==1
    nix = get_interp_number(mx);
    niy = get_interp_number(my);
    xlin = linspace(min(xin),max(xin),nix);
    ylin = linspace(min(yin),max(yin),niy);
    [x,y] = meshgrid(xlin,ylin);
    z = griddata(xin,yin,data,x,y);
    z = remove_values(v,x,y,z);
else
    error('Could not recognize input structure. Make sure matrix dimensions match.')
end
if strcmp(flag, 'contour')
    contourf(x,y,z);
elseif strcmp(flag, 'surface')
    surf(x,y,z);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = remove_values(v,x,y,data)
out = data;
in = inpolygon(x,y,v(:,1),v(:,2));
out(~in) = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_interp_number(in)
if in <= 1e2
    out = 20;
elseif in <= 1e3
    out = 30;
elseif in <= 1e4
    out = 40;
elseif in <= 1e5
    out = 50;
elseif in <= 1e6
    out = 60;
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%