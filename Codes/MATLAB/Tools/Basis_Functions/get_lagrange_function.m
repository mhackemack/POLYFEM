function out = get_lagrange_function(flag, ctype, dim)

if strcmp(flag, 'vals')
    if dim == 1
        out = @evaluate_1D_reference_values;
    elseif dim == 2
        if ctype == 1
            out = @evaluate_2D_triangle_reference_values;
        else
            out = @evaluate_2D_quad_reference_values;
        end
    else
        if ctype == 1
            out = @evaluate_3D_tet_reference_values;
        else
            out = @evaluate_3D_hex_reference_values;
        end
    end
elseif strcmp(flag, 'grads')
    if dim == 1
        out = @evaluate_1D_reference_gradients;
    elseif dim == 2
        if ctype == 1
            out = @evaluate_2D_triangle_reference_gradients;
        else
            out = @evaluate_2D_quad_reference_gradients;
        end
    else
        if ctype == 1
            out = @evaluate_3D_tet_reference_gradients;
        else
            out = @evaluate_3D_hex_reference_gradients;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Basis Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_1D_reference_values(deg, x)
if deg == 1
    out = [1-x, x];
elseif deg == 2
    out = [2*(1/2-x).*(1-x), -2*x.*(1/2-x), 4*x.*(1-x)];
elseif deg == 3
    out = [9/2*(1/3-x).*(2/3-x).*(1-x), 9/2*(1/3-x).*(2/3-x).*x, 27/2*(1-x).*(2/3-x).*x, -27/2*(1-x).*(1/3-x).*x];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_1D_reference_gradients(deg, x)
if deg == 1
    out = [-1, 1];
elseif deg == 2
    out = [4*x-3, 4*x-1, 4-8*x];
elseif deg == 3
    out = [-9/2*(3*x.^2 - 4*x + 11/9), 9/2*(3*x.^2 - 2*x + 2/9), 9/2*(9*x.^2 - 10*x + 2), -9/2*(9*x.^2 - 8*x + 1)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_triangle_reference_values(deg, x)
s = x(:,1); t = x(:,2);
if deg == 1
    out = [1-s-t,s,t];
elseif deg == 2
    out = [  2*s.^2+4*s.*t+2*t.^2-3*s-3*t+1,...
             2*s.^2-s,...
             2*t.^2-t,...
             -4*s.^2-4*s.*t+4*s,...
             4*s.*t,...
             -4*t.^2-4*s.*t+4*t...
          ];
elseif deg == 3
    out = [9/2*(1-s-t).*(1/3-s-t).*(2/3-s-t),...
           9/2*s.*(1/3-s).*(2/3-s),...
           9/2*t.*(1/3-t).*(2/3-t),...
           27/2*s.*(1-s-t).*(2/3-s-t),...
          -27/2*s.*(1/3-s).*(1-s-t),...
          -27/2*s.*t.*(1/3-s),...
          -27/2*s.*t.*(1/3-t),...
          -27/2*t.*(1-s-t).*(1/3-t),...
           27/2*t.*(1-s-t).*(2/3-s-t),...
           27*s.*t.*(1-s-t)];
elseif deg == 4
    out = [64/6*(1-s-t).*(1/4-s-t).*(2/4-s-t).*(3/4-s-t),...
          -64/6*s.*(1/4-s).*(2/4-s).*(3/4-s),...
          -64/6*t.*(1/4-t).*(2/4-t).*(3/4-t),...
           256/6*s.*(1-s-t).*(2/4-s-t).*(3/4-s-t),...
          -768/12*s.*(1-s-t).*(1/4-s).*(3/4-s-t),...
           256/6*s.*(1-s-t).*(1/4-s).*(2/4-s),...
           256/6*s.*t.*(1/4-s).*(2/4-s),...
           768/12*s.*t.*(1/4-s).*(1/4-t),...
           256/6*s.*t.*(1/4-t).*(2/4-t),...
           256/6*t.*(1-s-t).*(1/4-t).*(2/4-t),...
          -768/12*t.*(1-s-t).*(1/4-t).*(3/4-s-t),...
           256/6*t.*(1-s-t).*(2/4-s-t).*(3/4-s-t),...
           768/6*s.*t.*(1-s-t).*(3/4-s-t),...
          -768/6*s.*t.*(1-s-t).*(1/4-s),...
          -768/6*s.*t.*(1-s-t).*(1/4-t)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_quad_reference_values(deg, xx)
s = xx(:,1); t = xx(:,2);
if deg == 1
    out = [(1-s).*(1-t), s.*(1-t), s.*t, t.*(1-s)];
elseif deg == 2
    out = [ 4*(1/2-s).*(1-s).*(1/2-t).*(1-t),...
           -4*(1/2-s).*s.*(1/2-t).*(1-t),...
            4*(1/2-s).*s.*(1/2-t).*t,...
           -4*(1/2-s).*(1-s).*(1/2-t).*t,...
            8*s.*(1-s).*(1/2-t).*(1-t),...
           -8*s.*(1/2-s).*t.*(1-t),...
           -8*s.*(1-s).*t.*(1/2-t),...
            8*(1-s).*(1/2-s).*t.*(1-t),...
           16*s.*(1-s).*t.*(1-t)];
elseif deg == 3
    x=s;y=t;
    out = [    (x-1.0).*(x-1.0./3.0).*(x-2.0./3.0).*(y-1.0).*(y-1.0./3.0).*(y-2.0./3.0).*(8.1e1./4.0),...
               x.*(x-1.0./3.0).*(x-2.0./3.0).*(y-1.0).*(y-1.0./3.0).*(y-2.0./3.0).*(-8.1e1./4.0),...
               x.*y.*(x-1.0./3.0).*(x-2.0./3.0).*(y-1.0./3.0).*(y-2.0./3.0).*(8.1e1./4.0),...
               y.*(x-1.0).*(x-1.0./3.0).*(x-2.0./3.0).*(y-1.0./3.0).*(y-2.0./3.0).*(-8.1e1./4.0),...
               x.*(x-1.0).*(x-2.0./3.0).*(y-1.0).*(y-1.0./3.0).*(y-2.0./3.0).*(-2.43e2./4.0),...
               x.*(x-1.0).*(x-1.0./3.0).*(y-1.0).*(y-1.0./3.0).*(y-2.0./3.0).*(2.43e2./4.0),...
               x.*y.*(x-1.0./3.0).*(x-2.0./3.0).*(y-1.0).*(y-2.0./3.0).*(2.43e2./4.0),...
               x.*y.*(x-1.0./3.0).*(x-2.0./3.0).*(y-1.0).*(y-1.0./3.0).*(-2.43e2./4.0),...
               x.*y.*(x-1.0).*(x-1.0./3.0).*(y-1.0./3.0).*(y-2.0./3.0).*(-2.43e2./4.0),...
               x.*y.*(x-1.0).*(x-2.0./3.0).*(y-1.0./3.0).*(y-2.0./3.0).*(2.43e2./4.0),...
               y.*(x-1.0).*(x-1.0./3.0).*(x-2.0./3.0).*(y-1.0).*(y-1.0./3.0).*(2.43e2./4.0),...
               y.*(x-1.0).*(x-1.0./3.0).*(x-2.0./3.0).*(y-1.0).*(y-2.0./3.0).*(-2.43e2./4.0),...
               x.*y.*(x-1.0).*(x-2.0./3.0).*(y-1.0).*(y-2.0./3.0).*(7.29e2./4.0),...
               x.*y.*(x-1.0).*(x-1.0./3.0).*(y-1.0).*(y-2.0./3.0).*(-7.29e2./4.0),...
               x.*y.*(x-1.0).*(x-1.0./3.0).*(y-1.0).*(y-1.0./3.0).*(7.29e2./4.0),...
               x.*y.*(x-1.0).*(x-2.0./3.0).*(y-1.0).*(y-1.0./3.0).*(-7.29e2./4.0)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_triangle_reference_gradients(deg, xx)
if deg == 1
    out = [-1,-1;1,0;0,1];
elseif deg == 2
    s = xx(:,1); t = xx(:,2);
    out = [  4*s+4*t-3,  4*s+4*t-3;...
             4*s-1,         0;...
             0,                4*t-1;...
             -8*s-4*t+4,-4*s;...
             4*t,           4*s;...
             -4*t,    -8*t-4*s+4];
elseif deg == 3
    s = xx(:,1); t = xx(:,2);
    out = [ - ((9*s)/2 + (9*t)/2 - 9/2)*(s + t - 1/3) - ((9*s)/2 + (9*t)/2 - 9/2)*(s + t - 2/3) - (9*(s + t - 1/3)*(s + t - 2/3))/2, - ((9*s)/2 + (9*t)/2 - 9/2)*(s + t - 1/3) - ((9*s)/2 + (9*t)/2 - 9/2)*(s + t - 2/3) - (9*(s + t - 1/3)*(s + t - 2/3))/2;...
                                                                  (9*(s - 1/3)*(s - 2/3))/2 + (9*s*(s - 1/3))/2 + (9*s*(s - 2/3))/2, 0;...
                                                                                                                                  0, (9*(t - 1/3)*(t - 2/3))/2 + (9*t*(t - 1/3))/2 + (9*t*(t - 2/3))/2;...
                                                   (27*s*(s + t - 1))/2 + (27*s*(s + t - 2/3))/2 + (27*(s + t - 1)*(s + t - 2/3))/2, (27*s*(s + t - 1))/2 + (27*s*(s + t - 2/3))/2;...
                                                         - (27*s*(s + t - 1))/2 - (27*(s - 1/3)*(s + t - 1))/2 - (27*s*(s - 1/3))/2,-(27*s*(s - 1/3))/2;...
                                                                                                    (27*s*t)/2 + (27*t*(s - 1/3))/2, (27*s*(s - 1/3))/2;...
                                                                                                                 (27*t*(t - 1/3))/2, (27*s*t)/2 + (27*s*(t - 1/3))/2;...
                                                                                                                -(27*t*(t - 1/3))/2, - (27*t*(s + t - 1))/2 - (27*(t - 1/3)*(s + t - 1))/2 - (27*t*(t - 1/3))/2;...
                                                                                      (27*t*(s + t - 1))/2 + (27*t*(s + t - 2/3))/2, (27*t*(s + t - 1))/2 + (27*t*(s + t - 2/3))/2 + (27*(s + t - 1)*(s + t - 2/3))/2;...
                                                                                                        - 27*t*(s + t - 1) - 27*s*t, - 27*s*(s + t - 1) - 27*s*t];
%     out = [ - (3*s + 3*t - 1)*((3*s)/2 + (3*t)/2 - 1) - (3*(3*s + 3*t - 1)*(s + t - 1))/2 - 3*((3*s)/2 + (3*t)/2 - 1)*(s + t - 1), - (3*s + 3*t - 1)*((3*s)/2 + (3*t)/2 - 1) - (3*(3*s + 3*t - 1)*(s + t - 1))/2 - 3*((3*s)/2 + (3*t)/2 - 1)*(s + t - 1);...
%                                                                   (3*s*(3*s - 1))/2 + 3*s*((3*s)/2 - 1) + (3*s - 1)*((3*s)/2 - 1),                                                                                                                     0;...
%                                                                                                                                 0,                                                       (3*t*(3*t - 1))/2 + 3*t*((3*t)/2 - 1) + (3*t - 1)*((3*t)/2 - 1);...
%                                        (27*s*(s + t - 1))/2 + 9*((3*s)/2 + (3*t)/2 - 1)*(s + t - 1) + 9*s*((3*s)/2 + (3*t)/2 - 1),                                                                    (27*s*(s + t - 1))/2 + 9*s*((3*s)/2 + (3*t)/2 - 1);...
%                                                          - (9*s*(3*s - 1))/2 - (27*s*(s + t - 1))/2 - (9*(3*s - 1)*(s + t - 1))/2,                                                                                                    -(9*s*(3*s - 1))/2;...
%                                                                                                    (9*t*(3*s - 1))/2 + (27*s*t)/2,                                                                                                     (9*s*(3*s - 1))/2;...
%                                                                                                                 (9*t*(3*t - 1))/2,                                                                                        (9*s*(3*t - 1))/2 + (27*s*t)/2;...
%                                                                               - (9*t*(3*s - 1))/2 - 3*t*((9*s)/2 + (9*t)/2 - 9/2),                                                             - (9*t*(3*s - 1))/2 - (3*s - 1)*((9*s)/2 + (9*t)/2 - 9/2);...
%                                                                             9*t*((3*s)/2 + (3*t)/2 - 1) + (3*t*(9*s + 9*t - 9))/2,                       ((3*s)/2 + (3*t)/2 - 1)*(9*s + 9*t - 9) + 9*t*((3*s)/2 + (3*t)/2 - 1) + (3*t*(9*s + 9*t - 9))/2;...
%                                                                                                       - 27*t*(s + t - 1) - 27*s*t,                                                                                           - 27*s*(s + t - 1) - 27*s*t];
elseif deg == 4
    s = xx(:,1); t = xx(:,2);
    out = [ (32*(s + t - 1/2)*(s + t - 1/4)*(s + t - 3/4))/3 + ((32*s)/3 + (32*t)/3 - 32/3)*(s + t - 1/2)*(s + t - 1/4) + ((32*s)/3 + (32*t)/3 - 32/3)*(s + t - 1/2)*(s + t - 3/4) + ((32*s)/3 + (32*t)/3 - 32/3)*(s + t - 1/4)*(s + t - 3/4), (32*(s + t - 1/2)*(s + t - 1/4)*(s + t - 3/4))/3 + ((32*s)/3 + (32*t)/3 - 32/3)*(s + t - 1/2)*(s + t - 1/4) + ((32*s)/3 + (32*t)/3 - 32/3)*(s + t - 1/2)*(s + t - 3/4) + ((32*s)/3 + (32*t)/3 - 32/3)*(s + t - 1/4)*(s + t - 3/4);...
                                                                                                 (32*s*(s - 1/2)*(s - 1/4))/3 + (32*s*(s - 1/2)*(s - 3/4))/3 + (32*s*(s - 1/4)*(s - 3/4))/3 + (32*(s - 1/2)*(s - 1/4)*(s - 3/4))/3,                                                                                                                                                                                                                                 0;...
                                                                                                                                                                                                                                 0,                                                                                                 (32*t*(t - 1/2)*(t - 1/4))/3 + (32*t*(t - 1/2)*(t - 3/4))/3 + (32*t*(t - 1/4)*(t - 3/4))/3 + (32*(t - 1/2)*(t - 1/4)*(t - 3/4))/3;...
                                                             - (128*(s + t - 1)*(s + t - 1/2)*(s + t - 3/4))/3 - (128*s*(s + t - 1)*(s + t - 1/2))/3 - (128*s*(s + t - 1)*(s + t - 3/4))/3 - (128*s*(s + t - 1/2)*(s + t - 3/4))/3,                                                                                                               - (128*s*(s + t - 1)*(s + t - 1/2))/3 - (128*s*(s + t - 1)*(s + t - 3/4))/3 - (128*s*(s + t - 1/2)*(s + t - 3/4))/3;...
                                                                                               64*s*(s + t - 1)*(s + t - 3/4) + 64*(s - 1/4)*(s + t - 1)*(s + t - 3/4) + 64*s*(s - 1/4)*(s + t - 1) + 64*s*(s - 1/4)*(s + t - 3/4),                                                                                                                                                                         64*s*(s - 1/4)*(s + t - 1) + 64*s*(s - 1/4)*(s + t - 3/4);...
                                                                                     - (128*(s - 1/2)*(s - 1/4)*(s + t - 1))/3 - (128*s*(s - 1/2)*(s - 1/4))/3 - (128*s*(s - 1/2)*(s + t - 1))/3 - (128*s*(s - 1/4)*(s + t - 1))/3,                                                                                                                                                                                                    -(128*s*(s - 1/2)*(s - 1/4))/3;...
                                                                                                                                                     (128*s*t*(s - 1/2))/3 + (128*s*t*(s - 1/4))/3 + (128*t*(s - 1/2)*(s - 1/4))/3,                                                                                                                                                                                                     (128*s*(s - 1/2)*(s - 1/4))/3;...
                                                                                                                                                                                       64*s*t*(t - 1/4) + 64*t*(s - 1/4)*(t - 1/4),                                                                                                                                                                                       64*s*t*(s - 1/4) + 64*s*(s - 1/4)*(t - 1/4);...
                                                                                                                                                                                                     (128*t*(t - 1/2)*(t - 1/4))/3,                                                                                                                                                     (128*s*t*(t - 1/2))/3 + (128*s*t*(t - 1/4))/3 + (128*s*(t - 1/2)*(t - 1/4))/3;...
                                                                                                                                                                                                    -(128*t*(t - 1/2)*(t - 1/4))/3,                                                                                     - (128*(t - 1/2)*(t - 1/4)*(s + t - 1))/3 - (128*t*(t - 1/2)*(t - 1/4))/3 - (128*t*(t - 1/2)*(s + t - 1))/3 - (128*t*(t - 1/4)*(s + t - 1))/3;...
                                                                                                                                                                         64*t*(t - 1/4)*(s + t - 1) + 64*t*(t - 1/4)*(s + t - 3/4),                                                                                               64*t*(s + t - 1)*(s + t - 3/4) + 64*(t - 1/4)*(s + t - 1)*(s + t - 3/4) + 64*t*(t - 1/4)*(s + t - 1) + 64*t*(t - 1/4)*(s + t - 3/4);...
                                                                                                               - (128*t*(s + t - 1)*(s + t - 1/2))/3 - (128*t*(s + t - 1)*(s + t - 3/4))/3 - (128*t*(s + t - 1/2)*(s + t - 3/4))/3,                                                             - (128*(s + t - 1)*(s + t - 1/2)*(s + t - 3/4))/3 - (128*t*(s + t - 1)*(s + t - 1/2))/3 - (128*t*(s + t - 1)*(s + t - 3/4))/3 - (128*t*(s + t - 1/2)*(s + t - 3/4))/3;...
                                                                                                                                                     128*t*(s + t - 1)*(s + t - 3/4) + 128*s*t*(s + t - 1) + 128*s*t*(s + t - 3/4),                                                                                                                                                     128*s*(s + t - 1)*(s + t - 3/4) + 128*s*t*(s + t - 1) + 128*s*t*(s + t - 3/4);...
                                                                                                                                                           - 128*s*t*(s - 1/4) - 128*s*t*(s + t - 1) - 128*t*(s - 1/4)*(s + t - 1),                                                                                                                                                                                 - 128*s*t*(s - 1/4) - 128*s*(s - 1/4)*(s + t - 1);...
                                                                                                                                                                                 - 128*s*t*(t - 1/4) - 128*t*(t - 1/4)*(s + t - 1),                                                                                                                                                           - 128*s*t*(t - 1/4) - 128*s*t*(s + t - 1) - 128*s*(t - 1/4)*(s + t - 1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_2D_quad_reference_gradients(deg, xx)
s = xx(:,1); t = xx(:,2);
if deg == 1
    out = [t-1,s-1;1-t,-s;t,s;-t,1-s];
elseif deg == 2
    stt = 8*s*t*t; sst = 8*s*s*t;
    ss = s*s; st = s*t; tt = t*t;
    out = [  stt-12*st+4*s-6*tt+9*t-3,    sst-6*ss-12*st+9*s+4*t-3;...
             stt-12*st+4*s-2*tt+3*t-1,    sst-6*ss-4*st+3*s;...
             stt-4*st-2*tt+t,             sst-2*ss-4*st+s;...
             stt-4*st-6*tt+3*t,           sst-2*ss-12*st+3*s+4*t-1;...
           -2*stt+24*st-8*s+8*tt-12*t+4, -2*sst+12*ss+16*st-12*s;...
           -2*stt+16*st+4*tt-4*t,        -2*sst+8*ss+8*st-4*s;...
           -2*stt+8*st+8*tt-4*t,         -2*sst+4*ss+16*st-4*s;...
           -2*stt+16*st+12*tt-12*t,      -2*sst+8*ss+24*st-12*s-8*t+4;...
            4*stt-32*st-16*tt+16*t,       4*sst-16*ss-32*st+16*s];
elseif deg == 3
    out = [ (81*(s - 1)*(s - 1/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4 + (81*(s - 1)*(s - 2/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4 + (81*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4, (81*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 1/3))/4 + (81*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 2/3))/4 + (81*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4;...
                      - (81*s*(s - 1/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4 - (81*s*(s - 2/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4 - (81*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4,                 - (81*s*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 1/3))/4 - (81*s*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 2/3))/4 - (81*s*(s - 1/3)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4;...
                                          (81*t*(s - 1/3)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4 + (81*s*t*(s - 1/3)*(t - 1/3)*(t - 2/3))/4 + (81*s*t*(s - 2/3)*(t - 1/3)*(t - 2/3))/4,                               (81*s*(s - 1/3)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4 + (81*s*t*(s - 1/3)*(s - 2/3)*(t - 1/3))/4 + (81*s*t*(s - 1/3)*(s - 2/3)*(t - 2/3))/4;...
                            - (81*t*(s - 1)*(s - 1/3)*(t - 1/3)*(t - 2/3))/4 - (81*t*(s - 1)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4 - (81*t*(s - 1/3)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4,           - (81*t*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1/3))/4 - (81*t*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 2/3))/4 - (81*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4;...
                       - (243*s*(s - 1)*(t - 1)*(t - 1/3)*(t - 2/3))/4 - (243*s*(s - 2/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4 - (243*(s - 1)*(s - 2/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4,                    - (243*s*(s - 1)*(s - 2/3)*(t - 1)*(t - 1/3))/4 - (243*s*(s - 1)*(s - 2/3)*(t - 1)*(t - 2/3))/4 - (243*s*(s - 1)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4;...
                         (243*s*(s - 1)*(t - 1)*(t - 1/3)*(t - 2/3))/4 + (243*s*(s - 1/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4 + (243*(s - 1)*(s - 1/3)*(t - 1)*(t - 1/3)*(t - 2/3))/4,                      (243*s*(s - 1)*(s - 1/3)*(t - 1)*(t - 1/3))/4 + (243*s*(s - 1)*(s - 1/3)*(t - 1)*(t - 2/3))/4 + (243*s*(s - 1)*(s - 1/3)*(t - 1/3)*(t - 2/3))/4;...
                                             (243*t*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 2/3))/4 + (243*s*t*(s - 1/3)*(t - 1)*(t - 2/3))/4 + (243*s*t*(s - 2/3)*(t - 1)*(t - 2/3))/4,                                (243*s*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 2/3))/4 + (243*s*t*(s - 1/3)*(s - 2/3)*(t - 1))/4 + (243*s*t*(s - 1/3)*(s - 2/3)*(t - 2/3))/4;...
                                           - (243*t*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 1/3))/4 - (243*s*t*(s - 1/3)*(t - 1)*(t - 1/3))/4 - (243*s*t*(s - 2/3)*(t - 1)*(t - 1/3))/4,                              - (243*s*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 1/3))/4 - (243*s*t*(s - 1/3)*(s - 2/3)*(t - 1))/4 - (243*s*t*(s - 1/3)*(s - 2/3)*(t - 1/3))/4;...
                                         - (243*t*(s - 1)*(s - 1/3)*(t - 1/3)*(t - 2/3))/4 - (243*s*t*(s - 1)*(t - 1/3)*(t - 2/3))/4 - (243*s*t*(s - 1/3)*(t - 1/3)*(t - 2/3))/4,                                - (243*s*(s - 1)*(s - 1/3)*(t - 1/3)*(t - 2/3))/4 - (243*s*t*(s - 1)*(s - 1/3)*(t - 1/3))/4 - (243*s*t*(s - 1)*(s - 1/3)*(t - 2/3))/4;...
                                           (243*t*(s - 1)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4 + (243*s*t*(s - 1)*(t - 1/3)*(t - 2/3))/4 + (243*s*t*(s - 2/3)*(t - 1/3)*(t - 2/3))/4,                                  (243*s*(s - 1)*(s - 2/3)*(t - 1/3)*(t - 2/3))/4 + (243*s*t*(s - 1)*(s - 2/3)*(t - 1/3))/4 + (243*s*t*(s - 1)*(s - 2/3)*(t - 2/3))/4;...
                                 (243*t*(s - 1)*(s - 1/3)*(t - 1)*(t - 1/3))/4 + (243*t*(s - 1)*(s - 2/3)*(t - 1)*(t - 1/3))/4 + (243*t*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 1/3))/4,              (243*t*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1))/4 + (243*t*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1/3))/4 + (243*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 1/3))/4;...
                               - (243*t*(s - 1)*(s - 1/3)*(t - 1)*(t - 2/3))/4 - (243*t*(s - 1)*(s - 2/3)*(t - 1)*(t - 2/3))/4 - (243*t*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 2/3))/4,            - (243*t*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1))/4 - (243*t*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 2/3))/4 - (243*(s - 1)*(s - 1/3)*(s - 2/3)*(t - 1)*(t - 2/3))/4;...
                                                 (729*t*(s - 1)*(s - 2/3)*(t - 1)*(t - 2/3))/4 + (729*s*t*(s - 1)*(t - 1)*(t - 2/3))/4 + (729*s*t*(s - 2/3)*(t - 1)*(t - 2/3))/4,                                      (729*s*(s - 1)*(s - 2/3)*(t - 1)*(t - 2/3))/4 + (729*s*t*(s - 1)*(s - 2/3)*(t - 1))/4 + (729*s*t*(s - 1)*(s - 2/3)*(t - 2/3))/4;...
                                               - (729*t*(s - 1)*(s - 1/3)*(t - 1)*(t - 2/3))/4 - (729*s*t*(s - 1)*(t - 1)*(t - 2/3))/4 - (729*s*t*(s - 1/3)*(t - 1)*(t - 2/3))/4,                                    - (729*s*(s - 1)*(s - 1/3)*(t - 1)*(t - 2/3))/4 - (729*s*t*(s - 1)*(s - 1/3)*(t - 1))/4 - (729*s*t*(s - 1)*(s - 1/3)*(t - 2/3))/4;...
                                                 (729*t*(s - 1)*(s - 1/3)*(t - 1)*(t - 1/3))/4 + (729*s*t*(s - 1)*(t - 1)*(t - 1/3))/4 + (729*s*t*(s - 1/3)*(t - 1)*(t - 1/3))/4,                                      (729*s*(s - 1)*(s - 1/3)*(t - 1)*(t - 1/3))/4 + (729*s*t*(s - 1)*(s - 1/3)*(t - 1))/4 + (729*s*t*(s - 1)*(s - 1/3)*(t - 1/3))/4;...
                                               - (729*t*(s - 1)*(s - 2/3)*(t - 1)*(t - 1/3))/4 - (729*s*t*(s - 1)*(t - 1)*(t - 1/3))/4 - (729*s*t*(s - 2/3)*(t - 1)*(t - 1/3))/4,                                    - (729*s*(s - 1)*(s - 2/3)*(t - 1)*(t - 1/3))/4 - (729*s*t*(s - 1)*(s - 2/3)*(t - 1))/4 - (729*s*t*(s - 1)*(s - 2/3)*(t - 1/3))/4];
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_3D_tet_reference_values(deg, xx)
x = xx(1); y = xx(2); z = xx(3);
if deg == 1
    out = [1-x-y-z,x,y,z];
elseif deg == 2
    
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_3D_tet_reference_gradients(deg, xx)
% x = xx(1); y = xx(2); z = xx(3);
if deg == 1
    out = [-1,-1,-1;1,0,0;0,1,0;0,0,1];
elseif deg == 2
    
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_3D_hex_reference_values(deg, xx)
x = xx(:,1); y = xx(:,2); z = xx(:,3);
if deg == 1
    out = [(1-x).*(1-y).*(1-z),...
            x.*(1-y).*(1-z),...
            x.*y.*(1-z),...
            (1-x).*y.*(1-z),...
            (1-x).*(1-y).*z,...
            x.*(1-y).*z,...
            x.*y.*z,...
            (1-x).*y.*z];
elseif deg == 2
    
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = evaluate_3D_hex_reference_gradients(deg, xx)
x = xx(:,1); y = xx(:,2); z = xx(:,3);
if deg == 1
    out = [-y.*z+y+z-1,  -x.*z+x+z-1,  -x.*y+x+y-1;...
           (1-y).*(1-z),  x.*z-x,       x.*y-x;...
           y.*(1-z),      x.*(1-z),    -x.*y;...
           y.*z-y,       (1-x).*(1-z),  x.*y-y;...
           y.*z-z,       x.*z-z         (1-x).*(1-y);...
           (1-y).*z,     -x.*z          x.*(1-y);...
           y.*z,          x.*z,         x.*y;...
           -y.*z,         (1-x).*z,     (1-x).*y];
elseif deg == 2
    
else
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

