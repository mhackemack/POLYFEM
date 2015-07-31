function x = cleanup_small_vals(x)
global glob

v = x < glob.small;
x(v) = 0;