function val = test_func(input)
global G_var
G_var = G_var + 1;
a = G_var;
val = sin(input) + a;
end