eq_pts=[];
for i=8:12
    f=@(states)eq_hh(states,i)
    [vars,fval]=fsolve(f,[-60 0 0.6 0.3]);
    eq_pts=[eq_pts;vars];
end