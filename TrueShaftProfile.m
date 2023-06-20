function ORexp = TrueShaftProfile(xcoord, baseradius, steploc, stepval)

ORexp = sym(baseradius);

numsteps = length(steploc);

for i = 1:numsteps
    ORexp = ORexp + stepval(i)*heaviside(xcoord-steploc(i));
end

end