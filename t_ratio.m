function T = t_ratio(smallmatpc,arc)

[~,Vconv] = convhulln(smallmatpc); % volume of convex hull of the data.
[~,Vsimplex] = convhulln(arc'); % volume of the bound simplex defined by the archetypes.
T = Vsimplex / Vconv;