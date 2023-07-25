installPackage("A1BrouwerDegrees")

-- viewHelp A1BrouwerDegrees
-- check A1BrouwerDegrees

T4 = CC[v];
f5 = {v^4 + v^3 - v^2 - v};
f5GD = globalA1Degree(f5);
f5GD


ff = GF(17);
T3 = ff[y_1..y_3];
f3 = {y_1^2, y_2^2, y_3^2};
f4 = {y_2^2, y_3^2, y_1^2};
f3GD = globalA1Degree(f3)
f4GD = globalA1Degree(f4)


ff = GF(17);
T3 = ff[w_1..w_3];
f3 = {w_1^2, w_2^2, w_3^2};
globalA1Degree(f3)
