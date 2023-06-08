restart
needsPackage "Brackets"

G = gc(a .. f,3); 

ABC = (a*b*c)_G;
DEF = (d*e*f)_G;

aeLine = (a*e)_G;
bfLine = (b*f)_G;
cdLine = (c*d)_G;

afLine = (a*f)_G;
bdLine = (b*d)_G;
ceLine = (c*e)_G;

X = aeLine ^ bdLine;
Y = bfLine ^ ceLine;
Z = cdLine ^ afLine;

XYZ = X*Y*Z;


