newPackage(
	"QuadraticIdealExamplesByRoos",
	Version => "0.1",
	Date => "June, 2023",
	AuxiliaryFiles => false,
	Authors => {{Name => "David Eisenbud", Email => "de@msri.org"},
	    {Name => "Michael Perlman", Email => "mperlman@umn.edu"}, 
	    {Name => "Ritvik Ramkumar", Email => "ritvikr@cornell.edu"},
	    {Name => "Deepak Sireeshan"},
	    {Name => "Aleksandra Sobieska", Email => "asobieska@wisc.edu"},
	    {Name => "Teresa Yu", Email => "twyu@umich.edu"},
	    {Name => "Jacob Zoromski", Email => "jzoromsk@nd.edu"} },
	Headline => "Examples of Quadratic Ideals with Embedding Dimension Four by Jan-Erik Roos",
	PackageExports => {"Depth"},
	DebuggingMode => true)
export {
 "roosTable"
}



roosTable = (S = ZZ/101[x,y,z,u];
L= {
sub(ideal "0",S),
ideal "x2",
ideal "x2, y2",
ideal "x2, xy",
ideal "x2, y2, z2",
ideal "x2, y2 + xz, yz",
ideal "x2 + y2, z2 + u2, xz + yu",
ideal "x2, y2, xz",
ideal "x2, xy, y2",
ideal "x2, xy, xz",
ideal "x2, y2, z2, u2",
ideal "x2 + xy, y2 + xu, z2 + xu, u2 + zu",
ideal "x2 + z2 + u2, y2, xz, yu + zu",
ideal "xz, y2, z2 + u2, yu + zu",
ideal "xz, y2, yz + u2, yu + zu",
ideal "xy + z2 + yu, y2, yu + zu, xz",
ideal "xz, yz + xu, y2, yu + zu",
ideal "x2, y2, z2, yu",
ideal "xz, y2, yu + zu, u2",
ideal "xz,y2,yu+z2,yu+zu",
ideal "xz,y2,z2,yu+zu",
ideal "x2+xy,xu,xz+yu,y2",
ideal "xz,xu,y2,z2",
ideal "xz,y2,yz+z2,yu+zu",
ideal "x2,xy,xz,u2",
ideal "xz,y2,yu,zu",
ideal "xy,xz,y2,yz",
ideal "x2,xy,xz,xu",
ideal "x2+xy,y2+xu,z2+xu,zu+u2,yz",
ideal "xy+u2,xz,x2+z2+u2,y2,yu+zu",
ideal "x2-y2,y2-z2,z2-u2,xz+yu,-x2+xy-yz+xu",
ideal "x2+z2,xz,y2,yu+zu,u2",
ideal "x2+xy,y2+yz,y2+xu,z2+xu,zu+u2",
ideal "x2+xy+yu+u2,y2,xz,x2+z2+u2,yu+zu",
ideal "x2+z2+u2,y2,xz,xy+yz+yu,yu+zu",
ideal "x2+y2,z2,u2,yz-yu,xz+zu",
ideal "x2,y2,xy-zu,yz-xu,(x-y)(z-u)",
ideal "x2,y2,z2,zu,u2",
ideal "x2+yz+u2,xz+z2 +yu,xy, xu, zu",
ideal "x2-xu,xu-y2,y2-z2,z2-u2,xz+yu",
ideal "xy,y2,z2,zu,u2",
ideal "x2 +xy,zu,y2,xu,xz+yu",
ideal "x2,y2,yz,zu,u2",
ideal "xz,yz,y2,yu+zu,z2+u2",
ideal "xy+yz,xy+z2 +yu,yu + zu, y2, xz",
ideal "x2,xy,yz,zu,u2",
ideal "x2+xy,y2,xu,xz+yu,-x2+xz-yz",
ideal "xy,z2+yu,yu+zu,y2,xz",
ideal "xz, y2, z2, yu, zu",
ideal "x2, xy, xz, y2, z2",
ideal "xy,xz,yz+xu,z2,zu",
ideal "x2, xy, xz, y2, yz",
ideal "y2 -u2, xz, yz, z2, zu",
ideal "x2,xz,y2,z2,yu+zu,u2",
ideal "x2 +xy,xz+yu,xu,y2,z2,zu+u2",
ideal "x2+xz+u2,xy,xu,x2 -y2,z2,zu",
ideal "x2+yz+u2,xu,x2 +xy,xz+yu, zu+u2, y2+z2",
ideal "x2 + xy, x2 + zu, y2, z2, xz+yu, xu",
ideal "x2 - y2, xy, xu, z2, zu,xz + yu",
ideal "x2 + yz + u2, xz+yu, zu, xy, z2, xu",
ideal "x2-y2, xy, z2, xu, zu, u2",
ideal "x2 - y2, xy, xu, yz+yu, z2, zu",
ideal "x2, xy, xu, y2, z2, zu",
ideal "x2 - y2, xy, z2, xu, yu, zu",
ideal "x2, xy, xz, y2, yu+z2, yu+zu",
ideal "xz, y2, yu, z2, zu, u2",
ideal "xy, xz, y2, yu, z2, zu",
ideal "x2, xy, xz,y2, yz, z2",
ideal "x2, xz, xu, xy-zu, yz, z2",
ideal "x2, xy, xz xu, y2, yz",
ideal "x2, y2, z2, u2, xy, zu, yz+xu",
ideal "x2-y2, xy, yz, zu, z2, xz+yu, xu",
ideal "x2, y2, z2, u2, zu, yu, xu",
ideal "x2, xy+z2, yz, xu, yu, zu, u2",
ideal "x2,xy,xz,xu,y2,yz,u2",
ideal "x2, xy, xz, xu, z2, zu, yu",
ideal "x2, xy, xz, xu, y2, yz, yu",
ideal "x2,xy,y2,z2,zu,u2,xz+yu,yz-xu",
ideal "x2,xy,xz,xu,y2,yu,z2,zu",
ideal "x2,xy,xz,y2,yz,yu,z2,zu",
ideal "x2,y2,z2,u2,xy,xz,yz-xu,yu,zu",
ideal "x2,xy,xz,xu,y2,zu,u2,yz,yu",
ideal "x2,y2,z2,u2,xy,xz,xu,yz,yu,zu"
};
new HashTable from for i from 1 to 83 list i => L#(i-1)
)

---TO DO: write function that outputs the depth zero ones, as described in the tables

----need to export and document the following two functions
onedimIrrationalPoincare = (degs1 := {18,24,25,26,28,30,33}; --this example is from froberg-roos 2000, lofwall-lundqvist-roos
    ker map(QQ[t], QQ[w_1 .. w_7, Degrees => degs1], apply(degs1, a -> t^a))
)

twodimIrrationalPoincare = (degs2 := {{36,0}, {33,3}, {30,6}, {28,8}, {26,10}, {25,11}, {24,12}, {18,18}, {0,36}};
    ker map(QQ[t,s], QQ[w_1 .. w_9, Degrees => degs2], apply(degs2, a -> t^(a#0)*s^(a#1)))
    )



      -* Documentation section *-
      
beginDocumentation()

doc ///
Key
 "roosTable"
Headline
 Creates hashtable of Jan-Erik Roos' examples of quadratic ideals
Usage
 H = roosTable ()
Outputs
 H: HashTable
Description
  Text
    This is based on Main Theorem and Tables 3-8 in "Homological properties of the homology
    algebra of the Koszul complex of a local ring: Examples and questions" by Jan-Erik Roos, Journal of Algebra 
    465 (2016) 399-436. The ideals in this table exemplify 83 known cases of bi-graded Poincar\'e series of 
    quadratic ideals of embedding dimension four in characteristic zero.
  Example
    roosTable
///



-* Test section *-
TEST///
I=roosTable#7
S=ring I
assert(depth (S/I) == 0)
///

TEST///
I=roosTable#68
S=ring I
assert(depth (S/I) == 1)
///

end--

uninstallPackage "QuadraticIdealExamplesByRoos"
restart
installPackage "QuadraticIdealExamplesByRoos"

H0 = applyPairs(H, (i,I) -> if (depth(S/I) == 0 and not isBurch I and not isGolod (S/I)) then (i,I))

elapsedTime applyPairs(H0, (i,I) -> (i,socleSummandsSemigroup(I,7)))
--22 & 59 have semigroup {0,6} --when checking up to spot 6
H5 = new HashTable from {42 => H0#42, 34 => H0#34, 55 => H0#55, 60 => H0#60}
elapsedTime applyPairs(H5, (i,I) -> (i, socleSummandsSemigroup(I,9)))
