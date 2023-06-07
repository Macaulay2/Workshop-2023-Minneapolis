restart
uninstallPackage "MatrixSchubert"
restart
installPackage "MatrixSchubert"
restart
needsPackage "MatrixSchubert"
I = schubertDetIdeal {2,1,6,3,5,4}
schubs = schubertDecomposition I --should be {2,1,6,3,5,4}
isIntersectionSchubIdeals I --should be true
isASMIdeal I