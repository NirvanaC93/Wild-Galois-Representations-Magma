# Wild Galois Representations
 In this repository I include the Magma functions that compute the Galois representation attached to an elliptic curve with non abelian inertia image, or a hyperelliptic curve of the form y^2=f(x), where f is a polynomial of degree p over a p-adic field, with maximal inertia image.
 
**What's in this repository**
The file "Wild3.m" contains the function "GaloisRepresentationEllChar3inC3D4case", that computes the Galois representation on the étale cohomology of an elliptic curve over a 3-adic field K (i.e. a finite extension of Q3) of odd inertia degree over Q3, in the case that the image of inertia is non-abelian (equivalently, of order 12). This happens if and only if the residue characteristic of K is not a square, the discriminant has odd valuation and the Kodaira tipe of the curve is one of II, II*, IV, IV* (cf. [K90]). If this does not happen, the function "GaloisRepresentation" that is already in Magma computes the Galois representation.

The result is given in the form of the tensor product of an unramified character and a 2-dimensional irreducible representation of the finite group C3:D4 (cf. [GN]), as in the main theorem of [C20].

The file "Wild2.m" contains the function "GaloisRepresentationEllChar2NonAbelian", that similarly computes the Galois representation on the étale cohomology of an elliptic curve over a 2-adic field K of odd inertia degree over Q2, in the case that the image of inertia is non-abelian. See again K90 to see when this occurs. If the image of inertia is abelian, but K has even inertia degree over Q2, the Magma function "GaloisRepresentation" computes the correct result. This is based on [C20.2]. This files also contains the auxiliary function "ChangeModel", that is used in the main function.

Note that the function GaloisRepresentation that is in Magma does compute the Galois representation in each of these cases, but it is slow and it may not complete the computation.

*Edit: In the most recent version of Magma, the function GaloisRepresentation has been replaced with the one in this repository for the 3-adic case.


The file "WildP.m" contains the function "GaloisRepresentationHyperellipticLarge", which given a polynomial f of degree p over a p-adic field, computes the Galois representation on the étale cohomology of the hyperelliptic curve y^2=f(x), if this has maximal inertia degree, as explained in [C21+]. Note that, if p=3, this function coincides with that of "Wild3.m".

**How to use the functions**

1. Download the files "Wild3.m", "Wild2.m" and "WildP.m".
2. Open Magma, and run the command: Load "Wild3.m"; (respectively Wild2, WildP);
3. In order to use Wild3.m or Wild2.m, define an elliptic curve E over a 3-adic or 2-adic field; then run: GaloisRepresentationEllChar3inC3D4case(E); or GaloisRepresentationEllChar2NonAbelian(E); respectively.
3bis. In order to use WildP.m, define a polynomial f of degree p over a p-adic field; then run: GaloisRepresentationHyperellipticLarge(f);

**Examples**

E3:=EllipticCurve([pAdicField(3,40)|0,0,0,0,3]);
E2:=EllipticCurve([pAdicField(2,40)|0,0,0,2,2]);

for WildP.m:
p:=*choose a prime*;
R<x>:=PolynomialRing(pAdicField(p,40));
f:=x^p+p;
 
(Note: for p at most 7, the function work, even though slowly. For p larger than 7, it may exceed computation time).

**Acknowledgments**

I wish to thank T. Dokchitser for supervising and helping during the writing of this code.


References:
C20 N. Coppola "Wild Galois Representations: elliptic curves over a 3-adic field", Acta Arithmetica 195 No 3 (2020), pp. 289-303. https://doi.org/10.4064/aa190423-9-1
C20.2 N. Coppola "Wild Galois Representations: elliptic curves over a 2-adic field with non-abelian inertia action", International Journal of Number Theory 16 No 6 (2020), pp. 1199-1208. https://doi.org/10.1142/S179304212050061X
C21+ N. Coppola "Wild Galois representations: a family of hyperelliptic curves with large inertia image", arxiv preprint: https://arxiv.org/abs/2001.08287
GN T. Dokchister "Group Names" https://people.maths.bris.ac.uk/~matyd/GroupNames/
K90 A. Kraus "Sur le défaut de semi-stabilité des courbes elliptiques à réduction additive", Manuscripta Math 69, (1990) 353-385. https://doi.org/10.1007/BF02567933
