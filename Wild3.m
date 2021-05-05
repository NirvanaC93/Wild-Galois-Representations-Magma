function GaloisRepresentationEllChar3inC3D4case(E)   
// using Nirvana Coppola, arxiv: 1812.05651, doi: 10.4064/aa190423-9-1

  K:=BaseField(E); 
  error if Type(K) ne FldPad, "E should be a curve over a local field";
  
  E:=WeierstrassModel(E);  // Put to y^2=x^3+a4x+a6
  _,_,_,a4,a6:=Explode(aInvariants(E));
  
  OK:=Integers(K);                // Field data
  pi:=UniformizingElement(OK);
  k,m:=ResidueClassField(OK);
  q:=#k;
  p:=Characteristic(k);
  error if p ne 3, "Residue characteristic should be 3";
  
  D:=Discriminant(E);              
  v:=Valuation(D) div 4;
  D:=ChangePrecision((OK!D) div (pi^(4*v)),Precision(OK));    // Discriminant mod 4th powers
  
  loc,min:=LocalInformation(E);   // Local information at 3
  _,vdisc,vcond,cv,kod:=Explode(loc); kod:=Sprint(kod);
  
  // Check we are in the C3:D4 case
  error if IsSquare(q) or IsEven(vdisc) or (kod notin ["II","II*","IV","IV*"]), 
    "Not the C3:D4 case";

  R<x>:=PolynomialRing(K);  
  fourtorspoly:=(x^3+a4*x+a6)*(x^4-D);          // F = K(E[4])
  reps:=GaloisRepresentations(fourtorspoly);    // could be possibly made faster by specifying
  F:=Field(reps[1]);                            // an Eisenstein polynomial that defines the ramified part
  
  G:=Group(reps[1]);                       // G = C3:D4
  I:=reps[1]`I;                            // I = C3:C4
  
  a:=reps[1]`act;       // action of G on F
  
  alpha,beta,gamma:=Explode([r[1]: r in Roots(x^3+a4*x+a6,F)]);
  sigma:=[g: g in G | Order(g) eq 3][1];   // sigma = one of the two elements of order 3
  asigma:=a(sigma);
  if Valuation(asigma(alpha)-beta) lt Valuation(asigma(alpha)-gamma) then    
    temp:=beta; beta:=gamma; gamma:=temp;  // make sure that sigma: alpha->beta->gamma->alpha
  end if;
  
  sba:=Sqrt(beta-alpha);    // Sqrt(beta-alpha)
  sgb:=asigma(sba);         // Sqrt(gamma-beta)
  sag:=asigma(sgb);         // Sqrt(alpha-gamma)
  
  prec:=Precision(F) div 2;
  frobs:=[g: g in G | Order(g) eq 2 and g notin I and Valuation(a(g)(sba)-sba) gt prec];
  assert #frobs eq 1;
  phi:=frobs[1];            // phi = Frobenius element of K(sba,sgb,sag)
  assert sigma*phi eq phi*sigma;

  C<i>:=ComplexField();
  
  //  ansreps:=[r: r in reps | IsFaithful(Character(r)) and Imaginary(Character(r)(sigma*phi)) lt 0];
  //  assert #ansreps eq 1;
  //  rho:=UnramifiedCharacter(K,(i*Sqrt(C!q))) * ansreps[1];    
      // the representation on the Tate module as from the paper
 
  ansreps2:=[r: r in reps | IsFaithful(Character(r)) and Imaginary(Character(r)(sigma*phi)) gt 0];
  assert #ansreps2 eq 1;
  rho2:=UnramifiedCharacter(K,(i*Sqrt(p))^AbsoluteInertiaDegree(K)) * ansreps2[1];  // the other one
    // the representation on the etale cohomology (dual to rho up to twist by 1)

  return rho2;
end function;
