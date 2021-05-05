// Galois representations at 2 for elliptic curves with Q8 or SL_2(F3) inertia

function ChangeModel(E, r,s,t: Iso:=false)
   a1,a2,a3,a4,a6:=Explode(aInvariants(E));
   a6 := a6 + r*(a4 + r*(a2 + r)) - t*(a3 + r*a1 + t);
   a4 := a4 - s*a3 + 2*r*a2 - (t + r*s)*a1 + 3*r*r - 2*s*t;
   a3 := a3 + r*a1 +t+t;
   a2 := a2 - s*a1 + 3*r - s*s;
   a1 := a1 + s+s;
   F:=EllipticCurve([a1,a2,a3,a4,a6]);
   if not Iso then return F; end if;
   i:=Isomorphism(F,E,[r,s,t,1]);
   return F,i;
end function;


function GaloisRepresentationEllChar2NonAbelian(E)   
// using Nirvana Coppola

  K:=BaseField(E); 
  error if Type(K) ne FldPad, "E should be a curve over a local field";
  
  E:=WeierstrassModel(E);  // Put to y^2=x^3+a4x+a6
  _,_,_,a4,a6:=Explode(aInvariants(E));
  
  OK:=Integers(K);                // Field data
  pi:=UniformizingElement(OK);
  k,m:=ResidueClassField(OK);
  q:=#k;
  p:=Characteristic(k);
  error if p ne 2, "Residue characteristic should be 2";  
  
  D:=Discriminant(E);              
 
  error if Valuation(D) mod 3 eq 0 and not IsPower(K!D,3), 
    "The case where D^(1/3) gives a non-trivial unramified extension is not implemented";
  
  error if IsSquare(q), "K contains zeta_3 - wrong case";
  
  R<t>:=PolynomialRing(K);
  f:=t^8+18*a4*t^4+108*a6*t^2-27*a4^2; 
  reps:=GaloisRepresentations(f);    
  Fz:=Field(reps[1]);                      // F(zeta_3) 
  OFz:=Integers(Fz);                
  piFz:=UniformizingElement(OFz);
  kFz,m:=ResidueClassField(OFz);
  G:=Group(reps[1]);                       // G = GL_2(F3) or SD_16
  I:=reps[1]`I;                            // I = SL_2(F3) or Q8

  error if #G notin [16,48], "Not in the GL_2(F3) or SD16 case";
   
  slopes:=reps[1]`r;                       // 3-torsion point slopes
  error if Max([Valuation(Evaluate(f,a)):a in slopes]) lt Precision(Fz) div 2, 
    "slopes are not roots of f to the correct precision"; 

  xP:=[t^2/3: t in slopes];
  yP:=[(t^4+3*a4)/6/t: t in slopes];

  error if Valuation(Evaluate(t^3+a4*t+a6,xP[1])-yP[1]^2) lt Precision(Fz) div 2, 
    "point xP[1],yP[1] is not on E";

  a:=reps[1]`act;       // action of G on slopes
  
  E1:=ChangeModel(E,xP[1],slopes[1],yP[1]);    // shift xP,yP to (0,0), gradient to 0
  A1,A2,A3,A4,A6:=Explode(aInvariants(E1));
  
  ok,scale:=IsPower(Fz!A3,3); 
  assert ok;
  
  redxP:=[m((x-xP[1])/scale^2): x in xP];   
  redyP:=[m((yP[i]-yP[1]-slopes[1]*(xP[i]-xP[1]))/scale^3): i in [1..8]];   
  Qindex:=Position(redxP,1);
  error if Qindex eq 0, "Could not find a point Q which reduces to (1,z)";
  
  Ered:=EllipticCurve([kFz|0,0,1,0,0]);
  redpts:=[Ered![redxP[i],redyP[i]]: i in [1..8]];
  P:=redpts[1];
  Q:=redpts[Qindex];
  mPmQ:=Position(redpts,-P-Q);
  PmQ:=Position(redpts,P-Q);
  error if mPmQ eq 0 or PmQ eq 0, "Could not find -P-Q and P-Q";

  alist:=[g: g in G | (1^g eq mPmQ) and (Qindex^g eq PmQ)];
  error if #alist ne 1, "Could not lift [[2,1],[2,2]] to Gal(Fz/K)";
  a:=alist[1];
  
  ansreps:=[r: r in reps | Degree(r) eq 2 and IsFaithful(Character(r)) and 
    Imaginary(Character(r)(a)) lt 0];     // for etale cohomology want trace -sqrt(-2)
  assert #ansreps eq 1;

  C<i>:=ComplexField();
  rho:=UnramifiedCharacter(K,(i*Sqrt(C!q))) * ansreps[1];  
  return rho;    
end function;