function GaloisRepresentationHyperellipticLarge(f)
	X:=HyperellipticCurve(f);
	K:=BaseField(X);
	OK:=Integers(K); 
	pi:=UniformizingElement(OK);
	k,m:=ResidueClassField(OK);
	q:=#k;
	p:=Characteristic(k);
	D:=Discriminant(X);
	vD:=Valuation(D);
	if Degree(f) ne p or IsIrreducible(f) eq false or IsCoprime(vD,p-1) eq false then 
		error if RamificationDegree(SplittingField(f)) ne p*(p-1) or IsEven(vD), "Not the right case";
	end if;

		R<x>:=PolynomialRing(K);
		K0:=SplittingField(f);
		rts:=Roots(f,K0);
		alpha:=rts[1][1];
		beta:=rts[2][1];
		RK<t>:=PolynomialRing(K0);
		r:=ChangePrecision(alpha-beta,Precision(K0));
		v:=Valuation(r);
		assert IsOdd(v);
		r0:=r/UniformizingElement(K0)^(v-1);
		F0:=ext<K0|t^2-r0>;
		f2:=MinimalPolynomial(UniformizingElement(F0),K);
		
		reps:=GaloisRepresentations(f2);
		F:=Field(reps[1]);                            
		G:=Group(reps[1]); 
		I:=reps[1]`I; 
		a:=reps[1]`act;
		C<i>:=ComplexField();
		if IsSquare(q) then
			assert I eq G;
			ansreps:=[r: r in reps | Dimension(r) eq p-1 and IsFaithful(Character(r))];
			assert #ansreps eq 1;
			if p mod 4 eq 1 then
				rho:=UnramifiedCharacter(K,(Sqrt(p))^(AbsoluteInertiaDegree(K)))*ansreps[1];
			end if;
			if p mod 4 eq 3 then
				rho:=UnramifiedCharacter(K,(i*Sqrt(p))^(AbsoluteInertiaDegree(K)))*ansreps[1];
			end if;
			end if;

		if IsSquare(q) eq false then
			rts2:=Roots(f,F);
			alpha:=rts2[1][1];
			sigma:=[g: g in G | Order(g) eq p][1];  
			asigma:=a(sigma);
			beta:=asigma(alpha);
			sba:=Sqrt(beta-alpha);
			prec:=Precision(F) div 2;
			frobs:=[g: g in G | Order(g) eq 2 and g notin I and Valuation(a(g)(sba)-sba) gt prec];
			assert #frobs eq 1;
			phi:=frobs[1];  
			assert sigma*phi eq phi*sigma;
			if p mod 4 eq 1 then
				ansreps:=[r: r in reps | IsFaithful(Character(r)) and Real(Character(r)(sigma*phi)) gt 0];
				assert #ansreps eq 1;
				rho:=UnramifiedCharacter(K,(Sqrt(p))^(AbsoluteInertiaDegree(K)))*ansreps[1];
			end if;
			if p mod 4 eq 3 then
				ansreps:=[r: r in reps | IsFaithful(Character(r)) and Imaginary(Character(r)(sigma*phi)) gt 0];
				assert #ansreps eq 1;
				rho:=UnramifiedCharacter(K,(i*Sqrt(p))^(AbsoluteInertiaDegree(K)))*ansreps[1];
			end if;
		end if;
	return rho;
	
end function;