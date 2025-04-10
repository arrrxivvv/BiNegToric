
function runCalcCoeffBiNeg( numL::Int64; lambMax = 2, lambStep = 0.1 )
	# lambMax = 2;
	# lambStep = 0.1;
	lambALst = [lambStep:lambStep:lambMax;];
	lambBLst = [lambStep:lambStep:lambMax;];
	# beat = 1;
	helper = BiTransHelper();
	data = BiNegData( numL, helper );
	
	beta = 1;
	
	coeffBiLstLst = [ zeros( ntuple( x -> 2, 2*numL ) ) for iA = 1 : length(lambALst), iB = 1 : length(lambBLst) ];
	for (iLin, (lambA, lambB)) in Iterators.enumerate( Iterators.product( lambALst, lambBLst ) ) 
		calcBiNegCoeff( data, beta, lambA, lambB );
		coeffBiLstLst[iLin] .= data.coeffBiLst;
	end
	minBiLst = minimum.(coeffBiLstLst);
	
	valLst = Any[numL, lambMax, lambStep];
	fName = fNameFunc( fMainBiNegMinArr, attrLstBiNegMinArr, valLst, jld2Type );
	jldsave( fName; nSpin = numL, lambMax = lambMax, lambStep = lambStep, miniBiLst = minBiLst );
	
	# return coeffBiLstLst, minBiLst;
	return fName;
end


function runCalcCoeffBiNegTransHelper( numL::Int64 )
	lambMax = 2;
	lambStep = 0.1;
	lambALst = [lambStep:lambStep:lambMax;];
	lambBLst = [lambStep:lambStep:lambMax;];
	# beat = 1;
	helper = TransHelper();
	
	coeffBiLstLst = Matrix{Array{Float64,2*numL}}(undef, length(lambALst), length(lambBLst));
	for (iLin, (lambA, lambB)) in Iterators.enumerate( Iterators.product( lambALst, lambBLst ) ) 
		coeffBiLstLst[iLin] = calcCoeffBiPT( helper, numL, lambA, lambB );
	end
	
	return coeffBiLstLst;
end
