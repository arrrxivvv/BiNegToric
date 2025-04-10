#array of arrays initialization

module BiNegToric2d

using StaticArrays
using LinearAlgebra
using JLD2
using FilenameManip

using Infiltrator

const sigMatZ = @MMatrix [1 0; 0 -1];
const IsingJ = @MMatrix [1.0 -1; -1 1];
const IsingHRaw = @MMatrix [1.0 0; 0 -1];
const IsingH = ( 1 .- IsingHRaw ) ./ 2;
const idMat2d = @MMatrix [1 0; 0 1];

const bool2Lst = @MVector [false, true];
const sgn2Lst = @MVector [1, -1];
const it2Lst = @MVector [1, 2];

const fMainBiNegMinArr = "biNegMinArr";
const attrLstBiNegMinArr = ["nSpin", "lambMax", "lambStep"];

struct TransHelper
	transMat::MMatrix{2,2,Float64,4};
	matResult::MMatrix{2,2,Float64,4};
	matTmp::MMatrix{2,2,Float64,4};
end

function TransHelper()
	transMat = @MMatrix zeros(2,2);
	matResult = @MMatrix zeros(2,2);
	matTmp = @MMatrix zeros(2,2);
	
	TransHelper( transMat, matResult, matTmp );
end

struct BiTransHelper
	transMatLst::Vector{MMatrix{2,2,Float64,4}};
	matResult::MMatrix{2,2,Float64,4};
	matTmp::MMatrix{2,2,Float64,4};
end

function BiTransHelper()
	transMatLst = [ @MMatrix zeros(2,2) for ii = 1 : 2 ];
	matResult = @MMatrix zeros(2,2);
	matTmp = @MMatrix zeros(2,2);
	
	BiTransHelper( transMatLst, matResult, matTmp );
end

struct BiNegData{N,N2,Pow2N}
	numL::Int64;
	helper::BiTransHelper;
	iter2Prod::Iterators.ProductIterator{NTuple{N,MVector{2,Int64}}};
	# it2ProdLst::MArray{MVector{N,Int64},N};
	it2ProdLst::MArray{NTuple{N,2},NTuple{N,Int64},N,Pow2N};
	iter2ABProd::Iterators.ProductIterator{ NTuple{2, Iterators.ProductIterator{NTuple{N, MVector{2,Int64}}} } };
	coeff1Lst::Array{Float64,N2};
	coeff0Lst::Array{Float64,N2};
	coeffBiLst::Array{Float64,N2};
	spinALst::Vector{Bool};
	spinBLst::Vector{Bool};
	spinLst::Vector{Bool};
	iaALst::MVector{N,Int64};
	ibBLst::MVector{N,Int64};
	function BiNegData( numL::Int64, helper::BiTransHelper, iter2Prod::Iterators.ProductIterator{NTuple{N,MVector{2,Int64}}}, it2ProdLst::MArray{NTuple{N,2},NTuple{N,Int64},N,Pow2N}, iterABProd::Iterators.ProductIterator{NTuple{2,Iterators.ProductIterator{NTuple{N,MVector{2,Int64}}}}}, coeff1Lst::Array{Float64,N2}, coeff0Lst::Array{Float64,N2}, coeffBilst::Array{Float64,N2}, spinALst::Vector{Bool}, spinBLst::Vector{Bool}, spinLst::Vector{Bool}, iaALst, ibBLst ) where {N,N2,Pow2N}
		if N2 != 2*N || N != numL || Pow2N != 2^N
			error( "N2, 2N, Pow2N or numL mismatch" );
		end
		new{N,N2,Pow2N}( numL, helper, iter2Prod, it2ProdLst, iterABProd, coeff1Lst, coeff0Lst, coeffBilst, spinALst, spinBLst, spinLst, iaALst, ibBLst );
	end
end

function BiNegData( numL::Int64, transHelper::BiTransHelper )
	iter2Prod = Iterators.product( ntuple( x->it2Lst, numL )... );
	it2ProdLst = collect( iter2Prod );
	iter2ABProd = Iterators.product( ntuple( x->iter2Prod, 2 )... );
	coeff1Lst = zeros( ntuple( x->2, 2*numL ) );
	coeff0Lst = deepcopy( coeff1Lst );
	coeffBiLst = similar(coeff0Lst);
	spinALst = zeros(Bool, numL);
	spinBLst = similar(spinALst);
	spinLst = similar(spinALst);
	iaALst = @MVector ones(Int64,numL);
	ibBLst = similar(iaALst);
	
	data = BiNegData( numL, transHelper, iter2Prod, it2ProdLst, iter2ABProd, coeff1Lst, coeff0Lst, coeffBiLst, spinALst, spinBLst, spinLst, iaALst, ibBLst );
	calcPT0Coeff(data);
	
	return data;
end

function calcPT0Coeff( data::BiNegData )
	for (iAProd, iBProd) in data.iter2ABProd
		for ii = 1 : data.numL
			data.spinALst[ii] = bool2Lst[iAProd[ii]];
			data.spinBLst[ii] = bool2Lst[iBProd[ii]];
		end
		if reduce(xor, data.spinALst) && reduce(xor, data.spinBLst)
			data.spinLst[1] = false;
			for ii = 1 : data.numL-1
				data.spinLst[ii+1] = !xor( data.spinLst[ii], data.spinBLst[ii] );
			end
			ans = false;
			for ii = 1 : data.numL
				if data.spinALst[ii]
					ans = xor(ans, data.spinLst[ii]);
				end
			end
			data.coeff0Lst[iAProd..., iBProd...] = ans ? -2 : 2;
		end
	end
end

function calcPTCoeff( data::BiNegData, beta, lambA, lambB )
	matResult = data.helper.matResult;
	matTmp = data.helper.matTmp;
	
	kappaA = -log( tanh( beta * lambA ) );
	
	for ii = 1 : 2
		sgn = sgn2Lst[ii];
		data.helper.transMatLst[ii] .= exp.( -kappaA .* IsingH .+ beta .* lambB .* IsingJ .* sgn );
	end
	for ( iLin, ( iAProd, iBProd ) ) in Iterators.enumerate( data.iter2ABProd )
		matResult .= idMat2d;
		for ii = 1 : data.numL
			iASgn = iAProd[ii];
			iBSgn = iBProd[ii];
			if bool2Lst[iASgn]
				matResult .= mul!( matTmp, matResult, sigMatZ );
			end
			matResult .= mul!( matTmp, matResult, data.helper.transMatLst[iBSgn] );
		end
		data.coeff1Lst[iLin] = tr( matResult );
	end
end

function calcCoeffBiPTFromPT( data::BiNegData )
	# coeffBiLst = zeros( ntuple( x->2, 2*numL ) );
	itABProdWithLin = Iterators.enumerate( data.iter2ABProd );
	for (iLin, ( iAProd, iBProd ) ) in Iterators.enumerate( data.iter2ABProd )
		val = 0;
		# for (jLin, ( jAProd, jBProd ) ) in Iterators.enumerate( data.iter2ABProd )
		for (jLin, ( jAProd, jBProd ) ) in itABProdWithLin
			for ii = 1 : data.numL
				data.iaALst[ii] = iAProd[ii] == jAProd[ii] ? 1 : 2;
				data.ibBLst[ii] = iBProd[ii] == jBProd[ii] ? 1 : 2;
			end
			val += abs(data.coeff1Lst[data.iaALst..., data.ibBLst...]) * data.coeff0Lst[jLin];
		end
		data.coeffBiLst[iLin] = val;
	end
end

function calcBiNegCoeff( data::BiNegData, beta, lambA, lambB )
	calcPTCoeff( data, beta, lambA, lambB );
	calcCoeffBiPTFromPT( data );
end

function calcPT0Coeff( numL::Int64 )
	it2ProdLst = Iterators.product( ntuple( x -> it2Lst, numL )... );
	iter2ABProd = Iterators.product( ntuple( x -> it2ProdLst, 2 )... );
	
	coeffLst0 = zeros( ntuple(x->2, 2*numL) );
	spinALst = zeros(Bool, numL);
	spinBLst = similar(spinALst);
	spinLst = similar(spinALst);
	for (iAProd, iBProd) in iter2ABProd
		for ii = 1 : numL
			spinALst[ii] = bool2Lst[iAProd[ii]];
			spinBLst[ii] = bool2Lst[iBProd[ii]];
		end
		if reduce(xor, spinALst) && reduce(xor, spinBLst)
			spinLst[1] = false;
			for ii = 1 : numL-1
				spinLst[ii+1] = !xor( spinLst[ii], spinBLst[ii] );
			end
			ans = false;
			for ii = 1 : numL
				if spinALst[ii]
					ans = xor(ans, spinLst[ii]);
				end
			end
			coeffLst0[iAProd..., iBProd...] = ans ? -2 : 2;
		end
	end
	
	return coeffLst0;
end

function calcPTCoeff( helper::TransHelper, numL::Int64, beta, lambA, lambB )
	# matStart = @MMatrix [1 0; 0 1];
	matResult = helper.matResult
	matTmp = helper.matTmp

	it2ProdLst = Iterators.product( ntuple( x -> it2Lst, numL )... );
	iter2ABProd = Iterators.product( ntuple( x -> it2ProdLst, 2 )... );
	
	kappaA = -log( tanh( beta * lambA ) );
	
	coeffLst = zeros( ntuple( x->2, 2*numL ) );
	transMatPNLst = [ exp.( -kappaA .* IsingH .+ beta .* lambB .* IsingJ .* sgn ) for sgn in sgn2Lst ];
	for (iLin, ( iAProd, iBProd ) ) in Iterators.enumerate( iter2ABProd )
		matResult .= idMat2d;
		for ii = 1 : numL
			iASgn = iAProd[ii];
			iBSgn = iBProd[ii];
			if bool2Lst[iASgn]
				matResult .= mul!( matTmp, matResult, sigMatZ );
			end
			matResult .= mul!( matTmp, matResult, transMatPNLst[iBSgn] );
		end
		coeffLst[iLin] = tr( matResult );
	end
	
	for (iLin, ( iAProd, iBProd ) ) in Iterators.enumerate( iter2ABProd )
		
	end
	
	return coeffLst;
end

function calcPTCoeffTest( numL::Int64 = 3 )
	it2ProdLst = Iterators.product( ntuple( x -> it2Lst, numL )... );
	iter2ABProd = Iterators.product( ntuple( x -> it2ProdLst, 2 )... );
	
	lambA = 1;
	lambB = 1;
	beta = 1;
	kappaA = -log( tanh( beta * lambA ) );
	
	coeffLst = zeros( ntuple( x->2, 2*numL ) );
	ALst = zeros( Bool, numL );
	BLst = deepcopy(ALst);
	transMatPNLst = [ exp.( -kappaA .* IsingH .+ beta .* lambB .* IsingJ .* sgn ) for sgn in sgn2Lst ];
	matStart = @MMatrix [1 0; 0 1];
	matResult = @MMatrix zeros(2,2); 
	matTmp = @MMatrix zeros(2,2);
	for (iLin, ( iAProd, iBProd ) ) in Iterators.enumerate( iter2ABProd )
		matResult .= matStart;
		for ii = 1 : numL
			iASgn = iAProd[ii];
			iBSgn = iBProd[ii];
			if bool2Lst[iASgn]
				matResult .= mul!( matTmp, matResult, sigMatZ );
			end
			matResult .= mul!( matTmp, matResult, transMatPNLst[iBSgn] );
		end
		coeffLst[iLin] = tr( matResult );
	end
	
	for (iLin, ( iAProd, iBProd ) ) in Iterators.enumerate( iter2ABProd )
		
	end
	
	return coeffLst;
end

function calcCoeffBiPT( helper::TransHelper, numL::Int64, lambA, lambB )
	beta = 1;
	
	coeffLst1 = calcPTCoeff( helper, numL, beta, lambA, lambB );
	# coeffLst0 = calcPTCoeff( helper, numL, beta, lambA, lambB );
	coeffLst0 = calcPT0Coeff( numL );
	
	coeffBiLst = calcCoeffBiPTFromPT( numL, coeffLst1, coeffLst0 );
end

function calcCoeffBiPTFromPT( numL::Int64, coeffLst::Array, coeffLst0::Array )
	it2ProdLst = Iterators.product( ntuple( x -> it2Lst, numL )... );
	iter2ABProd = Iterators.product( ntuple( x -> it2ProdLst, 2 )... );
	
	coeffBiLst = zeros( ntuple( x->2, 2*numL ) );
	iaALst = @MVector zeros(Int64, numL);
	ibBLst = similar(iaALst);
	for (iLin, ( iAProd, iBProd ) ) in Iterators.enumerate( iter2ABProd )
		val = 0;
		for (jLin, ( jAProd, jBProd ) ) in Iterators.enumerate( iter2ABProd )
			for ii = 1 : numL
				iaALst[ii] = iAProd[ii] == jAProd[ii] ? 1 : 2;
				ibBLst[ii] = iBProd[ii] == jBProd[ii] ? 1 : 2;
			end
			val += abs(coeffLst[iaALst..., ibBLst...]) * coeffLst0[jLin];
		end
		coeffBiLst[iLin] = val;
	end
	
	return coeffBiLst;
end

include("biNegToric2d_IsingTest.jl");
include("biNegToric2d_runFunc.jl");

end
