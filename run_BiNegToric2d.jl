using BiNegToric2d
using JLD2
using FilenameManip
using SharedFNames
using DelimitedFiles

nSpinMin = 3;
nSpinMax = 6;

fMod = "";

fNameLst = Vector{String}(undef, nSpinMax - nSpinMin + 1);

lambMax = 2;
lambStep = 0.1;

for (iN,n) = Iterators.enumerate( nSpinMin : nSpinMax )
	@time fNameLst[iN] = BiNegToric2d.runCalcCoeffBiNeg( n; lambMax = lambMax, lambStep = lambStep; fMod = fMod );
end

writedlm( SharedFNames.dirLog * SharedFNames.fNameFileLstLst, fNameLst );
jldsave( SharedFNames.dirLog * SharedFNames.fNameFileLstJld2Lst; fNameArr = fNameLst );
jldsave( SharedFNames.dirLog * SharedFNames.fNameSaveParamsLst; nSpinMin, nSpinMax, lambMax, lambStep );
