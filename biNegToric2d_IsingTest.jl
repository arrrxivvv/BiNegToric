
function calc1dIsingGroup( numL::Int64 = 8 )
	J = 1;
	hMax = 2;
	hLst = [-hMax:0.01:hMax;]
	betaLst = [0.25,0.5,1,1.5];
	
	mArr = zeros( length(hLst), length(betaLst) );
	helper = TransHelper();
	
	for iBeta = 1 : length(betaLst), iH = 1 : length(hLst)
		beta = betaLst[iBeta];
		h = hLst[iH];
		mArr[iH, iBeta], _, _ = calc1dIsing( helper, numL, h, J, beta );
	end
	
	return mArr;
end

function calc1dIsing( helper::TransHelper, numL::Int64, h, J, beta )
	helper.transMat .= exp.( beta .* ( J .* IsingJ .+ h .* IsingHRaw ) );
	helper.matResult .= sigMatZ;
	
	for ii = 1 : numL
		helper.matResult .= mul!( helper.matTmp, helper.matResult, helper.transMat );
	end
	
	mRaw = tr(helper.matResult);
	
	helper.matResult .= idMat2d;
	for ii = 1 : numL
		helper.matResult .= mul!( helper.matTmp, helper.matResult, helper.transMat );
	end
	z = tr(helper.matResult);
	
	return mRaw / z, mRaw, z;
end
