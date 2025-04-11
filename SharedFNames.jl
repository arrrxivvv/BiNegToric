module SharedFNames

const dirLog = "./log/";
const dirSrc = "./src/";

const fNameFileLstLst = "fNameFileLstLst.txt";
const fNameFileLstJld2Lst = "fNameFileLstJld2Lst.txt";
const fNameSaveParamsLst = "fNameSaveParamsLst.txt";

const fNameTmpNameFileLst = "fNameTmpNameFileLst.txt";

const fNameShellVarScript = "sharedFNameVars.sh";

function defFNameShellVars()
	varNameLst = [ "fNameFileLstLst", "fNameTmpNameFileLst" ];
	fNameValBaseLst = [ fNameFileLstLst, fNameTmpNameFileLst ];
	fNameValLst = [ dirLog * fNameValBaseLst[ii] for ii = 1 : length(fNameValBaseLst) ];

	open( dirSrc * fNameShellVarScript, "w" ) do io
		println( io, "export" * " " * join( varNameLst, " " ) );
		for ii = 1 : length(varNameLst)
			println( io, varNameLst[ii] * "=" * "\"" * fNameValLst[ii] * "\"" );
		end
	end
end

end # module
