(* ::Package:: *)

(*#################################################################################################### Background functions*)
(***************************************************************************************************** redirect*)
exeDirectory=NotebookDirectory[];

redirect:=If[
	Directory[]<>"\\"!=exeDirectory,
	Module[{directory=exeDirectory},If[directory=!=$Canceled,SetDirectory[directory]]]
	];
	
redirect;

(***************************************************************************************************** filelists*)
filelist[directory_]:=
	Module[{filenames,joiner},
		SetDirectory[directory];
		filenames=FileNames[];
		redirect;
		joiner[logfile_]:=StringJoin[directory,"\\",logfile];
		Map[joiner,filenames]
	];
	
filenames[directory_]:=
	Module[{filenames,joiner,out},
		SetDirectory[directory];
		out=FileNames[];
		redirect;
		out
	];
	
(***************************************************************************************************** fast matrix export function*)
openwrite[path_,data_]:=
	With[
		{str=OpenWrite[path,PageWidth->Infinity],len=Length[data[[1]]]},
		Scan[Write[str,Sequence@@(Flatten[Table[{FortranForm[#[[i]]],OutputForm[","]},{i,len-1}]])~Join~{FortranForm[#[[len]]]}]&,data];
		Close[str];
	];



(*#################################################################################################### Installation Functions*)
(***************************************************************************************************** Find Anaconda/python.exe path*)
pythonExeDir:={
	redirect;
	Which[
		{(*check the search log to see if this system has been used before*)
			searchLog=Import[Directory[]<>"\\installation\\search_log.txt","Text"];
			candidateFiles=DeleteDuplicates[StringSplit[searchLog,"\n"]];
			logSearch=Flatten[Map[FileExistsQ,candidateFiles]];
			searchResults=candidateFiles[[Flatten[Position[logSearch,True]]]];
			searchResults=searchResults[[
				Flatten[
					Position[StringLength[searchResults],Sort[StringLength[searchResults]
					][[1]]]
				]
			]];
			logTrueQ=If[Count[logSearch,x]>0,AnyTrue[logSearch],False]
		}[[1]],
		searchResults,
		{(*check current search to see if the path works for this system*)
			searchResults=Import[Directory[]<>"\\search_results.txt","Text"];
			candidateFiles=StringSplit[
				searchResults,
				"\n"
				][[
					Flatten[
						Position[
							StringContainsQ[
								StringSplit[searchResults,"\n"],"Anaconda"~~x___~~"\\python.exe"
								],
							True
						]
					]
				]];
			AnyTrue[Map[FileExistsQ,candidateFiles]]
		}[[1]],
		searchResults,
		{(*if not search the c drive*)
			RunProcess[{"cmd","/c",Directory[]<>"\\installation\\find_anaconda_C.bat"}];
			searchResults=StringSplit[
				Import[Directory[]<>"\\search_results.txt","Text"],
				"\n"
				];
			Position[
				StringContainsQ[
					searchResults,
					"Anaconda"~~x___~~"\\python.exe"
					],
				True
			]
		}[[1]],
		searchResults,
		{(*if not search the f drive*)
			RunProcess[{"cmd","/c",Directory[]<>"\\installation\\find_anaconda_F.bat"}];
			searchResults=Import[Directory[]<>"\\search_results.txt","Text"];
			StringContainsQ[
				searchResults,
				"Anaconda"~~x___~~"\\python.exe"
			]
		}[[1]],
		searchResults
	];

	If[
		logTrueQ,
		anacondaPythonDir=searchResults,
		(
		candidateFiles=StringSplit[
			searchResults,
			"\n"
			][[
				Flatten[
					Position[
						StringContainsQ[
							StringSplit[searchResults,"\n"],
							"Anaconda"~~x___~~"\\python.exe"
						],
						True
					]
				]
			]];
		anacondaPythonDir=candidateFiles[[
			Flatten[
				Position[StringLength[candidateFiles],Sort[StringLength[candidateFiles]
				][[1]]]
			]
		]]
		)
	];

	newLine[entry_]:=entry<>"\n";
	candidateFiles=DeleteDuplicates[Append[candidateFiles,anacondaPythonDir]];
	searchLogText=StringJoin[Map[newLine,candidateFiles]];
	Export[Directory[]<>"\\installation\\search_log.txt",searchLogText,"Text"];
};

(*pythonExeDir;*)(*Only need this if you are running python functions - eg paramiko for ssh and sftp*)



(*#################################################################################################### Open a file in GaussView*)
GaussView[filePath_]:=Module[
	{},
	batFile=Export[
		Directory[]<>"\\GaussViewOpen.bat",
		"c:\\G16W\\gview.exe "<>filePath<>" %*",
		"Text"
	];
	RunProcess[{"cmd","/c","GaussViewOpen.bat"}];
];


OpenBabelPresets:={
	chargeOB=0,
	multiplicityOB=1
};
OpenBabelPresets;

(*#################################################################################################### ChemDraw/smiles to 3D coordinate conversion function*)
obabel[inputType_, input_]:=Module[
	{charge,multiplicity,runQ},

(***************************************************************************************************** make sure the input file exists and generate temp file*)
	Which[
		inputType=="smiles",
		{
			runQ=True;
			CreateDirectory[Directory[]<>"\\xyz_temp"];
			Export[
			Directory[]<>"\\smiles_temp.txt",
				temp="";
				If[
					Count[input,x_]==0,
					input,
					{					
						For[
							i=1,
							i<Count[input,x_]+1,
							i++,
							temp=temp<>input[[i]]<>"\n"
						];
						temp
					}[[1]]
				](*If*),
				"Text"
			];(*Export*)
		},
		Or[inputType=="chemdraw",inputType=="smiles -f"],
		If[
			FileExistsQ[input],
			{
				runQ=True;
				CreateDirectory[Directory[]<>"\\xyz_temp"];
			},
			{
				runQ=False;
				Print["File not found"];
			}
		];
	];

(***************************************************************************************************** generate coordinates*)
	If[
		runQ,
		Which[
(***************************************************************************************************** generate coordinates from ChemDraw input*)
			inputType=="smiles",
			{
			(*generate open babel command*)
			babelCommand="babel "<>
				"-ismiles "<>Directory[]<>"\\smiles_temp.txt"<>
				" -ogjf "<>Directory[]<>"\\xyz_temp\\temp.txt -m -h --gen3d";
			(*generate smiles*)
			RunProcess[{"cmd","/c",babelCommand}];
			},
		
(***************************************************************************************************** generate coordinates from ChemDraw input*)
			inputType=="chemdraw",
			{
			(*generate open babel command*)
			babelCommand="babel -icdx "<>input<>" -osmiles "<>Directory[]<>"\\smiles_temp.txt";
			(*generate smiles*)
			RunProcess[{"cmd","/c",babelCommand}];	
			(*generate open babel command*)
			babelCommand="babel "<>
				"-ismiles "<>Directory[]<>"\\smiles_temp.txt"<>
				" -ogjf "<>Directory[]<>"\\xyz_temp\\temp.txt -m -h --gen3d";
			(*generate coordinates*)
			RunProcess[{"cmd","/c",babelCommand}];
			},
		
(***************************************************************************************************** generate coordinates from smiles input*)
			inputType=="smiles -f",
			{
			(*generate open babel command*)
			babelCommand="babel"<>
				" -ismiles "<>input<>
				" -ogjf "<>Directory[]<>"\\xyz_temp\\temp.txt -m -h --gen3d";
			(*generate coordinates*)
			RunProcess[{"cmd","/c",babelCommand}];
			}
		],(*Which*)
		Null	
	];(*If*)

(***************************************************************************************************** convert open babel output to gaussian format*)
	gather={};
	If[
		runQ,	
		For[
			i=1,
			i<Count[filelist[Directory[]<>"\\xyz_temp"],x_]+1,
			i++,
			{
			text=Import[filelist[Directory[]<>"\\xyz_temp"][[i]],"Text"];
			
			If[
				chargeOB=="Automatic",
				charge=Flatten[StringCases[text,"\n\n"~~x:NumberString~~Whitespace~~y:NumberString->x]][[1]],
				charge=chargeOB
			];(*If*)
	
			If[
				multiplicityOB=="Automatic",
				multiplicity=Flatten[StringCases[text,"\n\n"~~x:NumberString~~Whitespace~~y:NumberString->y]][[1]],
				multiplicity=multiplicityOB
			];(*If*)
				
			Which[
				StringContainsQ[text,"Multiplicity"],
				{
				newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
				section1=StringReplace[text,"#Put Keywords Here, check Charge and Multiplicity."~~Whitespace...~~NumberString~~Whitespace...~~NumberString~~"\n"->""];
				section2=StringReplace[section1,Whitespace...~~"(Iso="~~x:NumberString~~")"->x~~"      "];
				gather=Append[gather,newhead<>section2];
				}
			];(*Which*)
						
			};
			](*For*)
		];(*If*)
		
(***************************************************************************************************** remove temp files*)		
	Quiet[
	{
	If[
		FileExistsQ[Directory[]<>"\\xyz_temp"],
		{
		DeleteFile[filelist[Directory[]<>"\\xyz_temp"]];
		DeleteDirectory[Directory[]<>"\\xyz_temp", DeleteContents->True];
		}
	];
	
	If[
		FileExistsQ[Directory[]<>"\\smiles_temp.txt"],
		DeleteFile[Directory[]<>"\\smiles_temp.txt"]
	];
	}
	];
(***************************************************************************************************** output*)	
	gather
];(*obabel*)


OpenBabelPresets:={
	chargeOB=0,
	multiplicityOB=1
};
OpenBabelPresets;

(*#################################################################################################### ChemDraw/smiles to 3D coordinate conversion function*)
OpenBabel[fileName_,fileType_]:=Module[
	{charge,multiplicity},

(***************************************************************************************************** set directories and indices*)
	Which[
		fileType=="chemdraw",
		inputPath=Directory[]<>"\\data\\0_computation_data\\0_chemdraw_input\\",
		fileType=="smiles",
		inputPath=Directory[]<>"\\data\\0_computation_data\\0_smiles_input\\"
	];

	outputPath=Directory[]<>"\\data\\0_computation_data\\0_coordinate_input\\"<>fileName;

(***************************************************************************************************** check to see if a file with the name already exists*)
	If[
		FileExistsQ[Directory[]<>"\\data\\0_computation_data\\0_coordinate_input\\"<>fileName],
		runBatchQ=False,
		runBatchQ=True
	];

(***************************************************************************************************** make sure the input file exists*)	
	If[
		FileExistsQ[inputPath<>fileName<>Which[fileType=="chemdraw",".cdx",fileType=="smiles",".txt"]],
		Null,
		runBatchQ=False
	];

(***************************************************************************************************** generate coordinates*)
	If[
		runBatchQ,
		Which[
	
(***************************************************************************************************** generate coordinates from ChemDraw input*)
			fileType=="chemdraw",
			{
			(*make an output directory*)
			CreateDirectory[Directory[]<>"\\data\\0_computation_data\\0_coordinate_input\\"<>fileName];
			(*generate open babel command*)
			babelCommand="babel -icdx "<>inputPath<>fileName<>".cdx -osmiles "<>Directory[]<>"\\data\\0_computation_data\\0_smiles_input\\"<>fileName<>".txt";
			(*generate smiles*)
			RunProcess[{"cmd","/c",babelCommand}];	
			(*generate open babel command*)
			babelCommand="babel "<>
				"-ismiles "<>Directory[]<>"\\data\\0_computation_data\\0_smiles_input\\"<>fileName<>".txt"<>
				" -ogjf "<>Directory[]<>"\\data\\0_computation_data\\0_coordinate_input\\"<>fileName<>"\\"<>fileName<>".txt -m -h --gen3d";
			(*generate coordinates*)
			RunProcess[{"cmd","/c",babelCommand}];
			},
		
(***************************************************************************************************** generate coordinates from smiles input*)
			fileType=="smiles",
			{
			(*make an output directory*)
			CreateDirectory[Directory[]<>"\\data\\0_computation_data\\0_coordinate_input\\"<>fileName];
			(*generate open babel command*)
			babelCommand="babel -ismiles "<>inputPath<>fileName<>".txt -ogjf "<>Directory[]<>"\\data\\0_computation_data\\0_coordinate_input\\"<>fileName<>"\\"<>fileName<>".txt -m -h --gen3d";
			(*generate coordinates*)
			RunProcess[{"cmd","/c",babelCommand}];
			}
		],(*Which*)
		Null	
	];(*If*)

(***************************************************************************************************** convert open babel output to gaussian format*)
	If[
		runBatchQ,	
		For[
			i=1,
			i<Count[filelist[outputPath],x_]+1,
			i++,
			{
			text=Import[filelist[outputPath][[i]],"Text"];
			
			If[
				chargeOB=="Automatic",
				charge=Flatten[StringCases[text,"\n\n"~~x:NumberString~~Whitespace~~y:NumberString->x]][[1]],
				charge=chargeOB
			];(*If*)
	
			If[
				multiplicityOB=="Automatic",
				multiplicity=Flatten[StringCases[text,"\n\n"~~x:NumberString~~Whitespace~~y:NumberString->y]][[1]],
				multiplicity=multiplicityOB
			];(*If*)
				
			Which[
				StringContainsQ[text,"Multiplicity"],
				{
				newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
				section1=StringReplace[text,"#Put Keywords Here, check Charge and Multiplicity."~~Whitespace...~~NumberString~~Whitespace...~~NumberString~~"\n"->""];
				section2=StringReplace[section1,Whitespace...~~"(Iso="~~x:NumberString~~")"->x~~"      "];
				Export[filelist[outputPath][[i]],newhead<>section2,"Text"];
				}
			](*Which*)
			};
			](*For*)
		];(*If*)
];(*OpenBabel*)


ConfGenPresets:={
	chargeOB=0,
	multiplicityOB=1
};
ConfGenPresets;

(*#################################################################################################### ChemDraw/smiles to 3D coordinate conversion function*)
ConfGen[cdxFilePath_]:=Module[
	{charge,multiplicity,runBatchQ,temp,outputPath},

(***************************************************************************************************** set directories and indices*)
	outputPath=Directory[]<>"\\data\\0_computation_data\\0_coordinate_input\\"<>StringReplace[cdxFilePath,{Longest[__~~"\\"]->"",".cdx"->""}];

(***************************************************************************************************** check to see if a file with the name already exists*)
	If[
		FileExistsQ[Directory[]<>"\\data\\0_computation_data\\0_coordinate_input\\"<>StringReplace[cdxFilePath,{Longest[__~~"\\"]->"",".cdx"->""}]],
		runBatchQ=False,
		runBatchQ=True
	];

(***************************************************************************************************** make sure the input file exists and generate temp file*)
	If[
		And[FileExistsQ[cdxFilePath],runBatchQ],
		{
			CreateDirectory[outputPath];
			CreateDirectory[Directory[]<>"\\xyz_temp"];
		},
		{
			runBatchQ=False;
			Print["File not found"];
		}
	];

(***************************************************************************************************** generate coordinates*)
	If[
		runBatchQ,

		(
(***************************************************************************************************** generate smiles from ChemDraw input*)
			(*generate open babel command*)
			babelCommand="babel "<>
				"-icdx "<>cdxFilePath<>
				" -osmiles "<>Directory[]<>"\\temp.txt";
			(*generate coordinates*)
			RunProcess[{"cmd","/c",babelCommand}];
			
(***************************************************************************************************** generate 3D coordinate (sdf) files from smiles list*)
			(*generate open babel command*)
			babelCommand="babel "<>
				"-ismiles "<>Directory[]<>"\\temp.txt"<>
				" -osdf "<>Directory[]<>"\\xyz_temp\\"<>StringReplace[cdxFilePath,{Longest[__~~"\\"]->"",".cdx"->""}]<>".sdf -m -h --gen3d";
			(*generate coordinates*)
			RunProcess[{"cmd","/c",babelCommand}];
			
(***************************************************************************************************** generate conformer coordiantes using genetic algorithm*)
			With[
				{xyztemp=filelist[Directory[]<>"\\xyz_temp"]},
				Table[
					(
					(*generate open babel command*)
					babelCommand="babel "<>
						"-isdf "<>xyztemp[[n]]<>
						" -ogjf "<>StringReplace[xyztemp[[n]],".sdf"->".txt"]<>" --conformer --writeconformers";
					(*generate coordinates*)
					RunProcess[{"cmd","/c",babelCommand}];
					),
				{n,1,Count[xyztemp,_],1}	
				]
			];(*With*)
			
(***************************************************************************************************** Break multiconformer files into individual conformers*)
			With[
				{allConfs=filelist[Directory[]<>"\\xyz_temp"][[Flatten[Position[StringContainsQ[filelist[Directory[]<>"\\xyz_temp"],".txt"],True]]]]},
				Table[
					(
						temp=StringSplit[Import[allConfs[[n]]],"#Put Keywords Here, check Charge and Multiplicity."];
						Table[
							Export[
								outputPath<>"\\"<>StringReplace[allConfs[[n]],{Longest[__~~"\\"]->"",".txt"->""}]<>"_conf"<>ToString[m]<>".txt",
								"#Put Keywords Here, check Charge and Multiplicity."<>temp[[m]],
								"Text"
							],
							{m,1,Count[temp,_],1}
						](*Table*)
					),
					{n,1,Count[allConfs,_],1}
				](*Table*)	
			];(*With*)

			
		),
		
		Null	
	];(*If*)
	
(***************************************************************************************************** convert open babel output to gaussian format*)

	If[
		runBatchQ,	
		For[
			i=1,
			i<Count[filelist[outputPath],x_]+1,
			i++,
			{
			text=Import[filelist[outputPath][[i]],"Text"];
			
			If[
				chargeOB=="Automatic",
				charge=Flatten[StringCases[text,"\n\n"~~x:NumberString~~Whitespace~~y:NumberString->x]][[1]],
				charge=chargeOB
			];(*If*)
	
			If[
				multiplicityOB=="Automatic",
				multiplicity=Flatten[StringCases[text,"\n\n"~~x:NumberString~~Whitespace~~y:NumberString->y]][[1]],
				multiplicity=multiplicityOB
			];(*If*)
				
			Which[
				StringContainsQ[text,"Multiplicity"],
				{
				newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
				section1=StringReplace[text,"#Put Keywords Here, check Charge and Multiplicity."~~Whitespace...~~NumberString~~Whitespace...~~NumberString~~"\n"->""];
				section2=StringReplace[section1,Whitespace...~~"(Iso="~~x:NumberString~~")"->x~~"      "];
				Export[filelist[outputPath][[i]],newhead<>section2,"Text"];
				}
			](*Which*)
			};
			](*For*)
		];(*If*)
		
(***************************************************************************************************** Remove temp files*)
	If[
		runBatchQ==True,
		(
			DeleteDirectory[Directory[]<>"\\xyz_temp",DeleteContents->True];
			DeleteFile[Directory[]<>"\\temp.txt"];
		)
	];		

];(*ConfGen*)


(*#################################################################################################### Gaussian input generation funciton*)
(***************************************************************************************************** preset g16 calculation parameters*)
ComputationalTheoryPresets:={
	atomsPerProcessor=6,
	maxProcessors=20,
	ramPerProcessor=2,
	theory={"B3LYP","6-31G*","LANL2DZ"},
	maxAtomicNumberSupported=36,
	chargeCT="Automatic",
	multiplicityCT="Automatic"
};
ComputationalTheoryPresets;

(***************************************************************************************************** processor and ram allocation*)
nProcMem[position1_]:={
		processors=If[
			Ceiling[Count[position1,x_]/atomsPerProcessor]>maxProcessors,
			maxProcessors,
			Ceiling[Count[position1,x_]/atomsPerProcessor]],
		ram=processors*ramPerProcessor
};

(***************************************************************************************************** inputGenerator*)
GenerateGaussianInputText[coordinateFilePath_,nproc_,mem_]:=Module[
	{(*coordinateFile,moleculeSpecs,charge,multiplicity*)},

(***************************************************************************************************** Process coordinate files*)
	(*import coodinte file as text and remove atom labels*)
	coordinateFile=Import[coordinateFilePath,"Text"];
	
	If[
		chargeCT=="Automatic",
		charge=Flatten[StringCases[coordinateFile,"\n\n"~~Whitespace...~~x:NumberString~~Whitespace~~y:NumberString->x]][[1]],
		charge=chargeCT
	];
	
	If[
		multiplicityCT=="Automatic",
		multiplicity=Flatten[StringCases[coordinateFile,"\n\n"~~Whitespace...~~x:NumberString~~Whitespace~~y:NumberString->y]][[1]],
		multiplicity=multiplicityCT
	];
	
	(*remove unwanted input from generation if present*)
	Which[
		StringContainsQ[coordinateFile,"Title Card Required"],
		{
		newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
		section1=StringReplace[coordinateFile,{x__~~"\n\n"~~Whitespace...~~NumberString~~Whitespace..~~NumberString->"","\n\n"~~LetterCharacter~~Whitespace...~~y__->""}];
		coordinateFile=newhead<>section1;
		},
		StringContainsQ[coordinateFile,"--Link1--"],
		{
		newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
		section1=StringReplace[coordinateFile,"--Link1--"~~x__->""];
		section2=StringReplace[coordinateFile,x__~~"\n\n"~~NumberString~~Whitespace~~NumberString->""];
		coordinateFile=newhead<>section2;
		},
		StringContainsQ[coordinateFile,"\n\n"~~NumberString~~Whitespace~~NumberString],
		{
		newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
		section1=StringReplace[coordinateFile,x__~~"\n\n"~~NumberString~~Whitespace~~NumberString->""];
		coordinateFile=newhead<>section1;
		}
	];(*Which*)
	
	(*split input into charge, multiplicity, and coodinates*)
	moleculeSpecs=Flatten[
		StringCases[
			coordinateFile,
			"Title Card Required\n"~~"\n"~~x:NumberString~~WhitespaceCharacter..~~y:NumberString~~"\n"~~z___->{x,y,z}
			]
		];
	(*split into individual strings*)
	section1=StringSplit[moleculeSpecs[[3]]];
	(*find the positions of atoms*)
	position1=Flatten[
			Position[
			StringMatchQ[section1,LetterCharacter..~~___],
			True
			]
		];
	(*list of atoms*)
	atomlist=section1[[position1]];
	
(***************************************************************************************************** Write coordinate string*)
	(*Find max number of characters to make aesthetically pleasing*)
	maxSpacing=1;
	For[
		i=1,
		i<Count[section1,x_]+1,
		i++,
		{
		count=Count[StringPartition[section1[[i]],1],x_];
		If[count>maxSpacing,maxSpacing=count];
		}
	];
	coordinateSpacer[string_]:={string,Table[" ",{n,1,maxSpacing-Count[StringPartition[string,1],x_]+5,1}]};
	(*Write coordinate text*)
	coordinates="";
	For[
		i=1,
		i<Count[position1,x_]+1,
		i++,
		{
		newLine=StringJoin[Flatten[{
			coordinateSpacer[ToString[section1[[position1[[i]]]]]],
			coordinateSpacer[ToString[section1[[position1[[i]]+1]]]],
			coordinateSpacer[ToString[section1[[position1[[i]]+2]]]],
			coordinateSpacer[ToString[section1[[position1[[i]]+3]]]],
			"\n"
			}]];
		coordinates=coordinates<>newLine;
		}
		];
	
(***************************************************************************************************** processor and ram allocation*)
	processorsMemory=If[
		{nproc,mem}=={"Automatic","Automatic"},
		nProcMem[position1],
		{nproc,mem}
	];
	
(***************************************************************************************************** Determine if an ECP is required*)
	(*elements and atomic numbers*)
	elements=Transpose[{{"H","1"},{"He","2"},{"Li","3"},{"Be","4"},{"B","5"},{"C","6"},{"N","7"},{"O","8"},{"F","9"},{"Ne","10"},{"Na","11"},{"Mg","12"},{"Al","13"},{"Si","14"},{"P","15"},{"S","16"},{"Cl","17"},{"Ar","18"},{"K","19"},{"Ca","20"},
		{"Sc","21"},{"Ti","22"},{"V","23"},{"Cr","24"},{"Mn","25"},{"Fe","26"},{"Co","27"},{"Ni","28"},{"Cu","29"},{"Zn","30"},{"Ga","31"},{"Ge","32"},{"As","33"},{"Se","34"},{"Br","35"},{"Kr","36"},{"Rb","37"},{"Sr","38"},{"Y","39"},{"Zr","40"},
		{"Nb","41"},{"Mo","42"},{"Tc","43"},{"Ru","44"},{"Rh","45"},{"Pd","46"},{"Ag","47"},{"Cd","48"},{"In","49"},{"Sn","50"},{"Sb","51"},{"Te","52"},{"I","53"},{"Xe","54"},{"Cs","55"},{"Ba","56"},{"La","57"},{"Ce","58"},{"Pr","59"},{"Nd","60"},
		{"Pm","61"},{"Sm","62"},{"Eu","63"},{"Gd","64"},{"Tb","65"},{"Dy","66"},{"Ho","67"},{"Er","68"},{"Tm","69"},{"Yb","70"},{"Lu","71"},{"Hf","72"},{"Ta","73"},{"W","74"},{"Re","75"},{"Os","76"},{"Ir","77"},{"Pt","78"},{"Au","79"},{"Hg","80"},
		{"Tl","81"},{"Pb","82"},{"Bi","83"},{"Po","84"},{"At","85"},{"Rn","86"},{"Fr","87"},{"Ra","88"},{"Ac","89"},{"Th","90"},{"Pa","91"},{"U","92"},{"Np","93"},{"Pu","94"},{"Am","95"},{"Cm","96"},{"Bk","97"},{"Cf","98"},{"Es","99"},{"Fm","100"},
		{"Md","101"},{"No","102"},{"Lr","103"},{"Rf","104"},{"Db","105"},{"Sg","106"},{"Bh","107"},{"Hs","108"},{"Mt","109"},{"Ds","110"},{"Rg","111"},{"Cn","112"},{"Uut","113"},{"Fl","114"},{"Uup","115"},{"Lv","116"},{"Uus","117"},{"Uuo","118"}}];
	(*unique atoms in the job*)
	atoms=DeleteDuplicates[StringReplace[atomlist,NumberString->""]];
	(*atomic numbers for atoms in the job*)
	atomicNumbers={};
	For[
		i=1,
		i<Count[atoms,x_]+1,
		i++,
		{
		entry=Position[elements[[1]],atoms[[i]]];
		atomicNumbers=Flatten[Append[atomicNumbers,entry]];
		}
	];(*For*)
	(*sepparate into atoms supported and not supported by the primary basis set*)
	lightatoms={};
	heavyatoms={};
	atomSupportedQ[atomicnumber_]:=If[
		atomicnumber<=maxAtomicNumberSupported,
		lightatoms=Append[lightatoms,elements[[1]][[atomicnumber]]],
		heavyatoms=Append[heavyatoms,elements[[1]][[atomicnumber]]]
	];
	Map[atomSupportedQ,atomicNumbers];
		
(***************************************************************************************************** Generate basis set specific input files*)
	jobName=StringReplace[coordinateFilePath,{x__~~"\\"->"","."~~x__->""}];
	If[
	
(***************************************************************************************************** If there are no heavy atoms generate a standard input file*)
		Count[DeleteDuplicates[heavyatoms],x_]==0,
		{
		
(***************************************************************************************************** Generate input for the first calculation*)
		chk1=jobName<>"1.chk";
		textJob1=StringJoin[
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Chk=",chk1,"\n",
			"# opt=CalcFc "<>theory[[1]]<>"/"<>theory[[2]]<>" scf=xqc","\n",
			"\n","Title Card Required","\n\n",moleculeSpecs[[1]]," ",moleculeSpecs[[2]],"\n",coordinates,"\n\n"
		];
		
(***************************************************************************************************** Generate input for the second calculation*)
		chk2=jobName<>"2.chk";
		textJob2=StringJoin[
			"\n","--Link1--","\n",
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Oldchk=",chk1,"\n",
			"%Chk=",chk2,"\n",
			"# freq "<>theory[[1]]<>"/"<>theory[[2]]<>" volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"
		];
		
(***************************************************************************************************** Generate input for the third calculation*)
		chk3=jobName<>"3.chk";
		textJob3=StringJoin[
			"\n\n\n\n","--Link1--","\n",
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Oldchk=",chk2,"\n",
			"%Chk=",chk3,"\n",
			"# TD(NStates=10, Root=1) "<>theory[[1]]<>"/"<>theory[[2]]<>" volume pop=NPA density=current Geom=AllCheck Guess=Read"
		];
		},
		
(***************************************************************************************************** If there are heavy atoms generate a genecp input file*)
		{
		
(***************************************************************************************************** Generate input for the first calculation*)
		chk1=jobName<>"1.chk";
		textJob1=StringJoin[Flatten[{
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Chk=",chk1,"\n",
			"# opt=CalcFc "<>theory[[1]]<>"/genecp scf=xqc","\n",
			"\n","Title Card Required","\n\n",moleculeSpecs[[1]]," ",moleculeSpecs[[2]],"\n",coordinates,
			Flatten[{"\n",Thread[{lightatoms,Table[" ",{n,1,Count[lightatoms,x_],1}]}]}],"0",
			"\n",theory[[2]],"\n****",
			Flatten[{"\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]],"\n****",
			Flatten[{"\n\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]],"\n\n\n"
		}]];
		
(***************************************************************************************************** Generate input for the second calculation*)
		chk2=jobName<>"2.chk";
		textJob2=StringJoin[Flatten[{
			"\n","--Link1--","\n",
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Oldchk=",chk1,"\n",
			"%Chk=",chk2,"\n",
			"# freq "<>theory[[1]]<>"/genecp volume NMR pop=NPA density=current Geom=AllCheck Guess=Read",
			Flatten[{"\n\n",Thread[{lightatoms,Table[" ",{n,1,Count[lightatoms,x_],1}]}]}],"0",
			"\n",theory[[2]],"\n****",
			Flatten[{"\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]],"\n****",
			Flatten[{"\n\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]]
		}]];
		
(***************************************************************************************************** Generate input for the third calculation*)
		chk3=jobName<>"3.chk";
		textJob3=StringJoin[Flatten[{
			"\n\n\n\n","--Link1--","\n",
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Oldchk=",chk2,"\n",
			"%Chk=",chk3,"\n",
			"# TD(NStates=10, Root=1) "<>theory[[1]]<>"/genecp volume pop=NPA density=current Geom=AllCheck Guess=Read",
			Flatten[{"\n\n",Thread[{lightatoms,Table[" ",{n,1,Count[lightatoms,x_],1}]}]}],"0",
			"\n",theory[[2]],"\n****",
			Flatten[{"\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]],"\n****",
			Flatten[{"\n\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]],"\n\n\n"
			}]];
		}
	];(*If*)

(***************************************************************************************************** generate final input text*)
	textJob1<>textJob2<>textJob3<>"\n\n\n"
];(*inputgenerator*)


(*#################################################################################################### Gaussian input generation funciton*)
GenerateGaussianInputTextTSPresets:={
jobTypeTS="opt/freq"(*opt/freq,stepwise,IRC_forward,IRC_reverse*)
};
GenerateGaussianInputTextTSPresets;

(***************************************************************************************************** inputGenerator*)
GenerateGaussianInputTextTS[coordinateFilePath_,nproc_,mem_]:=Module[
	{coordinateFile,moleculeSpecs,charge,multiplicity},
	
(***************************************************************************************************** Process coordinate files*)
	(*import coodinte file as text and remove atom labels*)
	coordinateFile=Import[coordinateFilePath,"Text"];
	
	If[
		chargeCT=="Automatic",
		charge=Flatten[StringCases[coordinateFile,"\n\n"~~Whitespace...~~x:NumberString~~Whitespace~~y:NumberString->x]][[1]],
		charge=chargeCT
	];
	
	If[
		multiplicityCT=="Automatic",
		multiplicity=Flatten[StringCases[coordinateFile,"\n\n"~~Whitespace...~~x:NumberString~~Whitespace~~y:NumberString->y]][[1]],
		multiplicity=multiplicityCT
	];
	
	(*remove unwanted input from generation if present*)
	Which[
		StringContainsQ[coordinateFile,"Title Card Required"],
		{
		newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
		section1=StringReplace[coordinateFile,{x__~~"\n\n"~~Whitespace...~~NumberString~~Whitespace..~~NumberString->"","\n\n"~~LetterCharacter~~Whitespace...~~y__->""}];
		coordinateFile=newhead<>section1;
		},
		StringContainsQ[coordinateFile,"--Link1--"],
		{
		newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
		section1=StringReplace[coordinateFile,"--Link1--"~~x__->""];
		section2=StringReplace[coordinateFile,x__~~"\n\n"~~NumberString~~Whitespace~~NumberString->""];
		coordinateFile=newhead<>section2;
		},
		StringContainsQ[coordinateFile,"\n\n"~~NumberString~~Whitespace~~NumberString],
		{
		newhead="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n"<>ToString[charge]<>"  "<>ToString[multiplicity]<>"\n";
		section1=StringReplace[coordinateFile,x__~~"\n\n"~~NumberString~~Whitespace~~NumberString->""];
		coordinateFile=newhead<>section1;
		}
	];(*Which*)
	
	(*split input into charge, multiplicity, and coodinates*)
	moleculeSpecs=Flatten[
		StringCases[
			coordinateFile,
			"Title Card Required\n"~~"\n"~~x:NumberString~~WhitespaceCharacter..~~y:NumberString~~"\n"~~z___->{x,y,z}
			]
		];
	(*split into individual strings*)
	section1=StringSplit[moleculeSpecs[[3]]];
	(*find the positions of atoms*)
	position1=Flatten[
			Position[
			StringMatchQ[section1,LetterCharacter..~~___],
			True
			]
		];
	(*list of atoms*)
	atomlist=section1[[position1]];
	
(***************************************************************************************************** Write coordinate string*)
	(*Find max number of characters to make aesthetically pleasing*)
	maxSpacing=1;
	For[
		i=1,
		i<Count[section1,x_]+1,
		i++,
		{
		count=Count[StringPartition[section1[[i]],1],x_];
		If[count>maxSpacing,maxSpacing=count];
		}
	];
	coordinateSpacer[string_]:={string,Table[" ",{n,1,maxSpacing-Count[StringPartition[string,1],x_]+5,1}]};
	(*Write coordinate text*)
	coordinates="";
	For[
		i=1,
		i<Count[position1,x_]+1,
		i++,
		{
		newLine=StringJoin[Flatten[{
			coordinateSpacer[ToString[section1[[position1[[i]]]]]],
			coordinateSpacer[ToString[section1[[position1[[i]]+1]]]],
			coordinateSpacer[ToString[section1[[position1[[i]]+2]]]],
			coordinateSpacer[ToString[section1[[position1[[i]]+3]]]],
			"\n"
			}]];
		coordinates=coordinates<>newLine;
		}
		];
	
(***************************************************************************************************** processor and ram allocation*)
	processorsMemory=If[
		{nproc,mem}=={"Automatic","Automatic"},
		nProcMem[position1],
		{nproc,mem}
	];
	
(***************************************************************************************************** Determine if an ECP is required*)
	(*elements and atomic numbers*)
	elements=Transpose[{{"H","1"},{"He","2"},{"Li","3"},{"Be","4"},{"B","5"},{"C","6"},{"N","7"},{"O","8"},{"F","9"},{"Ne","10"},{"Na","11"},{"Mg","12"},{"Al","13"},{"Si","14"},{"P","15"},{"S","16"},{"Cl","17"},{"Ar","18"},{"K","19"},{"Ca","20"},
		{"Sc","21"},{"Ti","22"},{"V","23"},{"Cr","24"},{"Mn","25"},{"Fe","26"},{"Co","27"},{"Ni","28"},{"Cu","29"},{"Zn","30"},{"Ga","31"},{"Ge","32"},{"As","33"},{"Se","34"},{"Br","35"},{"Kr","36"},{"Rb","37"},{"Sr","38"},{"Y","39"},{"Zr","40"},
		{"Nb","41"},{"Mo","42"},{"Tc","43"},{"Ru","44"},{"Rh","45"},{"Pd","46"},{"Ag","47"},{"Cd","48"},{"In","49"},{"Sn","50"},{"Sb","51"},{"Te","52"},{"I","53"},{"Xe","54"},{"Cs","55"},{"Ba","56"},{"La","57"},{"Ce","58"},{"Pr","59"},{"Nd","60"},
		{"Pm","61"},{"Sm","62"},{"Eu","63"},{"Gd","64"},{"Tb","65"},{"Dy","66"},{"Ho","67"},{"Er","68"},{"Tm","69"},{"Yb","70"},{"Lu","71"},{"Hf","72"},{"Ta","73"},{"W","74"},{"Re","75"},{"Os","76"},{"Ir","77"},{"Pt","78"},{"Au","79"},{"Hg","80"},
		{"Tl","81"},{"Pb","82"},{"Bi","83"},{"Po","84"},{"At","85"},{"Rn","86"},{"Fr","87"},{"Ra","88"},{"Ac","89"},{"Th","90"},{"Pa","91"},{"U","92"},{"Np","93"},{"Pu","94"},{"Am","95"},{"Cm","96"},{"Bk","97"},{"Cf","98"},{"Es","99"},{"Fm","100"},
		{"Md","101"},{"No","102"},{"Lr","103"},{"Rf","104"},{"Db","105"},{"Sg","106"},{"Bh","107"},{"Hs","108"},{"Mt","109"},{"Ds","110"},{"Rg","111"},{"Cn","112"},{"Uut","113"},{"Fl","114"},{"Uup","115"},{"Lv","116"},{"Uus","117"},{"Uuo","118"}}];
	(*unique atoms in the job*)
	atoms=DeleteDuplicates[StringReplace[atomlist,NumberString->""]];
	(*atomic numbers for atoms in the job*)
	atomicNumbers={};
	For[
		i=1,
		i<Count[atoms,x_]+1,
		i++,
		{
		entry=Position[elements[[1]],atoms[[i]]];
		atomicNumbers=Flatten[Append[atomicNumbers,entry]];
		}
	];(*For*)
	(*sepparate into atoms supported and not supported by the primary basis set*)
	lightatoms={};
	heavyatoms={};
	atomSupportedQ[atomicnumber_]:=If[
		atomicnumber<=maxAtomicNumberSupported,
		lightatoms=Append[lightatoms,elements[[1]][[atomicnumber]]],
		heavyatoms=Append[heavyatoms,elements[[1]][[atomicnumber]]]
	];
	Map[atomSupportedQ,atomicNumbers];
	
(***************************************************************************************************** Generate basis set specific input files*)
	jobName=StringReplace[coordinateFilePath,{x__~~"\\"->"","."~~x__->""}];
	If[
	
(***************************************************************************************************** If there are no heavy atoms generate a standard input file*)
		Count[DeleteDuplicates[heavyatoms],x_]==0,
		{
		
(***************************************************************************************************** Generate input for the first calculation*)
		chk1=jobName<>"1.chk";
		textJob1=StringJoin[
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Chk=",chk1,"\n",
			"# opt=(calcfc,ts,noeigentest) scf=xqc "<>theory[[1]]<>"/"<>theory[[2]]<>"\n",
			"\n","Title Card Required","\n\n",moleculeSpecs[[1]]," ",moleculeSpecs[[2]],"\n",coordinates,"\n\n"
		];
		
(***************************************************************************************************** Generate input for the second calculation*)
		chk2=jobName<>"2.chk";
		textJob2=StringJoin[
			"\n","--Link1--","\n",
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Oldchk=",chk1,"\n",
			"%Chk=",chk2,"\n",
			"# freq "<>theory[[1]]<>"/"<>theory[[2]]<>" volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"
		];
		
		},
		
(***************************************************************************************************** If there are heavy atoms generate a genecp input file*)
		{
		
(***************************************************************************************************** Generate input for the first calculation*)
		chk1=jobName<>"1.chk";
		textJob1=StringJoin[Flatten[{
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Chk=",chk1,"\n",
			"# opt=(calcfc,ts,noeigentest) scf=xqc "<>theory[[1]]<>"/genecp","\n",
			"\n","Title Card Required","\n\n",moleculeSpecs[[1]]," ",moleculeSpecs[[2]],"\n",coordinates,
			Flatten[{"\n",Thread[{lightatoms,Table[" ",{n,1,Count[lightatoms,x_],1}]}]}],"0",
			"\n",theory[[2]],"\n****",
			Flatten[{"\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]],"\n****",
			Flatten[{"\n\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]],"\n\n\n"
		}]];
		
(***************************************************************************************************** Generate input for the second calculation*)
		chk2=jobName<>"2.chk";
		textJob2=StringJoin[Flatten[{
			"\n","--Link1--","\n",
			"%nprocshared=",ToString[processorsMemory[[1]]],"\n",
			"%Mem=",ToString[processorsMemory[[2]]],"GB\n",
			"%Oldchk=",chk1,"\n",
			"%Chk=",chk2,"\n",
			"# freq "<>theory[[1]]<>"/genecp volume NMR pop=NPA density=current Geom=AllCheck Guess=Read",
			Flatten[{"\n\n",Thread[{lightatoms,Table[" ",{n,1,Count[lightatoms,x_],1}]}]}],"0",
			"\n",theory[[2]],"\n****",
			Flatten[{"\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]],"\n****",
			Flatten[{"\n\n",Thread[{heavyatoms,Table[" ",{n,1,Count[heavyatoms,x_],1}]}]}],"0",
			"\n",theory[[3]]
		}]];
		}
	];(*If*)
	
(***************************************************************************************************** generate final input text*)

(***************************************************************************************************** opt/freq versus stepwise control*)	
	Which[
		jobTypeTS=="stepwise",
		textJob1<>textJob2<>"\n",
		jobTypeTS=="opt/freq",
		StringReplace[textJob1,"opt=(calcfc,ts,noeigentest) scf=xqc"->"opt=(calcfc,ts,noeigentest) scf=xqc freq"],
		jobTypeTS=="IRC_forward",
		StringReplace[textJob1,"opt=(calcfc,ts,noeigentest) scf=xqc"->"IRC=(CalcFC,forward,maxpoints=10,stepsize=10)"],
		jobTypeTS=="IRC_reverse",
		StringReplace[textJob1,"opt=(calcfc,ts,noeigentest) scf=xqc"->"IRC=(CalcFC,reverse,maxpoints=10,stepsize=10)"]
	](*Which*)
];(*inputgenerator*)


(*#################################################################################################### Coordinate input generation from output files*)
(***************************************************************************************************** *)
ExtractCoordinates[logFilePath_]:=Module[
{},
	coordinates="";
	If[
(***************************************************************************************************** Run control*)
		FileExistsQ[logFilePath],
		(
		
(***************************************************************************************************** Find which output sections are present*)		
			logFileText=Import[logFilePath,"Text"];
			
			rootSections=StringCases[
				logFileText,
				Shortest["#"~~p1__~~EndOfLine~~p2__~~EndOfLine~~p3__~~EndOfLine],
				Overlaps->True
				];
			
			possibleSectionKeywords={
				{"opt","o"~~Whitespace...~~"p"~~Whitespace...~~"t"},
				{"freq","f"~~Whitespace...~~"r"~~Whitespace...~~"e"~~Whitespace...~~"q"},
				{"td","t"~~Whitespace...~~"d"~~Whitespace...~~"("}
			};
			
			outputSections={};
			
			For[
				i=1,
				i<Count[possibleSectionKeywords,x_]+1,
				i++,
				(
					If[
						Count[StringContainsQ[rootSections,possibleSectionKeywords[[i]][[2]],IgnoreCase->True],True]==2,
						outputSections=Append[outputSections,possibleSectionKeywords[[i]][[1]]]						
					](*If*)
				)
			];
(***************************************************************************************************** Extract atom labels from first input section*)
			With[
				{
				outputText=StringCases[
					logFileText,
					Shortest[rootSections[[1]]~~__~~"Normal termination"],
					Overlaps->True						
				][[1]](*StringCases*)
				},
				atomLabels=Transpose[
					Partition[
						StringSplit[
							StringReplace[
								StringCases[outputText,__~~Shortest["Multiplicity = "~~NumberString..~~p1__~~"\n\n"]~~__->p1][[1]],
								LetterCharacter..~~WhitespaceCharacter..~~LetterCharacter..~~__->""
							]
						],
						4
					]
				][[1]];	
			];(*With*)

(***************************************************************************************************** Extract optimized coordinates from prefered output sections*)

			elements=Transpose[{{"H","1"},{"He","2"},{"Li","3"},{"Be","4"},{"B","5"},{"C","6"},{"N","7"},{"O","8"},{"F","9"},{"Ne","10"},{"Na","11"},{"Mg","12"},{"Al","13"},{"Si","14"},{"P","15"},{"S","16"},{"Cl","17"},{"Ar","18"},{"K","19"},{"Ca","20"},
				{"Sc","21"},{"Ti","22"},{"V","23"},{"Cr","24"},{"Mn","25"},{"Fe","26"},{"Co","27"},{"Ni","28"},{"Cu","29"},{"Zn","30"},{"Ga","31"},{"Ge","32"},{"As","33"},{"Se","34"},{"Br","35"},{"Kr","36"},{"Rb","37"},{"Sr","38"},{"Y","39"},{"Zr","40"},
				{"Nb","41"},{"Mo","42"},{"Tc","43"},{"Ru","44"},{"Rh","45"},{"Pd","46"},{"Ag","47"},{"Cd","48"},{"In","49"},{"Sn","50"},{"Sb","51"},{"Te","52"},{"I","53"},{"Xe","54"},{"Cs","55"},{"Ba","56"},{"La","57"},{"Ce","58"},{"Pr","59"},{"Nd","60"},
				{"Pm","61"},{"Sm","62"},{"Eu","63"},{"Gd","64"},{"Tb","65"},{"Dy","66"},{"Ho","67"},{"Er","68"},{"Tm","69"},{"Yb","70"},{"Lu","71"},{"Hf","72"},{"Ta","73"},{"W","74"},{"Re","75"},{"Os","76"},{"Ir","77"},{"Pt","78"},{"Au","79"},{"Hg","80"},
				{"Tl","81"},{"Pb","82"},{"Bi","83"},{"Po","84"},{"At","85"},{"Rn","86"},{"Fr","87"},{"Ra","88"},{"Ac","89"},{"Th","90"},{"Pa","91"},{"U","92"},{"Np","93"},{"Pu","94"},{"Am","95"},{"Cm","96"},{"Bk","97"},{"Cf","98"},{"Es","99"},{"Fm","100"},
				{"Md","101"},{"No","102"},{"Lr","103"},{"Rf","104"},{"Db","105"},{"Sg","106"},{"Bh","107"},{"Hs","108"},{"Mt","109"},{"Ds","110"},{"Rg","111"},{"Cn","112"},{"Uut","113"},{"Fl","114"},{"Uup","115"},{"Lv","116"},{"Uus","117"},{"Uuo","118"}}
			];

			Which[
(***************************************************************************************************** *)
(***************************************************************************************************** If a freq calculation was carried out*)
(***************************************************************************************************** *)
				MemberQ[outputSections,"freq"],
				Catch[
					outputText=StringCases[
						logFileText,
						Shortest[
							rootSections[[
								Flatten[
									Position[
										StringContainsQ[rootSections,"f"~~Whitespace...~~"r"~~Whitespace...~~"e"~~Whitespace...~~"q",IgnoreCase->True],
									True]
								][[-2]]
							]]~~
							x__~~"Normal termination"
						],
						Overlaps->True						
					][[1]];(*StringCases*)
					
(***************************************************************************************************** First pass *)					
					templist=StringCases[
						outputText,
						Shortest["Standard orientation:"~~p1__~~"Rotational constants"]->p1
					];

					templist=If[
						Count[templist,_]>0,
						StringCases[
							templist[[-1]],
							NumberString~~Whitespace..~~
							p1:NumberString~~Whitespace..~~
							NumberString~~Whitespace..~~
							p2:NumberString~~Whitespace..~~
							p3:NumberString~~Whitespace..~~
							p4:NumberString~~Whitespace..
							->{p1,p2,p3,p4}
						](*StringCases*)
					];(*If*)

					coordinateList={};
					If[
						Count[templist[[1]],_]==4,
						For[
							i=1,
							i<Count[templist,_]+1,
							i++,
							coordinateList=Append[coordinateList,Join[elements[[1]][[Flatten[Position[elements[[2]],templist[[i]][[1]]]]]],templist[[i]][[2;;]]]]
						](*For*)
					];(*If*)

					coordinates="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n "<>
								ToString[DeleteDuplicates[IntegerPart[ToExpression[Flatten[StringCases[outputText,"Charge"~~Whitespace...~~"="~~Whitespace...~~p1:NumberString->p1]]]]][[1]]]<>
								" "<>
								ToString[DeleteDuplicates[IntegerPart[ToExpression[Flatten[StringCases[outputText,"Multiplicity"~~Whitespace...~~"="~~Whitespace...~~p1:NumberString->p1]]]]][[1]]]<>
								"\n";

(***************************************************************************************************** Second pass *)
					If[
						Count[DeleteDuplicates[Table[Count[coordinateList[[n]],_],{n,1,Count[coordinateList,_],1}]],_]>1,
						(
							templist=Flatten[
								StringCases[
									outputText,
									Shortest[
										"Charge"~~Whitespace...~~"="~~Whitespace...~~p1:NumberString..~~Whitespace...~~
										"Multiplicity"~~Whitespace...~~"="~~Whitespace...~~p2:NumberString..~~
										p3__~~"Grad"
									](*Shortest*)
									->{p1,p2,p3},
									Overlaps->True
								](*StringCases*)
							];(*Flatten*)
						
							coordinateList=If[
								Count[templist,_]==3,
								StringCases[
									templist[[3]],
									p1:WordCharacter~~Whitespace...~~","~~Whitespace...~~
									NumberString~~Whitespace...~~","~~Whitespace...~~
									p2:NumberString~~Whitespace...~~","~~Whitespace...~~
									p3:NumberString~~Whitespace...~~","~~Whitespace...~~
									p4:NumberString~~Whitespace...
									->{p1,p2,p3,p4}
								](*StringCases*)
							];(*If*)	
							
							coordinates="%chk=blank.chk\n# hf/3-21g\n\nTitle Card Required\n\n "<>templist[[1]]<>" "<>templist[[2]]<>"\n";
							
						)
					];(*If*)

(***************************************************************************************************** Add atom labels from the input section *)
					coordinateList=Transpose[Prepend[Transpose[coordinateList][[2;;]],atomLabels]];
					
(***************************************************************************************************** Collect *)
					If[
						Count[DeleteDuplicates[Table[Count[coordinateList[[n]],_],{n,1,Count[coordinateList,_],1}]],_]==1,
						For[
							i=1,
							i<Count[coordinateList,x_]+1,
							i++,
							coordinates=coordinates<>"  "<>coordinateList[[i]][[1]]<>"\t\t"<>coordinateList[[i]][[2]]<>"\t\t"<>coordinateList[[i]][[3]]<>"\t\t"<>coordinateList[[i]][[4]]<>"\n"
						],(*For*)
						coordinates="Error in coordinate extraction. Check .log file."
					];(*If*)	
						
				](*Catch*),

(***************************************************************************************************** *)				
(***************************************************************************************************** If a td calculation was carried out*)
(***************************************************************************************************** *)
				MemberQ[outputSections,"td"],
				(
					outputText=StringCases[
						logFileText,
						Shortest[
							rootSections[[
								Flatten[
									Position[
										StringContainsQ[rootSections,"t"~~Whitespace...~~"d"~~Whitespace...~~"(",IgnoreCase->True],
									True]
								][[-2]]
							]]~~
							x__~~"Normal termination"
						],
						Overlaps->True						
					][[1]];(*StringCases*)
					
					
				),

(***************************************************************************************************** *)				
(***************************************************************************************************** If an opt calculation was carried out*)
(***************************************************************************************************** *)
				MemberQ[outputSections,"opt"],
				(
					outputText=StringCases[
						logFileText,
						Shortest[
							rootSections[[
								Flatten[
									Position[
										StringContainsQ[rootSections,"o"~~Whitespace...~~"p"~~Whitespace...~~"t",IgnoreCase->True],
									True]
								][[-2]]
							]]~~
							x__~~"Normal termination"
						],
						Overlaps->True						
					][[1]];(*StringCases*)
				),
				True,
				Print["Error: calculation type not recognized."]
			];(*If*)
			
			
		)
	];(*If*)
	
	coordinates
];(*Module*)


(*#################################################################################################### Unix command file generator*)
(***************************************************************************************************** preset g16 calculation parameters*)
unixHostSpec=StringSplit[Import[Directory[]<>"\\control\\1_host_spec.txt","Text"],"\n"];
UnixCommandPresets:={
	clusterName=unixHostSpec[[1]],
	clusterPort=unixHostSpec[[2]],
	clusterPassword=unixHostSpec[[3]],
	userID=unixHostSpec[[4]],
	jobWallTime="23:59:00"
};
UnixCommandPresets;

(***************************************************************************************************** Generate unix command file text*)
GenerateCMDText[gaussianInputPath_,wallTime_]:=Module[
	{fileName,processorCount},
	
(***************************************************************************************************** Extract processor specification from input file*)
	processorCount=Flatten[
		StringCases[
			StringReplace[Import[gaussianInputPath,"Text"],"Title Card Required"~~y___->""],
			"%nprocshared="~~x__~~"\n%Mem="~~y__->x
			]
		][[1]];

(***************************************************************************************************** Generate command file text*)
	fileName=StringReplace[gaussianInputPath,{x__~~"\\"->"","."~~x__->""}]<>".com";
	StringJoin[
		"#!/bin/sh",
		"\n#SBATCH -N 1",
		"\n#SBATCH --ntasks-per-node=",
		processorCount,
		"\n#SBATCH -t ",
		wallTime,
		"\n#SBATCH -C haswell",
		"\nmkdir -p /scratch/",
		userID,
		"/job.$$",
		"\nexport GAUSS_SCRDIR=/scratch/",
		userID,
		"/job.$$",
		"\n\ng16 ",
		fileName,(*submit the Gaussian calculation*)
		"\n\nbash echo 'Job Complete' >",
		StringReplace[fileName,".com"->".done"],(*when the Gaussian calculation is done create a .done file*)
		"\n\nrm ",
		StringReplace[fileName,"123.com"->"1.chk "],(*then delete the .chk files because they take up a lot of space*)
		StringReplace[fileName,"123.com"->"2.chk "],
		StringReplace[fileName,"123.com"->"3.chk"]
		](*StringJoin*)
	](*Module*)



GenerateGaussianBatchPresets={
	jobType="equilibirum"
};
GenerateGaussianBatchPresets;

(*#################################################################################################### Job geneartion function*)
GenerateGaussianBatch[coordinateDirectoryPath_]:=Quiet[Module[
	{runBatchQ},
(***************************************************************************************************** set directories and indices*)
	coordinateDirectoryName=StringReplace[coordinateDirectoryPath,x__~~"\\"->""];
	outputDirectoryName=coordinateDirectoryName;

(***************************************************************************************************** make sure the input directory exists*)	
	If[
		FileExistsQ[coordinateDirectoryPath],
		runBatchQ=True,
		runBatchQ=False
	];

(***************************************************************************************************** make sure the submission directory doesn't already exist*)		
	If[
		FileExistsQ[Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName],
		(*If true give it a new name*)
		{
		name=outputDirectoryName;
		ChoiceDialog[
			Panel[Row[{"New name: ",InputField[Dynamic[name],String,FieldMasked->False]}],FrameMargins->8,ImageMargins->5],
			WindowSize->{300,400},
			WindowTitle->"Batch name is taken!"
			];
		outputDirectoryName=ToString[name];
		}
	];

(***************************************************************************************************** generate submission file*)	
	If[
		runBatchQ,
		{
	
(***************************************************************************************************** Make a working directory for submission and retreival of jobs*)
		CreateDirectory[Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName];
	
(***************************************************************************************************** gaussian .com file from coordinates*)
		jobList=filelist[coordinateDirectoryPath];
		Table[
			Export[
				Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName<>"\\"<>StringReplace[jobList[[i]],{x__~~"\\"->"","."~~x__->""}]<>"123.com",
				Which[
					jobType=="equilibrium",
					GenerateGaussianInputText[jobList[[i]],"Automatic","Automatic"],
					jobType=="transition_state",
					GenerateGaussianInputTextTS[jobList[[i]],"Automatic","Automatic"]
					](*Which*),
				"Text"
			],
			{i,1,Count[jobList,x_],1}
		];

(***************************************************************************************************** generate unix command files*)
		jobList=filelist[Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName];
		jobList[[Flatten[Position[StringContainsQ[jobList,".com"],True]]]];
			
		Table[
			Export[
				Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName<>"\\"<>StringReplace[jobList[[i]],{x__~~"\\"->"","."~~x__->""}]<>".cmd",
				GenerateCMDText[jobList[[i]],jobWallTime],
				"Text",
				DOSTextFormat->False
			],
			{i,1,Count[jobList,x_],1}
		];
	
(***************************************************************************************************** generate .sh submission script*)
		jobList=filelist[Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName];
		jobList=jobList[[Flatten[Position[StringContainsQ[jobList,".cmd"],True]]]];
			
		shScriptText="";
		For[
			i=1,
			i<Count[jobList,x_]+1,
			i++,
			shScriptText=shScriptText<>"sbatch "<>StringReplace[jobList[[i]],{x__~~"\\"->""}]<>"\n"
		];
		Export[
			Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName<>"\\"<>coordinateDirectoryName<>".sh",
			shScriptText,
			"Text",
			DOSTextFormat->False
		];
	
(***************************************************************************************************** generate SFTP commands - just in case of manual submission*)
		Export[
			Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName<>"\\"<>coordinateDirectoryName<>"_SFTP.txt",
			StringJoin[
				"put -r ",StringReplace[Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName,"\\"->"/"]," ",outputDirectoryName,
				"\n\n\n",
				"get -r ",outputDirectoryName," ",StringReplace[Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>outputDirectoryName,"\\"->"/"]
			],
			"Text"
		]
		};
	];(*If*)
]](*Module*)


(*#################################################################################################### Add encoded smiles strings to the input title cards*)
InsertTitleCode[batchFolderPath_,smilesListPath_]:=Module[
	{
		smiles=StringSplit[Import[smilesListPath,"Text"],"\n"],
		files=({(StringCases[#,Longest[x:NumberString..]~~"123.com"->x][[1]]//ToExpression),#}&/@Select[filelist[batchFolderPath],StringContainsQ[#,".com"]&]//Sort)[[All,2]]
	},
	
	If[
		Count[smiles,_]==Count[files,_],
		Table[
			(
				newtext=Import[files[[n]],"Text"]//StringReplace["Title Card Required"->"SmilesCodeStart "<>(({#//ToString," "}&/@(smiles[[n]]//ToCharacterCode))//Flatten//StringJoin)<>"SmilesCodeStop"];
				Export[files[[n]],newtext,"Text"];
			),
			{n,1,Count[files,_],1}
		],
		Print["Error: file and string lists not compatible"]
	];
	
];


(*#################################################################################################### Automate submission file generation*)
GenerateJobs[arguement_]:=Quiet[Module[
	{ChemDraw,Smiles,Gaussian},
(***************************************************************************************************** run any outstanding jobs for files in the input folders*)
	ChemDraw:={
		chemdrawJobs=StringReplace[filenames[Directory[]<>"\\data\\0_computation_data\\0_chemdraw_input\\"],".cdx"->""];
	
		Quiet[Table[
			OpenBabel[chemdrawJobs[[i]],"chemdraw"],
			{i,1,Count[chemdrawJobs,x_],1}
		]];
	};
	
	Smiles:={
		smilesJobs=StringReplace[filenames[Directory[]<>"\\data\\0_computation_data\\0_smiles_input\\"],".txt"->""];
	
		Quiet[Table[
			OpenBabel[smilesJobs[[i]],"smiles"],
			{i,1,Count[smilesJobs,x_],1}
		]];
	};

	Gaussian:={		
(***************************************************************************************************** get list of potential jobs*)
		coordinateDirectories=filelist[Directory[]<>"\\data\\0_computation_data\\0_coordinate_input"];
	
(***************************************************************************************************** see if a submission file with the same batch name already exists*)
		Table[
				{
				If[
					FileExistsQ[Directory[]<>"\\data\\0_computation_data\\1_submission\\"<>StringReplace[coordinateDirectories[[i]],x__~~"\\"->""]],
					Null,
					GenerateGaussianBatch[coordinateDirectories[[i]]]
				](*If*)
				},
				{i,1,Count[coordinateDirectories,x_],1}
			];(*Table*)
	};
	
	Which[
		arguement=="chemdraw",
		ChemDraw,
		arguement=="smiles",
		Smiles,
		arguement=="gaussian",
		Gaussian,
		arguement=="all",
		{ChemDraw, Smiles, Gaussian}
	];
(***************************************************************************************************** run arguemnt*)	
]];


(*#################################################################################################### SSH control function - run a list of commands*)
(***************************************************************************************************** presets for paramiko*)
unixHostSpec=StringSplit[Import[Directory[]<>"\\control\\1_host_spec.txt","Text"],"\n"];
ParamikoPresets:={
	paramikoAuthentication="key",
	clusterName=unixHostSpec[[1]],
	clusterPort=unixHostSpec[[2]],
	clusterAuthentication=unixHostSpec[[3]],
	userID=unixHostSpec[[4]],
	jobWallTime="23:59:00"
};
ParamikoPresets;

(***************************************************************************************************** SSH function*)
RunCommand[commands_]:=Module[
{},
(*commands are in the format {command 1, command 2, .... command n}*)
(***************************************************************************************************** check to see if paramiko is working*)

	(*If true continue*)
	
	(*If false break*)
	
(***************************************************************************************************** password or .pem key file*)
	Which[
	paramikoAuthentication=="key",
	{
		authentication=StringReplace["paramiko.RSAKey.from_private_key_file(\""<>Directory[]<>"\\aws\\.ssh\\"<>clusterAuthentication<>"\")","\\"->"/"];
		runQ=True;
	},
	paramikoAuthentication=="password",
	{
		authentication=clusterAuthentication;
		runQ=True;
	},
	True,
	runQ=False
	];
	
(***************************************************************************************************** Generate Paramiko Script*)
	If[
	runQ,
	paramikoScript=StringJoin[
		"import paramiko #sftp and ssh package",
		"\nimport os #operating system interface package",
		"\n\n################################### host information and file specificaiton",
		"\n\n#host information",
		"\nhost = ",
		clusterName,
		"\nport = ",
		clusterPort,
		"\nauth = ",
		authentication,
		"\nusername = ",
		userID,
		"\n\n################################### #Execute commands",
		"\ncommand = ",
		"'",StringJoin[Flatten[Table[{commands[[i]],"\\n"},{i,1,Count[commands,x_],1}]]],"'",
		"\nnbytes = 4096",
		"\n\nclient = paramiko.Transport((host, port))",
		"\nclient.connect(username=username,",If[paramikoAuthentication=="key"," pkey=auth)"," password=auth)"],
		"\n\nstdout_data = []",
		"\nstderr_data = []",
		"\nsession = client.open_channel(kind='session')",
		"\nsession.exec_command(command)     #commands submitted",
		"\nwhile True:",
		"\n    if session.recv_ready():",
		"\n        stdout_data.append(session.recv(nbytes))",
		"\n    if session.recv_stderr_ready():",
		"\n        stderr_data.append(session.recv_stderr(nbytes))",
		"\n    if session.exit_status_ready():",
		"\n        break",
		"\n\nsession.close()",
		"\nclient.close()"
	];
	
	];(*If*)

(***************************************************************************************************** Write paramiko script*)
	Export[
	Directory[]<>"\\submission\\cluster_command.py",
	paramikoScript,
	"Text"
	];

(***************************************************************************************************** Generate bat file for paramiko script*)
	Export[
	Directory[]<>"\\submission\\cluster_command.bat",
	anacondaPythonDir<>" "<>Directory[]<>"\\submission\\cluster_command.py %*",
	"Text"
	];

(***************************************************************************************************** Run bat file*)
	RunProcess[{"cmd","/c",Directory[]<>"\\submission\\cluster_command.bat"}];
]


(*#################################################################################################### Molecular Vibration Extraction Functions*)
(***************************************************************************************************** Extract frequency related data*)
FrequencyData[logFilePath_]:=Catch[
	If[
(***************************************************************************************************** Run control*)
	FileExistsQ[logFilePath],
	(

(***************************************************************************************************** Import log file text and extract freq section*)
	frequencyData={};
	logFileText=Import[logFilePath,"Text"];
	
	If[
		StringContainsQ[logFileText,"Harmonic frequencies"~~__~~"- Thermochemistry"],	
		freqSection=StringCases[
			logFileText,
			{p1__~~"Harmonic frequencies"~~p2__~~"- Thermochemistry -"~~p3__->p2}
			],
		freqSection={}
	];	
	
(***************************************************************************************************** Find number of atoms in molecule*)		
	nAtoms=ToExpression[
		DeleteDuplicates[
			StringCases[
				logFileText,
				"NAtoms="~~Whitespace...~~x:NumberString~~Whitespace->x
			]
		][[1]]
	];(*ToExpression*)

	If[
		nAtoms>1,
		(
(***************************************************************************************************** If there is more than one frequency section throw an error*)
		If[
			Count[freqSection,x_]==1,
			freqSection=StringCases[
				freqSection[[1]],
				p1__~~"normal coordinates:"~~p2__->p2
			][[1]],
			(
			Throw[Print["Error: log file has "<>ToString[Count[freqSection,x_]]<>" frequency sections."]]
			)
		];(*If*)

(***************************************************************************************************** Get a list of frequencies*)
		frequencyList=StringSplit[
			Flatten[
				StringCases[
					StringSplit[freqSection,"Atom"],
					"Frequencies --"~~x__~~"Red. masses"->x
				]
			],
			Whitespace..
		];(*StringSplit*)

(***************************************************************************************************** If there are frequencies carry on, otherwise throw an error*)
		If[
			Count[Flatten[{frequencyList}],x_]>1,
			(

(***************************************************************************************************** Extract data from possible headings*)		
			frequencyData1=
					Flatten[
						StringCases[
							StringSplit[freqSection,"Atom"],
								" Frequencies --"~~Whitespace...~~
								p1:NumberString~~Whitespace...~~
								p2:NumberString~~Whitespace...~~
								p3:NumberString~~Whitespace...~~
								rest__->
								{p1,p2,p3}
							](*StringCases*)
						];(*Flatten*)
					
			frequencyData2=
					Flatten[
						StringCases[
							StringSplit[freqSection,"Atom"],
								" Red. masses --"~~Whitespace...~~
								p1:NumberString~~Whitespace...~~
								p2:NumberString~~Whitespace...~~
								p3:NumberString~~Whitespace...~~
								rest__->
								{p1,p2,p3}
							](*StringCases*)
						];(*Flatten*)
					
			frequencyData3=
					Flatten[
						StringCases[
							StringSplit[freqSection,"Atom"],
								" Frc consts  --"~~Whitespace...~~
								p1:NumberString~~Whitespace...~~
								p2:NumberString~~Whitespace...~~
								p3:NumberString~~Whitespace...~~
								rest__->
								{p1,p2,p3}
							](*StringCases*)
						];(*Flatten*)

			frequencyData4=
					Flatten[
						StringCases[
							StringSplit[freqSection,"Atom"],
								" IR Inten    --"~~Whitespace...~~
								p1:NumberString~~Whitespace...~~
								p2:NumberString~~Whitespace...~~
								p3:NumberString~~Whitespace...~~
								rest__->
								{p1,p2,p3}
							](*StringCases*)
						];(*Flatten*)

			frequencyData5=
					Flatten[
						StringCases[
							StringSplit[freqSection,"Atom"],
								" Dip. str.   --"~~Whitespace...~~
								p1:NumberString~~Whitespace...~~
								p2:NumberString~~Whitespace...~~
								p3:NumberString~~Whitespace...~~
								rest__->
								{p1,p2,p3}
							](*StringCases*)
						];(*Flatten*)

			frequencyData6=
					Flatten[
						StringCases[
							StringSplit[freqSection,"Atom"],
								" Rot. str.   --"~~Whitespace...~~
								p1:NumberString~~Whitespace...~~
								p2:NumberString~~Whitespace...~~
								p3:NumberString~~Whitespace...~~
								rest__->
								{p1,p2,p3}
							](*StringCases*)
						];(*Flatten*)

			frequencyData7=
					Flatten[
						StringCases[
							StringSplit[freqSection,"Atom"],
								" E-M angle   --"~~Whitespace...~~
								p1:NumberString~~Whitespace...~~
								p2:NumberString~~Whitespace...~~
								p3:NumberString~~Whitespace...~~
								rest__->
								{p1,p2,p3}
							](*StringCases*)
						];(*Flatten*)
					
			frequencyDataProspects=Transpose[
				{
				{"Frequencies","Red. masses","Frc consts","IR Inten","Dip. str.","Rot. str.","E-M angle"},
				{frequencyData1,frequencyData2,frequencyData3,frequencyData4,frequencyData5,frequencyData6,frequencyData7}
				}
			];

(***************************************************************************************************** Select only headings which have data (if not present it will have 0 entries)*)								
			frequencyData={};
			frequencyDataLabels={};
			For[
				i=1,
				i<8,
				i++,
				(
					If[
						Count[frequencyDataProspects[[i]][[2]],x_]!=0,
						(
							frequencyData=Append[frequencyData,frequencyDataProspects[[i]][[2]]];
							frequencyDataLabels=Append[frequencyDataLabels,frequencyDataProspects[[i]][[1]]]
						)
					](*If*)
				)
			];(*For*)

(***************************************************************************************************** Make sure all components have the same numer of entries, otherwise throw an error*)		
			If[
				Count[DeleteDuplicates[Table[Count[frequencyData[[n]],x_],{n,1,Count[frequencyData,x_],1}]],x_]==1,
				(
					frequencyData=Transpose[frequencyData];
					frequencyTable=Prepend[
						frequencyData,
						frequencyDataLabels
					](*Prepend*)
				),
				Throw[Print["Error: output vectors unequal length."]]
			](*If*)
		),
		Throw[Print["Error: could not parse frequency table."]]
		](*If*)
	),
	Throw[Print["Error: 1 atom has no vibrations."]]
	]
	),
	Throw[Print["Error: file not found."]]
	](*If*)
](*Catch*);

(*#################################################################################################### Extract atomic movement vectors*)
MovementVectors[logFilePath_]:=Catch[
	If[
(***************************************************************************************************** Run control*)
	FileExistsQ[logFilePath],
	(
(***************************************************************************************************** Import log file text and extract freq section*)
	logFileText=Import[logFilePath,"Text"];
	Clear[out2];
	
	freqSection=StringCases[
		logFileText,
		{p1__~~"Harmonic frequencies"~~p2__~~"- Thermochemistry -"~~p3__->p2}
	];
	
(***************************************************************************************************** Find number of atoms in molecule*)		
	nAtoms=ToExpression[
		DeleteDuplicates[
			StringCases[
				logFileText,
				"NAtoms="~~Whitespace...~~x:NumberString~~Whitespace->x
			]
		][[1]]
	];(*ToExpression*)

	If[
	nAtoms>1,
	(
(***************************************************************************************************** If there is more than one frequency section throw an error*)
	If[
		Count[freqSection,x_]==1,
		freqSection=StringCases[
			freqSection[[1]],
			p1__~~"normal coordinates:"~~p2__->p2
			][[1]],
		Throw[Print["Error: log file has "<>ToString[Count[freqSection,x_]]<>" frequency sections."]]
	];(*If*)

(***************************************************************************************************** Parse data and find positions of Z coordinates*)
	textMatrix=StringSplit[freqSection,Whitespace..];

	positionIndex=Flatten[Position[textMatrix,"Z"]];

(***************************************************************************************************** Extract data from text matrix*)
	in=Transpose[Partition[positionIndex,3]][[3]];
	out2={};
	
	For[
		j=1,
		j<nAtoms+1,
		j++,
		(
			indices={};
			For[
				i=1,
				i<Count[in,x_]+1,
				i++,
				(
					indices=Flatten[Append[
					indices,
					Table[in[[i]]+n+11*(j-1),{n,3,11,1}]
					]];
					out1=Join[textMatrix[[{-2,-1}+indices[[1]]]],textMatrix[[indices]]]
				)
			];(*For*)
			out2=Append[out2,out1];
		)
	];(*For*)

	Prepend[
		out2,
		Flatten[
			{"Atom","Atomic_Number",
			Table[
				{"X"<>ToString[n],"Y"<>ToString[n],"Z"<>ToString[n]},
				{n,1,Quotient[Count[out2[[1]],x_],3],1}
			](*Table*)
			}
		](*Flatten*)
	](*Prepend*)
	),
	Throw[Print["Error: 1 atom has no vibrations."]]
	]
	),
	Throw[Print["Error: file not found."]]
	](*If*)

];(*Catch*)

(*#################################################################################################### Transform coordinates A to overlap with coordinates B*)
(*
Find shared atoms in a molecular structure set and compute vibrational similarity versus a standard labeled molecule
*)

(***************************************************************************************************** Molecular alignment function from Derek*)
rigidTransform[A_,B_]:=Module[
	{
(***************************************************************************************************** Compute the unweighted centroids*)
		centroidA=Map[Mean,Transpose[A]],
		centroidB=Map[Mean,Transpose[B]],
		H,U,s,V,R,t
	},
	
(***************************************************************************************************** Center points*)

(*not necessary for this application*)

(***************************************************************************************************** Find rotation marix*)
	H=Covariance[A,B];
	{U,s,V}=SingularValueDecomposition[H];
	R=V.Transpose[U];
	
(***************************************************************************************************** Special reflection case*)

(*not necessary for this application*)


(***************************************************************************************************** Compute translation vector*)
	t=((-R).centroidA)+centroidB;

	{R,t}
];(*rigidTransform*)

(*#################################################################################################### Compare vibrations for two labeled molecules*)

(*find correlations between mol1 and mol2 where mol1 is rotated and translated to overlap with mol2*)

(***************************************************************************************************** Preset search parameters*)
CompareVibrationsPresets:={getLabels={"auto"},corrThreshold=0.5,freqThreshold={500,500},noDuplicateVibrations=False};
CompareVibrationsPresets;

(***************************************************************************************************** Function*)
CompareVibrations[mol1path_,mol2path_]:=Module[
	{c1,c2,at1,at2,mv1,mv2,r2,R,t},

(***************************************************************************************************** Extract overlapping coordinates molecule 1*)
	ExtractCoordinates[mol1path];
	c1=coordinateList;(*Global variable from ExtractCoordinates*)
	
(***************************************************************************************************** Set of shared atoms which will be used to calculate vibrational correlations*)
	at1=Sort[
			With[
				{list={Range[Count[atomLabels,_]],StringReplace[atomLabels,__~~p1:NumberString->p1]},labels=Map[ToString,getLabels]},
				If[
					getLabels=={"auto"},
					(*get all numerically labeled atoms*)
					ToExpression[Select[Map[Flatten,Transpose[{list[[1]],StringCases[NumberString..][list[[2]]]}]],Count[#,_]==2&]],
					(*get only numerically labeled in getLabels vector*)
					ToExpression[Select[Map[Flatten,Table[If[MemberQ[labels,list[[2]][[n]]],Transpose[list][[n]],{}],{n,1,Count[list[[2]],_],1}]],Count[#,_]==2&]]
				](*If*)
			],
			#1[[2]]<#2[[2]]&
	];(*Sort*)
	
(***************************************************************************************************** Set of shared atoms which will be used align the molecules*)	
	alignment1=Sort[
		With[
			{list={Range[Count[atomLabels,_]],StringReplace[atomLabels,__~~p1:NumberString->p1]},labels=Map[ToString,getLabels]},
			ToExpression[Select[Map[Flatten,Transpose[{list[[1]],StringCases[NumberString..][list[[2]]]}]],Count[#,_]==2&]]
			],
		#1[[2]]<#2[[2]]&
	];
	
	Clear[atomLabels];
	
(***************************************************************************************************** Extract overlapping molecule 2*)
	ExtractCoordinates[mol2path];
	c2=coordinateList;(*Global variable from ExtractCoordinates*)
	
(***************************************************************************************************** Set of shared atoms which will be used to calculate vibrational correlations*)
	at2=Sort[
			With[
				{list={Range[Count[atomLabels,_]],StringReplace[atomLabels,__~~p1:NumberString->p1]},labels=Map[ToString,getLabels]},
				If[
					getLabels=={"auto"},
					(*get all numerically labeled atoms*)
					ToExpression[Select[Map[Flatten,Transpose[{list[[1]],StringCases[NumberString..][list[[2]]]}]],Count[#,_]==2&]],
					(*get only numerically labeled in getLabels vector*)
					ToExpression[Select[Map[Flatten,Table[If[MemberQ[labels,list[[2]][[n]]],Transpose[list][[n]],{}],{n,1,Count[list[[2]],_],1}]],Count[#,_]==2&]]
				](*If*)
			],
			#1[[2]]<#2[[2]]&
	];
(***************************************************************************************************** Set of shared atoms which will be used align the molecules*)	
	alignment2=Sort[
		With[
			{list={Range[Count[atomLabels,_]],StringReplace[atomLabels,__~~p1:NumberString->p1]},labels=Map[ToString,getLabels]},
			ToExpression[Select[Map[Flatten,Transpose[{list[[1]],StringCases[NumberString..][list[[2]]]}]],Count[#,_]==2&]]
			],
		#1[[2]]<#2[[2]]&
	];
	
	Clear[atomLabels];

(***************************************************************************************************** If molecules have labeled overlapping atoms*)
	If[
		Transpose[alignment1][[2]]==Transpose[alignment2][[2]],
		(
(***************************************************************************************************** Extract overlapping atom coordinates*)
		coords1=ToExpression[Transpose[Transpose[c1[[Transpose[alignment1][[1]]]]][[2;;]]]];
		coords2=ToExpression[Transpose[Transpose[c2[[Transpose[alignment2][[1]]]]][[2;;]]]];

(***************************************************************************************************** Find rotation matrix and translation vector*)
		{R,t}=rigidTransform[coords1,coords2];
		coords1=Transpose[R.Transpose[coords1]+t];

(***************************************************************************************************** Get vibrational movement vectors*)
		MovementVectors[mol1path];
		mv1=ToExpression[Transpose[Transpose[out2][[3;;]]][[Transpose[at1][[1]]]]];(*out2 is a global variable from MovementVectors*)
		mv1=Transpose[Table[Partition[mv1[[n]],3],{n,1,Count[mv1,_],1}]];
		mv1=Table[Flatten[Transpose[(R.Transpose[mv1[[n]]])]],{n,1,Count[mv1,_],1}];(*Translation is not necessary for vibraitonal correlation*)mv11=mv1;at11=at1;
		Clear[out2];
		
		MovementVectors[mol2path];
		mv2=ToExpression[Transpose[Transpose[out2][[3;;]]][[Transpose[at2][[1]]]]];mv21=mv2;
		mv2=Transpose[Table[Partition[mv2[[n]],3],{n,1,Count[mv2,_],1}]];
		mv2=Map[Flatten,mv2];mv23=mv2;
		Clear[out2];

(***************************************************************************************************** Calculate vibrational correlations and filter by threshold*)
		r2=Flatten[
			Table[
				If[And[StandardDeviation[mv1[[n]]]!=0.,StandardDeviation[mv2[[m]]]!=0.],{n,m,Correlation[mv1[[n]],mv2[[m]]]},{n,m,0.}],
				{n,1,Count[mv1,_],1},
				{m,1,Count[mv2,_],1}
			],
		1];

		r2=Complement[Table[If[Abs[r2[[n]][[3]]]>corrThreshold,r2[[n]]],{n,1,Count[r2,_],1}],{Null}];

(***************************************************************************************************** Filter by frequency value*)
		With[
			{fd1=FrequencyData[mol1path],fd2=FrequencyData[mol2path]},
			
			r2=Complement[
				Table[
					If[
						And[ToExpression[fd1[[r2[[n]][[1]]+1]][[1]]]>freqThreshold[[1]],
							ToExpression[fd2[[r2[[n]][[2]]+1]][[1]]]>freqThreshold[[1]],
							Abs[(ToExpression[fd1[[r2[[n]][[1]]+1]][[1]]]-ToExpression[fd2[[r2[[n]][[2]]+1]][[1]]])]<freqThreshold[[2]]
						],
						r2[[n]]
					],
					{n,1,Count[r2,_],1}
				],
				{Null}
			];

(***************************************************************************************************** For duplicate frequencies choose only the highest correlation*)
			If[
				noDuplicateVibrations==True,
				Null
			];(*If*)
			
(***************************************************************************************************** Package results*)
			Table[
				Join[r2[[n]],ToExpression[fd1[[r2[[n]][[1]]+1]]],ToExpression[fd2[[r2[[n]][[2]]+1]]]],
				{n,1,Count[r2,_],1}
			]
		]

	),
	Print["Error: labeled atoms don't match"]

	](*If*)

];

(*#################################################################################################### Find shared vibrations for a whole set of molecules*)
(*Find a set of shared vibrations between a set of molecules and a standard*)

SharedVibrations[setPaths_,standardPath_]:=Module[

{standardPosition,temp,list,vibrationList,allCorrelated,selected,correlatedSet,tempIterator},

(***************************************************************************************************** Compute vibration comparision array*)
	compareVibrations=ParallelTable[
		{StringReplace[setPaths[[n]],{Longest[__~~"\\"]->"",".log"->"","123.log"->""}],CompareVibrations[setPaths[[n]],standardPath]},
		{n,1,Count[setPaths,_],1}
	];
	
(***************************************************************************************************** Get list of correlated vibrations indexed by the standard*)
	list=Map[Transpose,Transpose[compareVibrations][[2]]];
	list=Table[list[[n]][[2]],{n,1,Count[list,_],1}];
	
(***************************************************************************************************** Find the intersection of all sets of indexed vibrations*)
	intersection=Intersection@@list;
	
(***************************************************************************************************** If the intersection is greater than zero list the vibrations in the index*)	
	If[
		Count[intersection,_]>0,
		(
(***************************************************************************************************** get the overlapping vibrations for each set in intersection*)	
			allCorrelated=Table[
				{
					compareVibrations[[n]][[1]],
					compareVibrations[[n]][[2]][[Flatten[Position[compareVibrations[[n]][[2]],Flatten[{__,intersection[[m]],Table[__,Range[15]]}]]]]]
				},
				{m,1,Count[intersection,_],1},
				{n,1,Count[compareVibrations,_],1}
			];		
(***************************************************************************************************** Sort the overlapping vibrations by correlation and choose the highest correlation*)				
			selected=Table[Table[Flatten[{allCorrelated[[m]][[n]][[1]],Sort[allCorrelated[[m]][[n]][[2]],Abs[#1[[3]]]>Abs[#2[[3]]]&][[1]]}],{n,1,Count[allCorrelated[[m]],_],1}],{m,1,Count[intersection,_],1}];

(***************************************************************************************************** Remove duplicate correlations if present*)	
			headers={"file_name","vibration","standard_vibration","correlation","frequency","reduced_mass","frc_const","IR_intensity","dip_strength","rot_strength","E-M_angle",
				"standard_frequency","standard_reduced_mass","standard_frc_const","standard_IR_intensity","standard_dip_strength","standard_rot_strength","standard_E-M_angle"	
				};
			
			headers=Table[
				("vib_"<>ToString[n]<>"_"<>#)&/@(headers[[2;;]]),
				{n,1,Count[selected,_],1}
			]//Flatten//Prepend[headers[[1]]];
									
			Table[
				Flatten[Table[selected[[n]][[m]][[2;;]],{n,2,Count[selected,_],1}]//Prepend[selected[[1]][[m]]]],
				{m,1,Count[selected[[1]],_],1}
			]//Prepend[headers]
			
		),
		Print["No similar vibrations with getLabels = "<>ToString[getLabels]<>", corrThreshold = "<>ToString[corrThreshold]<>", freqThreshold = "<>ToString[freqThreshold]]
	](*If*)
];(*Shared vibrations*)



(*#################################################################################################### %Buried Volume Calculator*)
(*
Takes gaussian log files or coordinates as input. The center atom can be specified as a unique element or as a 
numerically labeled atom.
*)
(***************************************************************************************************** Function presets*)
OccupiedVolumePresets:={
	setRadius=3.0,(*Angstroms*)
	meshCount="automatic",
	occupiedVolumeInput="log",
	selectMatrixPosition=False
};
OccupiedVolumePresets;

(***************************************************************************************************** Function*)
OccupiedVolume[input_,atom_]:=Module[
	{allCoordinates,atomPosition,radii,distances,selectedAtoms,mesh,gridSpacing,occupiedMesh},

(***************************************************************************************************** Van der Waals radii*)
	vanDerWaals={{"H",1.2`},{"Li",2.2`},{"Be",1.9`},{"B",1.8`},{"C",1.7`},{"N",1.6`},{"O",1.55`},{"F",1.5`},{"Na",2.4`},{"Mg",2.2`},
		{"Al ",2.1`},{"Si",2.1`},{"P",1.95`},{"S",1.8`},{"Cl",1.8`},{"K",2.8`},{"Ca",2.4`},{"Sc",2.3`},{"Ti",2.15`},{"V",2.05`},
		{"Cr",2.05`},{"Mn",2.05`},{"Fe",2.05`},{"Co",2.`},{"Ni",2.`},{"Cu",2.`},{"Zn",2.1`},{"Ga",2.1`},{"Ge",2.1`},{"As",2.05`},
		{"Se",1.9`},{"Br",1.9`},{"Rb",2.9`},{"Sr",2.55`},{"Y",2.4`},{"Zr",2.3`},{"Nb",2.15`},{"Mo",2.1`},{"Tc",2.05`},{"Ru",2.05`},
		{"Rh",2.`},{"Pd",2.05`},{"Ag",2.1`},{"Cd",2.2`},{"In",2.2`},{"Sn",2.25`},{"Sb",2.25`},{"Te",2.2`},{"I",2.1`},{"Cs",2.1`},
		{"Ba",2.1`},{"La",2.5`},{"Hf",2.25`},{"Ta",2.2`},{"W",2.1`},{"Re",2.05`},{"Os",2.2`},{"Ir",2.`},{"Pt",2.05`},{"Au",2.1`},
		{"Hg",2.05`},{"Tl",2.2`},{"Pb",2.3`},{"Bi",2.3`},{"Th",2.4`},{"U",2.3`}};
		
(***************************************************************************************************** Input type*)
	Which[
		occupiedVolumeInput=="log",
		(
			ExtractCoordinates[input];
			allCoordinates=ToExpression[Map[Rest,coordinateList]];
		),
		occupiedVolumeInput=="array",
		(
			coordinateList=input;
			allCoordinates=ToExpression[Map[Rest,coordinateList]];
		)
	];

(***************************************************************************************************** Find central atom*)
	Which[
		selectMatrixPosition==False,
		atomPosition=Flatten[Position[StringContainsQ[Transpose[coordinateList][[1]],ToString[atom]],True]],
		selectMatrixPosition==True,
		atomPosition={atom},
		True,
		atomPosition={}
	];

(***************************************************************************************************** If a unique atom is found carry out the calculation*)
	If[
		Count[atomPosition,_]==1,
		(

(***************************************************************************************************** Coordinate data for central atom and molecule*)
			atomPosition=allCoordinates[[atomPosition[[1]]]];

			radii=Transpose[vanDerWaals][[2]][[
				Flatten[
					Table[
						Position[Transpose[vanDerWaals][[1]],StringReplace[Transpose[coordinateList][[1]][[n]],NumberString..->""]],
						{n,1,Count[coordinateList,_],1}
					]
				]
			]];

(***************************************************************************************************** Remove atoms which are not within the sphere from calculation*)
			(*compute distances from the central atom to each other atom*)
			distances=Table[
				{n,radii[[n]],Sqrt[SquaredEuclideanDistance[atomPosition,allCoordinates[[n]]]]*1.},
				{n,1,Count[allCoordinates,_],1}
			];
			
			(*find all atoms to include in the calculation*)
			selectedAtoms=Transpose[Select[distances,#[[3]]<=#[[2]]+setRadius&]][[1]];
			
			(*set coordinates and radii*)
			allCoordinates=allCoordinates[[selectedAtoms]];
			radii=radii[[selectedAtoms]];

(***************************************************************************************************** Define a region and an integration grid*)
			(*set integration mesh if not specified*)
			If[meshCount=="automatic",gridSpacing=setRadius/10,gridSpacing=setRadius/meshCount];

			(*integration grid is a cube of points around the central atom of spacing gridSpacing*)
			mesh=Flatten[
				Table[
					{x,y,z},
					{x,atomPosition[[1]]-setRadius,atomPosition[[1]]+setRadius,gridSpacing},
					{y,atomPosition[[2]]-setRadius,atomPosition[[2]]+setRadius,gridSpacing},
					{z,atomPosition[[3]]-setRadius,atomPosition[[3]]+setRadius,gridSpacing}
				],
			2];

(***************************************************************************************************** Prune the grid to define the sphere*)
			(*compute distances from the central atom to each integration cube*)
			distances=Table[
				{n,Sqrt[SquaredEuclideanDistance[atomPosition,mesh[[n]]]]*1.},
				{n,1,Count[mesh,_],1}
			];
		
			(*redefine the mesh by the cubes inside the sphere*)
			mesh=mesh[[Transpose[Select[distances,#[[2]]<=setRadius&]][[1]]]];
			
(***************************************************************************************************** Find all cubes within the Van der Waals radii of each atom*)
			occupiedMesh={};
			For[
				i=1,
				i<Count[allCoordinates,_]+1,
				i++,
				(
					distances=Table[{n,Sqrt[SquaredEuclideanDistance[allCoordinates[[i]],mesh[[n]]]]},{n,1,Count[mesh,_],1}];
					output=Select[distances,#[[2]]<=radii[[i]]&];
					If[
						Count[output,_]>0,
						occupiedMesh=DeleteDuplicates[Join[occupiedMesh,Transpose[output][[1]]]]
					];
				)
			];

(***************************************************************************************************** Compute occupied volume*)
			volume=(gridSpacing^3*Count[occupiedMesh,_]);(*volume*)
			
			(gridSpacing^3*Count[occupiedMesh,_])/(gridSpacing^3*Count[mesh,_])
			
		),
		Print["Unique atom "<>ToString[atom]<>" not found."]
	]

];


(*#################################################################################################### Descriptor Extraction Function*)
(***************************************************************************************************** Function presets*)
ExtractDescriptorsPresets:={
	"ground_state",
	"excited_state",
	"atom_specific"
};
ExtractDescriptorsPresets;

(***************************************************************************************************** Function*)
ExtractDescriptors[logFilePath_]:=Catch[
	If[
(***************************************************************************************************** Run control*)
		FileExistsQ[logFilePath],
		(
		
(***************************************************************************************************** Clear descriptors*)			
			descriptorVector={};
			atomDescriptorMatrix={};
			atomDescriptorMatrixLabeled={};
			
(***************************************************************************************************** Find which output sections are present*)			
			logFileText=Import[logFilePath,"Text"];
			
			If[
			StringContainsQ[logFileText,"Normal"],
			
			(
			ExtractCoordinates[logFilePath];
			
			rootSections=StringCases[
				logFileText,
				Shortest["#"~~p1__~~EndOfLine~~p2__~~EndOfLine~~p3__~~EndOfLine],
				Overlaps->True
				];
			
			possibleSectionKeywords={
				{"opt","o"~~Whitespace...~~"p"~~Whitespace...~~"t"},
				{"freq","f"~~Whitespace...~~"r"~~Whitespace...~~"e"~~Whitespace...~~"q"},
				{"td","t"~~Whitespace...~~"d"~~Whitespace...~~"("}
			};
			
			outputSections={};
			
			For[
				i=1,
				i<Count[possibleSectionKeywords,x_]+1,
				i++,
				(
					If[
						Count[StringContainsQ[rootSections,possibleSectionKeywords[[i]][[2]],IgnoreCase->True],True]==2,
						outputSections=Append[outputSections,possibleSectionKeywords[[i]][[1]]]						
					](*If*)
				)
			];
			),
			
			outputSections={};
			];
			
(***************************************************************************************************** For each output section find which descriptors have been computed*)

(***************************************************************************************************** *)
(***************************************************************************************************** opt output*)
(***************************************************************************************************** *)
			If[
				MemberQ[outputSections,"opt"],
				(
					possibleOptKeywords={
						{"opt_coordinates","Stationary point found"}
					};
				)
			];(*If*)
			
(***************************************************************************************************** *)
(***************************************************************************************************** freq output*)
(***************************************************************************************************** *)
			If[
				MemberQ[outputSections,"freq"],
				
				(		
					freqSection=StringCases[
						logFileText,
						Shortest[
							rootSections[[
								Flatten[
									Position[
										StringContainsQ[rootSections,"f"~~Whitespace...~~"r"~~Whitespace...~~"e"~~Whitespace...~~"q",IgnoreCase->True],
									True]
								][[-2]]
							]]~~
							x__~~"Normal termination"
						],
						Overlaps->True						
					][[1]];(*StringCases*)
					
					If[
						AnyTrue[StringContainsQ[{freqSection},"Normal termination"],TrueQ],
						Null,
						(
							Throw[Print["Error: frequency calculation not sucessfully completed."]];
							freqSection="Calculation not complete.";
						)
					];
			
					possibleFreqKeywords={
						{"molar_mass","Molar Mass ="},
						{"molar_volume","Molar volume ="},
						{"fmo","Population analysis using the SCF density."},
						{"thermochem","Thermochemistry"},
						{"electronic_spatial_extent","Electronic spatial extent"},
						{"dipole","Dipole moment (field-independent basis, Debye)"}
					};
					
					
					If[
						MemberQ[
							ExtractDescriptorsPresets,
							"atom_specific"
						],
						possibleFreqKeywords=Join[
							possibleFreqKeywords,
							{
								{"Mulliken_population","Mulliken charges"},
								{"APT","APT charges:"},
								{"NPA","Summary of Natural Population Analysis:"},
								{"nmr","Magnetic shielding tensor (ppm):"}
							}
						]
					];(*If*)
				
					
					freqKeywords={};
					For[
						i=1,
						i<Count[possibleFreqKeywords,x_]+1,
						i++,
						(
							If[
								StringContainsQ[freqSection,possibleFreqKeywords[[i]][[2]]],
								freqKeywords=Append[freqKeywords,possibleFreqKeywords[[i]][[1]]]
							](*If*)
						)
					];(*For*)
					
					
(***************************************************************************************************** generic search and value calculation functions*)	
					valueSearchQ[descriptorName_,propertyValue_]:=
						If[
							MemberQ[freqKeywords,descriptorName],
							(
								descriptorVector=Append[
									descriptorVector,
									{
									descriptorName,
									If[
										Count[propertyValue,x_]==1,
										If[
											StringContainsQ[ToString[propertyValue[[1]]],LetterCharacter],
											propertyValue[[1]],
											ToExpression[propertyValue[[1]]]
										],
										"NA"	
									](*If*)
									}
								];(*Append*)
								
								descriptorVector=DeleteDuplicates[descriptorVector]
							)
						];(*If*)
						
					valueSearch[descriptorName_,propertyValue_]:=(
						descriptorVector=Append[
							descriptorVector,
							{
							descriptorName,
							If[
								Count[propertyValue,x_]==1,
								If[
									StringContainsQ[ToString[propertyValue[[1]]],LetterCharacter],
									propertyValue[[1]],
									ToExpression[propertyValue[[1]]]
								],
								"NA"	
							](*If*)
							}
						];(*Append*)
								
						descriptorVector=DeleteDuplicates[descriptorVector];
					);

					addComputedValue[descriptorName_,propertyValue_]:=(
						descriptorVector=Append[
							descriptorVector,
							{
							descriptorName,
							If[
								Count[propertyValue,x_]==1,
								propertyValue[[1]],
								"NA"	
							](*If*)
							}
						];(*Append*)
								
						descriptorVector=DeleteDuplicates[descriptorVector];
					);										
					
(***************************************************************************************************** stoichiometry*)					
					valueSearch[
						"stoichiometry",
						DeleteDuplicates[Flatten[StringCases[freqSection,Shortest["Stoichiometry"~~Longest[Whitespace...]~~p1__~~"\n"]->p1,Overlaps->True]]]
					];

(***************************************************************************************************** number_of_atoms*)					
					valueSearch[
						"number_of_atoms",
						DeleteDuplicates[Flatten[StringCases[freqSection,"NAtoms="~~Whitespace...~~p1:NumberString->p1]]]
					];
(***************************************************************************************************** charge*)					
					valueSearch[
						"charge",
						DeleteDuplicates[IntegerPart[ToExpression[Flatten[StringCases[freqSection,"Charge"~~Whitespace...~~"="~~Whitespace...~~p1:NumberString->p1]]]]]
					];

(***************************************************************************************************** multiplicity*)					
					valueSearch[
						"multiplicity",
						DeleteDuplicates[IntegerPart[ToExpression[Flatten[StringCases[freqSection,"Multiplicity"~~Whitespace...~~"="~~Whitespace...~~p1:NumberString->p1]]]]]
					];

(***************************************************************************************************** convergence*)					
					valueSearch[
						"convergence_criteria",
						If[
							Count[temp=Flatten[StringCases[freqSection,Shortest["Maximum"~~Whitespace...~~"Force"~~p1__~~Whitespace..~~"Predicted"]->p1]],_]==1,
									If[
										StringContainsQ[temp[[1]],"NO",IgnoreCase->True],
										(
											convTable=Map[StringSplit,Flatten[StringCases[temp,p1:NumberString~~Whitespace..~~p2:NumberString..~~Whitespace..~~WordCharacter..]]];
											If[
												AllTrue[
													{
													ToExpression[convTable[[1]][[1]]]<0.00001,
													ToExpression[convTable[[2]][[1]]]<0.00001
													},
													TrueQ
												],
												{"met"},
												{"failed"}
											](*If*)
										),
										{"met"}
									],
							{"NA"}
						](*If*)
					];

(***************************************************************************************************** dipole*)					
					valueSearchQ[
						"dipole",
						Flatten[StringCases[freqSection,Shortest["Dipole moment (field-independent basis, Debye):"~~__~~"Tot="~~Whitespace...~~p1:NumberString]->p1,Overlaps->True]]
					];
									
(***************************************************************************************************** molar mass*)					
					valueSearchQ[
						"molar_mass",
						Flatten[StringCases[freqSection,"Molar Mass ="~~WhitespaceCharacter..~~p1:NumberString->p1]]
					];
					
(***************************************************************************************************** molar volume*)											
					valueSearchQ[
						"molar_volume",
						Flatten[StringCases[freqSection,"Molar volume ="~~WhitespaceCharacter..~~x:NumberString->x]]
					];

(***************************************************************************************************** electronic_spatial_extent*)						
					valueSearchQ[
						"electronic_spatial_extent",
						Flatten[StringCases[freqSection,"Electronic spatial extent (au):  <R**2>="~~Whitespace...~~x:NumberString->x]]
					];
					
(***************************************************************************************************** homo, lumo, and related properties*)				
					templist=Flatten[
						StringSplit[
							Flatten[
								StringCases[freqSection,"Population analysis using the SCF density."~~x__~~"Condensed to atoms (all electrons):"->x]
							](*Flatten*)
						](*StringSplit*)
					];(*Flatten*)
					
					If[
						MemberQ[freqKeywords,"fmo"],
						valueSearch[
							"homo_energy",
							Flatten[{templist[[First[Flatten[Position[templist,"virt."]]]-2]]}]
						]
					];(*If*)
				
					If[
						MemberQ[freqKeywords,"fmo"],
						valueSearch[
							"lumo_energy",
							Flatten[{templist[[First[Flatten[Position[templist,"virt."]]]+3]]}]
						]
					];(*If*)	

								
					If[
						MemberQ[freqKeywords,"fmo"],
						addComputedValue[
							"electronegativity",
							{-0.5*(
									Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"lumo_energy"]]]][[1]]+
									Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"homo_energy"]]]][[1]]
									)
							}
						]
					];(*If*)
																
					If[
						MemberQ[freqKeywords,"fmo"],
						addComputedValue[
							"hardness",
							{0.5*(
									Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"lumo_energy"]]]][[1]]-
									Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"homo_energy"]]]][[1]]
									)
							}
						]
					];(*If*)
					
					If[
						MemberQ[freqKeywords,"fmo"],
						addComputedValue[
							"electrophilicity",
							{0.5*(
									Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"lumo_energy"]]]][[1]]-
									Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"homo_energy"]]]][[1]]
									)
							}
						]
					];(*If*)	
																																									
(***************************************************************************************************** electronic_energy*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"E_scf",
							Flatten[StringCases[freqSection,Shortest["SCF Done:  E"~~__~~"="~~Whitespace...~~p1:NumberString~~__~~"cycles"]->p1,Overlaps->True]]
						]
					];(*If*)	
														
(***************************************************************************************************** zero_point_correction*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"zero_point_correction",
							Flatten[StringCases[freqSection,"Zero-point correction="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)	
										
(***************************************************************************************************** E_thermal_correction*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"E_thermal_correction",
							Flatten[StringCases[freqSection,"Thermal correction to Energy="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)	
										
(***************************************************************************************************** H_thermal_correction*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"H_thermal_correction",
							Flatten[StringCases[freqSection,"Thermal correction to Enthalpy="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)
					
(***************************************************************************************************** G_thermal_correction*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"G_thermal_correction",
							Flatten[StringCases[freqSection,"Thermal correction to Gibbs Free Energy="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)
(***************************************************************************************************** E_zpe*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"E_zpe",
							Flatten[StringCases[freqSection,"Sum of electronic and zero-point Energies="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)	
										
(***************************************************************************************************** E*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"E",
							Flatten[StringCases[freqSection,"Sum of electronic and thermal Energies="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)	
										
(***************************************************************************************************** H*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"H",
							Flatten[StringCases[freqSection,"Sum of electronic and thermal Enthalpies="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)
					
(***************************************************************************************************** G*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"G",
							Flatten[StringCases[freqSection,"Sum of electronic and thermal Free Energies="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)

(***************************************************************************************************** G*)						
					If[
						MemberQ[freqKeywords,"thermochem"],
						valueSearch[
							"G",
							Flatten[StringCases[freqSection,"Sum of electronic and thermal Free Energies="~~WhitespaceCharacter..~~p1:NumberString->p1]]
						]
					];(*If*)	
																	
(***************************************************************************************************** Mulliken charges and spin densities*)
					If[
						MemberQ[freqKeywords,"Mulliken_population"],
						Which[
							StringContainsQ[freqSection,"Mulliken charges and spin densities:"],
							(
								templist=Flatten[
									StringSplit[
										Flatten[
											StringCases[freqSection,Shortest["Mulliken charges and spin densities:"~~p1__~~"Sum of Mulliken charges"]->p1]
										](*Flatten*)
									](*StringSplit*)
								];(*Flatten*)
							
								index=Flatten[Position[Map[LetterQ,templist],True]];
							
								output={};
								For[
									i=1,
									i<Count[index,x_]+1,
									i++,
									output=Append[output,templist[[index[[i]]+{-1,0,1,2}]]]
								];(*For*)
							
								For[
									i=1,
									i<5,
									i++,
									atomDescriptorMatrix=If[
										{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
										DeleteDuplicates[
											Append[
												atomDescriptorMatrix,
												Transpose[Prepend[output,{"atom_number","atom","Mulliken_charge","Mulliken_spin_density"}]][[i]]
											](*Append*)
										](*DeleteDuplicates*),
										atomDescriptorMatrix
									](*If*)
								];(*For*)
							),
						
							StringContainsQ[freqSection,"Mulliken charges:"],
							(
								templist=Flatten[
									StringSplit[
										Flatten[
											StringCases[freqSection,Shortest["Mulliken charges:"~~p1__~~"Sum of Mulliken charges"]->p1]
										](*Flatten*)
									](*StringSplit*)
								];(*Flatten*)
							
								index=Flatten[Position[Map[LetterQ,templist],True]];
							
								output={};
								For[
									i=1,
									i<Count[index,x_]+1,
									i++,
									output=Append[output,templist[[index[[i]]+{-1,0,1}]]]
								];(*For*)
							
								For[
									i=1,
									i<4,
									i++,
									atomDescriptorMatrix=If[
										{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
										DeleteDuplicates[
											Append[
												atomDescriptorMatrix,
												Transpose[Prepend[output,{"atom_number","atom","Mulliken_charge"}]][[i]]
											](*Append*)
										](*DeleteDuplicates*),
										atomDescriptorMatrix
									](*If*)
								];(*For*)
							),
							
							StringContainsQ[freqSection,"Mulliken atomic charges:"],
							(
								templist=Flatten[
									StringSplit[
										Flatten[
											StringCases[freqSection,Shortest["Mulliken atomic charges:"~~p1__~~"Sum of Mulliken atomic charges"]->p1]
										](*Flatten*)
									](*StringSplit*)
								];(*Flatten*)
							
								index=Flatten[Position[Map[LetterQ,templist],True]];
							
								output={};
								For[
									i=1,
									i<Count[index,x_]+1,
									i++,
									output=Append[output,templist[[index[[i]]+{-1,0,1}]]]
								];(*For*)
							
								For[
									i=1,
									i<4,
									i++,
									atomDescriptorMatrix=If[
										{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
										DeleteDuplicates[
											Append[
												atomDescriptorMatrix,
												Transpose[Prepend[output,{"atom_number","atom","Mulliken_charge"}]][[i]]
											](*Append*)
										](*DeleteDuplicates*),
										atomDescriptorMatrix
									](*If*)
								];(*For*)
							)
							
							
						];(*Which*)
					];(*If*)

(***************************************************************************************************** APT charges*)	
					If[
						MemberQ[freqKeywords,"APT"],
						(
							templist=Flatten[
								StringSplit[
									Flatten[
										StringCases[freqSection,Shortest["APT charges:"~~p1__~~"Sum of APT charges"]->p1]
									](*Flatten*)
								](*StringSplit*)
							];(*Flatten*)
							
							index=Flatten[Position[Map[LetterQ,templist],True]];
							
							output={};
							For[
								i=1,
								i<Count[index,x_]+1,
								i++,
								output=Append[output,templist[[index[[i]]+{-1,0,1}]]]
							];(*For*)
							
							For[
								i=1,
								i<4,
								i++,
								atomDescriptorMatrix=If[
									{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
									DeleteDuplicates[
										Append[
											atomDescriptorMatrix,
											Transpose[Prepend[output,{"atom_number","atom","APT_charge"}]][[i]]
										](*Append*)
									](*DeleteDuplicates*),
									atomDescriptorMatrix
								](*If*)
							];(*For*)
						)
					];(*If*)
					
(***************************************************************************************************** Natural population analysis*)
					If[
						MemberQ[freqKeywords,"NPA"],
						(
							templist=Flatten[
								StringSplit[
									Flatten[
										StringCases[freqSection,Shortest["Summary of Natural Population Analysis:"~~p1__~~"* Total *"]->p1]
									][[1]](*Flatten*)
								](*StringSplit*)
							];(*Flatten*)
							
							If[
								MemberQ[templist,"Total"],
								templist=templist[[Flatten[Position[templist,"Total"]][[1]]+2;;-2]]
							];(*If*)
							
							index=Flatten[Position[Map[LetterQ,templist],True]];
							
							output={};
							For[
								i=1,
								i<Count[index,x_]+1,
								i++,
								output=Append[output,templist[[index[[i]]+{0,1,2,3,4,5,6}]]]
							];(*For*)
							
							For[
								i=1,
								i<8,
								i++,
								atomDescriptorMatrix=If[
									{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
									DeleteDuplicates[
										Append[
											atomDescriptorMatrix,
											Transpose[Prepend[output,{"atom","atom_number","NPA_charge","NPA_core","NPA_valence","NPA_Rydberg","NPA_total"}]][[i]]
										](*Append*)
									](*DeleteDuplicates*),
									atomDescriptorMatrix
								](*If*)
							];(*For*)
						)
					];(*If*)
										
(***************************************************************************************************** NMR*)
					If[
						MemberQ[freqKeywords,"nmr"],
						(
							templist=Flatten[
								StringCases[
									freqSection,
									Longest["GIAO Magnetic shielding tensor (ppm):"~~p1__~~"Eigenvalues:"~~Whitespace...~~NumberString..~~Whitespace...~~NumberString...]->p1
								](*StringCases*)
							][[1]];
							
							output=StringCases[
								templist,
								p1:NumberString~~Whitespace...~~
								p2:LetterCharacter..~~Whitespace...~~"Isotropic ="~~Whitespace...~~
								p3:NumberString~~Whitespace...~~"Anisotropy ="~~Whitespace...~~
								p4:NumberString
								->{p1,p2,p3,p4}
								];
							
							For[
								i=1,
								i<5,
								i++,
								atomDescriptorMatrix=If[
									{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
									DeleteDuplicates[
										Append[
											atomDescriptorMatrix,
											Transpose[Prepend[output,{"atom_number","atom","NMR_shift","NMR_anisotropy"}]][[i]]
										](*Append*)
									](*DeleteDuplicates*),
									atomDescriptorMatrix
								](*If*)
							];(*For*)
						)
					];(*If*)
					
(***************************************************************************************************** *)
(***************************************************************************************************** td output*)
(***************************************************************************************************** *)	
			If[
				MemberQ[outputSections,"td"],
				
				(		
					tdSection=StringCases[
						logFileText,
						Shortest[
							rootSections[[
								Flatten[
									Position[
										StringContainsQ[rootSections,"td",IgnoreCase->True],
									True]
								][[1]]
							]]~~
							x__~~"Normal termination"
						],
						Overlaps->True						
					][[1]];(*StringCases*)
			
					possibleTDKeywords={
						{"standard_tddft","Excitation energies and oscillator strengths:"},
						{"ES_molar_volume","Molar volume ="},
						{"ES_electronic_spatial_extent","Electronic spatial extent"},
						{"ES_dipole","Dipole moment (field-independent basis, Debye)"}
					};
					
					If[
						MemberQ[
							ExtractDescriptorsPresets,
							"atom_specific"
						],
						possibleTDKeywords=Join[
							possibleTDKeywords,
							{
								{"ES_Mulliken_population","Mulliken charges"},
								{"ES_APT","APT charges:"},
								{"ES_NPA","Summary of Natural Population Analysis:"}
							}
						]
					];(*If*)
				
					
					tdKeywords={};
					For[
						i=1,
						i<Count[possibleTDKeywords,x_]+1,
						i++,
						(
							If[
								StringContainsQ[tdSection,possibleTDKeywords[[i]][[2]]],
								tdKeywords=Append[tdKeywords,possibleTDKeywords[[i]][[1]]]
							](*If*)
						)
					];(*For*)

(***************************************************************************************************** Root ES: dipole moment*)
					If[
						MemberQ[tdKeywords,"ES_dipole"],
						valueSearch[
							"ES_root_dipole",
							Flatten[StringCases[tdSection,Shortest["Dipole moment (field-independent basis, Debye):"~~__~~"Tot="~~Whitespace...~~p1:NumberString]->p1,Overlaps->True]]
						]
					];(*If*)
														
(***************************************************************************************************** Root ES: molar volume*)
					If[
						MemberQ[tdKeywords,"ES_molar_volume"],
						valueSearch[
							"ES_root_molar_volume",
							Flatten[StringCases[tdSection,"Molar volume ="~~WhitespaceCharacter..~~x:NumberString->x]]
						]
					];(*If*)
														
(***************************************************************************************************** Root ES: electronic spatial extent*)
					If[
						MemberQ[tdKeywords,"ES_electronic_spatial_extent"],
						valueSearch[
							"ES_root_electronic_spatial_extent",
							Flatten[StringCases[tdSection,"Electronic spatial extent (au):  <R**2>="~~Whitespace...~~x:NumberString->x]]
						]
					];(*If*)

(***************************************************************************************************** Standard TDDFT output*)
					With[
						{
							EStemp=StringCases[
								tdSection,
								Shortest[
									"Excited State"~~__~~p1:NumberString..~~Whitespace...~~"nm"~~Whitespace...~~
									"f="~~Whitespace...~~p2:NumberString..~~Whitespace...~~
									"<S**2>="~~Whitespace...~~p3:NumberString..~~Whitespace..
								]->{p1,p2,p3},
								Overlaps->False
							]
						},
						
						If[
							Count[EStemp,_]>0,
							(

								Table[
									(
									 addComputedValue["ES"<>ToString[n]<>"_transition",{EStemp[[n]][[1]]}];
									 addComputedValue["ES"<>ToString[n]<>"_osc_strength",{EStemp[[n]][[2]]}];
									 addComputedValue["ES"<>ToString[n]<>"_<S**2>",{EStemp[[n]][[3]]}];							 
									 ),
									{n,1,Count[EStemp,_],1}
								];

							)
						];(*If*)
					];(*With*)

		
(***************************************************************************************************** Root ES: Mulliken charges*)			
					If[
						MemberQ[tdKeywords,"ES_Mulliken_population"],
						Which[
							StringContainsQ[tdSection,"Mulliken charges and spin densities:"],
							(
								templist=Flatten[
									StringSplit[
										Flatten[
											StringCases[tdSection,Shortest["Mulliken charges and spin densities:"~~p1__~~"Sum of Mulliken charges"]->p1]
										](*Flatten*)
									](*StringSplit*)
								];(*Flatten*)
							
								index=Flatten[Position[Map[LetterQ,templist],True]];
							
								output={};
								For[
									i=1,
									i<Count[index,x_]+1,
									i++,
									output=Append[output,templist[[index[[i]]+{-1,0,1,2}]]]
								];(*For*)
							
								For[
									i=1,
									i<5,
									i++,
									atomDescriptorMatrix=If[
										{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
										DeleteDuplicates[
											Append[
												atomDescriptorMatrix,
												Transpose[Prepend[output,{"atom_number","atom","ES_root_Mulliken_charge","ES_root_Mulliken_spin_density"}]][[i]]
											](*Append*)
										](*DeleteDuplicates*),
										atomDescriptorMatrix
									](*If*)
								];(*For*)
							),
						
							StringContainsQ[tdSection,"Mulliken charges:"],
							(
								templist=Flatten[
									StringSplit[
										Flatten[
											StringCases[tdSection,Shortest["Mulliken charges:"~~p1__~~"Sum of Mulliken charges"]->p1]
										](*Flatten*)
									](*StringSplit*)
								];(*Flatten*)
							
								index=Flatten[Position[Map[LetterQ,templist],True]];
							
								output={};
								For[
									i=1,
									i<Count[index,x_]+1,
									i++,
									output=Append[output,templist[[index[[i]]+{-1,0,1}]]]
								];(*For*)
							
								For[
									i=1,
									i<4,
									i++,
									atomDescriptorMatrix=If[
										{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
										DeleteDuplicates[
											Append[
												atomDescriptorMatrix,
												Transpose[Prepend[output,{"atom_number","atom","ES_root_Mulliken_charge"}]][[i]]
											](*Append*)
										](*DeleteDuplicates*),
										atomDescriptorMatrix
									](*If*)
								];(*For*)
							),
							
							StringContainsQ[tdSection,"Mulliken atomic charges:"],
							(
								templist=Flatten[
									StringSplit[
										Flatten[
											StringCases[tdSection,Shortest["Mulliken atomic charges:"~~p1__~~"Sum of Mulliken atomic charges"]->p1]
										](*Flatten*)
									](*StringSplit*)
								];(*Flatten*)
							
								index=Flatten[Position[Map[LetterQ,templist],True]];
							
								output={};
								For[
									i=1,
									i<Count[index,x_]+1,
									i++,
									output=Append[output,templist[[index[[i]]+{-1,0,1}]]]
								];(*For*)
							
								For[
									i=1,
									i<4,
									i++,
									atomDescriptorMatrix=If[
										{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
										DeleteDuplicates[
											Append[
												atomDescriptorMatrix,
												Transpose[Prepend[output,{"atom_number","atom","ES_root_Mulliken_charge"}]][[i]]
											](*Append*)
										](*DeleteDuplicates*),
										atomDescriptorMatrix
									](*If*)
								];(*For*)
							)
							
							
						];(*Which*)
					];(*If*)
														
(***************************************************************************************************** Root ES: APT charges*)			
					If[
						MemberQ[tdKeywords,"ES_APT"],
						(
							templist=Flatten[
								StringSplit[
									Flatten[
										StringCases[tdSection,Shortest["APT charges:"~~p1__~~"Sum of APT charges"]->p1]
									](*Flatten*)
								](*StringSplit*)
							];(*Flatten*)
							
							index=Flatten[Position[Map[LetterQ,templist],True]];
							
							output={};
							For[
								i=1,
								i<Count[index,x_]+1,
								i++,
								output=Append[output,templist[[index[[i]]+{-1,0,1}]]]
							];(*For*)
							
							For[
								i=1,
								i<4,
								i++,
								atomDescriptorMatrix=If[
									{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
									DeleteDuplicates[
										Append[
											atomDescriptorMatrix,
											Transpose[Prepend[output,{"atom_number","atom","ES_root_APT_charge"}]][[i]]
										](*Append*)
									](*DeleteDuplicates*),
									atomDescriptorMatrix
								](*If*)
							];(*For*)
						)
					];(*If*)
														
(***************************************************************************************************** Root ES: NPA charge, core, valence, Rydberg, total*)											
					If[
						MemberQ[tdKeywords,"ES_NPA"],
						(
							templist=Flatten[
								StringSplit[
									Flatten[
										StringCases[tdSection,Shortest["Summary of Natural Population Analysis:"~~p1__~~"* Total *"]->p1]
									][[1]](*Flatten*)
								](*StringSplit*)
							];(*Flatten*)
							
							If[
								MemberQ[templist,"Total"],
								templist=templist[[Flatten[Position[templist,"Total"]][[1]]+2;;-2]]
							];(*If*)
							
							index=Flatten[Position[Map[LetterQ,templist],True]];
							
							output={};
							For[
								i=1,
								i<Count[index,x_]+1,
								i++,
								output=Append[output,templist[[index[[i]]+{0,1,2,3,4,5,6}]]]
							];(*For*)
							
							For[
								i=1,
								i<8,
								i++,
								atomDescriptorMatrix=If[
									{Count[output,x_]}==Transpose[descriptorVector][[2]][[Flatten[Position[Transpose[descriptorVector][[1]],"number_of_atoms"]]]],
									DeleteDuplicates[
										Append[
											atomDescriptorMatrix,
											Transpose[Prepend[output,{"atom","atom_number","ES_root_NPA_charge","ES_root_NPA_core","ES_root_NPA_valence","ES_root_NPA_Rydberg","ES_root_NPA_total"}]][[i]]
										](*Append*)
									](*DeleteDuplicates*),
									atomDescriptorMatrix
								](*If*)
							];(*For*)
						)
					];(*If*)
																
(***************************************************************************************************** Root ES: NPA core*)						
	
				)													
			];(*If*)
			
(***************************************************************************************************** *)			
(***************************************************************************************************** Final output sections*)	
(***************************************************************************************************** *)

(***************************************************************************************************** Data for labeled atoms*)					
				If[
					MemberQ[
						ExtractDescriptorsPresets,
						"atom_specific"
						],
					If[
						AnyTrue[StringContainsQ[Transpose[coordinateList][[1]],NumberString],TrueQ],
						(
							labelIndex=Sort[
								With[
									{
										list={
											Range[Count[Transpose[coordinateList][[1]],_]],
											StringReplace[Transpose[coordinateList][[1]],__~~p1:NumberString->p1]
										},
										labels=Map[ToString,getLabels]
									},
									If[
										getLabels=={"auto"},
										ToExpression[Select[Map[Flatten,Transpose[{list[[1]],StringCases[NumberString..][list[[2]]]}]],Count[#,_]==2&]],
										ToExpression[Select[Map[Flatten,Table[If[MemberQ[labels,list[[2]][[n]]],Transpose[list][[n]],{}],{n,1,Count[list[[2]],_],1}]],Count[#,_]==2&]]
									]
								],
								#1[[2]]<#2[[2]]&
							];(*Sort*)
							
							atomDescriptorMatrixLabeled=Transpose[ReplacePart[atomDescriptorMatrix,2->Join[{"atom"},Transpose[coordinateList][[1]]]]][[Flatten[{1,Transpose[labelIndex][[1]]+1}]]];
						),
						atomDescriptorMatrixLabeled={}
					];(*If*)
				](*If*)
																																															
				)
			];(*If*)

(***************************************************************************************************** Final tables*)																																																																												
		descriptorVector=descriptorVector;
		
		atomDescriptorMatrix=If[
			Count[atomDescriptorMatrix,_]>0,
			Transpose[atomDescriptorMatrix]
		];
			
		atomDescriptorMatrixLabeled=If[
			Count[atomDescriptorMatrixLabeled,_]>0,
			atomDescriptorMatrixLabeled
		];

		),
		Throw[Print["Error: file not found."]]
	](*If*)
](*Catch*)
