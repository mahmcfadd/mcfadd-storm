dirimp=getDirectory("Choose STORM stack Directory");

ALM=newArray("8-neighbourhood", "4-neighbourhood");
SLMmethod=newArray("PSF: Integrated Gaussian", "PSF: Gaussian", "PSF: Elliptical Gaussian (3D astigmatism)");
SLMfitmeth=newArray("Least squares", "Weighted Least squares", "Maximum likelihood");
Expform=newArray("CSV (comma separated)", "XLS (tab separated)");
Explabels=newArray("id", "frame", "x", "y", "z", "sigma", "intensity", "offset", "bkgstd", "uncertainty_xy", "uncertainty_z");
Expdef=newArray(true, true, true, true, true, true, true, true, true, true, true);

Dialog.create("TSBatchAnalysis");
Dialog.setInsets(0, 200, 0);
Dialog.addMessage("Analysis Parameters");

	Dialog.setInsets(10, 10, 0);
	Dialog.addMessage("Image filtering");
		Dialog.setInsets(0, 60, 0);
		Dialog.addMessage("Filter: Wavelet filter (B-Spline) \n");
		Dialog.setInsets(0, 60, 0);
		Dialog.addNumber("B-Spline order:", 3);
		Dialog.setInsets(0, 60, 0);
		Dialog.addNumber("B-Spline scale:", 2.0, 1, 5, "pixels");
		
	Dialog.setInsets(10, 10, 0);
	Dialog.addMessage("Approximate localization of molecules");
		Dialog.setInsets(0, 60, 0);
		Dialog.addMessage("Method: Local Maximum");
		Dialog.setInsets(0, 60, 0);
		Dialog.addString("Peak intensity threshold:", "std(Wave.F1)");
		Dialog.setInsets(0, 60, 0);
		Dialog.addRadioButtonGroup("Connectivity:", ALM, 2, 1, ALM[0]);
		
	Dialog.setInsets(10, 10, 0);
	Dialog.addMessage("\n \n Sub-pixel localization of molecules");
		Dialog.setInsets(0, 60, 0);
		Dialog.addChoice("Method:", SLMmethod, SLMmethod[2]);
		Dialog.setInsets(0, 60, 0);
		Dialog.addNumber("Fitting radius [px]:", 3);
		Dialog.setInsets(0, 60, 0);
		Dialog.addChoice("Fitting method:", SLMfitmeth, SLMfitmeth[2]);
		Dialog.setInsets(0, 60, 0);
		Dialog.addNumber("Initial sigma [px]:", 1.6);
		Dialog.setInsets(0, 60, 0);
		Dialog.addCheckbox("Select calibration file?", false);
		
	Dialog.setInsets(10, 10, 0);
	Dialog.addMessage("\n \n Reconstruction Export");
		Dialog.setInsets(0, 60, 0);
		Dialog.addString("Export folder name:", "Name");
		Dialog.setInsets(0, 30, 0);
		Dialog.addChoice("File format:", Expform, Expform[0]);
		Dialog.setInsets(0, 60, 0);
		Dialog.addCheckbox("Save measurement protocol", true);
		Dialog.setInsets(10, 60, 0);
		Dialog.addCheckboxGroup(13, 1, Explabels, Expdef);
Dialog.show();

Splineorder=Dialog.getNumber;
Splinescale=Dialog.getNumber;

ALMmeth="Local maximum";
ALMint=Dialog.getString;
ALMconn=Dialog.getRadioButton;

SLMmeth=Dialog.getChoice;
SLMfitr=Dialog.getNumber;
SLMfitmeth=Dialog.getChoice;
SLMinisig=Dialog.getNumber;
cali=Dialog.getCheckbox;

dirname=Dialog.getString;
expform=Dialog.getChoice;
proto=Dialog.getCheckbox;
expcheck=newArray(11);
for(i=0; i<11; i++){
	expcheck[i]=Dialog.getCheckbox;
};

for(i=0; i<11; i++) {
	if(expcheck[i]==0) {
	expcheck[i]="false";
	} else { 
	expcheck[i]="true";
	};
};
if(proto==1){
	proto="true";
} else { 
	proto="false";
};

if(cali==true) {
	cal=File.openDialog("Select calibration file");
	} else {
	cal="C:\\Users\\Maureen\\Documents\\ImageJMomo\\PSF STORM\\100xCalThunderstormNew";
	};

dir=File.getParent(dirimp); 
direxp=dir+File.separator+dirname;

File.makeDirectory(direxp);

STORMstks=getFileList(dirimp);
l=lengthOf(STORMstks);
nij=nImages;
for(i=0; i<l; i++) {
	stkname=dirimp+File.separator+STORMstks[i];
	//run("Bio-Formats Importer", "open=stkname color_mode=Default view=Hyperstack stack_order=XYCZT use_virtual_stack");
	run("Stack From List...", "open=stkname use");
	while(nij==0) {
	wait(100);
	nij=nImages;
	};
	analy="filter=[Wavelet filter (B-Spline)] scale=" + Splinescale + " order=" + Splineorder + " detector=[" + ALMmeth + "] connectivity=" + ALMconn;
	sis=" threshold=" + ALMint + " estimator=[" + SLMmeth + "] sigma=" + SLMinisig + " fitradius=" + SLMfitr + " method=[" + SLMfitmeth + "] calibrationpath=[" + cal + "] full_image_fitting=false mfaenabled=false renderer=[No Renderer]";
	analysis=analy + sis;
	//run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[ALMmeth] connectivity=ALMconn threshold=ALMint estimator=[SLMmeth] sigma=SLMinisig fitradius=SLMfitr method=[SLMfirmeth] calibrationpath=[cal] full_image_fitting=false mfaenabled=false renderer=[No Renderer]");
	run("Run analysis", analysis);
	if(expform=="CSV (comma separated)") {
		fileexp=direxp+File.separator+"Rec"+STORMstks[i]+".csv";
		//expform="[CSV (comma separated)]";
	} else {
		fileexp=direxp+File.separator+"Rec"+STORMstks[i]+".xls";
		//expform="[XLS (tab separated)]";
	};
	//export="filepath=[" + fileexp + "] fileformat=[CSV (comma separated)]";
	export="floatprecision=5 filepath=[" + fileexp + "] fileformat=[" + expform + "] sigma=" + expcheck[5] + " intensity=" + expcheck[6] + " offset=" + expcheck[7] + " saveprotocol=" + proto + " x=" + expcheck[2] + " y=" + expcheck[3] + " z=" + expcheck[4] + " bkgstd=" + expcheck[8] + " id=" + expcheck[0] + " uncertainty_xy=" + expcheck[9] + " frame=" + expcheck[1] + " uncertainty_z=" + expcheck[10];
	//run("Export results", "floatprecision=5 filepath=fileexp fileformat=[CSV (comma separated)] sigma=expcheck[5] intensity=expcheck[7] offset=expcheck[8] saveprotocol=proto x=expcheck[2] y=expcheck[3] bkgstd=expcheck[9] id=expcheck[0] uncertainty_xy=expcheck[9] frame=expcheck[1]");
	run("Export results", export);
	nij=nImages;
	if(nImages>0){
	close();
	};
};

waitForUser("Thank you for using Maureen's macro ;P");




