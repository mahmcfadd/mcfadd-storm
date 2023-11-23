
path=getDirectory("Choose STORM images directory");
close("*");

ROIcreate=newArray("Yes" , "No");

Dialog.create("");
filearray=getFileList(path);
Dialog.addString("Image name:", "ST001");
Dialog.addChoice("Select Channel 405:", filearray);
Dialog.addChoice("Select Channel Cy3:", filearray);
Dialog.addCheckbox("2-color?", false);
Dialog.addCheckbox("Create STORM ROIs", false);
Dialog.addCheckbox("Confirm ROIs", false);
Dialog.show();

In=Dialog.getString();
Ch1name="filepath=[" + path + Dialog.getChoice() + "] fileformat=[CSV (comma separated)] append=false startingframe=1";
Ch2name="filepath=[" + path + Dialog.getChoice() + "] fileformat=[CSV (comma separated)] append=false startingframe=1";
pa2th=path + In +"ROIs" + File.separator;
Ijcolor=Dialog.getCheckbox();
ROIcreate=Dialog.getCheckbox();
ROIname="ROI";
ROIconfirm=Dialog.getCheckbox();

ROIfilename=In + ROIname + "set";
ROIpath=pa2th + ROIfilename + ".zip";

if (ROIcreate==1) {

pathex=File.exists(pa2th);
if(pathex==1) {
exit(In + "ROIs already exist");
};
File.makeDirectory(pa2th);

run("Import results", Ch1name);
Dialog.create("STORM ROI creator");
Dialog.addMessage("Drift OK?");
Dialog.show();
close("Averaged*");
run("Visualization", "imleft=0.0 imtop=0.0 imwidth=256.0 imheight=256.0 renderer=[Averaged shifted histograms] magnification=16.0 colorizez=false shifts=3 threed=false");
rename("405");

if(Ijcolor==1){
	run("Import results", Ch2name);
	close("Averaged*");
	run("Visualization", "imleft=0.0 imtop=0.0 imwidth=256.0 imheight=256.0 renderer=[Averaged shifted histograms] magnification=16.0 colorizez=false shifts=3 threed=false");
	rename("Cy3");
	run("Merge Channels...", "c1=Cy3 c2=405 create");
};
run("ROI Manager...");
roiManager("Show All");
makeRectangle(1700, 1700, 200,200);
waitForUser("Select ROI set and press OK");

//CB1R ROIrename and export macro________________________________________________________________

i=roiManager("index");
c=roiManager("count");
ri=roiManager("count");

roiManager("Select All");
ROIfilename=In + ROIname + "set";
ROIpath=pa2th + ROIfilename + ".zip";
roiManager("Save", ROIpath);

roiManager("Select", 0);

//Cy3 ROI export__________________________________________________________________________________
if(Ijcolor==1){
for (i=0; i+1<=c; i++) {
	roiManager("Select", i);
	getSelectionBounds(x, y, width, height);
	roiManager("Rename", ROIname + i+1);

	roi="action=filter formula=[(x >" +  x*10 + "& x <" + x*10+2000 + "& y >" + y*10 + "& y <" + y*10+2000 + ")]";
	
	
	run("Show results table", roi);

	filename=In + ROIname + i+1 + "Cy3"+"_TS3D";
	ChROIname="filepath=[" + pa2th + filename + ".csv] fileformat=[CSV (comma separated)] sigma1=true sigma2=true intensity=true uncertainty_xy=true saveprotocol=false id=false frame=true bkgstd=true uncertainty_z=true offset=true detections=true z=true y=true x=true";
	
	run("Export results", ChROIname);
			};

close("Averaged*");

	
//405 ROI export__________________________________________________________________________________

run("Import results", Ch1name);
};

for (i=0; i+1<=c; i++) {
	
	roiManager("Select", i);
	getSelectionBounds(x, y, width, height);	
	roi="action=filter formula=[(x >" +  x*10 + "& x <" + x*10+2000 + "& y >" + y*10 + "& y <" + y*10+2000 + ")]";

	run("Show results table", roi);	
	file405name=In + ROIname + i+1 + "405"+"_TS3D";
	ChROIname="filepath=[" + pa2th + file405name + ".csv] fileformat=[CSV (comma separated)] sigma1=true sigma2=true intensity=true uncertainty_xy=true saveprotocol=false id=false frame=true bkgstd=true uncertainty_z=true offset=true detections=true z=true y=true x=true";
	run("Export results", ChROIname);	
	};

close("Averaged*");
close("*");

};

if (ROIconfirm==1) {

//ROI confirmation ______________________________________________________________________________

ri=roiManager("count");
i=roiManager("Index");
c=roiManager("count");

ROIarray=newArray("Yes" , "No");

k=1;

for (i=0; i+1<=c; i++) {
	Cy3filename=In + ROIname + i+1 + "Cy3"+"_TS3D";
	Cy3ROIname="filepath=[" + pa2th + Cy3filename + ".csv] fileformat=[CSV (comma separated)] sigma1=true sigma2=true intensity=true uncertainty_xy=true saveprotocol=false id=true frame=true bkgstd=true uncertainty_z=true offset=true detections=true z=true y=true x=true";
	file405name=In + ROIname + i+1 + "405"+"_TS3D";
	ChROIname="filepath=[" + pa2th + file405name + ".csv] fileformat=[CSV (comma separated)] sigma1=true sigma2=true intensity=true uncertainty_xy=true saveprotocol=false id=true frame=true bkgstd=true uncertainty_z=true offset=true detections=true z=true y=true x=true";
		

	run("Import results", Cy3ROIname);
	run("Visualization", "imleft=0.0 imtop=0.0 imwidth=256.0 imheight=256.0 renderer=[Scatter plot] magnification=16.0 colorizez=false threed=false");
	roiManager("Select", i);
	getSelectionBounds(x, y, width, height);

	open(pa2th+Cy3filename+".csv");
	wait(1000);
	updateResults();
	znoC=nResults;
	Cyz=newArray(znoC);
	Cyx=newArray(znoC);
	Cyy=newArray(znoC);
	for(d=0; d<znoC; d++) {
		Cyz[d]=getResult("z [nm]", d);
		Cyx[d]=getResult("x [nm]", d);
		Cyy[d]=getResult("y [nm]", d);
		};
	selectWindow("Results");
	run("Close");
	open(pa2th+file405name+".csv");
	wait(1000);
	znoF=nResults;
	Fourz=newArray(znoF);
	Fourx=newArray(znoF);
	Foury=newArray(znoF);
	for(d=0; d<znoF; d++) {
		Fourz[d]=getResult("z [nm]", d);
		Fourx[d]=getResult("x [nm]", d);
		Foury[d]=getResult("y [nm]", d);
		};
	selectWindow("Results");
	run("Close");
	Array.getStatistics(Cyz, Cyzmin, Cyzmax, Cyzmean, Cyzstdev);
	Array.getStatistics(Cyx, Cyxmin, Cyxmax, Cyxmean, Cyxstdev);
	Array.getStatistics(Cyy, Cyymin, Cyymax, Cyymean, Cyystdev);
	Array.getStatistics(Fourz, Fourzmin, Fourzmax, Fourzmean, Fourzstdev);
	Array.getStatistics(Fourx, Fourxmin, Fourxmax, Fourxmean, Fourxstdev);
	Array.getStatistics(Foury, Fourymin, Fourymax, Fourymean, Fourystdev);
	
	zmin=Fourzmin;
	zmax=Fourzmax;
	if(Cyzmin<zmin) {
		zmin=Cyzmin;
		};
	if(Cyzmax>zmax) {
		zmax=Cyzmax;
		};
	
	bounds="imleft=" + x/16 + " imtop=" + y/16 + " imwidth=12 imheight=12 renderer=[Scatter plot] magnification=16.0 colorize=false threed=true zrange=" + zmin + ":10.0:" + zmax;
	//print(bounds);
	close("Averaged*");
	run("Visualization", bounds);
	rename("Hmr");
	setMinAndMax(0, 1000);
	
	
	run("Import results", ChROIname);
	close("Averaged*");
	run("Visualization", bounds);
	rename("Bsn");

	run("Merge Channels...", "c1=Bsn c2=Hmr create");
	selectWindow("Composite");
	setMinAndMax(0, 1000);
	run("3D Project...", "projection=[Brightest Point] axis=Y-Axis slice=0.01 initial=0 total=360 rotation=10 lower=1 upper=255 opacity=0 surface=100 interior=0");
	selectWindow("Projections of Composite");
	rename(In + ROIname + i+1);

	Dialog.create("ROI state");
	Dialog.addCheckbox("ROI Ok? (check if yes)", true);
	waitForUser("When ready press OK");
	Dialog.show();
	
	ROIstate=Dialog.getCheckbox();
	
	if (ROIstate==0) {
		

		pa405th=pa2th + file405name + ".csv";
		paCy3th=pa2th + Cy3filename + ".csv";
			

		File.delete(pa405th);
		File.delete(paCy3th);
		
		c=roiManager("count");
		if (i+1<c) {
			for (f=k; f+1<=ri; f++) {
				pane405th=pa2th + In + ROIname + f+1 + "405_TS3D" + ".csv";
				paneCy3th=pa2th + In + ROIname + f+1 + "Cy3_TS3D" + ".csv";
				pa405th=pa2th + In + ROIname + f + "405_TS3D" + ".csv";
				paCy3th=pa2th + In + ROIname + f + "Cy3_TS3D" + ".csv";
				File.rename(pane405th, pa405th);
				File.rename(paneCy3th, paCy3th);
						};
			
			roiManager("Delete");
			c=roiManager("count");
			i=i-1;
			close("*");
			};
		if (i+1==c) {
			roiManager("Delete");
			c=roiManager("count");
				};
		}  else {
		roiManager("Rename", ROIname + i+1);		
		close("*");
			};
	
	};


c=roiManager("count");


roiManager("Select All");

roiManager("Save", ROIpath);

selectWindow("Angles");
saveAs("Results", path + In + "Angles.xls");
selectWindow(In + "Angles.xls");
//run("Close");

};
Dialog.create("Job Done!");
Dialog.addMessage("Thanks for using Maureen's macro!");
Dialog.show();

