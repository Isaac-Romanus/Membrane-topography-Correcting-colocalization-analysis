//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											--		Intro		--													//
//						This macro aims to simulate two images with varying intensities but							//
//						with similar membrane distrubution. A correlation calculation will							//
//						be done between the generated images before and after multiplicating						//
//						with the membrane-distribution-image. Printing a result table with							//
//						the differens.																				//

//SETTINGS------------------------------------------------------------------------------------------------------------
sz = 500; 		// Image size
prnt = 0;		// 1 - print during macro      				0 - no printing,			2 - Results printed in log
parts = 5		// 10 -10% intervals, standard				100 - 1% intervalls
runs = 10
batch = 1;		// 1 - hide images							0 - show images	
//------------------------------------------------------------------------------------------------------------SETTINGS

//CLEAR WORKSPACE-----------------------------------------------------------------------------------------------------
if (nImages>0) run("Close All"); 							// if there are 1 or more images - close all
if (isOpen("Results")) { 									// close results window, if open
	selectWindow("Results");
	close("Results");
}
print("\\Clear"); 											// empty log window
roiManager("reset"); 										// empty ROI manager
run("Options...", "iterations=1 count=1 black edm=32-bit"); // set Binary Options
//-----------------------------------------------------------------------------------------------------CLEAR WORKSPACE

//SETUP ROI AND RESULTS -----------------------------------------------------------------------------------------------
newImage("Base", "8-bit black", sz, sz, 1); // Creating an image with the size sz, all values 0
newImage("255", "8-bit white", sz, sz, 1); // Creating a image with size sz all values 255
run("Select All");
roiManager("Add"); // Creating a ROI the size of the selected image
ROI = roiManager("count")-1;
run("Select None");

// Creating arrays for result printing, paratanses for measurepoints
true_array = newArray(11); // The true correlation values between the images
corr_array=newArray(11);// Corelation values values
mem_correct_array = newArray(11); //Membrane corrected correlation values 
//-------------------------------------------------------------------------------------------------SETUP ROI AND RESULTS

setBatchMode(batch);

// Creating a new image, min max refering to amount of membrane
selectImage("Base");
run("RandomJ Uniform", "min=1 max=5 insertion=Additive");
rename("MembraneFraction");
// setMinAndMax(0.2, 0.8); // changeing display of brightness and contrast
getRawStatistics(nPixels, meanMemb, min, max, std, histogram);
if (prnt == 1) print(" meanMemb", meanMemb);

cell_nr = newArray(runs);

for (n = 0; n < runs; n++) {
	cell_nr[n] = n + 1;
		
	// Creating two opposite images with values between 0,2 and 0,8
		selectImage("Base");
		run("RandomJ Uniform", "min=51 max=204 insertion=Additive");
		// run("RandomJ Poisson", "mean=102 insertion=Additive");
		// run("RandomJ Gaussian", "mean=132 sigma=24 insertion=Additive");
		rename("Sub1"); // subdomain 1, fraction subdomain 1
		sub1 = getImageID();

		run("RandomJ Poisson", "mean=1 insertion=Modulatory");
		rename("Sub1Noise"); // probability of A
		sub1_n = getImageID();
		
		selectImage("255");
		imageCalculator("Subtract create 32-bit", "255","Sub1"); //Subdomain2 - remainder from Subdomain1   1.0 - Subdomain1
		rename("Sub2");// subdomain 2, fraction subdomain 2   sub1 + sub2 =1
		sub2 = getImageID();
		
		run("RandomJ Poisson", "mean=1 insertion=Modulatory");
		rename("Sub2Noise"); // probability of A
		sub2_n = getImageID();
	
	// Calculating the correlation for the created images with -1 correlation
		true_array[0] = correl(sub1, sub2, ROI, 0); //Getting the correlation value for the two opposite images and saving it in array.
		corr_array[0] = correl(sub1_n, sub2_n, ROI, 0); //Getting the correlation value for the two opposite images and saving it in array.
		if (prnt == 2) print("Nr", "0", "Pearson: ", corr_array[0]);
	
	// Creating new images, adjusting for membrane amount in pixel by multiplying with membrane-fraction-image   	
		imageCalculator("Multiply create 32-bit", "MembraneFraction","Sub1"); // fraction of each pixel occupied by Subdomain 1
		rename("SubFracSub1"); // probability of A
		run("RandomJ Poisson", "mean=1 insertion=Modulatory");
		rename("SubFracSub1Noise"); // probability of A
		sub_frac_sub1 = getImageID();
			
		imageCalculator("Multiply create 32-bit", "MembraneFraction","Sub2"); // fraction of each pixel occupied by Subdomain 2
		rename("SubFracSub2");// probability of B
		run("RandomJ Poisson", "mean=1 insertion=Modulatory");
		rename("SubFracSub2Noise"); // probability of A
		sub_frac_sub2 = getImageID();
	
	// Calculating correlation for the corrected images with -1 correlation and closing them
		mem_correct_array[0] = correl(sub_frac_sub1, sub_frac_sub2, ROI, 0);
		if (prnt == 2) print("Nr", "0", "Corrected_mem ", mem_correct_array[0]);
		
		close("Sub1Noise");
		close("Sub2Noise");
		close("SubFracSub1");
		close("SubFracSub2");
		close("SubFracSub1Noise");
		close("SubFracSub2Noise");
	 
	// Creating images with varying amounts of correlation, standard is in 10% intervals. Defined in settings section
	for (i = 1; i <= parts; i++) {
		division = sz / parts;  // part to be copied from one image to another, currently set to 10% steps
		if (prnt == 1) print("division", division);
		
		// Marking 10% of image 1 and pasting it onto image 2 to increase their correlation, done in steps from left to right
		selectImage("Sub1");
		sub1 = getImageID();
		makeRectangle(0, 0, round(division*i), sz); // Selecting slice with i*10% thickness
		roiManager("Add");
		division_ROI = roiManager("count")-1;
		run("Copy");
		
		selectImage(sub1);
		run("RandomJ Poisson", "mean=1 insertion=Modulatory");
		rename("Sub1Noise"); // probability of A
		sub1_n = getImageID();
		
		// Pasting the selected i*10% slice onto the second image
		selectImage("Sub2");
		sub2 = getImageID();
		roiManager("select", division_ROI);
		run("Paste");
		run("RandomJ Poisson", "mean=1 insertion=Modulatory");
		rename("Sub2Noise"); // probability of A
		sub2_n = getImageID();
		
		// Calculating correlatin for the generated images with varying correlation
		true_array[i] = correl(sub1, sub2, ROI, 0);		
		corr_array[i] = correl(sub1_n, sub2_n, ROI, 0);
		if (prnt == 2) print("nr", i, "Pearson: ", corr_array[i]);
		run("Select None");
		
		// Creating new images, adjusting for membrane amount in pixel by multiplying with membrane-fraction-image   	
		imageCalculator("Multiply create 32-bit", "MembraneFraction","Sub1"); // fraction of each pixel occupied by Subdomain 1
		rename("SubFracSub1"); // probability of A
		run("RandomJ Poisson", "mean=1 insertion=Modulatory");
		rename("SubFracSub1Noise"); // probability of A
		sub_frac_sub1 = getImageID();
			
		imageCalculator("Multiply create 32-bit", "MembraneFraction","Sub2"); // fraction of each pixel occupied by Subdomain 2
		rename("SubFracSub2");// probability of B
		run("RandomJ Poisson", "mean=1 insertion=Modulatory");
		rename("SubFracSub2Noise"); // probability of A
		sub_frac_sub2 = getImageID();
		
		// Calculating correlation for the corrected images and closing them
		mem_correct_array[i] = correl(sub_frac_sub1, sub_frac_sub2, ROI, 0);
		if (prnt == 2) print("Nr", i, "Corrected_mem ", mem_correct_array[i]);
		close("SubFracSub1");
		close("SubFracSub2");
		close("SubFracSub1Noise");
		close("SubFracSub2Noise");
		close("Sub1Noise");
		close("Sub2Noise");
	}
	
	//	Printing to results table
	for (i = 0; i <= parts; i++) {
		extend = nResults;
		setResult("True", extend, true_array[i]);
		setResult("Pearson", extend, corr_array[i]);
		setResult("Membrane_corrected", extend, mem_correct_array[i]);
		setResult("diff", extend, mem_correct_array[i] - corr_array[i]);
		setResult("cell", extend, cell_nr[n]);
		updateResults();
		selectWindow("Results");
	}
	close("Sub1");
	close("Sub2");
}

// saving an exelfile
	savedir = getDirectory("Save folder"); // selecting save folder
	savestr = getString("datafil", "First test");	
	run("Read and Write Excel", "dataset_label=[My new results data label]");
	run("Read and Write Excel", "no_count_column");
	run("Read and Write Excel", "file=[" + savedir + savestr + ".csv]");


function correl(img1_id, img2_id, ROI, prnt) { 
// This is a variant of the Pcorrel macro. It calculates Pearson between the two selected images
// Note that "ROI" has to refer to an exicting ROI. And prnt = 1 gives print, 0 no print
	
	run("Select None");// roiManager("select", ROI);
	getRawStatistics(n,mean,min,max,std);
    if (prnt == 1) print("stats for SRC     min ",min," max ",max,"  STD ",std);
	run("Select None");
	
	selectImage(img1_id);
	run("Select All");
	run("Duplicate...", " ");
	roiManager("select",ROI);
	getRawStatistics(n,mean,min,max,std);
	if (prnt == 1) print("stats for SRC     min ",min," max ",max,"  STD ",std);
	run("Select None");
	rename("RedPic");
	run("32-bit");

//Subtract red-mean
	RedSubmPicID=getImageID();
	roiManager("select",ROI);	
	getRawStatistics(nPix,Redmean);
	run("Subtract...", "value=Redmean");;	
	getRawStatistics(nPix,Redmean);
	
	selectImage(img2_id);
	run("Select All");
	run("Duplicate...", " ");
	rename("GrnPic");
	roiManager("select",ROI);
	getRawStatistics(n,mean,min,max,std);
	if (prnt == 1) print("stats for SRC     min ",min," max ",max,"  STD ",std);
	run("Select None");
	run("32-bit");
	
//Subtract green-mean
	GrnSubmPicID=getImageID();
	roiManager("select",ROI);
	getRawStatistics(nPix,Grnmean);
	run("Subtract...", "value=Grnmean");	
	getRawStatistics(nPix,GrnSubmean);
			
//Calculate using Pearson Eq
//Numerator
	run("Select None");
	imageCalculator("Multiply create 32-bit",RedSubmPicID,GrnSubmPicID);
	RmXGmID=getImageID();
	roiManager("select",ROI); // the specific ROI
	NumPicID=getImageID();
	rename("Numerator");
	getRawStatistics(nPix,MeanNumerator);
	numerator=nPix*MeanNumerator;

//Denominator
	selectImage(RedSubmPicID);
	roiManager("select",ROI); // the specific ROI
	getRawStatistics(nPix,Rmean);
	run("Square"); // square image
	getRawStatistics(nPix,R2Mean);
	ER2Mean=R2Mean*nPix;
	rootER2Mean=sqrt(ER2Mean);
	
	selectImage(GrnSubmPicID);
	roiManager("select",ROI); // the specific ROI
	getRawStatistics(nPix,GMean);
	run("Square");
	getRawStatistics(nPix,G2Mean);
	EG2Mean=G2Mean*nPix;
	rootEG2Mean=sqrt(EG2Mean);
	denom=(rootER2Mean*rootEG2Mean);
	
	corr = numerator/denom;

	
// tidy up - delete images used in calculation
	selectImage(NumPicID);
	close();
	selectImage(GrnSubmPicID);
	close();
	selectImage(RedSubmPicID);
	close();
//print('end correlation');
	return corr
} // end of function correl()	





	
