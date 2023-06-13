//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											--		Intro		--													//
//						This macro aims to simulate two images with pixels having varying							//
//						procentages of two different subdomains. Intensities will be gene-							//
//						rated depending on partitioning to the subdomains.											//
//																													//
//																													//

// Dependancies:
// -ImageJ RandomJ plugin

// ATT GÖRA
// - Add point spread function - before need to add fake membrane image (self avoiding random walk)
// - Add a third domain,(not any of the two domain) ( Can I just make this the inverse of the membrane image?)

// Görande
// - Thershold loop to select pixels that wil be acounted for in the correlation calculcation ish ( se rad 554)

// Gjort
// - Cartoon images of membrane and pixels
// - Loop to change subdomain fraction 
// - Large loop, get mean std
// - Add noise with randomJ (poission)
// - Increase mebrane range
// - Added to Spearman as an option
// - laga när det är 1-1



//SETTINGS------------------------------------------------------------------------------------------------------------
prnt = 0;		// 1 - print during macro      				0 - no printing,			2 - Results printed in log
spear = 0;		// 1 - Spearman correlation					0 - Pearson correlation
scat = 0;	 	// 1 - show no scatterplot					0 - no scatterplots 		*Does not work with batch!!
batch = 1;		// 1 - hide images							0 - show images

sz = 100; 						// Image size
noise = 0.05;					// Setting std for the rudimentary noise based on adding a gausian with mean 0 // Not in use
								

start_fraction = 0.1; 			// Fraction of membrane made up of subdomain 1, 0.8 will give 80%. This is start value
frac_step = 0.1;				// Change in fraction proportions, Currently set to not run, but works

mem_volume_start = newArray(1, 1);	// (low, high) start volume in multiple, eg 1,4 > pixel with most membrane = 4*lowest
mem_step = 0.25;					// Change in mabrane range each step of the loop

thr_setup = newArray(5, 0);	// (runs, step): reuns are amount of measurepoints and step is diffrens in intensity
								// between them. step is currently set to use 10% intensity steps - line 580

// Subdomains - start-partition 	1-raft  2-non raft   	3-non membrane		- assume fluoresecnce unaffected by medium
flr_a_partition	 = 		newArray(	10,		1,				0);	//				- starting partitioning for loop below
flr_b_partition	 =		newArray(	10,		1,				0);
partition_step = 0.1 	//															- steps to take when changing partition

// Subdomain distribution and distribution settings
sig = 16; 
lamb = 0.04; // ------------not implemented
clip = 0; // ------------not implemented - 
distribution = 1;	//	0 - Poisson,	1- half-Gaussian,	2 - Exponential,	3 - uniform (Only Poi / gau currently work)
//------------------------------------------------------------------------------------------------------------SETTINGS


// CLEAR WORKSPACE-----------------------------------------------------------------------------------------------------
if (nImages>0) run("Close All"); 							// if there are 1 or more images - close all
if (isOpen("Results")) { 									// close results window, if open
	selectWindow("Results");
	close("Results");
}
print("\\Clear"); 											// empty log window
roiManager("reset"); 										// empty ROI manager
run("Options...", "iterations=1 count=1 black edm=32-bit"); // set Binary Options
setOption("ScaleConversions", true);						// Set scaling for conversions
//-----------------------------------------------------------------------------------------------------CLEAR WORKSPACE

//SETUP ROI AND RESULTS -----------------------------------------------------------------------------------------------
newImage("Base", "8-bit black", sz, sz, 1); // Creating an image with the size sz

newImage("ones", "8-bit black", sz, sz, 1);
run("Add...", "value=1");

run("Select All");
roiManager("Add"); // Creating a ROI the size of the selected image
ROI = roiManager("count")-1;
run("Select None");

//-------------------------------------------------------------------------------------------------SETUP ROI AND RESULTS

setBatchMode(batch);

// add a loop here to run the simulation multiple times
for (n = 1; n < 20; n++) { // loop for running the simulation multiple times, multiple cells
	mem_volume = newArray(mem_volume_start[0],mem_volume_start[1]);
	
	for (i = 0; i < 6; i += 1) { // loop for varying range of different membrane amount in pixels
		if (i != 0) mem_volume[1] += mem_step; // the steps of increasing membrane raio in pixels
		
		if (mem_volume[0] == mem_volume[1]) { // RandomJ cannot create a uniform distribution with min=max
		selectImage("ones");
		run("Duplicate...", "title=membrane_fraction"); // the ones-image from before is thus used
		membrane_faction_id = getImageID();
		}
		else {
		selectImage("Base"); // Black image, all values 0
		run("RandomJ Uniform", "min="+mem_volume[0]+" max="+mem_volume[1]+" insertion=Additive"); // randomJ to create a uniform dist
		rename("membrane_fraction");
		membrane_faction_id = getImageID();
		}
		
		selectImage("ones");
		imageCalculator("Subtract create 32-bit", "ones","membrane_fraction");
		rename("non_membrane_fraction");
		
		for (j = 0; j < 5; j++) { // loop for adjusting the subdomain size
			 sub1_fraction = start_fraction + j*frac_step; // Increasing the procentage of subdomain 1 with the given value from settings each loop
			
			// Creating two opposite images with distribution around the percentage of domain 1 in entire image, 
			selectImage("Base"); // Image creaated previously in set up with only zeros.
			
			// Poission distribution
			if (distribution == 0) {
				run("RandomJ Poisson", "mean=" + sub1_fraction*1000 + " insertion=Additive");  	// fraction of the pixel thats made up of subdomain 1
				getRawStatistics(nPixels, mean, min, max, std, histogram);	// 800 will give 4:1 average ratio of subdomains
				changeValues(1000, max, 999); // all values between above 100% become 99%
				run("Divide...", "value=1000");
			}
			
			// Gaussian distribution
			if (distribution == 1) {
				run("RandomJ Gaussian", "mean="+ sub1_fraction*255+ " sigma=sig insertion=Additive");
				run("Divide...", "value=255");
				
				getRawStatistics(nPixels, mean, min, max, std, histogram);
				changeValues(min, 0, 0); // all values between above 100% become 99%
				changeValues(1, max, 1);
			}
			
			// Exponential distribution
			if (distribution == 2) {
				run("RandomJ Exponential", "lambda=lamb insertion=Additive");
				run("Enhance Contrast...", "saturated=0 normalize");
			}
			
			// Unifrom distribution
			if (distribution == 3) {
				run("RandomJ Uniform", "min=0 max=1 insertion=Additive");
			}
		
			// setMinAndMax(0, 0.99); //Adjustment in display
			rename("sub1"); // First domain is created using the randomJ plugin with the intensity distribution from "settings"
			sub1_id = getImageID();
			
			selectImage("ones"); // An image previously created with only ones.
			imageCalculator("Subtract create 32-bit", "ones","sub1"); // fraction of the pixel thats made up of subdomain 1
			rename("sub2");// subdomain 2, fraction subdomain 2;   sub1 + sub2 = 1
			sub2_id = getImageID();
			
			// The amount of each domain total -----------------------------------------Currently not in use
			imageCalculator("Multiply create 32-bit", "sub1","membrane_fraction");
			rename("amount_sub1");
			
			imageCalculator("Multiply create 32-bit", "sub2","membrane_fraction");
			rename("amount_sub2");
			// close("membrane_fraction"); //Closing the initial membrane_fraction image
			// The amount of each domain total -----------------------------------------Currently not in use				
			
			// Looping function for calculating correlation before and after multiplication with membrane amount
			b_part = newArray(flr_b_partition[0], flr_b_partition[1]); // Reseting the flr b partition
			varycorrelation_loop(flr_a_partition, b_part, partition_step, sub1_id, membrane_faction_id, sub1_fraction, ROI, n, spear, thr_setup, prnt);
			
			close("sub1");
			close("sub2");
			close("amount_sub1");
			close("amount_sub1");
		} // end of j-loop (adjustment of sudomain size)
		
		close("membrane_fraction");
		close("non_membrane_fraction");
		
	} // end of i-loop ( membrane fraction adjustment)
	
} //end of n-loop (amount of times simulation is run

print("Macro finished");
selectWindow("Log");
function scatter(size, img1_id, img2_id, ROI) { 
	// making a scatter plot from the intensities in two images on the x- and y-axis respectivley. Currently accepts squares
	// Size - size of the selected images, ex if 100x100: give size = 100
	
	// make scatterplot blank image
	newImage("Scatterplot", "16-bit black", 2*size, 2*size, 1);
	scatID = getImageID();
	
	// Getting the ROI as coordinate values
	roiManager("select", ROI);
	Roi.getContainedPoints(xArray, yArray);// all points within ROI
	run("Select None");
	
	// find intensity range for the two images - use to scale scatterplot
	selectImage(img1_id);
	setOption("ScaleConversions", true);
	resetMinAndMax(); // adjustments in brightness and contrast affect the 8-bit conversion
	run("16-bit");
	getRawStatistics(nPixels, meanX, minX, maxX);// have max x axis value
	
	selectImage(img2_id);
	setOption("ScaleConversions", true);
	resetMinAndMax(); // adjustments in brightness and contrast affect the 8-bit conversion
	run("16-bit");
	getRawStatistics(nPixels, meanY, minY, maxY);// have max y axis value
	setBatchMode(1);// stop updating screen - faster
	
	for(i = 0; i < xArray.length; i++) { // for every pixel in ROI
		selectImage(img1_id);
		vX = getPixel(xArray[i], yArray[i]);
		vXpixel = round(vX*(size*2)/maxX);
				
		selectImage(img2_id);
		vY = getPixel(xArray[i], yArray[i]);
		vYpixel=szY = 512 - round(vY*(size*2)/maxY);
			
		selectImage(scatID); // scatterPlot image
		vNew = 1 + getPixel(vXpixel, vYpixel); // read initial value
		setPixel(vXpixel, vYpixel, vNew); // update value in scatterplot
			
	}	// i loop
	setBatchMode(0);// stop updating screen - faster
	
	//display settings
	selectImage(scatID);// scatterplot
	getRawStatistics(nPixels, mean, minScat, maxScat);
	setMinAndMax(minScat, maxScat);
	
	return scatID
	corr
} // end of scatter-function

function correl(img1_id, img2_id, ROI, prnt) { 
// This is a variant of the Pcorrel macro. It calculates Pearson between the two selected images
// Note that "ROI" has to refer to an existing ROI-manager nr. And prnt = 1 gives print, 0 no print
	
	run("Select None");// roiManager("select", ROI);
	
	selectImage(img1_id);
	run("Select All");
	run("Duplicate...", " ");
	roiManager("select",ROI);
	getRawStatistics(n,mean,min,max,std);
	if (prnt == 1) print("stats for SRC-img1     min ",min," max ",max,"  STD ",std);
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
	if (prnt == 1) print("stats for SRC-img2     min ",min," max ",max,"  STD ",std);
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
	
	corr=numerator/denom;

// tidy up - delete images used in calculation
	selectImage(NumPicID);
	close();
	selectImage(GrnSubmPicID);
	close();
	selectImage(RedSubmPicID);
	close();
//print('end correlation');
	print("correlation",corr);
	return corr
} // end of function Pcorrel

function make_ranked(img1_id, img2_id, ROI, prnt) { //Makes image ranked in the ROI
	//set up ROI and Imagename
	ranked_img_array = newArray(2); // Creating an array for storing the newly created rank image-ids
	
	StackRankPicID = getImageID;// stack of 4 images
	Stackname = getTitle();
	
	roiManager("select",ROI); //selecting the chosen ROI-file
	setSelectionName("Area-Rank");
	roiManager("Add");// copy selection - check is there a selection?
	StackRankROIindex=roiManager("count")-1;
	run("Select None");// turn off selection

	// make a Binary image
	t0 = getTime();
	selectImage(img1_id);
	run("Duplicate...", "use");// one picture
	rename("Binary-rank");
	resetMinAndMax(); // adjustments in brightness and contrast affect the 8-bit conversion
	run("8-bit");
	roiManager("select", StackRankROIindex);
	getRawStatistics(nEntries);
	XposArray = newArray(nEntries);
	YposArray = newArray(nEntries);
	setColor(0);
	run("Clear Outside");
	setColor(255);
	run("Fill"); //run("Fill", "slice");
	
	// make arrays to hold XY coordinates for pixels within ROI
	getSelectionBounds(xst, yst, width, height);// smaller area to check than whole image
	CoorPosInArray = 0;
	for (xat = xst; xat < width + xst; xat = xat + 1) { // boundaries of ROI
		for (yat = yst; yat < height + yst; yat = yat + 1) {
			if (getPixel(xat, yat) == 255) {
			//if(selectionContains(xat, yat)) {  // important bit - slow
				XposArray[CoorPosInArray] = xat;
				YposArray[CoorPosInArray] = yat;
				   //print(XposArray[CoorPosInArray],YposArray[CoorPosInArray]);
					CoorPosInArray = CoorPosInArray + 1;// next location
					  //setPixel(xat,yat,999);// show which pixels accessed
			} // if inside
		} // yat
	} // xat

	//printing and time-check ???
	t1 = getTime();
	if (prnt==1) print("mk coordinate list from binary image",(t1-t0)/1000,"secs");
	if (prnt==1) print("N of points",XposArray.length);
	selectWindow("Log");
	close();// close binary image

	// have coordinate list in 2 arrays	XposArray & YposArray
	// Rank each image in the stack
	for (i = 1; i <= 2; i++) {   // for each of the two images given to the function
		
		if (i == 1){
			selectImage(img1_id);
			if (prnt == 1) print(" ran ranked for first image");
			run("Duplicate...", "title=ranked_img1");
			ranked_img_array[0] = getImageID();
		}
		else {
			selectImage(img2_id);
			run("Duplicate...", "title=ranked_img2");
			ranked_img_array[1] = getImageID();
			if (prnt == 1) print(" ran ranked for second image");
		}
		
		run("Select None");
		resetMinAndMax(); // adjustments in brightness and contrast affect the 8-bit conversion
		run("16-bit");// 16 bit works better for getRawStatistics-----------------------------------------????
		roiManager("select", StackRankROIindex);// ROI
		getRawStatistics(nEntries, mean, minI, maxI, std, histROIArray); //histArray = antal pixlar med den intensiteten
		
		// (1)Make an array covering the range of intensities with the replacement rank
		t0 = getTime();
		rankArray = newArray(maxI+1);// array filled with zeros
		if (prnt==1) print("rankArray.length",rankArray.length,"pixels ",nEntries);
		
		// for each intensity find the replacement rank
		Nconverted = 0;// number of values converted
		//print("rank array length",rankArray.length);
		for (n=minI;n<=maxI;n=n+1) {
			Nval=histROIArray[n];
			//print(n,Npixels);
			if(Nval>=1) { // non zero number of pixels
				mRank=Nconverted+Nval/2;// mean rank
				Nconverted=Nconverted+Nval;  // N pixels now ranked
				//print(n,"    Convert",Nval,"  rank",mRank,"     converted ",Nconverted);
				rankArray[n]=mRank + 0.5; // mean rank for each intensity
				//print("rank ",rankArray[n]);
			} // if Npixels
		} // n
		t1 = getTime();
		if (prnt == 1) print("make rank array",(t1-t0)/1000,"sec");
		
		// (2) using list of coordinates, replace value with rank
		t0 = getTime();
		run("32-bit");// some ranks maybe non integer
		for (j = 0; j < XposArray.length; j++) {
			x = XposArray[j];
			y = YposArray[j];
			//print(i,x,y);
				v = getPixel(x,y);
				vn = rankArray[v];
				setPixel(x, y,vn);
				//print(i,"pos ",x,y," val",v,"-",rankArray[v]);
		}
	t1 = getTime();
	if (prnt == 1) print("replaced with rank",(t1-t0)/1000,"sec");
	} // slc Loop
	if (prnt == 1) print("rank Images in make ranked:     rank_a: ",ranked_img_array[0], "rank_b: ", ranked_img_array[1]);
	return ranked_img_array;

}// end of function makeRanked

function varycorrelation_loop(flr_a_partition, flr_b_partition, partition_step, sub1_id, membrane_faction_id, sub1_size, roi, n, spear, thr_setup, prnt) { 
	// The main loop for creating images from subdomain and flourophore and measuring correlation coeffcient.
	// Takes:
	// - 	Start-partitioning in array for fluorphore a and b. and step in loop between partitions
	// -	ImageId for created subdomain 1 image.
	// -	ImageId for created membrane fraction image
	// -	id for the ROI in use
	// - 	n for the amount of times the macro has been run
	// - 	spear, binary option for Spearman or Pearson
	
	// Dependant on: "Scatter" and "Correl" functions.
	
	// Thresholding setup
	thr_n = thr_setup[0];
	thr_step = thr_setup[1];
	
	// Getting the title of the membrane fraction image for use in the image-calculator
	selectImage(membrane_faction_id);
	getRawStatistics(nPixels, mean, min_mem, max_mem, std, histogram);
	mem_fraction_title = getTitle();
	
	// Creating the sub2 image derived from the sub1, a complete opposite
	selectImage(sub1_id);
	getDimensions(width, height, channels, slices, frames);
	newImage("1s", "8-bit black", width, height, 1);
	run("Add...", "value=1");
	imageCalculator("Subtract create 32-bit", "1s","sub1"); // fraction of the pixel thats made up of subdomain 1
	rename("sub2");// subdomain 2, fraction subdomain 2;   sub1 + sub2 = 1
	sub2_id = getImageID();
	
	// Testing
	selectImage(sub2_id);
	run("Multiply...", "value=255");
	selectImage(sub1_id);
	run("Multiply...", "value=255");
	
	
	// Preparation for the measurment loop. flr_b_partition is derived from settings and reset with every new membrane trial
	// Creating arrays for storing the correlation values
	iterration = ((flr_a_partition[0]-1)/partition_step)*2 + 1; // every measuring step
	iterration = 45;


	value_nr = 0; // for seleting and printing values from the arrays
	corr_array = newArray(iterration*thr_n);
	corr_pre_array = newArray(iterration*thr_n);
	partition_array_a = newArray(iterration*thr_n);
	partition_array_b = newArray(iterration*thr_n);
	partition_array_diff = newArray(iterration*thr_n);
	pre_pixel_array = newArray(iterration*thr_n);
	mem_pixel_array = newArray(iterration*thr_n);
	
	part_k = flr_b_partition[1] / flr_b_partition[0];
	tot_flr = flr_b_partition[0] + flr_b_partition[1];
	
	for (j = 0; j < iterration; j += 1) {
	//	part_k -= 0.1
		
		
		
		if (j != 0) {
			if (flr_b_partition[0] == 1) {
				flr_b_partition[1] = flr_b_partition[1]*(1 + partition_step);
			}
			else {
				flr_b_partition[0] = flr_b_partition[0]*(1 - partition_step); // Adding to the partitioning to subdomain 2 for fluorophore B
				if (flr_b_partition[0] < 1) {
					flr_b_partition[0] = 1;
				}
			}
		}
	
		// The amount of the two fluorophores in each pixel. val is their partition to the specific domain
		// Fluorophore A
		selectImage(sub1_id);
		run("Duplicate...", "title=flr_a_sub1");
		val = flr_a_partition[0];
		run("Multiply...", "value="+val);
				
		selectImage(sub2_id);
		run("Duplicate...", "title=flr_a_sub2");
		val = flr_a_partition[1];
		if (prnt == 1) print("val flr a sub2:", val); 
		run("Multiply...", "value="+val);
						
		// Need to add the 3rd domain below here
		
		imageCalculator("Add create 32-bit", "flr_a_sub1","flr_a_sub2");
		rename("total_flr_a");
		fluor_a_tot = getImageID();
		
		imageCalculator("Multiply create 32-bit", "total_flr_a", mem_fraction_title);
		rename("total_flr_a_mem");
		fluor_a_tot_mem = getImageID();
				
		// Fluorophore B-------------------------------Currently where the alteration of fluorophore partitioning happens
		selectImage(sub1_id);
		run("Duplicate...", "title=flr_b_sub1");
		val = flr_b_partition[0];
		if (prnt == 1) print("val flr b sub1:", val);
		run("Multiply...", "value="+val);
					
		selectImage(sub2_id);
		run("Duplicate...", "title=flr_b_sub2");
		val = flr_b_partition[1];
		run("Multiply...", "value="+val);
				
		// need to add the 3rd domain below here
		
		imageCalculator("Add create 32-bit", "flr_b_sub1","flr_b_sub2");
		rename("total_flr_b");
		fluor_b_tot = getImageID();
				
		// Adding the membrane-factor created in loop i
		imageCalculator("Multiply create 32-bit", "total_flr_b", mem_fraction_title);
		rename("total_flr_b_mem");
		fluor_b_tot_mem = getImageID();
		
		// Adding poisson noise with the RandomJ plugins modulatory function
		selectImage(fluor_a_tot);
		resetMinAndMax(); // adjustments in brightness and contrast affect the 8-bit conversion
		run("8-bit");
		if (flr_a_partition[0] == flr_a_partition[1]) run("Add...", "value=1"); // will add 1 if all values are 0, this is to enable randomJ to create noise
		run("RandomJ Poisson", "mean=1 insertion=Modulatory"); // adding noise (modulatory takes random value for every pixel with that being the mean)
		rename("total_flr_a_noise");
		fluor_a_tot_noise = getImageID();
		
		selectImage(fluor_a_tot_mem);
		resetMinAndMax(); // adjustments in brightness and contrast affect the 8-bit conversion
		run("8-bit");
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		if (max == 0) run("Add...", "value=1"); // will add 1 if all values are 0, this is to enable randomJ to create noise
		run("RandomJ Poisson", "mean=1 insertion=Modulatory"); // adding noise (modulatory takes random value for every pixel with that being the mean)
		rename("total_flr_a_mem_noise");
		fluor_a_tot_mem_noise = getImageID();
	
		selectImage(fluor_b_tot);
		resetMinAndMax(); // adjustments in brightness and contrast affect the 8-bit conversion
		run("8-bit");
		if (flr_b_partition[0] == flr_b_partition[1]) run("Add...", "value=1"); // will add 1 if all values are 0, this is to enable randomJ to create noise
		run("RandomJ Poisson", "mean=1 insertion=Modulatory"); // adding noise (modulatory takes random value for every pixel with that being the mean)
		rename("total_flr_b_noise");
		fluor_b_tot_noise = getImageID();
		
		selectImage(fluor_b_tot_mem);
		resetMinAndMax(); // adjustments in brightness and contrast affect the 8-bit conversion
		run("8-bit");
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		if (max == 0) run("Add...", "value=1"); // will add 1 if all values are 0, this is to enable randomJ to create noise
		run("RandomJ Poisson", "mean=1 insertion=Modulatory"); // adding noise (modulatory takes random value for every pixel with that being the mean)
		rename("total_flr_b_mem_noise");
		fluor_b_tot_mem_noise = getImageID();
		
		// Setting a threshold with the option of looping to do it multiple times
		for (k = 0; k < 1; k++) {
			// Setting threshold for pre images k = 0; k < thr_n; k++
			run("Select None");
			resetThreshold;
			selectImage(fluor_b_tot_noise);
			getRawStatistics(nPixels, mean, b_min, b_max, std, histogram);
			
			if (b_min + k*thr_step < b_max) {// ------------------------------------------------- known problem here
				setThreshold(b_min + (b_max-min)*k*0.1, b_max); // setThreshold(b_min + k*std, b_max); // setThreshold(b_min + k*thr_step, b_max); setThreshold(b_min + (b_max-min)*k*0.1, b_max);
				run("Create Selection");
				roiManager("Add");
				a_tot_roi = roiManager("count")- 1;
				run("Select None");
				resetThreshold;
				
				selectImage(fluor_a_tot_noise);
				getRawStatistics(nPixels, mean, min, max, std, histogram);
				setThreshold(min + (max-min)*k*0.1, max); // setThreshold(min + k*std, max); // setThreshold(min + k*thr_step, max); setThreshold(min + (max-min)*k*0.1, max);
				run("Create Selection");
				roiManager("Add");
				b_tot_roi = roiManager("count")- 1;
				run("Select None");
				resetThreshold;
				
				roiManager("Select", newArray(a_tot_roi, b_tot_roi));
				roiManager("AND");
				roiManager("Add");
				pre_roi = roiManager("count")- 1;
				roiManager("Select", pre_roi);
				roiManager("rename", "pre_roi");
				
				// Setting threshold for mem-images
				selectImage(fluor_b_tot_mem_noise);
				getRawStatistics(nPixels, mean, min, max, std, histogram);
				setThreshold(min + (max-min)*k*0.1, max); // setThreshold(min + k*std, max); // setThreshold(min + k*thr_step, max); setThreshold(min + (max-min)*k*0.1, max);
				run("Create Selection");
				roiManager("Add");
				b_tot_mem_roi = roiManager("count")- 1;
				run("Select None");
				resetThreshold;
				
				selectImage(fluor_a_tot_mem_noise);
				getRawStatistics(nPixels, mean, min, max, std, histogram);
				setThreshold(min + (max-min)*k*0.1, max); // setThreshold(min + k*std, max); // setThreshold(min + k*thr_step, max); setThreshold(min + (max-min)*k*0.1, max);
				run("Create Selection");
				roiManager("Add");
				a_tot_mem_roi = roiManager("count")- 1;
				run("Select None");
				resetThreshold;
					
				roiManager("Select", newArray(a_tot_mem_roi, b_tot_mem_roi));
				roiManager("AND");
				roiManager("Add");
				mem_roi = roiManager("count")- 1;
				roiManager("Select", mem_roi);
				roiManager("rename", "mem_roi");

				if (prnt == 1) print("mem_roi: ", mem_roi, " pre_roi: ", pre_roi); // makeRectangle(0, 0, 1, 1);
			}
			else {
				pre_roi = roi;
				mem_roi = roi;
			}
			thr = 1;
			if (thr == 0) {
				pre_roi = roi;
				mem_roi = roi;
			}

			
			if (spear == 1) { // decided in setting weather spearman or Pearson correlation coeffcient should be used
				rank_array = make_ranked(fluor_a_tot_noise, fluor_b_tot_noise, pre_roi, prnt); // Using make_ranked-function to rank pixels
				rank_a = rank_array[0];
				rank_b = rank_array[1];
	
				rank_array = make_ranked(fluor_a_tot_mem_noise, fluor_b_tot_mem_noise, mem_roi, prnt); // Using make_ranked-function to rank pixels
				rank_a_mem = rank_array[0];
				rank_b_mem = rank_array[1];
				
				if (prnt == 1) print("rank Images in varycoloc:     rank_a: ",rank_array[0], "rank_b: ", rank_array[1]);
				
				// saving correlation data in array, using the "correl"-function to calculate Spearman correlation coefficient
				corr_pre_array[value_nr] = correl(rank_a, rank_b, pre_roi, prnt);
				corr_array[value_nr] = correl(rank_a_mem, rank_b_mem, mem_roi, prnt);	
			}
			else {
				// saving correlation data in array, using the "correl"-function to calculate Pearson correlation coefficient
				corr_pre_array[value_nr] = correl(fluor_a_tot_noise, fluor_b_tot_noise, pre_roi, prnt);
				corr_array[value_nr] = correl(fluor_a_tot_mem_noise, fluor_b_tot_mem_noise, mem_roi, prnt);
				
			}
					
			partition_array_a[value_nr] = toString(flr_a_partition[0] / flr_a_partition[1]);
			partition_array_b[value_nr] = toString(flr_b_partition[0] / flr_b_partition[1]);
			partition_array_diff[value_nr] = toString(flr_a_partition[1]+flr_b_partition[0]-flr_a_partition[0]-flr_b_partition[1]);
					
			// Create, display and close a scatterplot
			if (scat == 1) {
				scatID = scatter(sz, fluor_a_tot, fluor_b_tot, pre_roi);
				selectImage(scatID);// scatterplot
				wait(1000);
				close("Scatterplot");
			}
			
			// Getting ROI information
			roiManager("select", pre_roi);
			getRawStatistics(pre_roi_pixels);
			run("Select None");
			pre_pixel_array[value_nr] = pre_roi_pixels;
			
			roiManager("select", mem_roi);
			getRawStatistics(mem_roi_pixels);
			run("Select None");
			mem_pixel_array[value_nr] = mem_roi_pixels;
			
			// Clean-up in the ROI-manager
			if (b_min + k*thr_step < b_max) {
				roiManager("Select", newArray(a_tot_roi, b_tot_roi, a_tot_mem_roi, b_tot_mem_roi, pre_roi, mem_roi));
				roiManager("Delete");
			}
			
			if (prnt == 1) print("running the thresholding");
		
			// Uppdating the results
			extend = nResults;
			setResult("Pearson_pre_corr", extend, corr_pre_array[value_nr]);
			setResult("Pearson_mem_corr", extend, corr_array[value_nr]);
			setResult("Partitioning_fluorophore_A", extend, partition_array_a[value_nr]);
			setResult("Partitioning_fluorophore_B", extend, partition_array_b[value_nr]);
			// setResult("Partitioning_fluorophore_diff", extend, partition_array_diff[value_nr]);
			setResult("Membrane_range", extend, (max_mem / min_mem));
			setResult("Subdomain_1_size", extend, sub1_size);
			setResult("Pixels_in_threshold_pre", extend, pre_pixel_array[value_nr]);
			setResult("Pixels_in_threshold_mem", extend, mem_pixel_array[value_nr]);
			setResult("Threshold_run", extend, k+1);
			
			setResult("Cell_nr", extend, n);
			updateResults();   
			selectWindow("Results");
			value_nr += 1;
			
		}
	
		
		//Closing images
		close("flr_a_sub1");
		close("flr_a_sub2");
		close("total_flr_a");
		close("flr_b_sub1");
		close("flr_b_sub2");
		close("total_flr_b");
		close("total_flr_a_mem");
		close("total_flr_b_mem");
		close("1s");
		close("total_flr_a_noise");
		close("total_flr_a_mem_noise");
		close("total_flr_b_noise");
		close("total_flr_b_mem_noise");
		
		if (spear == 1) {
			close(rank_a);
			close(rank_b);
			close(rank_a_mem);
			close(rank_b_mem);
		}
	}// end of j
			
		/* Print the enairty of corr_array before it gets replaced
	for (j = 0; j < iterration; j++) {
		extend = nResults;
		setResult("Pearson_pre_corr", extend, corr_pre_array[j]);
		setResult("Pearson_mem_corr", extend, corr_array[j]);
		setResult("Partitioning_fluorophore_A", extend, partition_array_a[j]);
		setResult("Partitioning_fluorophore_B", extend, partition_array_b[j]);
		setResult("Partitioning_fluorophore_diff", extend, partition_array_diff[j]);
		setResult("Membrane_range", extend, (max_mem / min_mem));
		setResult("Subdomain_1_size", extend, sub1_size);
		setResult("Pixels_in_threshold", extend, pixel_array[j]);
		
		setResult("Cell_nr", extend, n);
		updateResults();   
		selectWindow("Results");
	}
		*/
		// close("membrane_fraction_id");
} // end function varycorrelation

