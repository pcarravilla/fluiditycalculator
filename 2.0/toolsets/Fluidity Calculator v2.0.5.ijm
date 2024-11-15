/* 
 *  ----------------------------------------------------------------------------------
 *  Author: Pablo Carravilla
 *  
 *  Name: Fluidity Calculator
 *  Version: 2.0.5
 *  Date: 15-11-2024
 *  
 *  Description: 
 *  This macro tool calculates GP images from multichannel microscopy images.
 *  
 *  For more information see: https://github.com/pcarravilla/fluiditycalculator
 *  Publication: Carravilla, Andronico, Schelegel et al (2024)
 *  
 *  Copyright (c) [2024] [Pablo Carravilla]
 *  
 *  License: GNU General Public License v3.0
 *  For more information, please refer to the LICENSE file provided with this project.
 *  
 *  ----------------------------------------------------------------------------------
 */
 
macro "Fluidity Calculator v2.0 Action Tool - icon:fluidity.png" {
		
	// Preparation
	setBatchMode(true);
	requires("1.54f");
	roiManager("reset");
		
	
	print(" ");
	print("Running Fluidity Calculator");
	
	
	// If there are no open images let the user open one
	imageList = getList("image.titles");
	if (imageList.length == 0) {
		openPath = File.openDialog("Open...");
		open(openPath);
	}
	
	
	
	// Get image name, location, bit depth and date/time
	imgName = getTitle();
	imgDir = getDirectory("image");
	imgDepth = bitDepth(); 
	dateTime = whatDateTime();
	nameDate = imgName + "_" + dateTime;
	defaultsPath = getDirectory("macros") + "toolsets/defaultsGP.csv";
	
	print("Image: " + imgName);
	print("Date/Time: " + dateTime);
	
	
	
	// Duplicate image to avoid modifying raw data and convert to 32-bit
	rawDataID = getImageID();
	run("Duplicate...", "title=dat duplicate");
	run("32-bit");
	getDimensions(dataWidth, dataHeight, dataChannels, dataSlices, dataFrames);
	dataID = getImageID();
	
	
	
	// Stop if the image does not have at least two channels
	if (dataChannels < 2) {
		exit("Images must have at least 2 spectral channels. \n" + 
			" \n" +
			"*Adjust Bio-Formats import settings or combine channels using: \n" +
			"      Image --> Color --> Merge Channels \n" +
			"            or \n" +
			"      Image --> Stack --> Images to Stack \n" +
			" \n" +
			"*If spectral channels are wrongly assigned as z slices or time points (Image --> Properties), re-assign them as channels using: \n" +
			"      Image --> Hyperstacks --> Re-order Hyperstack"
			
		);
		selectImage(dataID);
		close();}
	
	
	
	// Load defaults
	defaultSettingNumber = 13;
	defaultSettings = loadDefaults(defaultsPath, defaultSettingNumber); // Loads default settings to an array
	
	
	
	// Check if channel labels are wavelengths in spectral images (channels >= 5), if not relabel
	firstWavelengthDefault = defaultSettings[11];
	windowWavelengthDefault = defaultSettings[12];
	
	if (dataChannels >= 5) {newWavelengthDefaults = checkLabels(dataID, firstWavelengthDefault, windowWavelengthDefault);}
	else {newWavelengthDefaults = newArray(500, 8.9);}
	firstWavelengthUser = newWavelengthDefaults[0]; 
	windowWavelengthUser = newWavelengthDefaults[1];
	
	
	
	// Extrat wavelengths from channel labels
	wavelengthList = getChannelLabels(dataID);
	
	
	
	// Show GUI
	showGUI(defaultSettings, wavelengthList);
	
	
	
	// Get user input and store it in an array (userSettings)
	ch1LabelUser = Dialog.getChoice();
	ch2LabelUser = Dialog.getChoice();
	widthUser = Dialog.getNumber();
	doAveragingUser = Dialog.getCheckbox();
	averagingRadiusUser = Dialog.getNumber();
	thresholdChUser = Dialog.getChoice();
	satWarningPercUser = Dialog.getNumber();
	showSpectrumUser = Dialog.getCheckbox(); 
	saveResultsUser = Dialog.getCheckbox();
	saveAnalysisMaskUser = Dialog.getCheckbox();
	saveLocationUser = Dialog.getChoice();
	saveAsDefaultsUser = Dialog.getCheckbox();
	
	userSettings = newArray(defaultSettingNumber); // To store new defaults
	userSettings[0] = ch1LabelUser;
	userSettings[1] = ch2LabelUser;
	userSettings[2] = widthUser;
	userSettings[3] = doAveragingUser;
	userSettings[4] = averagingRadiusUser;
	userSettings[5] = thresholdChUser;
	userSettings[6] = satWarningPercUser;
	userSettings[7] = showSpectrumUser;
	userSettings[8] = saveResultsUser;
	userSettings[9] = saveAnalysisMaskUser;
	userSettings[10] = saveLocationUser;
	userSettings[11] = firstWavelengthUser;
	userSettings[12] = windowWavelengthUser;
	
	
	
	// Save settings as default
	if (saveAsDefaultsUser == 1) {saveAsDefaults(defaultsPath, userSettings);}
	
	
	
	// Get channel range
	channelRangeIndexesUser = defineChannelRange(dataID, wavelengthList, ch1LabelUser, ch2LabelUser, widthUser);
	ch1Low = channelRangeIndexesUser[0];
	ch1High = channelRangeIndexesUser[1];
	ch2Low = channelRangeIndexesUser[2];
	ch2High = channelRangeIndexesUser[3];
	
	
	
	// Define save location
	if (saveResultsUser == 1 || saveAnalysisMaskUser == 1) {
		if (saveLocationUser == "Image Directory") {
			saveDirUser = imgDir;}
		else if (saveLocationUser == "Image Directory/GP_Results") {
			saveDirUser = imgDir + "GP_Results" + File.separator;
			File.makeDirectory(saveDirUser);
		}
		else if (saveLocationUser == "Select") {
			saveDirUser = getDirectory("Select a directory");
		}
	}
	else saveDirUser = "";
	
	
	
	// Saturation
	
	// Create GP saturation mask (negative, i.e. saturated pixels = 0)
	// It will be used to mask the GP image
	newImage("GP Saturation Mask", "32-bit grayscale-mode", dataWidth, dataHeight, 1, dataSlices, dataFrames);
	gpSatMaskID = getImageID();
	
	for (z = 1; z <= dataSlices; z++) {
		for (f = 1; f <= dataFrames; f++) {
			generateSaturationMask(1, dataID, gpSatMaskID, imgDepth, f, z, ch1Low, ch1High, ch2Low, ch2High);
		}
	}
	
	
	// Create saturation mask (positive, i.e. saturated pixels = 1)
	// It will be used to calculate the percentage of saturated pixels
	newImage("Saturation Mask", "32-bit grayscale-mode", dataWidth, dataHeight, 1, dataSlices, dataFrames);
	percSatMaskID = getImageID();
	
	for (z = 1; z <= dataSlices; z++) {
		for (f = 1; f <= dataFrames; f++) {
			generateSaturationMask(0, dataID, percSatMaskID, imgDepth, f, z, ch1Low, ch1High, ch2Low, ch2High);
		}
	}
	
	
	
	//Threshold
	
	//Create threshold reference image
	if (thresholdChUser == "Channel 1+2") {
		extractChannels(dataID, ch1Low, ch1High); // Ch1
		rename("tempCh1");
		tempCh1ID = getImageID();
		extractChannels(dataID, ch2Low, ch2High); // Ch2
		rename("tempCh2");
		tempCh2ID = getImageID();
		
		imageCalculator("Add create 32-bit stack", tempCh1ID, tempCh2ID);
		rename("Set Threshold – Ch1+Ch2");
		refThrID = getImageID();
		
		selectImage(tempCh1ID); 
		close();	
		selectImage(tempCh2ID); 
		close();	
	}
	
	else if (thresholdChUser == "Channel 1") {
		extractChannels(dataID, ch1Low, ch1High);
		rename("Set Threshold – Ch1");
		refThrID = getImageID();
	}
	
	else if (thresholdChUser == "Channel 2") {
		extractChannels(dataID, ch2Low, ch2High);
		rename("Set Threshold – Ch2");
		refThrID = getImageID();
	}
	
	
	
	// Generate threshold mask and get threshold values
	thresholdValuesUser = userThreshold(refThrID, doAveragingUser, averagingRadiusUser);
	thrMaskID = getImageID();
	
	
	
	// Caculate saturated pixel percentage
	
	// Create mask of saturated pixels within the user thresholded pixels 
	imageCalculator("AND 32-bit create stack", percSatMaskID, thrMaskID);
	rename("Saturated Pixels Within User Treshold");
	satWithinThrMaskID = getImageID();
	getDimensions(maskWidth, maskHeight, maskChannels, maskSlices, maskFrames);
	
	selectImage(percSatMaskID);
	close();
	
	
	// Calculate satured pixels and show warning if above warning %
	
	// Warn user default values
	warnUserArray = newArray(2); // This array stores warnUser and showSatMask
	warnUser = 1; // Shows a warning if saturated pixel % above user defined value
	showSatMask = 0; // Shows saturation mask if saturated pixel % above user defined value in any frame/slice
	warnUserArray[0] = warnUser;
	warnUserArray[1] = showSatMask;
		
	for (z = 1; z <= maskSlices; z++) {
		for (f = 1; f <= maskFrames; f++) {
			warnUserArray = calculateSaturatedPixelPerc(thrMaskID, satWithinThrMaskID, f, z, satWarningPercUser, warnUserArray);
		}
	}
	
	// Get user warning/showSatMask choices 
	warnUser = warnUserArray[0]; 
	showSatMask = warnUserArray[1];
	
	
	// Show saturation/threshold mask (initial frame and z slice) if selected by user
	if (showSatMask == 1) {
		selectImage(satWithinThrMaskID);
		if (maskFrames > 1) {Stack.setFrame(1);}
		if (maskSlices > 1) {Stack.setSlice(1);}
		setBatchMode("show");
	}
	else if (showSatMask == 0) {
		selectImage(satWithinThrMaskID);
		close();
	}
	
	
	
	// Extract channels 1 and 2
	extractChannels(dataID, ch1Low, ch1High); // Ch1
	ch1ID = getImageID();
	rename("Ch1");
	extractChannels(dataID, ch2Low, ch2High); // Ch2
	ch2ID = getImageID();
	rename("Ch2");
	
	
	
	
	// Calculate GP from channels 1 and 2
	calculateGP(ch1ID, ch2ID);
	rename("GP_" + nameDate);
	gpID = getImageID();
	
	
	
	// Calculate and save (optional) GP analysis mask (pixels within threshold that are not saturated)
	imageCalculator("AND 32-bit create stack", thrMaskID, gpSatMaskID);
	analysisMaskID = getImageID();
	if (saveAnalysisMaskUser == 1) {saveAs("tiff", saveDirUser + "GP_Analysis_Mask_" + nameDate);}
	
	
	// Mask GP using analysis (threshold AND saturation) mask
	maskGP(gpID, analysisMaskID);
	
	
	
	// Show and save spectrum
	if (showSpectrumUser == 1 && dataChannels > 5) {
		showSpectrum(dataID, refThrID, analysisMaskID, nameDate);
		if (saveResultsUser == 1) {
			saveSpectrum(dataID, refThrID, analysisMaskID, nameDate, saveDirUser);
		}
	}
	
	else if (showSpectrumUser == 1&& dataChannels <= 5) {
		waitForUser("Emission spectrum cannot be shown. The image must have at least five channels");
	}
	
	
	
	// Generate and save (optional) results, parameters and histogram table
	generateResultsTable(gpID, saveResultsUser, saveDirUser, nameDate);
	if (saveResultsUser == 1) {
		numberOfBins = 100;
		saveHistogramTable(gpID, numberOfBins, saveDirUser, nameDate);
		saveParametersTable(userSettings, thresholdValuesUser, imgName, dateTime, saveDirUser, nameDate);
	}
	
	
	
	// Close intermediate images
	selectImage(thrMaskID);
	close();
	selectImage(gpSatMaskID);
	close();
	selectImage(dataID);
	close();
	selectImage(ch1ID);
	close();
	selectImage(ch2ID);
	close();
	selectImage(refThrID);
	close();
	selectImage(analysisMaskID);
	close();
	
	
	
	// GP LUT and Save
	selectImage(gpID);
	if (dataFrames > 1) {Stack.setFrame(1);}
	if (dataSlices > 1) {Stack.setSlice(1);}
	
	pathGPViridis = getDirectory("luts") + "gp-viridis.lut";
	pathMPLViridis = getDirectory("luts") + "mpl-viridis.lut";
	if (File.exists(pathGPViridis) == 1) {run("gp-viridis");} // Check if gp-viridis lut exists (Custom LUT filename)
	else if (File.exists(pathMPLViridis) == 1) {run("mpl-viridis");} // Check if mpl-viridis lut exists (Fiji default filename)
	
	if(saveResultsUser==1){saveAs("tiff", saveDirUser + "GP_" + nameDate);}
	
	// Histogram
	setMinAndMax(-1, 1); // Comment out to use automatic range
	run("Histogram", "bins=100 x_min=-1 x_max=1 y_max=Auto");
	
	if(saveResultsUser==1){saveAs("tiff", saveDirUser + "Histogram_" + nameDate);}
	
	// Tidy up
	selectImage(rawDataID);
	Stack.setChannel(ch1High);
	
	
	setBatchMode("exit and display");
	run("Tile");
	print("Finished running Fluidity Calculator");
	
	
	
	
	
	// ***************************************************** FUNCTIONS **************************************************************
	
	
	
	function generateSaturationMask(negativeMask, dataImg, satMaskImg, depth, satFrame, satSlice, lowCh1, highCh1, lowCh2, highCh2) { 
	// Checks pixel saturation in two channel ranges ([lowCh1-highCh1] and [lowCh2-highCh2]) and generates a mask
	// It can generate a positive mask, where saturated pixels = 1
	// or a negative mask, where saturated pixels = 0
	
		roiManager("reset");
		
		selectImage(dataImg);
		Stack.setFrame(satFrame);
		Stack.setSlice(satSlice);
		run("Select None"); // Unselect so that whole image is analysed
	
		getDimensions(dataImgWidth, dataImgHeight, dataImgChannels, dataImgSlices, dataImgFrames);
		
		// Check image depth and calculate maximum possible (saturated) pixel value
		saturatedPixelValue = Math.pow(2, depth)-1;
		
		// Iterate over the channels and store saturated ROIs
		for (c = 1; c <= dataImgChannels; c++) {
			// Check if the channel is in the range defined by the user for GP calculation
			if ((c >= lowCh1 && c <= highCh1) || (c >= lowCh2 && c <= highCh2)){
				Stack.setChannel(c);
				getStatistics(area, mean, min, max, std, histogram);
				if (saturatedPixelValue == max) { // Checks if any pixel is saturated
					setThreshold(max, max);
	
					// Measure saturated pixel number 
					satPixNumber = getValue("Area raw limit"); // raw to get pixel number (unscaled); limit to measure within threshold
									
					run("Create Selection");
					roiManager("add");
					run("Select None"); // Unselect ROI so that whole image is analysed in subsequent c iterations
					if (negativeMask == 1 && satPixNumber > 10) {print(satPixNumber + " saturated pixels in channel " + 
												c +	" frame " + satFrame + " of z slice " + satSlice + "."); 
											print("   The corresponding areas are excluded from the analysis.");
					}
				}
			}
		}
		
		// Create saturation mask
		roiNumber = RoiManager.size;
		
		// If there are saturated pixels, those get asigned as 0 in the mask
		if(roiNumber>0){
			//print("Saturation found in " + roiNumber + " channel(s) of frame " + satFrame + " slice " + satSlice);
			// Combine all ROIs from saturated pixels and NaN
			roiManager("Combine");
			roiManager("add");
			selectImage(satMaskImg); 
			// For some reason this has to be done before setting frame and slice. If not ROI Manager will change the frame (ImageJ bug?)
			roiManager("Select", roiNumber); // Selects the combination ROI (roiNumber is 1-based, index is 0-based)
			if (negativeMask == 1) {run("Make Inverse");}; // So that non saturated pixels = 1}
			if (dataImgFrames > 1) {Stack.setFrame(satFrame);}
			if (dataImgSlices > 1) {Stack.setSlice(satSlice);}
			run("Set...", "value=1");
		}
		
		// If there are no saturated pixels, the whole mask frame/slice will be 1
		else if (roiNumber == 0){
			selectImage(satMaskImg); 
			run("Select All");
			if (dataImgFrames > 1) {Stack.setFrame(satFrame);}
			if (dataImgSlices > 1) {Stack.setSlice(satSlice);}
			if (negativeMask == 1) {run("Set...", "value=1");}
			else if (negativeMask == 1) {run("Set...", "value=0");}
		}
		
		run("Select None"); // Unselect ROI so that whole image is analysed in subsequent f and z iterations
		
		roiManager("reset");
		
	}
	
	
	
	function extractChannels(dataImg, lowCh, highCh) { 
	// Extracts channels selected by user for GP calculation
	// It also removes saturated pixels 
	// lowCh and highCh are channel numbers
	
		
		// Duplicate dataID to re-oder
		// Channels are re-ordered to slices to sum them using the 'Z Project...' function
		selectImage(dataImg);
		run("Select None");
		run("Duplicate...", "title=reordered duplicate");
		
		run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
		reorderedImg = getImageID();
		run("Z Project...", "start=lowCh stop=highCh projection=[Sum Slices] all");
		// Correct order (z planes are now slices again)
		run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
		tempExtrImg = getImageID(); // There is a bug that prevents thresholding this image. It is solved by duplicating it.
		run("Duplicate...", "title=extractedChannels duplicate");
		extractedChannelsImg = getImageID();
		run("Grays");
			
		// Tidy up
		selectImage(reorderedImg);
		close();
		selectImage(tempExtrImg);
		close();
		
		selectImage(extractedChannelsImg); 
		// Output image is the sum of selected channels 
		
	}
	
	
	
	function calculateGP(ch1Img, ch2Img) {
	// Calculates GP image from ch1 and ch2 after thresholding
		
		imageCalculator("Subtract create 32-bit stack", ch1Img, ch2Img);
		numeratorImg = getImageID();
		rename("Numerator");
		
		imageCalculator("Add create 32-bit stack", ch1Img, ch2Img);
		denominatorImg = getImageID();
		rename("Denominator");
	        
		imageCalculator("Divide create 32-bit stack", numeratorImg, denominatorImg);
		gpImg = getImageID();
	
		// Tidy up
		selectImage(numeratorImg);
		close();
		selectImage(denominatorImg);
		close();
		
		selectImage(gpImg); 
		// Output image is GP image	
	
	}
	      
	
	
	function userThreshold(thrRefImg, averaging, radius) {
	// Create a mask (thrMaskImg) based on user set threshold on a reference image, after averaging (optional)
		
		selectImage(thrRefImg); 
		
		getDimensions(thrWidth, thrHeight, thrChannels, thrSlices, thrFrames);
		
		if(averaging==1){run("Mean...", "radius="+radius+" stack");} // Pixel averaging (conditional to user input)
		
		setBatchMode("show");
		if (thrFrames > 1) {Stack.setFrame(1);} // Select first frame if there are several frames
		if (thrSlices > 1) {Stack.setSlice(floor(thrSlices/2));} // Select slice in the middle of the stack
		run("Threshold...");
		waitForUser("Select threshold",
		            "Press OK after selecting a threshold value in the 'Threshold' window.\n"
		            + "Make sure to select a value significantly above the background level of your images.\n"
		            + " \n"
		            + "When processing a movie or a z-stack, select an adequate reference time point and \n"
		            + "slice, and check that the threshold is appropriate for all other frames and slices.\n"
	   	            + " \n"
		            + "Do not press Apply!");
	
		getThreshold(lower, upper); // Get the threshold values selected by the user
		close("Threshold");
		
		// Make pixel values outside the threshold range NaN
		selectImage(thrRefImg);
		setBatchMode("hide");
		
		// Create mask to store threshold values
		newImage("Threshold Mask", "32-bit grayscale-mode", thrWidth, thrHeight, thrChannels, thrSlices, thrFrames);
		thrMaskImg = getImageID();
		
		for (z = 1; z <= thrSlices; z++) {
			for (f = 1; f <= thrFrames; f++) {
				// Create threshold mask
				generateThresholdMask(thrRefImg, thrMaskImg, lower, upper, f, z);
			}
		}
		
	
		selectImage(thrMaskImg);	
		// Output image is threshold mask
		
		thresholdArray = newArray(lower, upper);
		return thresholdArray;
		
	}
	
	
	
	function generateThresholdMask(refImg, thrMaskImg, lowerThr, upperThr, thrFrame, thrSlice) {
	// thrMaskImg will show pixels within the threshold as 1, and those outside as 0 
	
		roiManager("reset");
		
		selectImage(refImg);
		getDimensions(refWidth, refHeight, refChannels, refSlices, refFrames);
		
		if (refFrames > 1) {Stack.setFrame(thrFrame);}
		if (refSlices > 1) {Stack.setSlice(thrSlice);}
		setThreshold(lowerThr, upperThr);
		run("Create Selection");
		selectionExists = getValue("selection.size");
		if (selectionExists != 0) {
			roiManager("Add");
			
			selectImage(thrMaskImg); // Make the thresholded pixels NaN in the target image
			// For some reason this has to be done before setting frame and slice. If not ROI Manager will change the frame (ImageJ bug?)
			roiManager("Select", 0); // Selects the first item in the ROI manager
			if (refFrames > 1) {Stack.setFrame(thrFrame);}
			if (refSlices > 1) {Stack.setSlice(thrSlice);}
			run("Set...", "value=1");
		}
		
		run("Select None");
		
		roiManager("reset");
		
	}
	
	
	
	function maskGP(targetImg, maskImg) { 
	// Apply a mask, e.g. threshold or saturation, to the an image to make pixels NaN
		
		selectImage(maskImg);
		setThreshold(1, 1e30);
		run("NaN Background", "stack");
		
		imageCalculator("Divide 32-bit stack", targetImg, maskImg); // Dividing by NaN results is NaN pixels
		
		selectImage(targetImg);
	}
	
	
	
	function calculateSaturatedPixelPerc(thrsMaskImg, satWithinThrMaskImg, maskFrame, maskSlice, warningPerc, warnArray) { 
	// It calculates the percentage of saturated pixels within thresholded pixels Using a saturated-within-threshold mask 
	// and a threshold mask as a reference, 
		
		selectImage(thrsMaskImg);
		getDimensions(thrMaskWidth, thrMaskHeight, thrMaskChannels, thrMaskSlices, thrMaskFrames);
		if (thrMaskFrames > 1) {Stack.setFrame(maskFrame);}
		if (thrMaskSlices > 1) {Stack.setSlice(maskSlice);}
		
		// Calculate area of thresholded ROI
		selectImage(thrsMaskImg);
		setThreshold(1, 1e30);
		thrPixelNumber = getValue("Area raw limit"); // raw to get pixel number (unscaled); limit to measure within threshold
		resetThreshold();
		
		// Count saturated pixels in the thresholded ROI
		selectImage(satWithinThrMaskImg);
		if (thrMaskFrames > 1) {Stack.setFrame(maskFrame);}
		if (thrMaskSlices > 1) {Stack.setSlice(maskSlice);}
		setThreshold(1, 1e30);
		List.setMeasurements("limit");
		satWithinThrPixelNumber = List.getValue("Area");
		resetThreshold();
		toUnscaled(satWithinThrPixelNumber);
		toUnscaled(satWithinThrPixelNumber); // Run twice because it is area, not length
		
		// Calculate percentage
		satPerc = satWithinThrPixelNumber / thrPixelNumber * 100;
		
		// Warnings
		showWarning = warnArray[0];
		
		if (satPerc >=  warningPerc && showWarning == 1) {
			Dialog.create("Saturated pixel warning!");
			Dialog.addMessage("Attention!\n"
		            + "" + satPerc + "% of the pixels are saturated in slice " + z + " of frame " +f + ".\n"
		            + "These pixels will not be used for GP quantification.");
			Dialog.addCheckbox("Do not warn me again for subsequent frames and z slices", 0);
			Dialog.addCheckbox("Show pixel saturation image", warnArray[1]);
			Dialog.show();
			doNotWarn = Dialog.getCheckbox();
			warnArray[1] = Dialog.getCheckbox(); // Show pixel saturation image
			if (doNotWarn == 1) {warnArray[0] = 0;} 
		}
		
		return warnArray; // First value is warnUser, second value is showSatMask
	}
	
	
	
	function whatDateTime() {
		 getDateAndTime(year, month, dayOfWeek, day, hour, minute, second, msec);
		 TimeString ="";
		 if (day<10) {TimeString = TimeString + "0";}
		 TimeString = TimeString + day + "-";
		 if (month<10) {TimeString = TimeString + "0";}
		 TimeString = TimeString + month + "-" + year + "_";
		 if (hour<10) {TimeString = TimeString + "0";}
		 TimeString = TimeString + hour + "-";
		 if (minute<10) {TimeString = TimeString + "0";}
		 TimeString = TimeString + minute + "-";
		 if (second<10) {TimeString = TimeString + "0";}
		 TimeString = TimeString + second;
		 return TimeString;
	}
	
	
	
	function generateResultsTable(gpImg, saveResults, saveDir, saveName) { 
	// Generates and saves (optional) a table containing GP and channel statistics
		
		selectImage(gpImg);
		getDimensions(gpImgWidth, gpImgHeight, gpImgChannels, gpImgSlices, gpImgFrames);
		
		// Create table
		gpTableName = "GP_Results_" + saveName;
		Table.create(gpTableName);
			
		// Add values to table
		// If the column does not exist it will be created in the first iteration
		row = 0;
		for (f = 1; f <= gpImgFrames; f++) { // First frames
			for (z = 1; z <= gpImgSlices; z++) {
				
				// GP values
				selectImage(gpImg);
				if (gpImgFrames > 1) {Stack.setFrame(f);}
				if (gpImgSlices > 1) {Stack.setSlice(z);}
				
				if (gpImgFrames > 1) {Table.set("Frame", row, f, gpTableName);}
				if (gpImgSlices > 1) {Table.set("Slice", row, z, gpTableName);}
				Table.set("Median", row, getValue("Median"), gpTableName);
				Table.set("Mean", row, getValue("Mean"), gpTableName);
				Table.set("StdDev", row, getValue("StdDev"), gpTableName);
				
				row++;
			}
		}
		
		Table.update(); // If not the table is not properly displayed/saved
		if(saveResults==1){Table.save(saveDir + gpTableName + ".csv", gpTableName);}
		
	}
	
	
	
	function saveHistogramTable(gpImg, binNumber, saveDir, saveName) { 
	// Saves GP histogram data
	
		selectImage(gpID);
		getDimensions(gpHistowidth, gpHistoHeight, gpHistoChannels, gpHistoSlices, gpHistoFrames);
		
		histoTableName = "Histogram_" + saveName;
		Table.create(histoTableName);
		
		for (f = 1; f <= gpHistoFrames; f++) { // First frames
			for (z = 1; z <= gpHistoSlices; z++) {
				if (gpHistoFrames > 1) {Stack.setFrame(f)};
				if (gpHistoSlices > 1) {Stack.setSlice(z)};
				getHistogram(gpValues, gpCounts, binNumber, -1, 1); // gpValues and gpCounts are arrays
				
				// Values (only first column)
				if (f == 1 && z == 1) {Table.setColumn("GP Values", gpValues);} // Bin values only the first time
				
				// Counts
				// For single images (no frames or slices)
				if (gpHistoSlices == 1 && gpHistoFrames == 1){Table.setColumn("Counts", gpCounts);}
				
				// For time lapse and/or z-stacks
				else {Table.setColumn("Frame " + f + " Slice " + z, gpCounts);}
	
			}
		}	
		
		Table.save(saveDir + histoTableName + ".csv", histoTableName);
		close(histoTableName);
		
	}
	
	
	
	function loadDefaults(path, defaultSettingNumber) { 
	// Returns an array with default settings loaded from the fluidityCalculatorDefaults.csv file
	// If this file does not exist it loads the default defaults (XD)
	
		path = getDirectory("macros") + "toolsets/defaultsGP.csv";
		
		// Check if the defaults.csv file exists, loads it and count the number of settings
		defaultsExist = File.exists(defaultsPath);
		csvCorrect = 0;
		
		if (defaultsExist == 1) {
			rawDefaultsFromCsv = File.openAsRawString(defaultsPath);
			defaults = split(rawDefaultsFromCsv, "\n"); // Splits the string into an array of substrings
			csvSettingNumber = lengthOf(defaults); 
			if (csvSettingNumber == defaultSettingNumber) {csvCorrect = 1;}
		}	
		
		// If the defaults.csv file does not exist or is wrong, load the default defaults
		if (defaultsExist == 0 || csvCorrect == 0) {
			defaults = newArray(defaultSettingNumber);
			defaults[0] = 0;					// ch1Label
			defaults[1] = 0;					// ch2Label
			defaults[2] = 1;					// width
			defaults[3] = 0;					// doAveraging
			defaults[4] = 1;					// radius
			defaults[5] = "Channel 1+2";		// thresholdCh
			defaults[6] = 10;					// satWarningPerc
			defaults[7] = 0;					// showSpectrum
			defaults[8] = 1;					// saveResults
			defaults[9] = 0;					// saveAnalysisMask
			defaults[10] = 0;					// saveLocation
			defaults[11] = 500;					// firstWavelength
			defaults[12] = 8.9;					// windowWavelength
		}
		
		return defaults;
	
	}
	
	
	function showGUI(defaultsArray, wavelengthArray) { 
	// Shows GUI with default settings
	
		// Default settings
		ch1LabelDefault = defaultsArray[0];
		ch2LabelDefault = defaultsArray[1];
		widthDefault = defaultsArray[2];
		doAveragingDefault = defaultsArray[3];
		averagingRadiusDefault = defaultsArray[4];
		thresholdChDefault = defaultsArray[5];
		satWarningPercDefault = defaultsArray[6];
		showSpectrumDefault = defaultsArray[7];
		saveResultsDefault = defaultsArray[8];
		saveAnalysisMaskDefault = defaultsArray[9];
		saveLocationDefault = defaultsArray[10];
		// firstWavelength and windowWavelength are not needed here
		
		// Create dialogue
		
		numberOfChannels = lengthOf(wavelengthArray);
		
		Dialog.create("Fluidity Calculator");
		Dialog.addChoice("Channel 1", wavelengthArray, ch1LabelDefault);
		Dialog.addChoice("Channel 2", wavelengthArray, ch2LabelDefault);
		
		Dialog.addMessage("");
		
		Dialog.addMessage("Advanced Options:");
		Dialog.addMessage("");
		Dialog.addSlider("Channel Window Range", 1, floor(numberOfChannels/2), widthDefault);
		Dialog.addCheckbox("Mean Filter (only for thresholding)", doAveragingDefault);
		Dialog.addNumber("Mean Filter Radius", averagingRadiusDefault);
			
		thresholdChList = newArray("Channel 1+2", "Channel 1", "Channel 2");
		Dialog.addChoice("Thresholding Channel", thresholdChList, thresholdChDefault);
		Dialog.addNumber("Saturated Pixel Warning %", satWarningPercDefault); 
		Dialog.addMessage("");
		Dialog.addCheckbox("Calculate Emission Spectrum", showSpectrumDefault);
		
		Dialog.addMessage("");
		
		Dialog.addMessage("Save Options:");
		Dialog.addMessage("");
		Dialog.addCheckbox("Save Results", saveResultsDefault);
		Dialog.addCheckbox("Save Analysis Mask", saveAnalysisMaskDefault);
		saveLocationList = newArray("Image Directory", "Image Directory/GP_Results", "Select");
		Dialog.addChoice("Save to", saveLocationList, saveLocationDefault);
		Dialog.addCheckbox("Save as Default Settings", 1);
		Dialog.addHelp("https://github.com/pcarravilla/fluiditycalculator/wiki/Documentation");
		
		Dialog.show();
	
	}
	
	
	
	function saveAsDefaults(path, newDefaultsArray) {
		// Save settings as default for next use
		// Pass array
		file = File.open(path);
		File.saveString(newDefaultsArray[0] + "\n" +	// ch1Label
						newDefaultsArray[1] + "\n" +	// ch2Label
						newDefaultsArray[2] + "\n" +	// width
						newDefaultsArray[3] + "\n" +	// doAveraging
						newDefaultsArray[4] + "\n" +	// radius
						newDefaultsArray[5] + "\n" +	// thresholdCh
						newDefaultsArray[6] + "\n" +	// satWarningPerc
						newDefaultsArray[7] + "\n" +	// showSpectrum
						newDefaultsArray[8] + "\n" +	// saveResults
						newDefaultsArray[9] + "\n" +	// saveAnalysisMask
						newDefaultsArray[10] + "\n" +	// saveLocation
						newDefaultsArray[11] + "\n" +	// firstWavelength
						newDefaultsArray[12] + "\n",	// windowWavelength
						path);
	}
	
	
	
	function checkLabels(dataImg, firstDefault, windowDefault) {
		// Check if labels are wavelengths and if not writes them
		
		newWavelengthDefaultsArray = newArray(firstDefault, windowDefault);
		selectImage(dataImg);
		getDimensions(width, height, channels, slices, frames);
		
		Stack.setChannel(2); // Check the second channel, as the first one might be a reference signal
		label = getInfo("slice.label");
		numberPattern = "\\d+(\\.\\d+)?"; // Regular expression to match numbers
		
		// If channel labels are not wavelengths, relabel
		
		if (!matches(label, numberPattern)) {
		    Dialog.create("Add Emission Wavelengths");
			Dialog.addMessage("Update your channel labels. \n" +
								"Attention! All channels will be assigned a wavelength. \n" +
								"Make sure all channels are spectral channels.");
			Dialog.addNumber("First Channel Wavelength (nm)", firstDefault);
			Dialog.addNumber("Channel Wavelegth Range (nm)", windowDefault);
			Dialog.addCheckbox("Do not relabel channels", 0);
			Dialog.show();
			
			firstWavelength = Dialog.getNumber();
			windowWavelength = Dialog.getNumber();
			doNotRelabel = Dialog.getCheckbox();
			
			if (doNotRelabel != 1) {
				selectImage(dataImg);
				for (c = 1; c <= channels; c++) {
					Stack.setChannel(c);
					run("Set Label...", "label=" + firstWavelength + windowWavelength * (c-1));
				}
				
				// Store values to save as defaults	
				newWavelengthDefaultsArray[0] = firstWavelength;
				newWavelengthDefaultsArray[1] = windowWavelength;
			}
		}
		
		return newWavelengthDefaultsArray;
		
	}
	
	
	
	function getChannelLabels(dataImg) { 
	// Returns the channel labels
	
		selectImage(dataImg);
		getDimensions(width, height, channels, slices, frames);
		
		wavelengthArray = newArray(channels);
		
		// Get channel labels
		for (c = 1; c <= channels; c++) {
				Stack.setChannel(c);
				label = getInfo("slice.label");
				// Relabel if labels are empty
				if (lengthOf(label) == 0){
					label = "Channel " + c;
				}
				wavelengthArray[c - 1] = label;
		}
	
		return wavelengthArray;
		
	}
	
	
	
	function defineChannelRange(dataImg, wavelengthArray, ch1Label, ch2Label, channelWidth) {
	// This will take channel labels and width (number of channels to combine and define the range for Ch1 and Ch2
	// It will give an error if the channel ranges overlap, or if the indexes are out of range
	
		getDimensions(width, height, channels, slices, frames);
		
		// Find the index of the selected channel labels
		// Note that indexes are 0-based but number of channels are 1-based
		// Channel 1
		ch1Index = -1;
		for (idx = 0; idx <= channels - 1; idx++) {
			if (wavelengthArray[idx] == ch1Label) {
				ch1Index = idx + 1; // Channel numbers are handled with 1-based indexes
				break;
			}
		}
		
		// Channel 2
		ch2Index = -1;
		for (idx = 0; idx <= channels - 1; idx++) {
			if (wavelengthArray[idx] == ch2Label) {
				ch2Index = idx + 1; // Channel numbers are handled with 1-based indexes
				break;
			}
		}
		
		// Calculate low and high wavelength channels according to the width
		// If width is even, take one more channel for chHigh
		ch1LowIndex = ch1Index - floor((channelWidth - 1) / 2); 
		ch1HighIndex = ch1Index + floor(channelWidth / 2);
		ch2LowIndex = ch2Index - floor((channelWidth - 1) / 2);
		ch2HighIndex = ch2Index + floor(channelWidth / 2);
		
		// Errors
		// Ch1 channel number is lower
		if (ch1Index < ch2Index && ch1HighIndex >= ch2LowIndex) {
			exit("Channel 1 and Channel 2 ranges overlap! Set a lower 'Channel Window' value.");
		}
		// Ch2 channel number is lower
		else if (ch2Index < ch1Index && ch2HighIndex >= ch1LowIndex) {
			exit("Channel 1 and Channel 2 ranges overlap! Set a lower 'Channel Window' value.");
		}
		
		if (ch1LowIndex < 1 || ch1HighIndex > channels) {exit("Channel 1 outside of the user-defined range! Set a lower 'Channel Window Range' value.");}
		if (ch2LowIndex < 1 || ch2HighIndex > channels) {exit("Channel 2 outside of the user-defined range! Set a lower 'Channel Window Range' value.");}
	 	
	 	channelRangeIndexes = newArray(ch1LowIndex, ch1HighIndex, ch2LowIndex, ch2HighIndex);
	 	
		return channelRangeIndexes;
	
	}
	
	
	
	function saveParametersTable (parametersArray, thresholdValues, fileName, analysisTime, saveDir, saveName) {
	// Saves user parameters (reproducibility)
		
	
		parameters = newArray(
						"Image name", 
						"Date/Time", 
						"Channel 1", 
						"Channel 2", 
						"Channel Window Range", 
						"Mean filter", 
						"Mean filter radius", 
						"Thresholding channel", 
						"Lower threshold", 
						"Upper threshold");
		
		doAveraging = parametersArray[3];
		if (doAveraging == 1) {
			averagingPerformed = "Yes";
			filterRadius = parametersArray[4];
		}
		
		else if (doAveraging == 0) {
			averagingPerformed = "No";
			filterRadius = "–";
		}
	
		values = newArray(
						fileName,
						analysisTime,
						parametersArray[0],								// ch1LabelUser
						parametersArray[1],								// ch2LabelUser
						parametersArray[2],								// widthUser
						averagingPerformed,								// doAveraging
						filterRadius,									// radius
						parametersArray[5],								// thresholdChUser
						thresholdValues [0],							// lower	
						thresholdValues [1]);							// upper
		
		// Create and save table
		parametersTableName = "Parameters_" + saveName;
		Table.create (parametersTableName); 
		Table.setColumn("Parameters", parameters);
		Table.setColumn("Values", values);
		Table.save(saveDir + parametersTableName + ".csv", parametersTableName);
		//Table.setLocationAndSize(700, 500, 600, 275);
		close(parametersTableName);
	
	}
	
	
	
	function showSpectrum(dataImg, thrRefImg, analysisMaskImg, saveName) {
	// Shows and saves the spectrum
	// dataImg is where the values are taken from
	// thrRefImg is used to let the user select a frame/slice
	// analysisMaskImg is used to discard thresholded and saturated pixels
	
		run("Select None");
		roiManager("reset");
	
		selectImage(dataImg);
		getDimensions(dataImgWidth, dataImgHeight, dataImgChannels, dataImgSlices, dataImgFrames);
	
		// If multiple frames/slices, let the user select the one to be analysed
		if (dataImgFrames > 1 || dataImgSlices > 1) {
			selectImage(thrRefImg);
			run("Select None");
			resetThreshold;
			rename("Show Spectrum – Select Frame/Slice");
			if(dataImgSlices > 1) {Stack.setSlice(1);} 
			if(dataImgFrames > 1) {Stack.setFrame(1);}
			setBatchMode("show");
			waitForUser("Select Slice/Frame, 'Show Spectrum' can only show a single frame or slice. \n" +
						"Please select a time frame and z slice to display its spectrum. \n" +
						"Data saved in the Spectrum csv file will contain the spectrum of all frames and z slices.");
			Stack.getPosition(userChannel, userSlice, userFrame);
			setBatchMode("hide");
			print("Showing emission spectrum of slice " + userSlice + "/frame " + userFrame);
		}
		
		
		// Find areas within threshold and not saturated in analysisMaskID
		wavelengthsArray = getChannelLabels(dataID);
		intensitiesArray = newArray(dataImgChannels);
		
		selectImage(analysisMaskID);
		if(dataImgSlices > 1) {Stack.setSlice(userSlice);}
		if(dataImgFrames > 1) {Stack.setFrame(userFrame);}
		setThreshold(1, 1e30); // Select all pixels in the mask
		run("Create Selection");
		selectionExists = getValue("selection.size"); // To avoid no selection errors
		
		
		// Measure mean intensity values in all channels
		if (selectionExists != 0) {
			if (dataImgFrames > 1 || dataImgSlices > 1) {
				Roi.setPosition(1, userSlice, userFrame); // This is super important! Set the Slice and Frame we want to check to the Roi that is added
			}
			roiManager("Add");
			resetThreshold; // In analysisMaskID
			run("Select None"); // In analysisMaskID
			selectImage(dataID);
			for (c = 1; c <= dataImgChannels; c++) {
				roiManager("select", 0);
				Stack.setChannel(c); // After selecting the ROI
				if(dataImgFrames > 1) {Stack.setSlice(userFrame);}
				if(dataImgSlices > 1) {Stack.setSlice(userSlice);}
				
				intensitiesArray[c-1] = getValue("Mean");
			}
		}
		
		else {waitForUser("Emission spectrum cannot be shown.", "No pixels within user-defined threshold values.");}
		
		
		resetThreshold;
		run("Select None");
		roiManager("reset"); 
		
		
		// Plot spectrum
		Plot.create("Emission spectrum_" + saveName, "Wavelength (nm)", "Intensity (mean pixel value)", wavelengthsArray, intensitiesArray);
		Plot.setLineWidth(2);
		Plot.freeze(true);
		Plot.show();
	    
	}
	
	
	
	function saveSpectrum(dataImg, thrRefImg, analysisMaskImg, saveName, saveDir) {
	// Shows and saves the spectrum
	// dataImg is where the values are taken from
	// thrRefImg is used to let the user select a frame/slice
	// analysisMaskImg is used to discard thresholded and saturated pixels
	
		run("Select None");
		roiManager("reset");
	
		selectImage(dataImg);
		getDimensions(dataImgWidth, dataImgHeight, dataImgChannels, dataImgSlices, dataImgFrames);
	
		
		// Find areas within threshold and not saturated in analysisMaskID
		wavelengthsArray = getChannelLabels(dataID);
		
		// Create table
		spectrumTableName = "Spectrum_" + saveName;
		Table.create(spectrumTableName);
		Table.setColumn("Channel", wavelengthsArray);
		
		// Get mean intensity values
		for (z = 1; z <= dataSlices; z++) {
			for (f = 1; f <= dataFrames; f++) {	
				
				intensitiesArray = newArray(dataImgChannels);			
				
				selectImage(analysisMaskID);
				if(dataImgSlices > 1) {Stack.setSlice(z);}
				if(dataImgFrames > 1) {Stack.setFrame(f);}
				setThreshold(1, 1e30); // Select all pixels in the mask
				run("Create Selection");
				selectionExists = getValue("selection.size"); // To avoid no selection errors
				
				
				// Measure mean intensity values in all channels
				if (selectionExists != 0) {
					if (dataImgFrames > 1 || dataImgSlices > 1) {
						Roi.setPosition(1, z, f); // This is super important! Set the Slice and Frame we want to check to the Roi that is added
					}
					roiManager("Add");
					resetThreshold; // In analysisMaskID
					run("Select None"); // In analysisMaskID
					selectImage(dataID);
					for (c = 1; c <= dataImgChannels; c++) {
						roiManager("select", 0);
						Stack.setChannel(c); // After selecting the ROI
						if(dataImgFrames > 1) {Stack.setSlice(f);}
						if(dataImgSlices > 1) {Stack.setSlice(z);}
						
						intensitiesArray[c-1] = getValue("Mean");
					}
				}
				
				// If there is no selection
				else {
					for (c = 1; c <= dataImgChannels; c++) {
						intensitiesArray[c-1] = "–";
					}
				}
				
				if (dataImgSlices == 1 && dataImgFrames == 1){Table.setColumn("Counts", intensitiesArray);}
				
				// For time lapse and/or z-stacks
				else {Table.setColumn("Frame " + f + " Slice " + z, intensitiesArray);}
					
				
				resetThreshold;
				run("Select None");
				roiManager("reset"); 
			}
		}
		
		Table.save(saveDir + spectrumTableName + ".csv", spectrumTableName);
		close(spectrumTableName);
		
		
	}
	
	

}