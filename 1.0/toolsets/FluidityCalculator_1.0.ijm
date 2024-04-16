// Look for PCNotes *************
    
macro "Fluidity Calculator Action Tool - icon:fluidity.png" {
	
// 1. PREPARATION

	setBatchMode(true); //Start batch mode
	
	imageList = getList("image.titles");
	if (imageList.length == 0) {exit("Error: No images open for analysis. Please open an image file before proceeding.");}
	// Get image name and time
	name = getTitle();
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	date = "_" + dayOfMonth + "_" + month + "_" + year;
	time = "_" + hour + "-" + minute + "-" + second;

	// Clear previous log, ROIs and selections
	close("Results");
	print("\\Clear");	
	run("Select None");
	roiManager("reset");
	
	// Duplicate original to avoid modifications and get dimensions
	run("Duplicate...", "title=dat duplicate");
	getDimensions(dataWidth, dataHeight, dataChannels, dataSlices, dataFrames);
	
	// Stop if the image does not have at least two channels
	if (dataChannels < 2) {exit("The image must have at least 2 channels.");}
	
	// Change if removing saturated pixels is an option
	removeSaturated = 1; 

	dataID = getImageID();

	// Store wavelength list
	wavelengthList = newArray(dataChannels);
	nChannels = wavelengthList.length;
	for (c = 1; c <= dataChannels; c++) {
		Stack.setChannel(c);
		label = getInfo("slice.label");
		wavelengthList[c - 1] = label;
	}

	// Load default settings
	defaultsPath = getDirectory("macros") + 'toolsets/defaults.csv';
	defaultsExist = File.exists(defaultsPath);
	// Check if the file contains the exact number of settings
	defaultsCorrect=0;
	if (defaultsExist==1) {
		csvDefaults = File.openAsRawString(defaultsPath);
		lines = split(csvDefaults, "\n");
		numLines = lengthOf(lines);
		if (numLines == 10) {defaultsCorrect=1;}
	}
	// If the file is corrrect, extract settings 
	if (defaultsExist==1 && defaultsCorrect==1) {
		csvDefaults = File.openAsRawString(defaultsPath);
		lines = split(csvDefaults, "\n");
		numLines = lengthOf(lines);
		defaults = split(csvDefaults, "\n");	
		ch1Default = parseFloat(defaults[0]);
		ch2Default = parseFloat(defaults[1]);
		widthDefault = parseFloat(defaults[2]);
		averagingDefault = parseFloat(defaults[3]);
		radiusDefault = parseFloat(defaults[4]);
		thresholdingChannelDefault = defaults[5];
		saturatedWarningDefault = parseFloat(defaults[6]);
		showSpectrumDefault = parseFloat(defaults[7]);
		showTablesDefault = parseFloat(defaults[8]);
		saveResultsDefault = parseFloat(defaults[9]);
		saveSettingsDefault = 0;
	}
	// If not, these are the default settings
	else {
		ch1Default = 0; // 1
		ch2Default = 0;
		widthDefault = 1;
		averagingDefault = 0;
		radiusDefault = 1; // 5
		thresholdingChannelDefault = "Channel 1+2";
		saturatedWarningDefault = 5;
		showSpectrumDefault = 0;
		showTablesDefault = 0; 
		saveResultsDefault = 0; //10
		saveSettingsDefault = 1;
		
	}

// 2. USER INPUT
	
	// Create dialogue
	Dialog.create("Fluidity Calculator");
	Dialog.addChoice("Channel 1", wavelengthList, ch1Default);
	Dialog.addChoice("Channel 2", wavelengthList, ch2Default);
	Dialog.addSlider("Number of Channels", 1, dataChannels, widthDefault);
	Dialog.addCheckbox("Mean Filter", averagingDefault);
	Dialog.addNumber("Mean Filter Radius", radiusDefault);
	thresholdingChannelList = newArray("Channel 1+2", "Channel 1", "Channel 2");
	Dialog.addChoice("Thresholding Channel", thresholdingChannelList, thresholdingChannelDefault);
	Dialog.addNumber("Saturated Pixel Warning %", saturatedWarningDefault); 
	Dialog.addMessage("");
	Dialog.addCheckbox("Show Emission Spectrum", showSpectrumDefault);
	Dialog.addCheckbox("Show Tables", showTablesDefault);
	Dialog.addMessage("");
	Dialog.addCheckbox("Save Results", saveResultsDefault);
	Dialog.addCheckbox("Save as Default Settings", saveSettingsDefault);
	Dialog.addHelp('https://github.com/pcarravilla');
	
	Dialog.show();

	// Get user input	
	ch1Label = Dialog.getChoice();
	ch2Label = Dialog.getChoice();
	width = Dialog.getNumber();
	averaging = Dialog.getCheckbox();
	radius = Dialog.getNumber();
	thresholdingChannel = Dialog.getChoice();
	saturatedWarning = Dialog.getNumber();
	showSpectrum = Dialog.getCheckbox();
	showTables = Dialog.getCheckbox();
	saveResults = Dialog.getCheckbox();
	saveAsDefault = Dialog.getCheckbox();
	
	if (saveResults == 1) {saveDir = getDirectory("Select a directory");}
	
	// Save settings as default for next use
	if (saveAsDefault == 1) {
		file = File.open(defaultsPath);
		File.saveString(ch1Label + '\n' +
						ch2Label + '\n' +
						width + '\n' +
						averaging + '\n' +
						radius + '\n' +
						thresholdingChannel + '\n' +
						saturatedWarning + '\n' +
						showSpectrum + '\n' +
						showTables + '\n' +
						saveResults
						, defaultsPath);
	}

	// Find the index of the selected channel labels
	// Channel 1
	ch1Index = -1;
	for (i = 0; i < nChannels; i++) {
		if (wavelengthList[i] == ch1Label) {
			ch1Index = i+1; //channels are 1-based
			break;
		}
	}
	// Channel 2
	ch2Index = -1;
	for (i = 0; i < nChannels; i++) {
		if (wavelengthList[i] == ch2Label) {
			ch2Index = i+1; //channels are 1-based
			break;
		}
	}
	
	// Calculate low and high wavelength channels according to the width
	ch1Low = ch1Index - Math.floor((width - 1) / 2);
	ch1High = ch1Index + Math.floor(width / 2);
	ch2Low = ch2Index - Math.floor((width - 1) / 2);
	ch2High = ch2Index + Math.floor(width / 2);
	
	// Exit if channels 1 and 2 overlap
	if (ch1High >= ch2Low) {exit("Channel 1 and Channel 2 overlap. Set lower Number of Channels");}

// 3. GP CALCULATION (PART 1)

	// Extract Lo and Ld channels from user-defined channels
	selectImage(dataID);
	run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
	reorderedID = getImageID();
	run("Z Project...", "start=ch1Low stop=ch1High projection=[Sum Slices] all");
	rename("ch1");
	ch1ID = getImageID(); 
	selectImage(reorderedID);
	run("Z Project...", "start=ch2Low stop=ch2High projection=[Sum Slices] all");
	rename("ch2");
	ch2ID = getImageID();
	selectImage(reorderedID);
	run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
	dataID = getImageID();
	
	// Calculate numerator (ch1 - ch2) and denominator (ch1 + ch2)

	// Numerator
	imageCalculator("Subtract create 32-bit stack", ch1ID, ch2ID);
	numeratorID = getImageID();
	rename("Numerator");

	// Denominator
	imageCalculator("Add create 32-bit stack", ch1ID, ch2ID);
	denominatorID = getImageID();
	rename("Denominator");


// 4. THRESHOLDING

	if (thresholdingChannel == "Channel 1+2") {thresholdDuplicateID = denominatorID;}
	else if (thresholdingChannel == "Channel 1") {thresholdDuplicateID = ch1ID;}
	else if (thresholdingChannel == "Channel 2") {thresholdDuplicateID = ch2ID;}

	// Threshold image based on the selected channel
	selectImage(thresholdDuplicateID);
	run("Duplicate...", "title=Thresholding duplicate");
	thresholdChannelID = getImageID();
	getDimensions(thresWidth, thresHeight, thresChannels, thresSlices, thresFrames);
	// Pixel averaging (conditional to user input)
	if(averaging==1){run("Mean...", "radius="+radius+" stack");}
	setBatchMode("show");
	if (thresFrames > 1) {Stack.setFrame(1);} // Select first frame if there are several frames
	run("Threshold...");
	waitForUser("Select threshold",
	            "Press OK after selecting a threshold.\n"
	            + "Make sure to select a value significantly above the noise level of your images.\n"
	            + "If data is a movie, check that the threshold is appropriate for all frames.\n"
	            + "Do not press Apply!");

	// Get the threshold values selected by the user
	getThreshold(lower, upper);
	close("Threshold");

	// Make pixel values outside the threshold range NaN
	selectImage(thresholdChannelID);
	run("Set Scale...", "distance=1 known=1 unit=pixel"); // To measure number of pixels
	setBatchMode("hide");
	for (f = 1; f <= thresFrames; f++) {
		if (thresFrames > 1) {Stack.setFrame(f)};
		setThreshold(lower, upper);
		run("Create Selection");
		run("Make Inverse");
		roiManager("Add");
		// Make the thresholded pixels NaN in the denominator to avoid having them in the GP calculation
		selectImage(denominatorID);
		if (thresFrames > 1) {Stack.setFrame(f)};
		roiManager("Select", 0); // Selects the first item in the ROI manager
		run("Set...", "value=NaN");
		roiManager("reset");
		selectImage(thresholdChannelID);
	}


// 5. REMOVE AND COUNT SATURATED PIXELS

	// Creates an array that stores the percentage of saturated pixels for each frame
	saturationArray = newArray(dataFrames); 
	
	if (removeSaturated==1) {

		// Check image depth and store maximum possible pixel value
		selectImage(dataID);
		maxPossible=Math.pow(2, bitDepth())-1;
		
		warnUser = 1; // To stop warning the user about saturated pixels
		
		// Iterate over the channels and NaN saturated pixels in all channels
		for (f = 1; f <= dataFrames; f++) {

			// NaN saturated pixels in this frame
			if (dataFrames > 1) {Stack.setFrame(f)};
			roiCounter = 0;
			for (c = 1; c <= dataChannels; c++) {
				if ((c >= ch1Low && c <= ch1High) || (c >= ch2Low && c <= ch2High)){
					Stack.setChannel(c);
					getStatistics(area, mean, min, max, std, histogram);
					if (maxPossible==max) { // Checks if any pixel is saturated
						setThreshold(max, max);
						run("Create Selection");
						roiManager("add");
						selectImage(denominatorID);
						if (dataFrames > 1) {Stack.setFrame(f)};
						roiManager("Select", roiCounter);
						run("Set...", "value=NaN stack");
						run("Select None");
						selectImage(dataID);
						roiCounter++;
					}
				}
			}
	
			// Check number of saturated pixels for this frame
			// This is done by combining the ROIs generated in the previous step
			if(RoiManager.size>0){
				// Combine the saturated pixel region of all areas
				selectImage(dataID);
				run("Set Scale...", "distance=1 known=1 unit=pixel");
				if (dataFrames > 1) {Stack.setFrame(f)};
				roiManager("Combine"); // ROI1: Pixels saturated in Ch1 or Ch2
				run("Create Mask");
				if (dataFrames > 1){rename("Saturated pixels_frame_" + f);}
				else {rename("Saturated pixels");}
				saturatedID = getImageID();

				// Check number of pixels within threshold
				selectImage(thresholdChannelID);
				if (dataFrames > 1) {Stack.setFrame(f)};
				run("Set Scale...", "distance=1 known=1 unit=pixel");
				setThreshold(lower, upper); // Threshold set by the user before
				run("Create Selection");
				List.setMeasurements; 
				pixelsWithinThreshold = List.getValue("Area");
				run("Create Selection"); 
				roiManager("add"); // ROI2: Area within threshold
				roiNumberArray = Array.getSequence(2); // To select both ROIS
				roiManager("select", roiNumberArray);
				roiManager("AND"); //Pixels in ROI1 AND ROI2 (saturated pixels within the threshold)
				selectionExists = selectionType(); // Returns -1 if there is no selection
				saturatedPixels = 0;
				if (selectionExists != -1) {
					List.setMeasurements; 
					saturatedPixels = List.getValue("Area");
				}
				percSaturated = 100 * saturatedPixels / pixelsWithinThreshold;
				saturationArray[f-1] = percSaturated;
				// Warning if above threshold set by user
				if (percSaturated > saturatedWarning && warnUser == 1){
					// Show saturation mask
					selectImage(saturatedID);
					setBatchMode("show");
					
					// Show warning
					Dialog.create("Saturated pixel warning!");
					Dialog.addMessage("Attention!\n"
				            + "" + percSaturated + "% of the pixels are saturated in frame " +f + ".\n"
				            + "These pixels cannot be used for GP quantification.");
					Dialog.addCheckbox("Warn me again for subsequent frames.", 1);	   
					Dialog.show();
					warnUser = Dialog.getCheckbox();
				}
			}
		}
	}
	roiManager("reset");
	

// 6. GP CALCULATION (PART II)

	// Caclulate GP Image
	imageCalculator("Divide create 32-bit stack", numeratorID,denominatorID);
	rename("GP_" + name + date + time);
	gpID = getImageID();
	if(saveResults==1){saveAs("tiff", saveDir + name + time + "_gp");}


	// Export histogram
	if(saveResults==1){
		nBins = 100; //Number of bins
		Table.create("Histogram" + name + date + time);
		selectImage(gpID);
		for (f = 1; f <= dataFrames; f++) {
			if (dataFrames > 1) {Stack.setFrame(f)};
			getHistogram(values, counts, nBins);
			if (f==1) {Table.setColumn("GP", values);} //Only the first frame
			Table.setColumn("Counts"+f, counts);
		}
		saveAs("text", saveDir + "Histogram_" + name + date + time + ".csv");
		close("Histogram" + name + date + time);
	}

	// Save results to a table and export
	resultsTableName = "Results_" + name + date + time;
	Table.create(resultsTableName);
	for (f = 1; f <= dataFrames; f++) {
		if (dataFrames > 1) {Stack.setFrame(f)};
		List.setMeasurements;
		medianGP = List.getValue("Median");
		meanGP = List.getValue("Mean");
		sdGP = List.getValue("StdDev");
		Table.set("Frame", f-1, f);
		Table.set("Median", f-1, medianGP);
		Table.set("Mean", f-1, meanGP);
		Table.set("SD", f-1, sdGP);
		// Store percentage of saturated pixels if these are removed. If not, that value is not calculated
		if (removeSaturated==1){
			sat = saturationArray[f-1];
			Table.set("Saturated_Pixel_%", f-1, sat);
			}
		}
	if(saveResults==1){saveAs("text", saveDir + "Results_" + name + date + time + ".csv");}
	Table.setLocationAndSize(500, 800, 500, 200);


	// Viridis LUT
	// Check if viridis lut exists
	lutViridis = getDirectory("luts") + 'gp-viridis.lut';
	if (dataFrames > 1) {Stack.setFrame(1)};
	if (File.exists(lutViridis) == 1) {run("gp-viridis");} 
	else {run("mpl-viridis")}; // In case our custom gp-viridis LUT does not exist
	run("Enhance Contrast", "saturated=0.35");
	
	// Scale bar
	// run("Scale Bar...", "width=5 height=6 font=30 color=White background=None location=[Lower Right] bold overlay");
	
	// Plot histogram
	selectImage(gpID);
	setBatchMode("show");
	run("Histogram", "bins=100 use y_max=Auto");
	if(saveResults==1){saveAs("tiff", saveDir + name + time + "_histo");}
	setBatchMode("show");


// 7. PARAMETERS TABLE

	parametersTableName = "Parameters_" + name + date + time;
	Table.create (parametersTableName); // add name and time
	parameters = newArray("Image name", "Date/Time", "Channel 1", "Channel 2", "Number of channels", "Mean filter", "Mean filter radius", "Thresholding channel", "Lower threshold", "Upper threshold");
	values = newArray(name, date+time, ch1Label, ch2Label, width, averaging, radius, thresholdingChannel, lower, upper);
	Table.setColumn("Parameters", parameters);
	Table.setColumn("Values", values);
	if(saveResults==1){saveAs("text", saveDir + "Parameters_" + name + date + time + ".csv");}
	Table.setLocationAndSize(700, 500, 600, 275);
	

// 8. SHOW SPECTRUM
	
	close("Results");
	if (showSpectrum==1) {	
		selectImage(denominatorID);
		if (dataFrames > 1) {Stack.setFrame(1)}; // Only for first frame
		setThreshold(-1e30, 1e30); // Select all pixels except NaNs
		run("Create Selection");
		roiManager("Add"); 
		roiManager("select", 0);
		selectImage(dataID);
		roiManager("Multi Measure");
		
		// Store values in table
		selectWindow("Results");
		meanIntensity = Table.getColumn("Mean1");
	    sdIntensity = Table.getColumn("StdDev1");
	    spectrumTableName = "Spectral data_" + name + date + time;
		Table.create (spectrumTableName);
		Table.setColumn("Wavelenght", wavelengthList);
	    Table.setColumn("Mean_Intensity", meanIntensity);
		Table.setColumn("SD_Intensity", sdIntensity);
		if(saveResults==1){saveAs("text", saveDir + "Spectrum_" + name + date + time + ".csv");}
		
		// Plot
	    Plot.create("Emission spectrum_"+name+time, "Wavelength_nm", "Intensity", wavelengthList, meanIntensity);
	    Plot.setLineWidth(2);
	    Plot.freeze(true);
	    Plot.show();
	    if (showTables == 0) {close(spectrumTableName);}
	    close("Results");
	}

// 9. FINISH

	setBatchMode(false); //Stop batch mode
	close("Log");
	if(showTables == 0) {close(parametersTableName);}
	run("Tile");
