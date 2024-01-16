//Fixed parameters:

viralChannel = 2;
nucleiChannel = 1;

//Infected cell thresholding style
viralThresholdMethod = "Triangle";

//Preferred thresholding method for the nuclei
nucleiThresholdMethod = "Otsu";

setBatchMode(true);


backgroundSubtractionBallPixelRadius = 100;


#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".dv") suffix

//How circular can a cell before being considered dead?
#@ Double (label="Maximum allowed circularity", style="slider", min=0, max=1.0, stepSize=0.01, value=0.75) maxCircularity

//Nucleus size parameters
#@ Integer (label="Minimum nucleus size (micron^2)", style="slider", min=0, max=200, value=40) nucleiMin
#@ Integer (label="Maximum nucleus size (micron^2)", style="slider", min=0, max=2000, value=1000) nucleiMax

//Cell size parameters
#@ Integer (label="Minimum cell size (micron^2)", style="slider", min=0, max=500, value=200) cellMin
#@ Integer (label="Maximum cell size (micron^2)", style="slider", min=0, max=10000, value=10000) cellMax


processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	IJ.renameResults("Summary","allResults");
	 selectWindow("allResults");
	 Table.deleteColumn("Total Area");
	 Table.deleteColumn("Average Size");
	 Table.deleteColumn("%Area");
	 Table.deleteColumn("Mean");
	 Table.deleteColumn("IntDen");
	 Table.deleteColumn("Median");
	 
	 Table.renameColumn("Slice", "Infected_cell");
	 Table.renameColumn("Count", "Nuclei_count");
	 
	 saveResults("allResults", "Results_infected", "csv");
	
}

function processFile(input, output, file) {
	run("Bio-Formats Macro Extensions");
	Ext.setId(input + File.separator + file);
	Ext.getSeriesCount(seriesCount);
	for (i = 1; i <= seriesCount; i++)
	{
		processSeries(input,output,file, i);
	}
}

function processSeries(input,output,file, series)
{
	importString = "open=[" + input + File.separator + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+series;
	run("Bio-Formats Importer", importString);
	
	currentImage = getTitle();
	
	numInfectedCells = segmentInfectedCells(currentImage);
	nucleiSegmentation = segmentNuclei(currentImage);
	
	countNucleiInInfectedCells(nucleiSegmentation, numInfectedCells);
	

	roiManager("reset");
	
	close(currentImage);
	close(nucleiSegmentation);
}


function segmentNuclei(imageTitle)
{
	duplicateChannelAndProjectStack(imageTitle, nucleiChannel,"Median");
	
	nucleiSegmentation = getTitle();
	
	//Subtract the background
	backgroundSubtractionString = "rolling=" + backgroundSubtractionBallPixelRadius;
	run("Subtract Background...", backgroundSubtractionString);
	thresholdString = nucleiThresholdMethod + " dark no-reset";
	
	setAutoThreshold(thresholdString);
	
	run("Convert to Mask");

	run("Close-");
	
	run("Watershed");
	
	return nucleiSegmentation;
}

function countNucleiInInfectedCells(nuclei, numCells) {
	
	selectWindow(nuclei);
	
	for(i = 0; i < numCells; i++)
	{
		roiManager("Select",i)
		
		anaPartString = "size=[" + nucleiMin + " - " + nucleiMax + "] show=Nothing summarize";
		run("Analyze Particles...", anaPartString);
		
	}

}

function segmentInfectedCells(imageTitle)
{
	duplicateChannelAndFocus(imageTitle,viralChannel);
	
	viralDuplicate = getTitle();
	
	//Subtract the background
	backgroundSubtractionString = "rolling=" + backgroundSubtractionBallPixelRadius;
	run("Subtract Background...", backgroundSubtractionString);
	
	
	thresholdString = viralThresholdMethod + " dark no-reset";
	
	setAutoThreshold(thresholdString);
	
	run("Convert to Mask");
	cellMin = 198;
	cellMax = 10000;
	
	anaPartString = "size=[" + cellMin + " - " + cellMax + "] circularity=[0.00-" + maxCircularity + "] exclude add";
	run("Analyze Particles...", anaPartString);
	numCells=roiManager("count");
	close(viralDuplicate);
	return numCells;
}

function duplicateChannel(imageTitle, channel)
{
	selectWindow(imageTitle);
	duplicationString = "duplicate channels=" + channel;
	run("Duplicate...", duplicationString);
}

//Duplicate single Z-slice from an image to make new image
function duplicateSlice(imageTitle, slice)
{
	selectWindow(imageTitle);
	duplicationString = "slices=" + slice;
	run("Make Substack...", duplicationString);

}


//Extract a single channel from a multi-channel image and then select most in focus Z-slice
function duplicateChannelAndFocus(imageTitle, channel)
{
	//Duplicate specific channel from image to extract it as a new image
	duplicateChannel(imageTitle, channel);
	getDimensions(width, height, channels, slices, frames);
	
	//Measure standard deviations to select most in focus Z-slice and duplicate this
	//as new image
	if (slices > 1)
	{
		stackTitle = getTitle();
		
		//Highest standard deviation observed so far
		maxStDev = 0;
		
		//Initial slice assumed to be most in focus before beginning loop
		focussedSlice = 1;
		
		//Cycle through slices
		for (i = 1; i <= slices; i++)
		{
			setSlice(i);
			stDev = getValue("StdDev");
			
			//If this slice's standard deviation is higher than the current
			//maximum, then update the maximum and the most focussed slice
			if (stDev > maxStDev)
			{
				maxStDev = stDev;
				focussedSlice = i;
			}
		}
		
		//Duplicate the most focussed slice
		duplicateSlice(stackTitle, focussedSlice);
		focussedTitle = getTitle();
		
		close(stackTitle);
		
		selectWindow(focussedTitle);
	}
}


//Extract a single channel from a multi-channel image and then do a Z project
function duplicateChannelAndProjectStack(imageTitle, channel, projectPreference)
{
	//Duplicate specific channel from image to extract it as a new image
	duplicateChannel(imageTitle, channel);
	
	//Do Z project to get single frame if the image is a Z-stack
	getDimensions(width, height, channels, slices, frames);
	if (slices > 1)
	{
		projectString = "projection=[" + projectPreference + "]";
		stackTitle = getTitle();
		run("Z Project...", projectString);
		zProjectTitle = getTitle();
		close(stackTitle);
		
		//Select the new image before finishing
		selectWindow(zProjectTitle);
	}
}

//Save results as CSV, with image name, appended with "_foci_per_nucleus"
function saveResults(table, name, extension)
{
	outCSVName = name + "." + extension;
	fullPath = output + "/" + outCSVName;
	
	extendNumber = 1;
	
	while(File.exists(fullPath))
	{
		outCSVName = name + extendNumber + "." + extension;
		fullPath = output + "/" + outCSVName;
		extendNumber ++;
	}
	
	selectWindow(table);
	saveAs(table, fullPath);
	close(outCSVName);
}