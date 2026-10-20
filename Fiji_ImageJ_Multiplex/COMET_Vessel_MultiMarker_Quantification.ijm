/*
 * COMET vessel-profile multi-marker intensity analysis
 * Version 2.0.2 (2026-10-04)
 * Requires Fiji/ImageJ 1.54a or later.
 *
 * PURPOSE
 *   1. Segment vascular profiles from one user-selected vessel-marker image
 *      (CD31 by default).
 *   2. Assign stable IDs (V0001, V0002, ...).
 *   3. Measure the vessel marker plus 1-5 additional markers in the identical
 *      ROIs (2-6 separate channel images in total).
 *   4. Preserve traceability with an ROI ZIP, binary mask, numbered map,
 *      full-image pixel coordinates and an analysis-settings record.
 *
 * BEFORE RUNNING
 *   - Open every required channel as a separate, matched, single-plane
 *     grayscale image in Fiji.
 *   - All selected images must represent exactly the same tissue crop and
 *     have identical width, height, orientation and spatial calibration.
 *   - Use native grayscale OME-TIFF data (8/16/32 bit), not RGB composites,
 *     JPEGs or screenshots.
 *   - Enter the full-resolution crop X/Y origin when prompted if full-image
 *     coordinates are needed.
 *
 * IMPORTANT ANALYTICAL NOTES
 *   - The macro intentionally pauses for threshold review and vessel QC.
 *   - Touching/branching marker-positive objects are not automatically
 *     guaranteed to be separate biological vessels. Correct the ROI list at
 *     the QC pause.
 *   - Lumen holes are not filled. Measurements therefore describe
 *     marker-positive vascular profiles in a 2D section.
 *   - Background-corrected mean = raw ROI mean - pooled mean of the shared
 *     user-selected background ROIs for that channel.
 *   - Compare a given marker across samples acquired and processed under
 *     identical conditions. Do not compare absolute intensities between
 *     different fluorophores as if they were on the same scale.
 *
 * OUTPUTS
 *   <ID>_Vessel_measurements.csv
 *   <ID>_Vessel_ROIs.zip
 *   <ID>_<vessel-marker>_binary_mask.tif
 *   <ID>_Vessel_ID_map.png
 *   <ID>_Analysis_settings.txt
 */

requires("1.54a");

if (nImages < 2)
    exit("Open at least two matched, separate grayscale channel images before running this macro.");

titles = getList("image.titles");
targetCountOptions = newArray("1", "2", "3", "4", "5");

Dialog.create("Number of target channels");
Dialog.addMessage("Choose how many additional markers will be measured inside the vessel-marker-positive ROIs.\n\nTotal selected images = 1 vessel marker + 1-5 target markers.");
Dialog.addChoice("Additional target channels", targetCountOptions, "2");
Dialog.show();

nTargets = parseInt(Dialog.getChoice());
nChannels = nTargets + 1;

if (nImages < nChannels)
    exit("You selected " + nChannels + " total channels, but only " + nImages + " image windows are open.");

markerNames = newArray("", "", "", "", "", "");
safeMarkers = newArray("", "", "", "", "", "");
imageTitles = newArray("", "", "", "", "", "");
bitDepths = newArray(0, 0, 0, 0, 0, 0);
backgroundPixels = newArray(0, 0, 0, 0, 0, 0);
backgroundSums = newArray(0, 0, 0, 0, 0, 0);
backgroundMeans = newArray(0, 0, 0, 0, 0, 0);

Dialog.create("Assign matched channel images");
Dialog.addMessage("Enter the biological marker name and select its matching image window.\nAll images must be separate grayscale exports from exactly the same crop.");
Dialog.addString("Vessel marker name", "CD31", 18);
Dialog.addChoice("Vessel marker image", titles, titles[0]);
for (k = 1; k < nChannels; k++) {
    Dialog.addString("Target " + k + " marker name", "", 18);
    Dialog.addChoice("Target " + k + " image", titles, titles[k]);
}
Dialog.show();

markerNames[0] = String.trim(Dialog.getString());
imageTitles[0] = Dialog.getChoice();
for (k = 1; k < nChannels; k++) {
    markerNames[k] = String.trim(Dialog.getString());
    imageTitles[k] = Dialog.getChoice();
}

// Validate names and prevent ambiguous channel assignments.
for (k = 0; k < nChannels; k++) {
    if (markerNames[k] == "")
        exit("Every selected channel requires a marker name.");
    safeMarkers[k] = safeLabel(markerNames[k]);
    if (safeMarkers[k] == "")
        exit("Marker name " + (k + 1) + " does not contain a usable character.");

    for (j = 0; j < k; j++) {
        if (imageTitles[k] == imageTitles[j])
            exit("Each channel must use a different image window. '" + imageTitles[k] + "' was selected more than once.");
        if (toLowerCase(safeMarkers[k]) == toLowerCase(safeMarkers[j]))
            exit("Marker names must be unique after filename/column-name cleanup. Please distinguish '" + markerNames[j] + "' and '" + markerNames[k] + "'.");
    }
}

// Validate dimensions, planes, bit depth and calibration for every channel.
calTol = 1e-9;
for (k = 0; k < nChannels; k++) {
    selectWindow(imageTitles[k]);
    getDimensions(tmpW, tmpH, tmpC, tmpZ, tmpT);
    tmpBD = bitDepth();
    getPixelSize(tmpUnit, tmpPW, tmpPH);
    bitDepths[k] = tmpBD;

    if (tmpC != 1 || tmpZ != 1 || tmpT != 1)
        exit("'" + imageTitles[k] + "' is not a single-channel, single-Z, single-time-point image.");
    if (tmpBD == 24)
        exit("'" + imageTitles[k] + "' is flattened RGB (24-bit). Use the original individual grayscale channel image.");

    if (k == 0) {
        imageWidth = tmpW;
        imageHeight = tmpH;
        lengthUnit = tmpUnit;
        pixelWidth = tmpPW;
        pixelHeight = tmpPH;
    } else {
        if (tmpW != imageWidth || tmpH != imageHeight)
            exit("The selected images do not have identical pixel dimensions. Re-import all channels using the same X, Y, width and height.");
        if (tmpUnit != lengthUnit || abs(tmpPW - pixelWidth) > calTol || abs(tmpPH - pixelHeight) > calTol)
            exit("Spatial calibration differs among the selected channels. Confirm that all files are matching exports.");
    }
}

defaultMinArea = 20;
if (indexOf(toLowerCase(lengthUnit), "pixel") >= 0)
    defaultMinArea = 200;

thresholdModes = newArray("Automatic + review", "Fixed lower threshold + review");
thresholdMethods = getList("threshold.methods");

Dialog.create("Vessel-analysis settings");
Dialog.addString("Sample/ROI ID", "Sample01_ROI01", 24);
Dialog.addNumber("Full-resolution crop origin X (pixels)", 0);
Dialog.addNumber("Full-resolution crop origin Y (pixels)", 0);
Dialog.addChoice(markerNames[0] + " threshold mode", thresholdModes, "Automatic + review");
Dialog.addChoice("Automatic threshold method", thresholdMethods, "Triangle");
Dialog.addNumber("Fixed lower threshold (raw units; ignored in automatic mode)", 0);
Dialog.addNumber("Gaussian sigma for segmentation (pixels; 0 = off)", 1.0);
Dialog.addNumber("Rolling-ball radius for segmentation (pixels; 0 = off)", 0);
Dialog.addNumber("Minimum vascular-profile area (" + lengthUnit + "^2)", defaultMinArea);
Dialog.addCheckbox("Exclude profiles touching the crop edge", true);
Dialog.addCheckbox("Pause for visual threshold review", true);
Dialog.addCheckbox("Measure background from shared user-drawn ROIs", true);
Dialog.addNumber("Number of background ROIs", 3);
Dialog.addMessage("Background correction affects the measurements only; the original channel images are not modified.");
Dialog.show();

sampleID = String.trim(Dialog.getString());
cropX = Dialog.getNumber();
cropY = Dialog.getNumber();
thresholdMode = Dialog.getChoice();
thresholdMethod = Dialog.getChoice();
fixedLower = Dialog.getNumber();
gaussianSigma = Dialog.getNumber();
rollingRadius = Dialog.getNumber();
minArea = Dialog.getNumber();
excludeEdges = Dialog.getCheckbox();
reviewThreshold = Dialog.getCheckbox();
measureBackground = Dialog.getCheckbox();
nBackgroundROIs = round(Dialog.getNumber());

if (sampleID == "")
    exit("A Sample/ROI ID is required.");
if (cropX < 0 || cropY < 0)
    exit("Crop origin X and Y cannot be negative.");
if (gaussianSigma < 0 || rollingRadius < 0 || minArea < 0)
    exit("Sigma, rolling-ball radius and minimum area cannot be negative.");
if (measureBackground && nBackgroundROIs < 1)
    exit("Use at least one background ROI, or disable interactive background measurement.");
if (thresholdMode == "Fixed lower threshold + review" && fixedLower <= 0)
    exit("Enter a fixed lower threshold greater than zero, or select automatic thresholding.");

safeID = safeLabel(sampleID);
outputDir = getDirectory("Choose output folder");
if (outputDir == "")
    exit("No output folder was selected.");

resultsPath = outputDir + safeID + "_Vessel_measurements.csv";
roiPath = outputDir + safeID + "_Vessel_ROIs.zip";
maskPath = outputDir + safeID + "_" + safeMarkers[0] + "_binary_mask.tif";
mapPath = outputDir + safeID + "_Vessel_ID_map.png";
settingsPath = outputDir + safeID + "_Analysis_settings.txt";

// Use sequential checks for compatibility with ImageJ builds that reject a
// logical-OR expression when File.exists() is assigned to a variable.
outputExists = false;
if (File.exists(resultsPath))
    outputExists = true;
if (File.exists(roiPath))
    outputExists = true;
if (File.exists(maskPath))
    outputExists = true;
if (File.exists(mapPath))
    outputExists = true;
if (File.exists(settingsPath))
    outputExists = true;
if (outputExists) {
    overwriteOK = getBoolean("One or more output files for '" + safeID + "' already exist. Overwrite them?", "Overwrite", "Cancel");
    if (!overwriteOK)
        exit("Analysis cancelled. Use a different Sample/ROI ID or output folder.");
}

// Warn before clearing an existing ROI Manager list.
existingROIs = roiManager("count");
if (existingROIs > 0) {
    clearOK = getBoolean("ROI Manager contains " + existingROIs + " ROI(s). They must be cleared for this analysis.", "Clear and continue", "Cancel");
    if (!clearOK)
        exit("Analysis cancelled without changing the existing ROI list.");
}
roiManager("reset");

// Optional shared spatial background ROIs, measured separately in every channel.
if (measureBackground) {
    for (b = 0; b < nBackgroundROIs; b++) {
        selectWindow(imageTitles[0]);
        run("Select None");
        setTool("rectangle");
        waitForUser("Select background ROI " + (b + 1) + " of " + nBackgroundROIs,
            "On the " + markerNames[0] + " image, draw a reasonably large tissue-background ROI\n" +
            "containing no vessel or marker-specific staining.\n\n" +
            "Use non-overlapping regions distributed across the crop.\n" +
            "The identical spatial ROI will be measured in every selected channel.\n" +
            "Click OK here after drawing it.");

        st = selectionType;
        if (st < 0 || (st >= 5 && st <= 8) || st == 10)
            exit("Each background selection must be an area ROI (rectangle, oval, polygon, freehand, traced or composite).");

        roiManager("add");
        bgIndex = roiManager("count") - 1;

        for (k = 0; k < nChannels; k++) {
            selectWindow(imageTitles[k]);
            roiManager("select", bgIndex);
            getRawStatistics(bgN, bgMean);
            if (bgN <= 0)
                exit("A background ROI contained no pixels in channel '" + markerNames[k] + "'.");
            backgroundPixels[k] = backgroundPixels[k] + bgN;
            backgroundSums[k] = backgroundSums[k] + bgN * bgMean;
        }
    }

    for (k = 0; k < nChannels; k++) {
        if (backgroundPixels[k] <= 0)
            exit("Background could not be calculated for channel '" + markerNames[k] + "'.");
        backgroundMeans[k] = backgroundSums[k] / backgroundPixels[k];
    }
    roiManager("reset");
}

for (k = 0; k < nChannels; k++) {
    selectWindow(imageTitles[k]);
    run("Select None");
}

// Create a processed duplicate for vessel segmentation; originals remain unchanged.
maskWorkingTitle = safeID + "_" + safeMarkers[0] + "_segmentation_working";
selectWindow(imageTitles[0]);
run("Duplicate...", "title=[" + maskWorkingTitle + "]");
selectWindow(maskWorkingTitle);

if (rollingRadius > 0)
    run("Subtract Background...", "rolling=" + rollingRadius);
if (gaussianSigma > 0)
    run("Gaussian Blur...", "sigma=" + gaussianSigma);

getRawStatistics(allN, allMean, allMin, allMax, allSD);
maskBitDepth = bitDepth();

if (thresholdMode == "Automatic + review") {
    thresholdOptions = thresholdMethod + " dark";
    if (maskBitDepth == 16)
        thresholdOptions = thresholdOptions + " 16-bit";
    setAutoThreshold(thresholdOptions);
} else {
    setThreshold(fixedLower, allMax, "raw");
}

if (reviewThreshold) {
    waitForUser("Review " + markerNames[0] + " threshold",
        "Inspect the red threshold overlay on the segmentation working image.\n" +
        "If necessary, use Image > Adjust > Threshold to refine it.\n\n" +
        "Include genuine " + markerNames[0] + "-positive vascular profiles and exclude background.\n" +
        "Click OK here only when the threshold is satisfactory.");
}

getThreshold(thresholdLow, thresholdHigh);
if (thresholdLow < 0 || thresholdHigh < 0)
    exit("No active threshold was found on the vessel-marker working image.");

saveSettings();
setOption("BlackBackground", true);
run("Convert to Mask");
restoreSettings();
rename(safeID + "_" + safeMarkers[0] + "_threshold_mask_working");
thresholdMaskID = getImageID();

roiManager("reset");
particleOptions = "size=" + minArea + "-Infinity circularity=0.00-1.00 show=Nothing add";
if (excludeEdges)
    particleOptions = particleOptions + " exclude";
run("Analyze Particles...", particleOptions);

nVessels = roiManager("count");
if (nVessels == 0)
    exit("No " + markerNames[0] + "-positive profiles were detected. Recheck the threshold and minimum-area setting.");

// Analyze Particles can add an outer particle outline to ROI Manager even
// when its size measurement excludes an interior hole. Intersect every
// candidate with the exact binary foreground mask so lumens/background holes
// are also excluded from subsequent intensity measurements.
selectImage(thresholdMaskID);
setThreshold(255, 255, "raw");
run("Create Selection");
if (selectionType < 0)
    exit("The binary foreground mask could not be converted to a selection.");
roiManager("add");
globalMaskIndex = nVessels;

for (i = 0; i < nVessels; i++) {
    pair = newArray(i, globalMaskIndex);
    roiManager("select", pair);
    roiManager("and");
    if (selectionType < 0)
        exit("A vessel candidate had an empty intersection with the binary foreground mask.");
    roiManager("add");
}

deleteIndexes = newArray(nVessels + 1);
for (i = 0; i <= nVessels; i++)
    deleteIndexes[i] = i;
roiManager("select", deleteIndexes);
roiManager("delete");

if (roiManager("count") != nVessels)
    exit("The vessel ROI count changed unexpectedly while preserving interior holes.");

// Mandatory visual QC against the original vessel-marker pixels.
selectWindow(imageTitles[0]);
roiManager("show all without labels");
waitForUser("Quality-control the vessel ROIs",
    "Inspect every ROI on the ORIGINAL " + markerNames[0] + " image. In ROI Manager:\n\n" +
    "- delete false-positive ROIs;\n" +
    "- add missed vascular profiles;\n" +
    "- manually correct fragmented or incorrectly joined profiles.\n\n" +
    "If a profile contains a lumen, keep the lumen excluded rather than drawing it as a separate vessel.\n" +
    "Click OK when the ROI list is final.");

nVessels = roiManager("count");
if (nVessels == 0)
    exit("The final ROI list is empty.");

// Confirm that every retained item is an area ROI before measuring it.
for (i = 0; i < nVessels; i++) {
    selectWindow(imageTitles[0]);
    roiManager("select", i);
    st = selectionType;
    if (st < 0 || (st >= 5 && st <= 8) || st == 10)
        exit("ROI " + (i + 1) + " is not an area ROI. Correct or delete it in ROI Manager, then rerun the macro.");
}

// Assign stable IDs after manual QC.
for (i = 0; i < nVessels; i++) {
    roiManager("select", i);
    roiManager("rename", "V" + IJ.pad(i + 1, 4));
}
roiManager("remove slice info");
RoiManager.useNamesAsLabels(true);
roiManager("save", roiPath);

// Rebuild the saved binary mask from the final, manually QC-approved ROIs.
// This keeps the mask, ROI ZIP, measurements and numbered map synchronized.
finalMaskTitle = safeID + "_" + safeMarkers[0] + "_binary_mask";
newImage(finalMaskTitle, "8-bit black", imageWidth, imageHeight, 1);
setVoxelSize(pixelWidth, pixelHeight, 1, lengthUnit);
finalMaskID = getImageID();
saveSettings();
setForegroundColor(255, 255, 255);
for (i = 0; i < nVessels; i++) {
    selectImage(finalMaskID);
    roiManager("select", i);
    run("Fill");
}
restoreSettings();
selectImage(finalMaskID);
run("Select None");
saveAs("Tiff", maskPath);

// Build a dedicated results table without changing any pre-existing Results table.
resultsTitle = safeID + "_Vessel_Results";
Table.create(resultsTitle);

for (i = 0; i < nVessels; i++) {
    vesselID = RoiManager.getName(i);

    // Geometry is calculated once from the vessel-marker image.
    selectWindow(imageTitles[0]);
    roiManager("select", i);
    getRawStatistics(profilePixels, profileMean, profileMin, profileMax, profileSD);
    centroidX = getValue("X raw");
    centroidY = getValue("Y raw");
    perimeterPx = getValue("Perim. raw");
    feretPx = getValue("Feret raw");
    circularity = getValue("Circ.");
    getSelectionBounds(boxX, boxY, boxW, boxH);

    areaCalibrated = profilePixels * pixelWidth * pixelHeight;

    Table.set("Sample_ID", i, sampleID, resultsTitle);
    Table.set("Vessel_ID", i, vesselID, resultsTitle);
    Table.set("Crop_origin_X_full_px", i, cropX, resultsTitle);
    Table.set("Crop_origin_Y_full_px", i, cropY, resultsTitle);
    Table.set("Centroid_X_local_px", i, centroidX, resultsTitle);
    Table.set("Centroid_Y_local_px", i, centroidY, resultsTitle);
    Table.set("Centroid_X_full_px", i, cropX + centroidX, resultsTitle);
    Table.set("Centroid_Y_full_px", i, cropY + centroidY, resultsTitle);
    Table.set("Bounding_X_local_px", i, boxX, resultsTitle);
    Table.set("Bounding_Y_local_px", i, boxY, resultsTitle);
    Table.set("Bounding_width_px", i, boxW, resultsTitle);
    Table.set("Bounding_height_px", i, boxH, resultsTitle);
    Table.set("Area_px2", i, profilePixels, resultsTitle);
    Table.set("Area_calibrated", i, areaCalibrated, resultsTitle);
    Table.set("Length_unit", i, lengthUnit, resultsTitle);
    Table.set("Perimeter_px", i, perimeterPx, resultsTitle);
    Table.set("Feret_px", i, feretPx, resultsTitle);
    Table.set("Circularity", i, circularity, resultsTitle);

    // Measure every channel inside exactly the same vessel ROI.
    for (k = 0; k < nChannels; k++) {
        selectWindow(imageTitles[k]);
        roiManager("select", i);
        getRawStatistics(channelN, channelMean, channelMin, channelMax, channelSD);

        rawIntDen = channelN * channelMean;
        correctedMean = channelMean - backgroundMeans[k];
        correctedIntDen = channelN * correctedMean;
        prefix = safeMarkers[k];

        Table.set(prefix + "_Mean_raw", i, channelMean, resultsTitle);
        Table.set(prefix + "_SD_raw", i, channelSD, resultsTitle);
        Table.set(prefix + "_Min_raw", i, channelMin, resultsTitle);
        Table.set(prefix + "_Max_raw", i, channelMax, resultsTitle);
        Table.set(prefix + "_RawIntDen", i, rawIntDen, resultsTitle);
        Table.set(prefix + "_Background_mean", i, backgroundMeans[k], resultsTitle);
        Table.set(prefix + "_Mean_bg_corrected", i, correctedMean, resultsTitle);
        Table.set(prefix + "_CorrectedIntDen", i, correctedIntDen, resultsTitle);
    }

    Table.set(safeMarkers[0] + "_threshold_low_raw", i, thresholdLow, resultsTitle);
    Table.set(safeMarkers[0] + "_threshold_high_raw", i, thresholdHigh, resultsTitle);
}

selectWindow(resultsTitle);
Table.update;
Table.save(resultsPath);

// Create a permanently numbered reference image; never use this PNG for measurement.
mapWorkingTitle = safeID + "_Vessel_ID_map_working";
selectWindow(imageTitles[0]);
run("Select None");
run("Duplicate...", "title=[" + mapWorkingTitle + "]");
mapWorkingID = getImageID();
roiManager("show all with labels");
run("Flatten");
rename(safeID + "_Vessel_ID_map");
saveAs("PNG", mapPath);

// Close temporary images only if they still exist. The final QC-derived mask
// remains open for inspection. The ID checks prevent the prior end-of-run
// "No window with the title ... binary_mask found" error.
if (isOpen(mapWorkingID)) {
    selectImage(mapWorkingID);
    close();
}
if (isOpen(thresholdMaskID)) {
    selectImage(thresholdMaskID);
    close();
}

// Restore labeled ROIs on the original vessel-marker image for lookup.
selectWindow(imageTitles[0]);
roiManager("show all with labels");

settingsText =
    "Macro_version=2.0.2\n" +
    "ImageJ_version=" + getVersion() + "\n" +
    "Sample_ID=" + sampleID + "\n" +
    "Number_of_channels=" + nChannels + "\n" +
    "Number_of_target_channels=" + nTargets + "\n" +
    "Crop_width_px=" + imageWidth + "\n" +
    "Crop_height_px=" + imageHeight + "\n" +
    "Crop_origin_X_full_px=" + cropX + "\n" +
    "Crop_origin_Y_full_px=" + cropY + "\n" +
    "Pixel_width=" + pixelWidth + "\n" +
    "Pixel_height=" + pixelHeight + "\n" +
    "Length_unit=" + lengthUnit + "\n";

for (k = 0; k < nChannels; k++) {
    channelRole = "Target";
    if (k == 0)
        channelRole = "Vessel_segmentation";
    settingsText = settingsText +
        "Channel_" + (k + 1) + "_role=" + channelRole + "\n" +
        "Channel_" + (k + 1) + "_marker=" + markerNames[k] + "\n" +
        "Channel_" + (k + 1) + "_column_prefix=" + safeMarkers[k] + "\n" +
        "Channel_" + (k + 1) + "_image=" + imageTitles[k] + "\n" +
        "Channel_" + (k + 1) + "_bit_depth=" + bitDepths[k] + "\n" +
        "Channel_" + (k + 1) + "_background_mean_raw=" + backgroundMeans[k] + "\n";
}

settingsText = settingsText +
    "Threshold_marker=" + markerNames[0] + "\n" +
    "Threshold_mode=" + thresholdMode + "\n" +
    "Automatic_threshold_method=" + thresholdMethod + "\n" +
    "Final_threshold_low_raw=" + thresholdLow + "\n" +
    "Final_threshold_high_raw=" + thresholdHigh + "\n" +
    "Gaussian_sigma_px=" + gaussianSigma + "\n" +
    "Rolling_ball_radius_px=" + rollingRadius + "\n" +
    "Minimum_profile_area_" + safeLabel(lengthUnit) + "2=" + minArea + "\n" +
    "Exclude_edge_profiles=" + excludeEdges + "\n" +
    "Background_ROI_used=" + measureBackground + "\n" +
    "Number_of_background_ROIs=" + nBackgroundROIs + "\n" +
    "Final_profile_count=" + nVessels + "\n" +
    "Automatic_ROI_definition=Particle candidates intersected with the exact binary foreground mask to preserve interior holes.\n" +
    "Binary_mask_definition=Union of final manually QC-approved vessel-profile ROIs.\n" +
    "Measurement_definition=Intensity measured within manually QC-verified " + markerNames[0] + "-positive profile ROIs; interior holes excluded unless altered during manual ROI editing.\n";

File.saveString(settingsText, settingsPath);

showMessage("Analysis complete",
    nVessels + " " + markerNames[0] + "-positive vascular profiles were measured across " + nChannels + " channels.\n\n" +
    "Saved:\n" +
    "- vessel measurements CSV\n" +
    "- editable ROI ZIP\n" +
    "- " + markerNames[0] + " binary mask\n" +
    "- numbered vessel map\n" +
    "- analysis settings record\n\n" +
    "Select a Vessel_ID in ROI Manager to relocate that profile on the original " + markerNames[0] + " image.");

function safeLabel(s) {
    s = replace(s, " ", "_");
    s = replace(s, "/", "_");
    s = replace(s, "\\", "_");
    s = replace(s, ":", "_");
    s = replace(s, ";", "_");
    s = replace(s, ",", "_");
    s = replace(s, "*", "_");
    s = replace(s, "?", "_");
    s = replace(s, "\"", "_");
    s = replace(s, "'", "_");
    s = replace(s, "<", "_");
    s = replace(s, ">", "_");
    s = replace(s, "|", "_");
    s = replace(s, "[", "_");
    s = replace(s, "]", "_");
    s = replace(s, "(", "_");
    s = replace(s, ")", "_");
    s = replace(s, "{", "_");
    s = replace(s, "}", "_");
    s = replace(s, "=", "_");
    s = replace(s, "+", "_");
    s = replace(s, "-", "_");
    s = replace(s, ".", "_");
    s = replace(s, "\t", "_");
    s = replace(s, "\r", "_");
    s = replace(s, "\n", "_");
    return s;
}
