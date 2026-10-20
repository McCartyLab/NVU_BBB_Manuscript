/*
 * COMET DAPI-defined tumor WNT7A quantification
 * Version 1.1.1 (2026-10-12)
 * Requires Fiji/ImageJ 1.54a or later.
 *
 * PURPOSE
 *   Quantify WNT7A immunofluorescence in a DAPI-defined tumor ROI and an
 *   independently drawn adjacent non-tumor brain ROI on the same section.
 *
 * INTERPRETATION
 *   This macro measures WNT7A fluorescence/protein signal. WNT7A staining alone
 *   does not establish canonical WNT pathway activity. A signaling claim needs
 *   an appropriate downstream readout (for example, nuclear beta-catenin or a
 *   validated target such as AXIN2/LEF1) and suitable controls.
 *
 * REQUIRED INPUTS
 *   1. Original single-channel grayscale DAPI image.
 *   2. Exactly aligned original single-channel grayscale WNT7A image.
 *
 * IMPORTANT
 *   - Do not use RGB composites, JPEG screenshots, flattened images, or images
 *     with unequal acquisition settings for quantitative comparisons.
 *   - Define tumor without viewing the WNT7A result whenever possible.
 *   - DAPI contrast is not tumor-specific. Validate the tumor boundary against
 *     H&E or a tumor-specific marker when available.
 *   - Adjacent brain is a biological reference, not technical background.
 *   - A positive-area threshold must be established from matched negative
 *     controls and applied unchanged to the complete cohort.
 *   - The animal, not each ROI or section, is normally the biological replicate.
 *
 * OUTPUTS
 *   <ID>_WNT7A_region_measurements.csv
 *   <ID>_WNT7A_tumor_vs_nontumor_summary.csv
 *   <ID>_WNT7A_analysis_ROIs.zip
 *   <ID>_Tumor_mask.tif
 *   <ID>_Adjacent_nontumor_mask.tif
 *   <ID>_WNT7A_positive_mask.tif       (fixed-threshold mode only)
 *   <ID>_DAPI_ROI_QC_map.png           (visualization only)
 *   <ID>_WNT7A_ROI_QC_map.png          (visualization only)
 *   <ID>_WNT7A_settings.txt
 */

requires("1.54a");

if (nImages < 2)
    exit("Open the matched DAPI and WNT7A grayscale images before running this macro.");

titles = getList("image.titles");
backgroundModes = newArray(
    "No scalar background subtraction",
    "Measure pooled technical background from user-drawn ROIs"
);

Dialog.create("DAPI-defined tumor WNT7A analysis");
Dialog.addChoice("DAPI image for tumor delineation", titles, titles[0]);
Dialog.addChoice("WNT7A image for measurement", titles, titles[1]);
Dialog.addString("Animal/Sample ID", "Sample01", 20);
Dialog.addString("Section/ROI ID", "ROI01", 20);
Dialog.addString("Experimental group", "", 20);
Dialog.addNumber("Full-resolution crop origin X (pixels)", 0);
Dialog.addNumber("Full-resolution crop origin Y (pixels)", 0);
Dialog.show();

dapiTitle = Dialog.getChoice();
wntTitle = Dialog.getChoice();
sampleID = String.trim(Dialog.getString());
roiID = String.trim(Dialog.getString());
groupName = String.trim(Dialog.getString());
cropX = Dialog.getNumber();
cropY = Dialog.getNumber();

if (dapiTitle == wntTitle)
    exit("Select two different image windows: DAPI and WNT7A.");
if (sampleID == "" || roiID == "")
    exit("Animal/Sample ID and Section/ROI ID are required.");
if (cropX < 0 || cropY < 0)
    exit("Crop origin X and Y cannot be negative.");
if (cropX != round(cropX) || cropY != round(cropY))
    exit("Crop origin X and Y must be whole-pixel coordinates.");

selectedTitles = newArray(dapiTitle, wntTitle);
calRelativeTolerance = 0.00001;
for (k = 0; k < 2; k++) {
    selectWindow(selectedTitles[k]);
    getDimensions(tmpW, tmpH, tmpC, tmpZ, tmpT);
    tmpBD = bitDepth();
    getPixelSize(tmpUnit, tmpPW, tmpPH);
    tmpNormalizedUnit = normalizeLengthUnit(tmpUnit);

    if (tmpC != 1 || tmpZ != 1 || tmpT != 1)
        exit("'" + selectedTitles[k] + "' is not a single-channel, single-Z, single-time-point image.");
    if (tmpBD == 24)
        exit("'" + selectedTitles[k] + "' is flattened RGB (24-bit). Use the original individual grayscale channel image.");

    if (k == 0) {
        imageWidth = tmpW;
        imageHeight = tmpH;
        lengthUnit = tmpUnit;
        normalizedLengthUnit = tmpNormalizedUnit;
        pixelWidth = tmpPW;
        pixelHeight = tmpPH;
        dapiBitDepth = tmpBD;
    } else {
        if (tmpW != imageWidth || tmpH != imageHeight)
            exit("DAPI and WNT7A do not have identical pixel dimensions.");
        widthTolerance = calRelativeTolerance * maxOf(abs(pixelWidth), abs(tmpPW));
        heightTolerance = calRelativeTolerance * maxOf(abs(pixelHeight), abs(tmpPH));
        if (widthTolerance < 1e-9)
            widthTolerance = 1e-9;
        if (heightTolerance < 1e-9)
            heightTolerance = 1e-9;
        if (tmpNormalizedUnit != normalizedLengthUnit || abs(tmpPW - pixelWidth) > widthTolerance || abs(tmpPH - pixelHeight) > heightTolerance)
            exit("DAPI and WNT7A spatial calibrations differ. Assign the exact same pixel width, pixel height and unit without resizing either image.");
        wntBitDepth = tmpBD;
    }
}

if (pixelWidth <= 0 || pixelHeight <= 0)
    exit("Invalid spatial calibration: pixel width and height must be greater than zero.");

Dialog.create("WNT7A measurement settings");
Dialog.addChoice("Background correction", backgroundModes, backgroundModes[0]);
Dialog.addNumber("Number of technical-background ROIs", 3);
Dialog.addCheckbox("Calculate WNT7A-positive area using a fixed raw threshold", false);
Dialog.addNumber("Fixed WNT7A-positive threshold (raw intensity)", 0);
Dialog.addMessage("Use positive-area analysis only with a threshold established from matched negative controls and applied unchanged to all samples.");
Dialog.show();

backgroundMode = Dialog.getChoice();
nBackgroundROIs = round(Dialog.getNumber());
calculatePositiveArea = Dialog.getCheckbox();
positiveThreshold = Dialog.getNumber();

measureBackground = false;
if (backgroundMode == backgroundModes[1])
    measureBackground = true;
if (measureBackground && nBackgroundROIs < 1)
    exit("Use at least one technical-background ROI, or select no background subtraction.");
if (calculatePositiveArea && positiveThreshold <= 0)
    exit("Positive-area analysis requires a fixed raw threshold greater than zero.");

safeID = safeLabel(sampleID + "_" + roiID);
outputDir = getDirectory("Choose output folder");
if (outputDir == "")
    exit("No output folder was selected.");

resultsPath = outputDir + safeID + "_WNT7A_region_measurements.csv";
summaryPath = outputDir + safeID + "_WNT7A_tumor_vs_nontumor_summary.csv";
roiZipPath = outputDir + safeID + "_WNT7A_analysis_ROIs.zip";
tumorMaskPath = outputDir + safeID + "_Tumor_mask.tif";
nonTumorMaskPath = outputDir + safeID + "_Adjacent_nontumor_mask.tif";
positiveMaskPath = outputDir + safeID + "_WNT7A_positive_mask.tif";
dapiQcMapPath = outputDir + safeID + "_DAPI_ROI_QC_map.png";
wntQcMapPath = outputDir + safeID + "_WNT7A_ROI_QC_map.png";
settingsPath = outputDir + safeID + "_WNT7A_settings.txt";

outputExists = false;
if (File.exists(resultsPath) || File.exists(summaryPath) || File.exists(roiZipPath) ||
    File.exists(tumorMaskPath) || File.exists(nonTumorMaskPath) ||
    File.exists(dapiQcMapPath) || File.exists(wntQcMapPath) ||
    File.exists(settingsPath))
    outputExists = true;
if (calculatePositiveArea && File.exists(positiveMaskPath))
    outputExists = true;
if (outputExists) {
    overwriteOK = getBoolean("One or more WNT7A outputs for '" + safeID + "' already exist. Overwrite them?", "Overwrite", "Cancel");
    if (!overwriteOK)
        exit("Analysis cancelled. Use a different Sample/ROI ID or output folder.");
}

existingROIs = roiManager("count");
if (existingROIs > 0) {
    clearOK = getBoolean("ROI Manager contains " + existingROIs + " ROI(s). Clear them for this WNT7A analysis?", "Clear and continue", "Cancel");
    if (!clearOK)
        exit("Analysis cancelled without changing the existing ROI list.");
}
roiManager("reset");

// Define tumor independently on DAPI.
selectWindow(dapiTitle);
run("Select None");
setTool("freehand");
waitForUser("Define DAPI tumor ROI",
    "On the ORIGINAL DAPI image, draw one area ROI around the tumor.\n\n" +
    "Use DAPI morphology/contrast only; do not use the WNT7A pattern to choose the boundary.\n" +
    "Exclude tissue-free space, folds, tears and clearly necrotic cavities.\n" +
    "A composite area selection is allowed.\n\n" +
    "Click OK here after the tumor outline is complete.");

st = selectionType;
if (st < 0 || (st >= 5 && st <= 8) || st == 10)
    exit("The tumor selection must be an area ROI, not a line or point selection.");
getRawStatistics(tumorPixelsCheck);
if (tumorPixelsCheck <= 0)
    exit("The tumor ROI contains no pixels.");
roiManager("add");
tumorIndex = 0;
roiManager("select", tumorIndex);
roiManager("rename", "Tumor_DAPI_defined");

// Define adjacent non-tumor brain independently on DAPI.
selectWindow(dapiTitle);
roiManager("show all without labels");
run("Select None");
setTool("freehand");
waitForUser("Define adjacent non-tumor brain ROI",
    "The stored tumor outline is visible. On the ORIGINAL DAPI image, draw an\n" +
    "adjacent non-tumor mouse-brain tissue ROI.\n\n" +
    "Exclude tumor, tissue-free space, ventricles, folds, tears and damaged tissue.\n" +
    "Choose the reference region using DAPI morphology, not WNT7A intensity.\n\n" +
    "Click OK here after the non-tumor outline is complete.");

st = selectionType;
if (st < 0 || (st >= 5 && st <= 8) || st == 10)
    exit("The adjacent non-tumor selection must be an area ROI.");
getRawStatistics(nonTumorPixelsCheck);
if (nonTumorPixelsCheck <= 0)
    exit("The adjacent non-tumor ROI contains no pixels.");
roiManager("add");
nonTumorIndex = 1;
roiManager("select", nonTumorIndex);
roiManager("rename", "Adjacent_non_tumor_brain");
roiManager("show none");

// The two biological regions must not share pixels. A temporary binary image
// makes this check deterministic even when the ROI Manager AND operation would
// otherwise return an empty selection.
overlapCheckTitle = safeID + "_ROI_overlap_check";
newImage(overlapCheckTitle, "8-bit black", imageWidth, imageHeight, 1);
saveSettings();
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0);
roiManager("select", tumorIndex);
run("Fill");
roiManager("select", nonTumorIndex);
getRawStatistics(overlapTestN, overlapTestMean);
restoreSettings();
overlapPixels = round(overlapTestN * overlapTestMean / 255);
close();
if (overlapPixels > 0)
    exit("Tumor and adjacent non-tumor ROIs overlap by " + overlapPixels + " pixels. Rerun and draw non-overlapping regions.");

// Optional pooled technical background. This is kept distinct from adjacent
// brain, which may contain genuine WNT7A and is therefore not 'background'.
backgroundPixels = 0;
backgroundSum = 0;
if (measureBackground) {
    selectWindow(wntTitle);
    roiManager("show all without labels");
    run("Select None");

    for (b = 0; b < nBackgroundROIs; b++) {
        selectWindow(wntTitle);
        run("Select None");
        setTool("rectangle");
        waitForUser("Select WNT7A technical background ROI " + (b + 1) + " of " + nBackgroundROIs,
            "On the ORIGINAL WNT7A image, draw a representative technical-background\n" +
            "rectangle. Avoid tissue signal, staining deposits, folds and saturated pixels.\n\n" +
            "Do not use adjacent normal brain as background because it may contain genuine WNT7A.\n\n" +
            "Click OK here after drawing the rectangle.");

        if (selectionType != 0)
            exit("Each technical-background selection must be a rectangular area ROI.");
        roiManager("add");
        backgroundIndex = roiManager("count") - 1;
        roiManager("select", backgroundIndex);
        roiManager("rename", "Technical_background_" + IJ.pad(b + 1, 2));

        selectWindow(wntTitle);
        roiManager("select", backgroundIndex);
        getRawStatistics(bgN, bgMean);
        if (bgN <= 0)
            exit("Technical-background ROI " + (b + 1) + " contained no pixels.");
        backgroundPixels = backgroundPixels + bgN;
        backgroundSum = backgroundSum + bgN * bgMean;
    }
    backgroundMean = backgroundSum / backgroundPixels;
    roiManager("show none");
} else {
    backgroundMean = 0;
    nBackgroundROIs = 0;
}

roiManager("remove slice info");
RoiManager.useNamesAsLabels(true);
roiManager("save", roiZipPath);

selectWindow(dapiTitle);
dapiImageID = getImageID();
selectWindow(wntTitle);
wntImageID = getImageID();
getRawStatistics(wntAllN, wntAllMean, wntAllMin, wntAllMax, wntAllSD);
if (calculatePositiveArea && positiveThreshold > wntAllMax)
    exit("The fixed WNT7A-positive threshold is greater than the maximum raw WNT7A value (" + wntAllMax + ").");

// Save exact binary masks for audit.
tumorMaskTitle = safeID + "_Tumor_mask";
newImage(tumorMaskTitle, "8-bit black", imageWidth, imageHeight, 1);
setVoxelSize(pixelWidth, pixelHeight, 1, lengthUnit);
tumorMaskID = getImageID();
saveSettings();
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0);
roiManager("select", tumorIndex);
run("Fill");
restoreSettings();
run("Select None");
saveAs("Tiff", tumorMaskPath);

nonTumorMaskTitle = safeID + "_Adjacent_nontumor_mask";
newImage(nonTumorMaskTitle, "8-bit black", imageWidth, imageHeight, 1);
setVoxelSize(pixelWidth, pixelHeight, 1, lengthUnit);
nonTumorMaskID = getImageID();
saveSettings();
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0);
roiManager("select", nonTumorIndex);
run("Fill");
restoreSettings();
run("Select None");
saveAs("Tiff", nonTumorMaskPath);

// Create a fixed raw-intensity positive mask only when requested.
positiveMaskID = 0;
if (calculatePositiveArea) {
    selectImage(wntImageID);
    run("Duplicate...", "title=[" + safeID + "_WNT7A_positive_mask_working]");
    positiveMaskID = getImageID();
    setThreshold(positiveThreshold, wntAllMax, "raw");
    saveSettings();
    setOption("BlackBackground", true);
    run("Convert to Mask");
    restoreSettings();
    rename(safeID + "_WNT7A_positive_mask");
    saveAs("Tiff", positiveMaskPath);
}

regionNames = newArray("Tumor", "Adjacent_non_tumor_brain");
regionIndexes = newArray(tumorIndex, nonTumorIndex);
regionAreas = newArray(2);
regionMeans = newArray(2);
regionCorrectedMeans = newArray(2);
regionPositiveFractions = newArray(2);

resultsTitle = safeID + "_WNT7A_Region_Results";
Table.create(resultsTitle);

for (r = 0; r < 2; r++) {
    regionName = regionNames[r];
    regionIndex = regionIndexes[r];

    selectImage(wntImageID);
    roiManager("select", regionIndex);
    getRawStatistics(regionPixels, wntMean, wntMin, wntMax, wntSD);
    if (regionPixels <= 0)
        exit(regionName + " contained no measurable WNT7A pixels.");

    selectImage(dapiImageID);
    roiManager("select", regionIndex);
    getRawStatistics(dapiPixels, dapiMean, dapiMin, dapiMax, dapiSD);
    if (dapiPixels != regionPixels)
        exit("DAPI and WNT7A pixel counts differ inside " + regionName + ". Confirm exact image alignment and dimensions.");

    areaCalibrated = regionPixels * pixelWidth * pixelHeight;
    rawIntDen = regionPixels * wntMean;
    correctedMean = wntMean - backgroundMean;
    correctedIntDen = regionPixels * correctedMean;

    thresholdOutput = "NA_not_enabled";
    positivePixels = 0;
    positiveArea = 0;
    positiveFraction = 0;
    if (calculatePositiveArea) {
        thresholdOutput = positiveThreshold;
        selectImage(positiveMaskID);
        roiManager("select", regionIndex);
        getRawStatistics(posN, posMean);
        positivePixels = round(posN * posMean / 255);
        positiveArea = positivePixels * pixelWidth * pixelHeight;
        positiveFraction = 100 * positivePixels / regionPixels;
    }

    regionAreas[r] = areaCalibrated;
    regionMeans[r] = wntMean;
    regionCorrectedMeans[r] = correctedMean;
    regionPositiveFractions[r] = positiveFraction;

    Table.set("Animal_Sample_ID", r, sampleID, resultsTitle);
    Table.set("Section_ROI_ID", r, roiID, resultsTitle);
    Table.set("Experimental_group", r, groupName, resultsTitle);
    Table.set("Region", r, regionName, resultsTitle);
    Table.set("Crop_origin_X_full_px", r, cropX, resultsTitle);
    Table.set("Crop_origin_Y_full_px", r, cropY, resultsTitle);
    Table.set("Area_px2", r, regionPixels, resultsTitle);
    Table.set("Area_calibrated", r, areaCalibrated, resultsTitle);
    Table.set("Area_unit", r, lengthUnit + "^2", resultsTitle);
    Table.set("DAPI_Mean_raw", r, dapiMean, resultsTitle);
    Table.set("DAPI_SD_raw", r, dapiSD, resultsTitle);
    Table.set("WNT7A_Mean_raw", r, wntMean, resultsTitle);
    Table.set("WNT7A_SD_raw", r, wntSD, resultsTitle);
    Table.set("WNT7A_Min_raw", r, wntMin, resultsTitle);
    Table.set("WNT7A_Max_raw", r, wntMax, resultsTitle);
    Table.set("WNT7A_RawIntDen", r, rawIntDen, resultsTitle);
    Table.set("WNT7A_Background_mean", r, backgroundMean, resultsTitle);
    Table.set("WNT7A_Mean_bg_corrected", r, correctedMean, resultsTitle);
    Table.set("WNT7A_CorrectedIntDen", r, correctedIntDen, resultsTitle);
    if (calculatePositiveArea) {
        Table.set("WNT7A_positive_threshold_raw", r, thresholdOutput, resultsTitle);
        Table.set("WNT7A_positive_pixels", r, positivePixels, resultsTitle);
        Table.set("WNT7A_positive_area_calibrated", r, positiveArea, resultsTitle);
        Table.set("WNT7A_positive_area_fraction_pct", r, positiveFraction, resultsTitle);
    }
}

selectWindow(resultsTitle);
Table.update;
Table.save(resultsPath);

tumorMinusNonTumorRaw = regionMeans[0] - regionMeans[1];
tumorMinusNonTumorCorrected = regionCorrectedMeans[0] - regionCorrectedMeans[1];
tumorToNonTumorRawRatio = "NA_nonpositive_denominator";
tumorToNonTumorCorrectedRatio = "NA_nonpositive_value";
tumorVsNonTumorRawPctChange = "NA_nonpositive_denominator";
if (regionMeans[1] > 0) {
    tumorToNonTumorRawRatio = regionMeans[0] / regionMeans[1];
    tumorVsNonTumorRawPctChange = 100 * tumorMinusNonTumorRaw / regionMeans[1];
}
if (regionCorrectedMeans[0] > 0 && regionCorrectedMeans[1] > 0)
    tumorToNonTumorCorrectedRatio = regionCorrectedMeans[0] / regionCorrectedMeans[1];

summaryTitle = safeID + "_WNT7A_Tumor_vs_NonTumor";
Table.create(summaryTitle);
Table.set("Animal_Sample_ID", 0, sampleID, summaryTitle);
Table.set("Section_ROI_ID", 0, roiID, summaryTitle);
Table.set("Experimental_group", 0, groupName, summaryTitle);
Table.set("Tumor_area_calibrated", 0, regionAreas[0], summaryTitle);
Table.set("Adjacent_nontumor_area_calibrated", 0, regionAreas[1], summaryTitle);
Table.set("Tumor_WNT7A_Mean_raw", 0, regionMeans[0], summaryTitle);
Table.set("Adjacent_nontumor_WNT7A_Mean_raw", 0, regionMeans[1], summaryTitle);
Table.set("Tumor_minus_nontumor_Mean_raw", 0, tumorMinusNonTumorRaw, summaryTitle);
Table.set("Tumor_to_nontumor_Mean_raw_ratio", 0, tumorToNonTumorRawRatio, summaryTitle);
Table.set("Tumor_vs_nontumor_Mean_raw_pct_change", 0, tumorVsNonTumorRawPctChange, summaryTitle);
Table.set("WNT7A_Background_mean", 0, backgroundMean, summaryTitle);
Table.set("Tumor_WNT7A_Mean_bg_corrected", 0, regionCorrectedMeans[0], summaryTitle);
Table.set("Adjacent_nontumor_WNT7A_Mean_bg_corrected", 0, regionCorrectedMeans[1], summaryTitle);
Table.set("Tumor_minus_nontumor_Mean_bg_corrected", 0, tumorMinusNonTumorCorrected, summaryTitle);
Table.set("Tumor_to_nontumor_Mean_bg_corrected_ratio", 0, tumorToNonTumorCorrectedRatio, summaryTitle);
if (calculatePositiveArea) {
    Table.set("WNT7A_positive_threshold_raw", 0, thresholdOutput, summaryTitle);
    Table.set("Tumor_WNT7A_positive_area_fraction_pct", 0, regionPositiveFractions[0], summaryTitle);
    Table.set("Adjacent_nontumor_WNT7A_positive_area_fraction_pct", 0, regionPositiveFractions[1], summaryTitle);
    Table.set("Tumor_minus_nontumor_positive_area_fraction_pct", 0, regionPositiveFractions[0] - regionPositiveFractions[1], summaryTitle);
}
selectWindow(summaryTitle);
Table.update;
Table.save(summaryPath);

// Create labeled DAPI and WNT7A reference maps for audit. Never measure PNGs.
selectImage(dapiImageID);
run("Select None");
run("Duplicate...", "title=[" + safeID + "_DAPI_ROI_QC_map_working]");
dapiQcWorkingID = getImageID();
roiManager("show all with labels");
run("Flatten");
rename(safeID + "_DAPI_ROI_QC_map");
saveAs("PNG", dapiQcMapPath);
if (isOpen(dapiQcWorkingID)) {
    selectImage(dapiQcWorkingID);
    close();
}

selectImage(wntImageID);
run("Select None");
run("Duplicate...", "title=[" + safeID + "_WNT7A_ROI_QC_map_working]");
wntQcWorkingID = getImageID();
roiManager("show all with labels");
run("Flatten");
rename(safeID + "_WNT7A_ROI_QC_map");
saveAs("PNG", wntQcMapPath);
if (isOpen(wntQcWorkingID)) {
    selectImage(wntQcWorkingID);
    close();
}

settingsText =
    "Macro_name=COMET_DAPI_Defined_Tumor_WNT7A_Quantification\n" +
    "Macro_version=1.1.1\n" +
    "ImageJ_version=" + getVersion() + "\n" +
    "Animal_Sample_ID=" + sampleID + "\n" +
    "Section_ROI_ID=" + roiID + "\n" +
    "Experimental_group=" + groupName + "\n" +
    "DAPI_image=" + dapiTitle + "\n" +
    "WNT7A_image=" + wntTitle + "\n" +
    "Image_width_px=" + imageWidth + "\n" +
    "Image_height_px=" + imageHeight + "\n" +
    "Crop_origin_X_full_px=" + cropX + "\n" +
    "Crop_origin_Y_full_px=" + cropY + "\n" +
    "Pixel_width=" + pixelWidth + "\n" +
    "Pixel_height=" + pixelHeight + "\n" +
    "Length_unit=" + lengthUnit + "\n" +
    "DAPI_bit_depth=" + dapiBitDepth + "\n" +
    "WNT7A_bit_depth=" + wntBitDepth + "\n" +
    "Background_mode=" + backgroundMode + "\n" +
    "Number_of_background_ROIs=" + nBackgroundROIs + "\n" +
    "WNT7A_background_mean_raw=" + backgroundMean + "\n" +
    "Positive_area_enabled=" + calculatePositiveArea + "\n" +
    "WNT7A_positive_threshold_raw=" + thresholdOutput + "\n" +
    "Tumor_definition=DAPI morphology/contrast; user-drawn ROI independent of WNT7A.\n" +
    "Adjacent_brain_definition=User-drawn non-overlapping adjacent non-tumor brain ROI selected using DAPI morphology.\n" +
    "Interpretation=WNT7A immunofluorescence/protein signal; not a direct measurement of WNT pathway activity.\n" +
    "Positive_threshold_rule=Must be established from matched negative controls and applied unchanged to the cohort.\n" +
    "Biological_replicate_rule=Summarize sections/ROIs within each animal; do not treat ROIs or pixels as independent animals.\n";
File.saveString(settingsText, settingsPath);

if (calculatePositiveArea) {
    if (isOpen(positiveMaskID)) {
        selectImage(positiveMaskID);
        close();
    }
}

selectImage(dapiImageID);
roiManager("show all with labels");
showProgress(1, 1);
showStatus("DAPI-defined tumor WNT7A analysis complete");

completionText =
    "WNT7A fluorescence was quantified in the DAPI-defined tumor and adjacent non-tumor brain ROIs.\n\n" +
    "Primary comparisons:\n" +
    "- Tumor_WNT7A_Mean_raw versus Adjacent_nontumor_WNT7A_Mean_raw\n" +
    "- Tumor_minus_nontumor_Mean_raw\n" +
    "- WNT7A-positive area fraction only when a control-derived fixed threshold was enabled\n\n" +
    "Use the QC map and saved masks to verify region placement. Report WNT7A fluorescence/abundance, not signaling activity, unless an appropriate downstream pathway readout was also measured.";
showMessage("WNT7A analysis complete", completionText);

function normalizeLengthUnit(unitText) {
    normalized = toLowerCase(String.trim(unitText));
    normalized = replace(normalized, "µ", "u");
    normalized = replace(normalized, "μ", "u");
    canonicalUnit = normalized;
    if (normalized == "um" || normalized == "micron" || normalized == "microns" ||
        normalized == "micrometer" || normalized == "micrometers")
        canonicalUnit = "micron";
    if (normalized == "pixel" || normalized == "pixels" || normalized == "pix")
        canonicalUnit = "pixel";
    if (normalized == "nm" || normalized == "nanometer" || normalized == "nanometers")
        canonicalUnit = "nanometer";
    if (normalized == "mm" || normalized == "millimeter" || normalized == "millimeters")
        canonicalUnit = "millimeter";
    return canonicalUnit;
}

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
    s = replace(s, "\t", "_");
    s = replace(s, "\r", "_");
    s = replace(s, "\n", "_");
    return s;
}
