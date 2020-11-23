# This script train and classify feeding calls using a random forest classifier.
#
# 
#
# set working directory in terminal: cd /Users/imor384/Desktop/vocal_signatures/vocal_signature_random_forest
#
# To run this script, type:
#	Rscript Random_forest.R
# in the terminal.
#
##GitHub repo by nick Young: https://github.com/UoA-eResearch/bird_recognition
## Nick Young's code was modified by Stilianos Louca and Ines Moran
## Load the necessary packages


library(bioacoustics)
library(tuneR)
library(soundgen)
library(seewave)
library(dplyr)
library(tools)
library(randomForest)
library(stringr)
library(keras)
library(pbapply)
library(doParallel)
library(foreach)
library(caret)
library(parallel)
library(e1071)

# settings
DEFAULT_PLOT_WIDTH=4
DEFAULT_PLOT_HEIGHT=3
DEFAULT_NTHREADS=4
MFCC_NTHREADS=2

#RFID_bird_ID_full = "../vocal_signature_random_forest/input/sanitized_clips/RFID_bird_ID_full"
#Nest_ID= "../vocal_signature_random_forest/input/sanitized_clips/Nest_ID"
# 

# correct names of non-normalized clips
wavfiles = list.files("../vocal_signature_random_forest/input/sanitized_not_combined", "*\\.\\.wav", full.names=TRUE, include.dirs=FALSE, recursive = TRUE)
if(length(wavfiles)>0){
	cat(sprintf("Correcting %d wav file names..\n",length(wavfiles)))
	for(wavfile in wavfiles){
		file.rename(wavfile, paste0(substr(wavfile,1,nchar(wavfile)-nchar("..wav")),".wav"))
	}
}

#Normalizing calls to -1dB from soundgen
cat(sprintf("Normalizing calls (by RMS)..\n"))
for(grouping_normalized in c("Individual","Nest")){
	input_subdirs=list.files(sprintf("../vocal_signature_random_forest/input/sanitized_not_combined/%s", grouping_normalized), full.names=TRUE, include.dirs=TRUE, recursive = FALSE)
	for(input_subdir in input_subdirs){
		if(!dir.exists(input_subdir)) next
		output_folder = sprintf("../vocal_signature_random_forest/input/sanitized_normalized_clips/%s/%s", grouping_normalized,basename(input_subdir))
		if(dir.exists(output_folder)){
			cat(sprintf("  Note: Skipping subdir '%s', since the normalized output dir already exists\n",basename(output_folder)))
		}else{
			dir.create(output_folder, showWarnings = FALSE, recursive=TRUE)
			wavfiles = list.files(input_subdir, "*\\.wav", full.names=TRUE, include.dirs=FALSE, recursive = FALSE)
			for(wavfile in wavfiles){
				system(sprintf("ffmpeg-normalize '%s' -of '%s' -nt rms -f -ext wav",wavfile,output_folder))
			}
# 			soundgen::normalizeFolder(
# 				myfolder = input_subdir,
# 				type = "peak",
# 				maxAmp = -1,
# 				summaryFun = "mean",
# 				savepath = output_folder
# 			)
		}
	}
}


## Load the wave files. Render some spectrograms

for(normalization in c("sanitized_clips", "sanitized_normalized_combined_clips")){
	for(grouping in c("Individual","Nest")){
		cat(sprintf("Analyzing signature at level of %s..\n",grouping))
	
		# read wav files from a directory. Wav files should be named according to group, and each wav file includes all calls recorded from a particular group. "Group" can refer to an individual, to a nest, etc.
		files = list.files(sprintf("input/%s/%s",normalization,grouping), "*.wav", full.names=TRUE)
		files_without_extension = basename(file_path_sans_ext(files))
		files = files[files_without_extension != "golf"]
		cat(sprintf("  Note: Found %d wav files\n",length(files)))
		files_without_extension = files_without_extension[files_without_extension != "golf"]
		wavs = setNames(
		  mclapply(
			files,
			read_audio,
			mc.cores=DEFAULT_NTHREADS
		  ),
		  files_without_extension
		)
		samples = setNames(
		  mclapply(
			files,
			read_audio,
			from = 0,
			to = 4,
			mc.cores=DEFAULT_NTHREADS
		  ),
		  files_without_extension
		)

		#invisible(suppressWarnings(lapply(names(samples), function(birdid) {
		#  spectro(samples[[birdid]], main = birdid, osc = TRUE, flim=c(5,15), wl=300)
		#})))


		## Calculate Mel-frequency Cepstral Coefficients (MFCC)
		cat(sprintf("  Calculating Mel-frequency Cepstral Coefficients (MFCC) for %s..\n",grouping))
		wintime=0.25 # duration of each time frame used as a "sample". A longer wintime means more calls are included in one "sample", thus making classification easier, however it also means from a particular wav file you can extract fewer "samples". Note that the hoptime should be the same as wintime, so that successive time frames do not overlap and thus the samples are truly independent.
		mels = pbsapply(wavs, melfcc, wintime=wintime, hoptime=wintime, numcep = 50, cl=MFCC_NTHREADS) # Calculate all MFCCs
		head(mels[[1]]) # Result is a matrix of 20 coefficients across nrow time frames

		# plot MFCCs
		cat(sprintf("  Plotting MFCCs for %s..\n",grouping))
		margins = c(1,1,1,0.5)
		plotfile= sprintf("output/plots/MFCCs %s %s.pdf", grouping,normalization)
		dir.create(dirname(plotfile), showWarnings = FALSE, recursive=TRUE)
		pdf(file=plotfile, width=DEFAULT_PLOT_WIDTH+margins[2]+margins[4], height=DEFAULT_PLOT_HEIGHT+margins[1]+margins[3])
		image(mels[[1]])
		invisible(dev.off());


		# randomize the order of time frames (e.g., in case bird clips have not been concatenated randomly)
		cat(sprintf("  Randomizing the order of time frame for %s..\n",grouping))
		for (b in 1:length(mels)) {
		 mels[[b]] = mels[[b]][sample.int(n=nrow(mels[[b]]), size=nrow(mels[[b]]), replace=FALSE),]
		}

		# save MFCCs
		cat(sprintf("  Saving MFCCs for %s..\n",grouping))
		for (birdid in names(mels)) {
		  filename = paste0(sprintf("output/MELFCC/%s/%s",grouping,normalization), birdid, ".csv")
		  dir.create(dirname(filename), showWarnings = FALSE, recursive=TRUE)
		  write.csv(mels[[birdid]], filename, row.names = FALSE)
		}
		margins = c(2,2,2,0.5)
		plotfile=sprintf("output/plots/Ntime_frames_per_%s_%s.pdf", grouping,normalization)
		dir.create(dirname(plotfile), showWarnings = FALSE, recursive=TRUE)
		pdf(file=plotfile, width=DEFAULT_PLOT_WIDTH+margins[2]+margins[4], height=DEFAULT_PLOT_HEIGHT+margins[1]+margins[3])
		barplot(sapply(mels, nrow), las=2)
		invisible(dev.off())


		# The data is quite imbalanced

		## Train a random forest classifier
		cat(sprintf("  Training a random forest classifier for %s..\n",grouping))
		set.seed(1337)
		Ntrain=50 # number of samples per group (e.g. individual or nest) to use for training
		min_Ntest=10 # minimum number of samples per group used for testing; any group with insufficient samples will be omitted from the analysis
		sampled_mels = mels[lapply(mels, nrow) >= Ntrain+min_Ntest] # only keep groups with a sufficient number of samples (i.e. time frames)
		get_balanced_slice = function(offset = 0, size = Ntrain) {
		  sampled_mels = lapply(sampled_mels, function(x) x[offset:min(nrow(x),offset+size),])
		  trainM = do.call(rbind, sampled_mels) # melfcc gives a matrix - rbind to cast from 3D to 2D across the whole training dataset
		  trainM.labels = rep(factor(names(sampled_mels)), lapply(sampled_mels, nrow)) # Create bird labels for each MFCC
		  rownames(trainM) = trainM.labels
		  trainM
		}

		# train the random forest classifier
		trainM = get_balanced_slice(offset=0, size=Ntrain)
		rf = randomForest(trainM, factor(rownames(trainM)))
		sink(file=sprintf("output/random_forest_%s_%s.txt",grouping,normalization)); print(rf); sink();

		# estimate the accuracy of the classifier with the remaining samples from each group
		cat(sprintf("  Estimating the accuracy of the classifier with the remaining sample from each %s..\n",grouping))
		testM = get_balanced_slice(offset=Ntrain, size=Inf)
		predictions = predict(rf, testM, type="response")
		test_results = confusionMatrix(predictions, factor(rownames(testM)))
		sink(file=sprintf("output/accuracy_%s_%s.txt",grouping,normalization)); print(test_results); sink();
	}
}
