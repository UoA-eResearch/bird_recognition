# This script detects audio signals using cross-correlations using monitoR.
#
#
# To run this script, type:
#	Rscript Signal_detection_monitoR.R
# in the terminal.
# This script was written by S. Louca

########################
# OPTIONS

#set working directory
setwd("~/Desktop/Signal_detection/Sound_templates") 
path_to_files = "/Users/imor384/Desktop/Signal_detection/Sound_templates"

DETECTION_SCORE_THRESHOLD=0.3

#########################
# PREPARATIONS

# load packages
suppressWarnings(suppressMessages(library(warbleR,quiet=TRUE)))
library(warbleR)
library(monitoR)



#########################
# AUXILIARY FUNCTIONS

# modified from: monitoR::Methods-plot (plotting template scores), version 1.0.7
# corrected a minor sprintf formatting bug in the original package
plot_peak_search=function(	x,  
							flim=c(0, 12),                 # Frequency limits for the spectrogram.
							scorelim,                     # Plot limits for scores
							which.one=names(x@templates), # Name(s) of templates to plot
							box=TRUE,                     # Set to FALSE to surpress boxes in spectrogram showing hits
							spec.col=gray.2(),            # Color palette for spectrogram
							t.each=30,                    # Time shown for each individual plot (s)
							hit.marker='lines',           # Markers for hits in score plot
							color=c('red', 'blue', 'green', 'orange', 'purple', 'pink', 'darkgreen', 'turquoise', 'royalblue', 'orchid4', 'brown', 'salmon2'), # Colors for individual templates
							legend=TRUE,                  # Set to FALSE to surpress legend
							all.peaks=FALSE,              # Set to TRUE to indicate locations of all peaks
							ask=if(dev.list() == 2) TRUE else FALSE){
    survey <- x@survey
    t.survey <- length(survey@left)/survey@samp.rate 
    n.plots <- ceiling(t.survey/t.each)
    t.start <- 1:n.plots*t.each - t.each
    if(n.plots == 1) t.each <- t.survey
    t.end <- t.start + t.each
    t.end[t.end>t.survey] <- t.survey
    t.start[n.plots] <- t.end[n.plots] - t.each # Adjust start of last plot back so it has the same length as the others

    # Pull out spectrogram data from scores object
    # Based on first template
    amp <- x@survey.data[[1]]$amp
    t.bins <- x@survey.data[[1]]$t.bins
    frq.bins <- x@survey.data[[1]]$frq.bins
 
    # Sort out colors for lines and boxes
    names.t <- names(x@templates)
    n.templates <- length(names.t)
    color <- c(rep(color, n.templates %/% length(color)), color[1:n.templates%%length(color)])
    names(color) <- names.t

    # Get scorelim
    if(missing(scorelim)) {
      upr <- 0
      for(i in seq(length(x@scores))) {
        upr <- max(upr, x@scores[[i]]$score)
      }
      scorelim <- c(0, upr)
    }

    oldpar <- par(mar=c(1, 4,1, 1), oma=c(6, 0,0, 0), mfrow=c(2, 1))
    oldask <- par(ask=par('ask'))
    on.exit(par(c(oldpar, oldask)))

    # Loop through time windows, plotting a spectrogram for each time
    for(i in 1:length(t.start)) {

      message(paste(t.start[i], 'to', t.end[i], 'seconds'))
    
      times <- t.bins[t.bins>=t.start[i] & t.bins<=t.end[i]]
      amp.clip <- amp[, t.bins %in% times]
      image(x=times, y=frq.bins, t(amp.clip), ylim=flim, col=spec.col, xlab='', ylab='Frequency (kHz)', xaxt='n', las=1)
 
      # Loop through templates and add boxes around detections
      for(j in which.one) {
        template <- x@templates[[j]]
        if(all.peaks) pks <- x@peaks[[j]] else pks <- x@detections[[j]]
        pks.clip <- pks[pks$time + template@duration >= t.start[i] & pks$time - template@duration <= t.end[i], ]

        if(box & nrow(pks.clip)>0) {
          for(k in 1:nrow(pks.clip)) { 
            xleft <- pks.clip$time[k] - template@duration/2
            xright <- pks.clip$time[k] + template@duration/2
            ylwr <- template@frq.lim[1]
            yupr <- template@frq.lim[2]
            polygon(x=c(xleft, xleft, xright, xright), y=c(ylwr, yupr, yupr, ylwr), border=color[j], lwd=1)
          }
        }
      }

      # Make plot of scores. Can't sort out xlab for some reason.
      plot(NULL, xlim=c(t.start[i], t.end[i]), ylim=scorelim, xlab='', ylab='Score', type='n', xaxs='i', las=1, mgp=c(3, 1,0))
      mtext("Time (s or min:sec)", 1,2.5, outer=TRUE)

      # Add x axis as mm:ss
      xaxp.sec <- par('xaxp')
      labs.sec <- seq(xaxp.sec[1], xaxp.sec[2], length.out=xaxp.sec[3]+1)
      labs.mmss <- paste(sprintf('%.2g', labs.sec%/%60), ':', sprintf('%.2g', labs.sec%%60), sep='')
      axis(1, at=labs.sec, labels=labs.mmss, mgp=c(3, 1.9, 0))

      if(legend) legend('topright', which.one, lty=1, col=color[which.one], cex=0.7)

      # Loop through templates 
      for(j in which.one) {
        template <- x@templates[[j]]
        score <- x@scores[[j]]         # score output from sccDetect. The correlation coefficients within which hits were found.
        if(all.peaks) pks <- x@peaks[[j]] else pks <- x@detections[[j]]
        cutoff <- template@score.cutoff      # If given, will plot a horizontal line at the correlation coefficient cutoff.

        score.clip <- score[score$time>=t.start[i] & score$time<=t.end[i], ]
        pks.clip <- pks[pks$time + template@duration >= t.start[i] & pks$time - template@duration <= t.end[i], ]

        lines(score.clip$time, score.clip$score, col=color[j])
        if(hit.marker == 'points') points(pks.clip$time, pks.clip$score, col=color[j]) else
          if(hit.marker == 'lines') abline(v=pks.clip$time, col=color[j]) 
        if(is.vector(cutoff)) abline(h=cutoff, lty=2, col=color[j])
      }
      par(ask=ask)
    }
}


getClip <- function(
  clip, 
  name="clip", 
  output="file", 
  write.wav=FALSE
) {

  if(class(clip) == "list" | (class(clip) == "character" && length(clip)>1)) {
    clist <- list()
    if(grepl(", ",name)) {
      name <- gsub(".*\\(", "", name)
      name <- gsub("\\)", "", name)
      name <- strsplit(name, ",")[[1]]
    }
    for(i in seq(length(clip))) {
      clist[[i]] <- getOneClip(clip[[i]], paste0(name[i], i), output, write.wav) 
    }
    return(clist)
  } 

  return(getOneClip(clip, name, output, write.wav))

}

getOneClip <- function(
  clip, 
  name, 
  output, 
  write.wav
) {

  if(output == "file") {
    if(class(clip) == "Wave") {
      fname <- paste0(name, ".wav")
      if(!write.wav) {
	stop("output argument is \"file\" but write.wav argument is FALSE so this function will not create a file. Set write.wav=TRUE to create a file, or else specify a wav file instead of a Wave object.")
      }
      if(file.exists(fname)) stop("Will not create a wav file from this clip because a file with name ", fname, " already exists.")
      else tuneR::writeWave(clip, fname) 
      return(fname)
    } else 
    if(class(clip) == "character") {
      if(!file.exists(clip)) stop("clip argument seems to be a file name but no file with the name ", clip, " exists.")
      return(clip)
    } else 
    stop("Can\'t figure out what to do with this clip:", clip, "with class:", class(clip))
  }

  if(output == "Wave") {
    if(class(clip) == "Wave") {
      return(clip)
    } else 
    if(class(clip) == "character") {
      if(!file.exists(clip)) stop("clip argument seems to be a file name but no file with the name ", clip, " exists!")
      file.ext <- tolower(gsub(".*\\.", "", clip))
      if(file.ext == "wav") 
        clip <- tuneR::readWave(clip) else 
      if(file.ext == "mp3") 
        clip <- readMP3(clip) else stop("File extension must be wav or mp3, but got ", file.ext)
      return(clip)
    } else 
    stop("Can\'t figure out what to do with this clip:", clip, "with class:", class(clip))
  }

}   

# Reads a single wav or mp3 file 
readClip <- function(clip) {

  if(class(clip) != "character" | length(clip) != 1) stop("Expected a length-one character vector for clip, but got a length ", length(clip), " ", class(clip), " object.")
  if(!file.exists(clip)) stop("clip argument seems to be a file name but no file with the name ", clip, " exists!")
 
  file.ext <- tolower(gsub(".*\\.", "", clip))
  if(file.ext == "wav") return(tuneR::readWave(filename=clip))
  if(file.ext == "mp3") return(readMP3(filename=clip)) 
  stop("File extension must be wav or mp3, but got ", file.ext)

}


spectro2 = function (
   wave, 
   wl = 512, 
   wn = "hanning", 
   zp = 0, 
   ovlp = 0, 
   fftw = FALSE, 
   dB = "max0", 
   dBref = NULL,
   ... # dots permit users to keep arg 'warn' in calls to binMatch/corMatch
   ) 
{
    if (!is.null(dB) && all(dB != c("max0", "A", "B", "C", "D"))) 
        stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")
    f <- wave@samp.rate
    wave <- as.matrix(wave@left)
    n <- nrow(wave)
    #step <- seq(1, n - wl, wl - (ovlp*wl/100)) # This is original, and will just leave out any partial stuff at the end
    step <- seq(1, n - wl + 1, wl - (ovlp*wl/100)) # New version squeezes out one more time bin in some cases
    # Two lines below drop parts of wave that go beyond step vector. I added them, but are they needed? I think so, at least for n, so the X is correct
    # New calculation of n is a bit tricky because of tails with no overlap at beginning and end
    n <- length(step)*wl  - (length(step)-1)*ovlp*wl/100 
    wave <- wave[1:n, ,drop=FALSE]
    z <- stft(wave = wave, f = f, wl = wl, zp = zp, step = step, wn = wn, fftw = fftw)
    X <- seq(0, (n - wl)/f, length.out = length(step)) # X is time, and here is for the left edge of time bins (could be right too, or center even). Note that the left edge of the last time bin is not affected by ovlp.
    Y <- seq((f/1000)/(wl + zp), f/2000, length.out = nrow(z))

    if (!is.null(dB)) {
        if (is.null(dBref)) 
            z <- 20 * log10(z)
        else z <- 20 * log10(z/dBref)
        if (dB == "max0") 
            z <- z
        if (dB == "A") 
            z <- dBweight(Y * 1000, dBref = z)$A
        if (dB == "B") 
            z <- dBweight(Y * 1000, dBref = z)$B
        if (dB == "C") 
            z <- dBweight(Y * 1000, dBref = z)$C
        if (dB == "D") 
            z <- dBweight(Y * 1000, dBref = z)$D
    }

    return(list(time = X, freq = Y, amp = z))
}


corMatch2 = function(
  survey,                 # Complete survey which is to be analyzed for calls. Wave object or vector.
  templates,              # Template list, made with makeTemplates
  parallel=FALSE,         # If TRUE, mclapply is used for correlation calculations, for parallel processing (Linux or Mac OS X only). If FALSE lapply is used.
  show.prog=FALSE,        # If TRUE, progress is displayed during correlation calculations 
  cor.method='pearson',   # Method used by cor function (see ?cor)
#  warn=TRUE,              # Set to FALSE to surpress warnings
  time.source='filename', # 'filename' or 'fileinfo' as the mtime source
  rec.tz=NA,              # Time zone setting for recorders 
  write.wav=FALSE,        # Set to TRUE to allow creation of file of survey in working directory
  quiet=FALSE,            # TRUE will suppress calls to cat()
  ...                     # Additional arguments to the spectro function
) {

  # Check arguments
  if(missing(survey)) stop('Required argument survey is missing.')
  if(missing(templates)) stop('Required argument templates is missing.')

  # Packages
  if(parallel) {
    lapplyfun <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores=parallel::detectCores())
  } else lapplyfun <- lapply

  # Quiet mode
  if(quiet) {
    catfun <- quietcat
  } else {
    catfun <- cat
  }

  # Start tracking time (after loading packages)
  t.start <- Sys.time()

  # Work with survey outside template loop
  # Creates a wav file for survey if it isn't already a file
  survey <- getClip(survey, name=deparse(substitute(survey)), write.wav=write.wav)

  if(time.source == 'fileinfo') {
     file.time <- file.info(survey)$mtime
     if(is.na(rec.tz)) rec.tz <- format(file.time, format='%Z')
     file.time <- as.POSIXct(format(file.time, tz=rec.tz), tz=rec.tz)
  } else if(time.source == 'filename') {
     survey.short <- strsplit(survey, '/')[[1]][length(strsplit(survey, '/')[[1]])]
     date.time.info <- regmatches(survey.short, regexpr('[0-9]{4}-[0-9]{2}-[0-9]{2}[ _][0-9]{6}[ _][A-Z0-9]{1,7}', survey.short))
     date.time.info <- gsub("_", " ", date.time.info)
     if(length(date.time.info) == 1 && nchar(date.time.info) %in% 19:23)
       file.time <- as.POSIXct(substr(date.time.info, start=1, stop=17), tz=substr(date.time.info, start=19, stop=length(date.time.info)), format='%Y-%m-%d %H%M%S') 
     else {
       warning('time.source was set to \"filename\" but file name does not have date and time info, so using \"fileinfo\" instead')
       file.time <- file.info(survey)$mtime
     }
  } else stop('time.source argument, ', time.source, ' not recognized.')

  survey.path <- survey
  survey <- readClip(survey)

  # score.L is a list for storing results
  score.L <- list()
  survey.data <- list()

  # Loop through templates
  for(i in names(templates@templates)) {
    catfun('\nStarting ', i,'. . .')

    # Working with a single template
    template <- templates@templates[[i]]

    if(i == names(templates@templates)[1] || any(template@wl != wl, template@ovlp != ovlp, template@wn != wn)) {
      catfun('\n\tFourier transform on survey . . .')
      wl <- template@wl
      ovlp <- template@ovlp
      wn <- template@wn
      # Perform Fourier transform on survey
      survey.spec <- spectro2(wave=survey, wl=wl, ovlp=ovlp, wn=wn, ...)
      # NTS arbitrary adjustment to eliminate -Inf
      survey.spec$amp[is.infinite(survey.spec$amp)] <- min(survey.spec$amp[!is.infinite(survey.spec$amp)]) - 10
      frq.bins <- survey.spec$freq
      t.bins <- survey.spec$time
      t.survey <- length(survey@left)/survey@samp.rate
      t.step <- t.bins[2] - t.bins[1]
      frq.step <- frq.bins[2] - frq.bins[1]
      catfun('\n\tContinuing. . .\n')
    }

    # Switch the order of columns in pt.on and pt.off to use them directly for indexing
    pts <- template@pts[, c(2:1, 3)]

    # Throws error if step sizes do not match
    if(!all.equal(template@t.step, t.step, tolerance=t.step/1E4)) {
      stop('For ', i,' time step doesn\'t match survey time step: ', t.step, ' != ', template@t.step, '. Adjust one or the other and try again.')
      ##pts[, 't'] <- round(pts[, 't']*template@t.step/t.step)
      ##if(warn) warning('For ', i,' time step doesn\'t match survey time step: ', t.step, ' != ', template@t.step)
    }
    if(!all.equal(template@frq.step, frq.step, tolerance=frq.step/1E6)) {
      stop(i, ' frequency step does\'t match survey frequency step, ', frq.step, ' != ', template@frq.step, '. Adjust one or the other and try again.')
      #pts[, 'frq'] <- round(pts[, 'frq']*template@frq.step/frq.step)
      #if(warn) warning(i, ' frequency step does\'t match survey frequency step, ', frq.step, ' != ', template@frq.step)
    }

    # Determine the frequency limits from the template points
    frq.lim <- frq.bins[range(pts[, 'frq'])] 

    # Get number of time windows/bins in frequency domain data
    n.t.survey <- length(survey.spec$time)
  
    # Pare down amplitude matrix based on filter frequencies 
    which.frq.bins <- which(survey.spec$freq >= frq.lim[1] & survey.spec$freq <= frq.lim[2])
    amp.survey <- survey.spec$amp[which.frq.bins, ]

    # Shift frq indices in pts. The t indices already start at 1.
    pts[, 'frq'] <- pts[, 'frq'] - min(which.frq.bins) + 1
    n.t.template <- max(pts[, 't'])
    n.frq.template <- max(pts[, 'frq'])

    # Translate pts matrix of indices into a vector index so indexing is faster within the lapplyfun call
    pts.v <- (pts[, 't'] - 1)*n.frq.template + pts[, 'frq']
    amp.template <- pts[, 'amp']
    amp.survey.v <- c(amp.survey)  

    # Create progress bar object if requested
    if(show.prog & !parallel) pb <- txtProgressBar(max=(n.t.survey-n.t.template)*n.frq.template + 1, char='.', width=0, style=3)

    # Perform analysis for each time value (bin) of survey 
    # Starting time value (bin) of correlation window, set up as a list to use mclapply

    c.win.start <- as.list(1:(n.t.survey-n.t.template)*n.frq.template) # Starting position of subset of each survey amp matrix  
    score.survey <- unlist(
      lapplyfun(X=c.win.start, FUN=function(x) 
        {
        if(!parallel && show.prog) setTxtProgressBar(pb, x)
        # Unpack columns of survey amplitude matrix for correlation analysis
        cor(amp.template, amp.survey.v[x + pts.v], method=cor.method, use='complete.obs')  
        }
      )
    )

    # Collect score results and time (center of time bins) in data frame
    score.L[[i]] <- data.frame(
      date.time=file.time + survey.spec$time[1:(n.t.survey-n.t.template)+n.t.template/2] - t.survey, 
      time=survey.spec$time[1:(n.t.survey-n.t.template)+n.t.template/2], 
      score=score.survey
    )
    survey.data[[i]] <- list(amp=survey.spec$amp, t.bins=t.bins, frq.bins=frq.bins)
    catfun('\n\tDone.\n')
  }

  # Calculate total run time
  t.run <- signif(as.numeric(difftime(Sys.time(), t.start, units='secs')), 4)
  time.info <- c(t.exe=as.character(t.run), RTfactor=paste(signif(t.survey/as.numeric(t.run, units='secs'), 4), 'x', sep=''))

  # Return results 
  scores.obj <- new('templateScores', survey.name=survey.path, survey=survey, survey.data=survey.data, templates=templates@templates, scores=score.L, time=time.info)
  return(scores.obj)
}


#########################
# MAIN SCRIPT BODY

# load sound files and selection tables
wavnames = c(	"ZOOM0049_NNK6_LNVX_F_FH_sel.01.ch01.181124.094641.12."
				"ZOOM0049_NNK6_LNVX_F_FH_sel.01.ch01.181124.094641.12..Table.1.selections.txt")

# write files to disk
#writeWave(ZOOM0049_NNK6_LNVX_F_FH_sel.01.ch01.181124.094641.12.,"ZOOM0049_NNK6_LNVX_F_FH_sel.01.ch01.181124.094641.12..wav")
#writeWave(ZZOOM0049_NNK6_LNVX_F_P_sel.03.ch01.181124.094642.05.,"ZOOM0049_NNK6_LNVX_F_P_sel.03.ch01.181124.094642.05..wav")
#writeWave(ZZOOM0049_NNK6_PXYW_M_FH_sel.04.ch01.181124.095747.94.,"ZOOM0049_NNK6_PXYW_M_FH_sel.04.ch01.181124.095747.94..wav")
# writeWave(ZZOOM0049_NNK6_PXYW_M_FH_sel.06.ch01.181124.100005.50.,"ZOOM0049_NNK6_PXYW_M_FH_sel.06.ch01.181124.100005.50..wav")
# writeWave(ZOOM0049_NNK6_PXYW_M_P_sel.05.ch01.181124.095750.38.,"ZOOM0049_NNK6_PXYW_M_P_sel.05.ch01.181124.095750.38..wav")


# use a single wav file as template
wavname="ZOOM0049_NNK6_LNVX_F_FH_sel.01.ch01.181124.094641.12."
selection_table=read.table(	file=sprintf("%s/%s_NNK6_LNVX_F_FH_sel.01.ch01.181124.094641.12..Table.1.selections.txt",path_to_files,substr(wavname,1,8)),
							header=TRUE,
							sep="\t",
							quote="",
							strip.white=TRUE,
							stringsAsFactors=FALSE,
							row.names=NULL)
selection=1
freq_lim = c(selection_table$Low.Freq..Hz.[selection],selection_table$High.Freq..Hz.[selection])/1000
soundt1 = monitoR::makeCorTemplate(	clip=sprintf("%s/%s.wav",path_to_files,wavname), 
									wl = 300, 
									ovl=90, 
									score.cutoff=DETECTION_SCORE_THRESHOLD,
									frq.lim=freq_lim,
									dens=1, 
									name=sprintf("%s_template",wavname))
							

# use Pearson correlation to detect signals
cat(sprintf("Calculating correlations between template and target..\n"))
cscoreszip = corMatch(	survey = sprintf("%s/ZOOM0015_NNT7_NXBW_M_FH_sel.01.ch01.181121.110558.109.wav",path_to_files),
						templates = soundt1, 
						parrallel = T, 
						show.prog = F, 
						time.source = "fileinfo", 
						cor.method = "pearson",
						warn = F, 
						write.wav = T)
						
cat(sprintf("Looking for peaks..\n"))
cdetects_soundt1 = monitoR::findPeaks(cscoreszip, parallel = FALSE)
print(cdetects_soundt1@peaks)

#View results
cat(sprintf("Plotting results..\n"))
plot_peak_search(	cdetects_soundt1, 
					scorelim = c(0,1),
					hit.marker = "lines", 
					flim=c(0,18), 
					legend=FALSE)



cat("Done")