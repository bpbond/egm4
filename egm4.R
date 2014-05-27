# egm4.R
# R script (written using 3.0.3) to process EGM-4 data
# BBL April 2014

# To use:
# 	1. Set SYSTEM_VOLUME and CHAMBER_AREA variables below
#	2. Set MEAS_INTERVAL variable
#	3. If using plot-specific ancillary data, set PLOTDATA variable
#			and make sure file is in correct format
#	4. Set INPUT_DIR to location of EGM4 output files
#	5. Check the flux calculations in compute_flux() for your use case
#	6. source() this file.
#			start with a simple test case, and carefully check log file and outputs!

# Important variable definitions, esp. data source & destination
SCRIPTNAME		<- "egm4.R"
INPUT_DIR		<- "sampledata"  # directory names shouldn't end with / (Windows)
OUTPUT_DIR		<- "outputs"
LOG_DIR			<- "logs"

# The optional PLOTDATA file *must* have a 'Plot' field
# It *may* have 'Mass' and/or 'Area' fields, which will be divided into the computed flux
#	- if present, the 'Area' field will override CHAMBER_AREA below
# It *may* have a 'Volume' field, cm3, which will be added to SYSTEM_VOLUME for each plot
PLOTDATA		<- "sampledata/plotdata.csv"

SYSTEM_VOLUME	<- 15.15 / 100^3		# m^3
CHAMBER_AREA	<- 5.5					# cm^2

SEPARATOR		<- "-------------------"

# This is a critical variable, since seconds don't appear in the EGM4 file
# Time between successive measurements (s)
MEAS_INTERVAL	<- 10


# -----------------------------------------------------------------------------
# Time-stamped output function
printlog <- function( msg="", ..., ts=TRUE, cr=TRUE ) {
	if( ts ) cat( date(), " " )
	cat( msg, ... )
	if( cr ) cat( "\n")
} # printlog

# -----------------------------------------------------------------------------
# Print dimensions of data frame
printdims <- function( d, dname=deparse( substitute( d ) ) ) {
	stopifnot( is.data.frame( d ) )
	printlog( dname, "rows =", nrow( d ), "cols =", ncol( d ) )
} # printdims

# -----------------------------------------------------------------------------
# Save a ggplot figure
saveplot <- function( pname, p=last_plot(), ptype=".pdf" ) {
	stopifnot( file.exists( OUTPUT_DIR ) )
	fn <- paste0( OUTPUT_DIR, "/", pname, ptype )
	printlog( "Saving", fn )
	ggsave( fn, p )
} # saveplot

# -----------------------------------------------------------------------------
# Save a data frame
savedata <- function( df, extension=".csv" ) {
	stopifnot( file.exists( OUTPUT_DIR ) )
	fn <- paste0( OUTPUT_DIR, "/", deparse( substitute( df ) ), extension )
	printlog( "Saving", fn )
	write.csv( df, fn, row.names=F )
} # savedata

# -----------------------------------------------------------------------------
# Open a csv file and return data
read_csv <- function( fn, datadir="." ) {
	fqfn <- paste( datadir, fn, sep="/" )
	printlog( "Opening", fqfn )
	stopifnot( file.exists( fqfn ) )
	read.csv( fqfn, stringsAsFactors=F )
} # read_csv

# -----------------------------------------------------------------------------
# Load requested libraries
loadlibs <- function( liblist ) {
	printlog( "Loading libraries..." )
	loadedlibs <- vector()
	for( lib in liblist ) {
		printlog( "Loading", lib )
		loadedlibs[ lib ] <- require( lib, character.only=T )
		if( !loadedlibs[ lib ] )
			warning( "this package is not installed!" )
	}
	invisible( loadedlibs )
} # loadlibs

# -----------------------------------------------------------------------------
# read a process a single EGM4 output file, returning data frame
read_egmfile <- function( fn ) {
	fqfn <- paste0( INPUT_DIR, "/", fn )
	printlog( "Reading", fqfn )
	stopifnot( file.exists( fqfn ) )
	d <- read.table( fqfn, comment.char=";", sep="\t" )
	printdims( d )
	# TODO: check that the file terminates with a "Logging stopped after X records" line
	
	# TODO: get these columns directly from file
	egm4names <- c( "Plot", "RecNo", "Day", "Month", "Hour", "Min", "CO2_Ref", "mb_Ref",
		 "mbR_Temp", "Input_A", "Input_B", "Input_C", "Input_D", "Input_E", "Input_F", 
		 "Input_G", "Input_H", "ATMP", "Probe Type" )
	stopifnot( length( names( d ) )==length( egm4names ) )
	names( d ) <- egm4names
	
	# Add ancillary data
	d$filename <- fn
	
	# Seconds don't appear in the EGM-4 output, which is weird and sucks. We fill them in
	# based on MEAS_INTERVAL (e.g. this is 5 if measurements were made every 5 secs)
	printlog( "Adding seconds (interval =", MEAS_INTERVAL, ")" )
	d <- ddply( d, .( Plot ), mutate, Sec=seq( from=0, length.out=length( Plot ), by=MEAS_INTERVAL ) )
	
	# Quality control
	# TODO: improve this
	printlog( "Computing CO2~Time R2 values for quality control..." )
	mods <- dlply( d, .( Plot ), lm, formula = CO2_Ref ~ Sec )
	r2 <- ldply( mods, .fun=function( x ){ round( summary( x )$r.squared, 2 ) } )
	names( r2 ) <- c( "Plot", "R2" )
	r2 <- r2[ order( r2$R2 ), ]
	print( r2 )
    
#	p <- qplot( factor( Plot ), R2, geom="bar", data=r2, stat="identity" )
#    saveplot( "r2_values", p ) 
    
	return( d )
} # read_egmfile


# -----------------------------------------------------------------------------
# read plot data, if it exists, and merge with EGM4 data
read_plotdata <- function( fn=PLOTDATA ) {
	d <- NULL
	if( file.exists( fn ) ) {
		d <- read_csv( fn )
		if( any( names( d )=="Plot" ) ) {
			printlog( "Plot data read OK" )
		} else {
			printlog( "Plot data file read, but no 'Plot' field!" )
			warning( "No plot field!" )
		}
	} else {
		printlog( "Plot data file", fn, "not found" )
	}
	
	return( d )
} # read_plotdata

# -----------------------------------------------------------------------------
# compute fluxes
compute_flux <- function( d ) {

	m <- lm( CO2_Ref ~ Sec, data=d )
	resp_raw <- as.numeric( coef( m )[ 2 ] )	# i.e. the slope
	
	# We want to convert raw respiration (d[CO2]/dt) to a flux using
	# A = dC/dt * V/S * Pa/RT (e.g. Steduto et al. 2002), where
	# 	A is CO2 flux (umol/m2/s)
	#	dC/dt is raw respiration as above (mole fraction/s)
	# 	V is total chamber volume (m3)
	#		...correcting for varying headspaces in the cores, if applicable
	#	S is ground surface area (m2), if applicable
	# 	M is sample dry mass (g), if applicable
	#	Pa is atmospheric pressure (kPa)
	#	R is universal gas constant (8.3 x 10-3 m-3 kPa mol-1 K-1)
	#	T is air temperature (K)

	S 			<- CHAMBER_AREA		# note cm2, not m2!
	if( any( names( d )=="Area" ) ) {
		 S <- mean( d$Area )
	}
	V 			<- SYSTEM_VOLUME
	if( any( names( d )=="Volume" ) ) {
		V <- V + mean( d$Volume )
	}
	M 			<- 1.0
	if( any( names( d )=="Mass" ) ) {
		 M <- mean( d$Mass )
	}
	
	Pa 			<- 101						# kPa
	R 			<- 8.3e-3					# m-3 kPa mol-1 K-1

	Tair <- mean( d$Input_C )		# assumes EGM temperature probe connected
	
	# Calculate mass- (or area-) corrected respiration, umol/g soil/s or umol/cm2/s
	resp_corrected <- resp_raw * V/S/M * Pa/( R*( 273.1+Tair ) )

	# Convert from umol/g soil/s to mgC/kg soil/day or whatever
	# NOTE: you probably want to change this line for your specific setup
	flux <- resp_corrected / 1e6 * 12 * 1000 * 1000 * 60 * 60 * 24
printlog(nrow(d),length(flux))
	return( c( Tair=Tair, V=V, S=S, Day=mean( d$Day ), Month=mean( d$Month ), N=nrow( d ), flux=flux ) )
}

# ==============================================================================
# Main

if( !file.exists( OUTPUT_DIR ) ) {
	printlog( "Creating", OUTPUT_DIR )
	dir.create( OUTPUT_DIR )
}
if( !file.exists( LOG_DIR ) ) {
	printlog( "Creating", LOG_DIR )
	dir.create( LOG_DIR )
}

sink( paste( LOG_DIR, paste0( SCRIPTNAME, ".txt" ), sep="/" ), split=T )

printlog( "Welcome to", SCRIPTNAME )

loadlibs( c( "ggplot2", "reshape2", "plyr" ) )
theme_set( theme_bw() )

alldata <- data.frame()
filelist <- list.files( path=INPUT_DIR, pattern="dat$", recursive=T )
printlog( "We have", length( filelist ), "files to process" )
for( fn in filelist ) {
	printlog( SEPARATOR )
	alldata <- rbind( alldata, read_egmfile( fn ) )
}

printlog( SEPARATOR )
printlog( "All done reading data." )
printdims( alldata )

plotdata <- read_plotdata()
if( !is.null( plotdata ) ) {
	printlog( "Merging respiration data with dry mass data..." )
	alldata <- merge( plotdata, alldata )
}

printlog( "Computing fluxes..." )
fluxes <- ddply( alldata, .( filename, Plot ), .fun=compute_flux )
fluxes$Plot <- as.factor( fluxes$Plot )

print( summary( fluxes ) )

p <- ggplot( fluxes, aes( Day, flux, group=Plot, colour=Plot ) ) + geom_point() + geom_line()
print( p )
saveplot( "flux_summary" )

p <- ggplot( fluxes, aes( Day, flux ) ) + geom_line() + facet_wrap( ~Plot, nrow=6, ncol=6 )
print( p )
saveplot( "flux_summary2" )

printlog( SEPARATOR )
printlog( "Saving flux data..." )
savedata( alldata )
printlog( "Saving flux data..." )
savedata( fluxes )

printlog( "All done with", SCRIPTNAME )
print( sessionInfo() )
sink()
