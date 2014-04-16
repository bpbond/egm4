# egm4.R
# R script (written using 3.0.3) to process EGM-4 data
# BBL April 2014

# Important variable definitions, esp. data source & destination

SCRIPTNAME		<- "egm4.R"
INPUT_DIR		<- "sampledata/"
OUTPUT_DIR		<- "outputs/"
LOG_DIR			<- "logs/"
SEPARATOR		<- "-------------------"

# This is a critical variable, since seconds don't appear in the EGM4 file
MEAS_INTERVAL	<- 10

# Support functions and common definitions

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
# Return matrix of memory consumption
object.sizes <- function()
{
    return( rev( sort( sapply( ls( envir=.GlobalEnv ), function( object.name ) 
        object.size( get( object.name ) ) ) ) ) )
}
 
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
} # saveplot

# -----------------------------------------------------------------------------
# Open a csv file and return data
read_csv <- function( fn, datadir="." ) {
	fqfn <- paste( datadir, fn, sep="/" )
	printlog( "Opening", fqfn )
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
	fqfn <- paste0( INPUT_DIR, fn )
	printlog( "Reading", fqfn )
	stopifnot( file.exists( fqfn ) )
	d <- read.table( fqfn, comment.char=";", sep="\t" )
	printdims( d )
	names( d ) <- c( "Plot", "RecNo", "Day", "Month", "Hour", "Min", "CO2_Ref", "mb_Ref",
		 "mbR_Temp", "Input_A", "Input_B", "Input_C", "Input_D", "Input_E", "Input_F", 
		 "Input_G", "Input_H", "ATMP", "Probe Type" )
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

	return( d )
} # read_egmfile

# -----------------------------------------------------------------------------
# compute fluxes
compute_flux <- function( d ) {

	m <- lm( CO2_Ref ~ Sec, data=d )
	Resp_raw <- as.numeric( coef( m )[ 2 ] )	# i.e. the slope
	
	# TODO: height?
	Height <- 1
	# TODO: dry mass?
	Dry.mass <- 1
	warning( "Currently using constant mass/area for all plots." )
	
	# We want to convert raw respiration (d[CO2]/dt) to a flux using
	# A = dC/dt * V/S * Pa/RT (e.g. Steduto et al. 2002), where
	# 	A is CO2 flux (umol/m2/s)
	#	dC/dt is raw respiration as above (mole fraction/s)
	# 	V is total chamber volume (m3)
	#		...we are correcting for varying headspaces in the cores
	#	S is ground surface area (m2)
	#		...but we're computing per kg of soil, so using dry mass instead
	#	Pa is atmospheric pressure (kPa)
	#	R is universal gas constant (8.3 x 10-3 m-3 kPa mol-1 K-1)
	#	T is air temperature (K)

	# Note this is currently written for a lab incubation, computing mass-specific
	# respiration. TODO: change to mass or area basis, user's choice.
	
	sleeve_diam <- 3.5			# diameter, cm
	sleeve_ht	<- 15.2 + 2		# height, cm; extra 2 is for cap
	egm4_vol	<- 9			# internal system volume, cm3
	lines_vol 	<- ( 1/8 * 2.54 / 2 ) ^2 * 122 * 2	# two 122-cm 1/8" lines, cm3
	S 			<- ( sleeve_diam / 2 ) * pi	# note cm2, not m2!
	sleeve_vol <- ( sleeve_ht - Height ) * S
	V	<- ( egm4_vol + lines_vol + sleeve_vol ) / 100^3		# m3
	Pa 			<- 101						# kPa				(Richland is ~120 m asl)
	R 			<- 8.3e-3					# m-3 kPa mol-1 K-1

	Tair <- mean( d$Input_C )		# assumes EGM temperature probe connected
	
	# Calculate mass-corrected respiration, umol/g soil/s
	Resp_mass <- Resp_raw * V/Dry.mass * Pa/( R*( 273.1+Tair ) )

	# Convert from umol/g soil/s to mgC/kg soil/day
	flux <- Resp_mass / 1e6 * 12 * 1000 * 1000 * 60 * 60 * 24

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

sink( paste0( LOG_DIR, SCRIPTNAME, ".txt" ), split=T )

printlog( "Welcome to", SCRIPTNAME )

loadlibs( c( "ggplot2", "reshape2", "plyr" ) )
theme_set( theme_bw() )

alldata <- data.frame()
filelist <- list.files( path=INPUT_DIR, pattern="*.dat" )
for( fn in filelist ) {
	printlog( SEPARATOR )
	alldata <- rbind( alldata, read_egmfile( fn ) )
}

printlog( SEPARATOR )
printlog( "All done reading data." )
printdims( alldata )

printlog( "Merging respiration data with dry mass data..." )
# TODO

printlog( "Computing fluxes..." )
fluxes <- ddply( alldata, .( filename, Plot ), .fun=compute_flux )

print( summary( fluxes ) )

p <- ggplot( fluxes, aes( Plot, flux ) ) + geom_point()
print( p )
saveplot( "flux_summary" )

printlog( SEPARATOR )
printlog( "Saving flux data..." )
savedata( alldata )
printlog( "Saving flux data..." )
savedata( fluxes )

printlog( "All done with", SCRIPTNAME )
print( sessionInfo() )
sink()
