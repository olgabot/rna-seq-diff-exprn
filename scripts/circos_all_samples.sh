#!/bin/sh

# Author: Olga Botvinnik (olga.botvinnik@gmail.com)
# Date: 21 June 2012
# Purpose of this script: Finds genome-wide counts
# of a BED file using both bedtools and HTSeq
# (PS: in my experience, genome-wide coverage is EXACTLY
# the same regardless of whether you use bedtools or HTSeq)


# In this file, COVERAGE = GENOME COVERAGE
# TODO: Need to consider male/female about whether to \
# include ChrY or not

# OUT_DIR=$1
# BAM=$1
# SAM_SORTED=$2
# # BED=$4
# GENDER=$3  # male or female
# ID=$4
COMMON_VARS=$1

source $COMMON_VARS

####### BEGIN create circos files & execute circos #######
####### create circos files #######

THIS_CIRCOS_OUT=$CIRCOS_OUT_DIR/all_samples



if [[ ! -d $THIS_CIRCOS_OUT ]]; then
    # If this directory doesn't yet exist, create it
    mkdir -p $THIS_CIRCOS_OUT
fi

# --- Removing this directory for debugging for now --- #
rm -r $THIS_CIRCOS_OUT/*

# CIRCOS_TEMPLATE_DIR=$SCRIPTS_DIR/circos-templates
cp -r $CIRCOS_TEMPLATE_DIR/all_samples/* \
    $THIS_CIRCOS_OUT

# TODO: Need to check if ALL genders are female. If ANY SINGLE
# sample is male, then need to show the Y chromosomes
if [[ $GENDER == female ]] ;
# Gender is female --> remove ChrY from circos plot
# by uncommenting the 'chromosomes = -chrY line
then
    sed -i '' -e's/^#chromosomes/chromosomes/' \
	$THIS_CIRCOS_OUT/etc/circos.conf
fi

# ----- Replace keywords in template file with variables ----- #

# --- Begin setting global parameters --- #
# Make the upper limit of the plot be the mean + (2 * std dev)
UPPER_LIMIT=`echo $BEDTOOLS_STATS | tr '=' " " | \\
    awk -F' ' '{ mu=$5 ; sigma=$8 ; print mu+2*sigma}'`

# Divide UPPER_LIMIT by 100 to get minimum value change
# Last awk command is a floor function
HISTOGRAM_MASTER=$THIS_CIRCOS_OUT/etc/histogram.master.conf
MIN_VALUE_CHANGE=`echo $UPPER_LIMIT | awk '{print $1/100}' | \\
    awk -F. '{print $1}'`

# Average all the upper limits to get a universal one
UPPER_LIMIT_STATS=`echo $UPPER_LIMITS | sed 's/^,//' | \\
    tr , '\n' | awk -F' ' \\
    ' {sumsq+=$1*$1; sum+=$1<0?-$1:$1} END \\
    {print int(sqrt(sumsq/NR - (sum/NR)**2)); print int(sum/NR) } '`

# Second field is the average upper limit/max value
UPPER_LIMIT=`echo $UPPER_LIMIT_STATS | \\
    awk -F' ' ' { print $2 } '`

# First field is standard deviation of upper limit/max values
MIN_VALUE_CHANGE=`echo $UPPER_LIMIT_STATS | \\
    awk -F' ' ' { print $1 } '`


sed -i '' -e"s:UPPER_LIMIT:${UPPER_LIMIT}:" \
    $HISTOGRAM_MASTER

sed -i '' -e"s:KARYOTYPE:${KARYOTYPE}:" \
    $THIS_CIRCOS_OUT/etc/circos.conf

sed -i '' -e"s:MIN_VALUE_CHANGE:${MIN_VALUE_CHANGE}:" \
    $HISTOGRAM_MASTER

# TODO: Need to check if GENE_DENSITY and GC_CONTENT
# were actually provided
sed -i '' -e"s:GENE_DENSITY:${GENE_DENSITY}:" \
    $HISTOGRAM_MASTER

sed -i '' -e"s:GC_CONTENT:${GC_CONTENT}:" \
    $HISTOGRAM_MASTER
# --- END setting global parameters --- #

HISTOGRAM_TEMPLATE=$THIS_CIRCOS_OUT/etc/histogram.single.plot.conf

# 0.99 = where the plot starts, just off of the karyotype
PREV_INNER_RADIUS=0.99

NUM_PLOTS=$END

if [[ $GENE_DENSITY_FILE != '' ]]; then
    NUM_PLOTS=`echo "$NUM_PLOTS+1" | bc`
fi

if [[ $GC_CONTENT_FILE != '' ]]; then
    NUM_PLOTS=`echo "$NUM_PLOTS+1" | bc`
fi

# RING_THICKNESS is how thick each plotted ring histogram is
RING_THICKNESS=`echo "scale=3; 0.8 / $NUM_PLOTS" | bc`

REVERSE_CIRCOS_ORDER=TRUE

if [[ $REVERSE_CIRCOS_ORDER == 'TRUE' ]]; then
    FIRST_VALUE=`echo "$END-1" | bc`
    BEGIN_LOOP="i = $FIRST_VALUE"
    END_LOOP='i >= 0'
    INCREMENT=--i
else
    BEGIN_LOOP='i = 0'
    END_LOOP="i < $END"
    INCREMENT=++i
fi

for (( $BEGIN_LOOP ; $END_LOOP ; $INCREMENT )); do
    ID=${ID_ARRAY[$i]}
    GROUP=`grep $ID $COND_WITHOUT_COMMENTS | cut -f3`
    HISTOGRAM_GROUP=$THIS_CIRCOS_OUT/etc/histogram.$GROUP.conf

    SAMPLE_FILE=$BEDTOOLS_DIR/$ID/$COVERAGE_BEDTOOLS_PREFIX.circos

    OUTER_RADIUS=`echo "scale=3; $PREV_INNER_RADIUS - 0.005" | bc`
    INNER_RADIUS=`echo "scale=3; \\
        $OUTER_RADIUS - $RING_THICKNESS" | bc`

    SAMPLE_COLOR=${CIRCOS_ID_COLOR_ARRAY[$i]}

    # --- BEGIN setting sample-specific parameters --- #
    sed -e"s:SAMPLE_FILE:$SAMPLE_FILE:" \
        $HISTOGRAM_TEMPLATE \
        >> $HISTOGRAM_GROUP
    sed -i '' -e"s:OUTER_RADIUS:$OUTER_RADIUS:"\
        $HISTOGRAM_GROUP
    sed -i '' -e"s:INNER_RADIUS:$INNER_RADIUS:" $HISTOGRAM_GROUP
    sed -i '' -e"s:SAMPLE_COLOR:$SAMPLE_COLOR:" $HISTOGRAM_GROUP
    # --- END setting sample-specific parameters --- #

    PREV_INNER_RADIUS=$INNER_RADIUS
done

if [[ $GENE_DENSITY_FILE != '' ]]; then
    OUTER_RADIUS=`echo "scale=3; $PREV_INNER_RADIUS - 0.005" | bc`
    INNER_RADIUS=`echo "scale=3; \\
        $OUTER_RADIUS - $RING_THICKNESS" | bc`

    HISTOGRAM_GENE_DENSITY=$THIS_CIRCOS_OUT/etc/histogram.gene.density.conf
    
    sed -i '' -e"s:GENE_DENSITY_FILE:$GENE_DENSITY_FILE:" \
        $HISTOGRAM_GENE_DENSITY
    sed -i '' -e"s:OUTER_RADIUS:$OUTER_RADIUS:"\
        $HISTOGRAM_GENE_DENSITY
    sed -i '' -e"s:INNER_RADIUS:$INNER_RADIUS:" $HISTOGRAM_GENE_DENSITY
    sed -i '' -e"s:GENE_DENSITY_COLOR:$GENE_DENSITY_COLOR:" \
        $HISTOGRAM_GENE_DENSITY
    sed -i '' -e"s:MAX_GENE_DENSITY:$MAX_GENE_DENSITY:" \
        $HISTOGRAM_GENE_DENSITY
    sed -i '' -e"s:GENE_DENSITY_MIN_VALUE_CHANGE:$GENE_DENSITY_MIN_VALUE_CHANGE:" \
        $HISTOGRAM_GENE_DENSITY

    PREV_INNER_RADIUS=$INNER_RADIUS

    # Add this file into the master histogram file
    sed -i '' "
/max = $UPPER_LIMIT/ a\\
<<include histogram.gene.density.conf>>
" \
        $HISTOGRAM_MASTER
fi

if [[ $GC_CONTENT_FILE != '' ]]; then
    OUTER_RADIUS=`echo "scale=3; $PREV_INNER_RADIUS - 0.005" | bc`
    INNER_RADIUS=`echo "scale=3; \\
        $OUTER_RADIUS - $RING_THICKNESS" | bc`

    HISTOGRAM_GC_CONTENT=$THIS_CIRCOS_OUT/etc/histogram.gc.content.conf
    
    sed -i '' -e"s:GC_CONTENT_FILE:$GC_CONTENT_FILE:" \
        $HISTOGRAM_GC_CONTENT
    sed -i '' -e"s:OUTER_RADIUS:$OUTER_RADIUS:"\
        $HISTOGRAM_GC_CONTENT
    sed -i '' -e"s:INNER_RADIUS:$INNER_RADIUS:" $HISTOGRAM_GC_CONTENT
    sed -i '' -e"s:GC_CONTENT_COLOR:$GC_CONTENT_COLOR:" \
        $HISTOGRAM_GC_CONTENT

    # Add this file into the master histogram file
    sed -i '' "
/max = $UPPER_LIMIT/ a\\
<<include histogram.gc.content.conf>>
" \
        $HISTOGRAM_MASTER
fi

# Initialize the histogram file for all the groups and
# add the include
for GROUP in `echo $TREATMENT_GROUPS | tr , ' '`; do
    HISTOGRAM_GROUP=$THIS_CIRCOS_OUT/etc/histogram.$GROUP.conf

    # Add an informative message at the top of this histogram file
    echo "# This is the Circos histogram configuration file for
# the samples of group type: $GROUP" \
        | cat - $HISTOGRAM_GROUP > $HISTOGRAM_GROUP.temp
    mv $HISTOGRAM_GROUP.temp $HISTOGRAM_GROUP

    # Include the histogram file for this group in the master
    # histogram file
    sed -i '' "
/max = $UPPER_LIMIT/ a\\
<<include histogram.$GROUP.conf>>
" \
        $HISTOGRAM_MASTER
done


####### execute circos #######
$CIRCOS_BIN -conf $THIS_CIRCOS_OUT/etc/circos.conf \
    -outputdir $THIS_CIRCOS_OUT
####### end create circos files & execute circos #######
