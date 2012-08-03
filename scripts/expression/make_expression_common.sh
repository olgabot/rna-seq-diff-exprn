#!/bin/sh -x

# make_expression_common.sh
# Created by: Olga Botvinnik on 07-31-2012
# Purpose of this script: Make variables needed for expression
# analysis handily available in an R script.
# Input(s): common_vars.sh
# Output(s): creates expression_common.R

COMMON_VARS=$1
EXPRESSION_COMMON=$2

source $COMMON_VARS

# Replace this text and create a new file
sed -e's/CONDITIONS\_FILE/"$COND_WITHOUT_COMMENTS"' \
	<$EXPRN_SCRIPTS_DIR/expression_common.R \
	>$EXPRESSION_COMMON


# --- Assign TREATMENT_COLORS_SELECTED ---#
# This gives the number of total group colors, ie the number
# of total treatment groups. For example, if you have an
# untreated and treated group, you will have two colors.

# Replace the text "in-place," do not create a new file via
# `sed -i`
MAKE_PINKOGRAM_R=$EXPRN_SCRIPTS_DIR/make.pinkogram.R
sed -i '' -e's/MAKE_PINKOGRAM_R/"$MAKE_PINKOGRAM_R"' \
	$EXPRESSION_COMMON
TREATMENT_COLORS=( "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666" )

# The number of colors we need to use to differentiate the
# different treatment groups is the number of groups in the
# groups array.
TREATMENT_COLORS_SELECTED_ARRAY=( ${TREATMENT_COLORS[@]:1:${#GROUPS_ARRAY[@]}} )
TREATMENT_COLORS_SELECTED=`echo \
	${TREATMENT_COLORS_SELECTED_ARRAY[@]} | 
	sed 's/^/"/' | sed 's/$/"/' | sed 's/ /", "/'`
	# Maybe try: sed 's/ /",\\n"/' for the last one

sed -i '' 's/TREATMENT_COLORS_SELECTED/"$TREATMENT_COLORS_SELECTED"' \
	$EXPRESSION_COMMON


# --- Assign INDIVIDUAL_TREATMENT_{COLORS,NAMES} --- #
# This will assign a treatment group color and name to each 
# individual sample, for example if you have two untreated samples 
# and four treated samples, then this will first create an array
# with two of one color, then four of the second color, where the
# colors are specified in the TREATMENT_COLORS_SELECTED_ARRAY.

# Now we need to distribute those treatment colors to the
# individuals in the treatment groups
for (( i = 0; i < ${#TREATMENT_GROUPS_ALL[@]}; i++ )); do
    # Get the name of the group for this sample index
    GROUP=${TREATMENT_GROUPS_ALL[$i]}

    # Find the index of this group in the array of unique group names
    IND=`IndexOf $GROUP ${GROUPS_ARRAY[@]}`

    # Assign this ID the IND'th color in the list of possible colors
    INDIVIDUAL_TREATMENT_COLORS_ARRAY[$i]=${TREATMENT_COLORS_SELECTED_ARRAY[$IND]}
    INDIVIDUAL_TREATMENT_NAMES_ARRAY[$i]=${GROUPS_ARRAY[$IND]}
done

INDIVIDUAL_TREATMENT_COLORS=`echo ${INDIVIDUAL_TREATMENT_COLORS_ARRAY[@]}\
	| sed 's/^/"/' | sed 's/$/"/' | sed 's/ /", "/'`

sed -i '' -e's/INDIVIDUAL_TREATMENT_COLORS/$INDIVIDUAL_TREATMENT_COLORS/' \
	$EXPRESSION_COMMON

# --- Assign INDIVIDUAL_TREATMENT_NAMES --- #
INDIVIDUAL_TREATMENT_NAMES=`echo ${INDIVIDUAL_TREATMENT_NAMES_ARRAY[@]}\
	| sed 's/^/"/' | sed 's/$/"/' | sed 's/ /", "/'`

sed -i '' -e's/INDIVIDUAL_TREATMENT_NAMES/$INDIVIDUAL_TREATMENT_NAMES/' \
	$EXPRESSION_COMMON

# --- Assign GROUP_TREATMENT_COLORS --- #
# If the number of groups per treatment group is greater than 0,
# then assign colors to each of the groups within the treatment
# groups. For example, if you had 6 untreated and 8 treated samples,
# and you specified NUM_GROUPS=2 then you would have two groups of
# 3 for the untreated samples, and two groups of 4 for the treated
# samples, so a total of four groups, two untreated and two treated.
# So then this would create an array of four colors, two for the
# untreated samples and two for the treated.

if [[ $NUM_GROUPS > 0 ]]; then
	GROUP_IDS_ARRAY=( echo $GROUP_IDS | tr , ' ' )
	for i in (( i = 0; i < ${#GROUP_IDS_ARRAY[@]}; i++ )); do
		# Get this group ID
		G=${GROUP_IDS_ARRAY[$i]}

		# Remove the "_groupNofN," e.g. "_group1of2" suffix
		# so we have just the raw group type and we can match
		# it to the group names

		# This `-E` specifies extended (modern) regular expressions
		# which may have backwards compatibility issues.
		echo $G | sed -E 's/_group[[:digit:]]+of[[:digit:]]+//'

		# Find the index of this group in the array of unique 
		# group names
		IND=`IndexOf $G ${GROUPS_ARRAY[@]}`

		# Assign this ID the IND'th color in the list of possible colors
		GROUP_TREATMENT_COLORS_ARRAY[$i]=${TREATMENT_COLORS_SELECTED_ARRAY[$IND]}
		GROUP_TREATMENT_NAMES_ARRAY[$i]=${GROUP_IDS_ARRAY[$IND]}
	done
	GROUP_TREATMENT_COLORS=`echo ${GROUP_TREATMENT_COLORS_ARRAY[@]}
		| sed 's/^/"/' | sed 's/$/"/' | sed 's/ /", "/'`
	sed -i '' 's/"GROUP_TREATMENT_COLORS"/$GROUP_TREATMENT_COLORS/' \
		$EXPRESSION_COMMON 

	GROUP_TREATMENT_NAMES=`echo ${GROUP_TREATMENT_NAMES_ARRAY[@]}
		| sed 's/^/"/' | sed 's/$/"/' | sed 's/ /", "/'`
	sed -i '' 's/"GROUP_TREATMENT_NAMES"/$GROUP_TREATMENT_NAMES/' \
		$EXPRESSION_COMMON 
fi