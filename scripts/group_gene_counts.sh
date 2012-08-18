#!/bin/sh -x

COMMON_VARS=$1
source $COMMON_VARS

echo "\n# Variables from calculating group gene counts" | \
	cat - >> $COMMON_VARS

GROUP_IDS=''

for GROUP in `echo $TREATMENT_GROUPS | tr , ' '`; do
	echo "finding gene counts for:" $GROUP
	# THIS_GROUP_DIR=$TREATMENT_GROUPS_DIR/$GROUP
	# echo "THIS_GROUP_DIR:" $THIS_GROUP_DIR

	echo "finding samples...."
	SAMPLES=`grep $GROUP $COND_WITHOUT_COMMENTS | cut -f2 | tr "\n" " "`
	echo "samples:" $SAMPLES
	declare -a SAMPLE_ARRAY=( `echo $SAMPLES` )

	# Find the index of the first time this treatment
	# group appears in TREATMENT_GROUPS_ALL
	FIRST_INDEX=`IndexOf $GROUP ${TREATMENT_GROUPS_ALL[@]}`

	# for (( i = 0 ; i < ${#GROUPS_ARRAY[@]} ; i++ ));  do
	  
	# done

	MIN_SAMPLES_PER_GROUP=`echo ${#SAMPLE_ARRAY[@]} $NUM_GROUPS | \\
		awk -F' ' '{ print $1/$2 }' | awk -F. '{print $1}'`

	for ((i=1;i<=$NUM_GROUPS;++i)); do
		GROUP_i=$GROUP\_group$i\of$NUM_GROUPS
		# GROUP_IDS=$GROUP_IDS,$GROUP_i

		GROUP_i_BAM_DIR=$TREATMENT_GROUPS_DIR/$GROUP_i

		# If the directory holding the BAM files doesn't exist yet,
		# create it.
		if [[ ! -d $GROUP_i_BAM_DIR ]] ; then
			mkdir -p $GROUP_i_BAM_DIR
		fi

		# Get the beginning index of the range of samples we're using
		ind0=`echo "$FIRST_INDEX+($i-1)*$MIN_SAMPLES_PER_GROUP" | \\
			bc`

		# Get the number of items for the range of samples we're using
		if [[ $i -lt $NUM_GROUPS ]] ; then
			ind1=$MIN_SAMPLES_PER_GROUP
		else
			# This only happens at the end of making groups for
			# a particular treatment type, i.e. only once.
			ind1=`echo "${#SAMPLE_ARRAY[@]}-($i-1)*$MIN_SAMPLES_PER_GROUP" | bc`
		fi


		# Get the sample names
		GROUP_i_SAMPLES=`echo ${ID_ARRAY[@]:ind0:ind1} | \\
			sed 's/ $//'`
		BAM_IN=`echo ${BAM_PREFIX_ARRAY[@]:ind0:ind1} | \\
			sed 's/ /.bam /g' | sed 's/$/.bam/' `
		BAM_IN_NEWLINE=`echo $BAM_IN | tr ' ' "\n"`

		# GROUP_i_NUM_SAMPLES=`awk -v i=ind0 -v j=ind1 ' { print j-i }'`

		if [[ $ind1 -eq 1 ]]; then
			echo "WARNING: $GROUP_i contains only one sample, \
			$GROUP_i_SAMPLES." #\
			# Not creating a group from this, since you can look \
			# up the gene expression counts of this individual sample, \
			# located in $BEDTOOLS_DIR/$GROUP_i_SAMPLES and \
			# $HTSEQ_DIR/$GROUP_i_SAMPLES."
			# continue
		fi

		GROUP_IDS=$GROUP_IDS,$GROUP_i

		# Make a readme file into the folder where the merged BAM 
		# file will go.
		echo "This group of samples, $GROUP_i, contains these samples: \\
		$GROUP_i_SAMPLES\nCreated from these BAM files:\\
		$BAM_IN_NEWLINE" > $GROUP_i_BAM_DIR/README.txt

		# Use the first Gender/Strandedness in the array:
		GROUP_i_GENDER=`echo ${GENDER_ARRAY[ind0]}`
		GROUP_i_STRAND=`echo ${STRAND_ARRAY[ind0]}`
		# for (( i = ind0+1 ; i <= ind1; i++ )); do
		# 	if [[  ]]; then
		# 		#statements
		# 	fi
		# done

		BAM_OUT_PREFIX=$GROUP_i_BAM_DIR/$GROUP_i
		BAM_OUT=$BAM_OUT_PREFIX.bam

		if [[ $ind1 -eq 1 ]]; then
			cp $BAM_IN $BAM_OUT
			cp -R $BEDTOOLS_DIR/$GROUP_i_SAMPLES/ $BEDTOOLS_DIR/$GROUP_i
			cp -R $HTSEQ_DIR/$GROUP_i_SAMPLES/ $HTSEQ_DIR/$GROUP_i			
		else
			if [[ ! -e $BAM_OUT ]] ; then
				# echo "merging:" $BAM1_IN " >" $BAM1_OUT
				samtools merge $BAM_OUT $BAM_IN
			fi
			# $SCRIPTS_DIR/gene_counts.sh \
			# 	$BAM_OUT_PREFIX \
			# 	$GROUP_i_GENDER \
			# 	$GROUP_i \
			# 	$GROUP_i_STRAND \
			# 	$COMMON_VARS
		fi

		# exit

			# 2>$GROUP_i_BAM_DIR/gene-counts.sh.err
	done
done

GROUP_IDS=`echo $GROUP_IDS | sed 's/^,//'`
echo "GROUP_IDS='$GROUP_IDS'" | cat - >> $COMMON_VARS

