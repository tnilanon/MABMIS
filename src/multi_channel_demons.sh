#!/bin/bash

if [[ ${#} -ne 12 ]]; then
	echo 'Not enough input arguments'
	echo 'Usage: multi_channel_demons.sh'
	echo '	[fixed_intensity_image]'
	echo '	[fixed_label_image]'
	echo '	[moving_intensity_image]'
	echo '	[moving_label_image]'
	echo '	[input_deformation_field] (input "NULL" if not exist)'
	echo '	[output_deformation_field]'
	echo '	[output_intensity_image]'
	echo '	[output_label_image]'
	echo '	[sigma_deformation_field]'
	echo '	[iteration_in_each_of_the_three_levels]'
	exit
fi

fixed_intensity=${1}
fixed_label=${2}
unaligned_moving_intensity=${3}
unaligned_moving_label=${4}
input_deformation_field=${5}
output_deformation_field=${6}
output_intensity=${7}
output_label=${8}
sigma_deformation_field=${9}
iteration_level_1=${10}
iteration_level_2=${11}
iteration_level_3=${12}

