#!/bin/bash

appendToConfigFile() {
	configFile=${1}
	variableName=${2}
	variableValue=${3}
	
	echo "${variableName}=\"${variableValue}\"" >> ${configFile}
}
		
configFile=${1}

source ${configFile}

# set up job/task ID variables
jobDir=${reduceScratchDir}/cWorld__stage1-rdc__${UUID}__${reduceID}__${LSB_JOBID}
mkdir -p ${jobDir}

appendToConfigFile ${configFile} "jobDir" ${jobDir}

# set up perl/python/shell paths
determineQVEncoding=$cMapping/perl/determineQVEncoding.pl
drawIterativePlots=$cMapping/R/Draw_CF_plots.R
reduce=$cMapping/mapReduce/tyler_reduceHiC.sh

#initialize isilon result folders
mkdir -p ${outputFolder}

# get name of files
side1FileName=`basename "${side1File}"`
side2FileName=`basename "${side2File}"`

#handle the file copy
cp ${side1File} ${jobDir}/.
cp ${side2File} ${jobDir}/.

if [[ ${zippedFlag} = 1 ]]
then
	side1FileName=s_${laneNum}_1_sequence.txt.gz
	side2FileName=s_${laneNum}_2_sequence.txt.gz
	cat `ls ${jobDir}/*_R1_[0-9][0-9][0-9].fastq.gz | sort -t _ -k3,3n -k4,4n` > ${jobDir}/${side1FileName}
	cat `ls ${jobDir}/*_R2_[0-9][0-9][0-9].fastq.gz | sort -t _ -k3,3n -k4,4n` > ${jobDir}/${side2FileName}
	rm ${jobDir}/*_R1_[0-9][0-9][0-9].fastq.gz
	rm ${jobDir}/*_R2_[0-9][0-9][0-9].fastq.gz
	
else
	gzip ${jobDir}/${side1FileName}
	gzip ${jobDir}/${side2FileName}
	side1FileName=${side1FileName}.gz
	side2FileName=${side2FileName}.gz
fi

appendToConfigFile ${configFile} "side1FileName" ${side1FileName}
appendToConfigFile ${configFile} "side2FileName" ${side2FileName}

side1ShortFileName=${side1FileName}
side1ShortFileName=${side1ShortFileName//.gz/}

side2ShortFileName=${side2FileName}
side2ShortFileName=${side2ShortFileName//.gz/}

appendToConfigFile ${configFile} "side1ShortFileName" ${side1ShortFileName}
appendToConfigFile ${configFile} "side2ShortFileName" ${side2ShortFileName}

nReads1=`zcat ${jobDir}/${side1FileName} | wc -l | awk -F" " '{ print $1 }'`
let "nReads1 = ($nReads1/4)"
nReads2=`zcat ${jobDir}/${side2FileName} | wc -l | awk -F" " '{ print $1 }'`
let "nReads2 = ($nReads2/4)"	
	
nReads=0
if [[ ${nReads1} = ${nReads2} ]] 
then
	nReads=${nReads1}
else 
	echo "error - input files do not contain equal number of lines"
	echo "${side1FileName} = ${nReads1}"
	echo "${side2FileName} = ${nReads2}"
	exit;
fi

appendToConfigFile ${configFile} "nReads" ${nReads}

qvEncoding1=`perl ${determineQVEncoding} -i ${jobDir}/${side1FileName} -aln ${aligner}`
qvEncoding2=`perl ${determineQVEncoding} -i ${jobDir}/${side2FileName} -aln ${aligner}`

# check for valid QV encoding on both sides
if [[ ${qvEncoding1} = ${qvEncoding2} ]] && [[ ${qvEncoding1} != "error" ]] && [[ ${qvEncoding2} != "error" ]]
then
	qvEncoding=${qvEncoding1}
else
	echo "error with QV format - exiting."
	exit;
fi

appendToConfigFile ${configFile} "qvEncoding" ${qvEncoding}

#run the map reduce on the reads 
${reduce} ${configFile} ${LSB_JOBID} ${jobDir}

# copy results back to output dir
cp ${jobDir}/${jobName}.tar.gz ${outputFolder}/.
cp ${jobDir}/${jobName}.validPair.txt.gz ${outputFolder}/.
cp ${jobDir}/${jobName}.validPair.itx.gz ${outputFolder}/.
cp ${configFile} ${outputFolder}/${jobName}.cfg

if [ $keepSAM = 1 ]
then
	cp ${jobDir}/${jobName}__side1.sam.gz ${outputFolder}/.
	cp ${jobDir}/${jobName}__side2.sam.gz ${outputFolder}/.
fi
