#!/bin/bash

appendToConfigFile() {
	configFile=${1}
	variableName=${2}
	variableValue=${3}
	
	echo "${variableName}=\"${variableValue}\"" >> ${configFile}
}		

checkForErrors() {
	if [ ! -z $1 ]
	then
		verbose=1
	else
		verbose=0
	fi
		
	if [ $verbose = 1 ]
	then
		echo -e "mapReduce (${jobID}) - checking for errors...";
	fi
	
	for ((  i = 0;  i < ${nMaps};  i++  )) 
	do
		errorFile=${mapReduceDir}/error/${maps[${i}]}.error
		if [ -s ${errorFile} ]
		then
			
			echo -e "mapReduce reported error! ( ${maps[${i}]}.error )\n----------"
			cat ${errorFile}
			echo -e "----------\n"
			
			bkill ${jobID}
			exit
		fi
	done
	
	if [ $verbose = 1 ]
	then
		echo -e "\tAOK\n";
	fi
	
}
		
#take a input file to work on
configFile=${1}
jobID=${2}
jobDir=${3}

source ${configFile}

# load R module
module load R/3.0.2 &> /dev/null
module load python/2.7.5 &> /dev/null

# set up perl/python/shell paths
map=${cMapping}/mapReduce/mapHiC.sh
mapReduceHome=${cMapping}/mapReduce
splitFile=${cMapping}/perl/splitFile.pl
collapseStatFile=${cMapping}/perl/collapse_stats_file.pl
plotMoleculeSizeHistogram=${cMapping}/perl/plotMoleculeSizeHistogram.pl
plotStrandBias=${cMapping}/perl/plotStrandBias.pl
strandBiasPlot=${cMapping}/R/plotStrandBias.R
moleculeSizeHistogram=${cMapping}/R/moleculeSizeHistogram.R
aggregrateInteractionLog=${cMapping}/perl/aggregrateInteractionLog.pl
plotInteractionLog=${cMapping}/R/plotInteractionLog.R
aggregrateMappingLog=${cMapping}/perl/aggregrateMappingLog.pl
drawIterativePlots=${cMapping}/R/Draw_CF_plots.R
createHiCMappingLog=${cMapping}/perl/createHiCMappingLog.pl
aggregrateAlleleLog=${cMapping}/perl/aggregrateAlleleLog.pl
initialEmail=${cMapping}/perl/emailInitial.pl
completionEmail=${cMapping}/perl/emailResults.pl
collapseHiCValidPair=${cMapping}/perl/collapseHiCValidPairs.pl

# set up perl/python/shell paths

jobDir=${jobDir/%\//} #strip file / from jobDir
mapReduceDir=${jobDir}/mapReduce

appendToConfigFile ${configFile} "mapReduceDir" ${mapReduceDir}

#create scratch space
mkdir -p ${mapReduceDir}
mkdir -p ${mapReduceDir}/state
mkdir -p ${mapReduceDir}/chunks
mkdir -p ${mapReduceDir}/stats
mkdir -p ${mapReduceDir}/sam
mkdir -p ${mapReduceDir}/intervals
mkdir -p ${mapReduceDir}/allele-log
mkdir -p ${mapReduceDir}/moleculeSize-log
mkdir -p ${mapReduceDir}/strandBias-log
mkdir -p ${mapReduceDir}/validPairs
mkdir -p ${mapReduceDir}/interaction-log
mkdir -p ${mapReduceDir}/mapping-log
mkdir -p ${mapReduceDir}/aligner-log
mkdir -p ${mapReduceDir}/log
mkdir -p ${mapReduceDir}/error
mkdir -p ${mapReduceDir}/plots

# send initial email
perl ${initialEmail} -j ${jobID} -jn ${mapReduceDir}/log/${jobName} -cf ${configFile}

# split the input file into N chunks
nSide1Chunks=`perl ${splitFile} -i ${jobDir}/${side1FileName} -s ${splitSize} -g 4 -o ${mapReduceDir}/chunks`
nSide2Chunks=`perl ${splitFile} -i ${jobDir}/${side2FileName} -s ${splitSize} -g 4 -o ${mapReduceDir}/chunks`

restrictionFragmentFileName=`basename ${restrictionFragmentPath}`
sort -k1,1 -k2,2n ${restrictionFragmentPath} -o ${mapReduceDir}/intervals/${restrictionFragmentFileName}
appendToConfigFile ${configFile} "mapReduceRestrictionFragmentFile" ${mapReduceDir}/intervals/${restrictionFragmentFileName}

nChunks=0
nJobs=0
if [ $nSide1Chunks -ne $nSide2Chunks ]
then
	echo "ERROR - files are not of equal size $nSide1Chunks / $nSide2Chunks"
	echo "exiting..."
	bkill ${jobID}
	exit
else 
	nChunks=$nSide1Chunks
	let "nJobs = (($nChunks/8)+1)";
fi

# now clean up input fastq files (to save space)
if [ $debugModeFlag = 0 ]; then rm ${jobDir}/${side1FileName}; fi
if [ $debugModeFlag = 0 ]; then rm ${jobDir}/${side2FileName}; fi 

#submit all the map segments.
nMaps=0
for ((  i = 0;  i < ${nJobs};  i++  )) 
do
	let "chunkStart = $i * 8";
	let "chunkEnd = ((($i+1) * 8)-1)";
	
	
	#if chunkEnd > #chunks - correct.
	if [ $chunkEnd -gt $nChunks ]
	then
		chunkEnd=${nChunks}
	fi
	let "nTasks = (($chunkEnd - $chunkStart)+1)"
		
	for ((  i2 = ${chunkStart};  i2 <= ${chunkEnd};  i2++  )) 
	do	
		maps[${nMaps}]=${jobName}.c${i2}
		let nMaps++
	done
	
	mapID=`uuidgen | rev | cut -d '-' -f 1`
	let adjustedMapMemoryNeededMegabyte=($mapMemoryNeededMegabyte+$nTasks-1)/$nTasks; # adjust memory usage by nCPU requested
	let adjustedMapMemoryNeededMegabyte=$adjustedMapMemoryNeededMegabyte+1024 # add 1024 of working memory per process
	bsub -n ${nTasks} -q $mapQueue -R "span[hosts=1]" -R "rusage[mem=$adjustedMapMemoryNeededMegabyte:tmp=$mapScratchSize]" -W $mapTimeNeeded -N -u $userEmail -J mapHiC -o $userHomeDirectory/lsf_jobs/LSB_%J.log -e $userHomeDirectory/lsf_jobs/LSB_%J.err -Ep "${cMapping}/scripts/garbageCollectTmp.sh ${UUID} ${mapID} $userHomeDirectory/lsf_jobs" ${map} ${configFile} ${chunkStart} ${chunkEnd} ${mapID} &> /dev/null
	
	sleep 5
done

#now look for maps to report back that they are complete.
completeMaps=0
while [ ${completeMaps} -lt ${nMaps} ]
do
	# check for any MAP errors reported
	checkForErrors

	completeMaps=0
	for ((  i = 0;  i < ${nMaps};  i++  )) 
	do		
		mapFile=${mapReduceDir}/state/${maps[${i}]}
		
		if [ -f ${mapFile} ]
		then
			let completeMaps++
		fi
	done	   
	sleep 30
done

# check for any MAP errors reported
checkForErrors 1

# combine all stat files
i=${iterativeMappingStart}
lastIteration=0
while [ $i -le $iterativeMappingEnd ]
do 
	
	cat ${mapReduceDir}/stats/i${i}/${side1ShortFileName}.c*.i${i}.stats > ${mapReduceDir}/stats/${side1ShortFileName}.i${i}.stats
	perl ${collapseStatFile} -i ${mapReduceDir}/stats/${side1ShortFileName}.i${i}.stats -s 1
	cat ${mapReduceDir}/stats/i${i}/${side2ShortFileName}.c*.i${i}.stats > ${mapReduceDir}/stats/${side2ShortFileName}.i${i}.stats
	perl ${collapseStatFile} -i ${mapReduceDir}/stats/${side2ShortFileName}.i${i}.stats -s 2
	
	i=$(( $i + $iterativeMappingStep ))
	if [ $i -gt $iterativeMappingEnd ] && [ $(( $i - $iterativeMappingEnd )) -lt $iterativeMappingStep ]
	then
		i=$readLength
	fi
	
done

# correct command to use to prevent argument list too long
# find ${mapReduceDir}/validPairs/ -name "${jobName}.c*.validPair.txt.gz" -print0 | xargs -0 gunzip
# find ${mapReduceDir}/validPairs/ -name "${jobName}.c*.validPair.txt" -print0 | xargs -0 sort -m -k6,6n -k12,12n -k3,3n -k9,9n -o ${mapReduceDir}/validPairs/${jobName}.validPair.txt

# [optional] - collapse all sam files
if [ $keepSAM = 1 ]
then
	cat ${mapReduceDir}/sam/${side1ShortFileName}.c*.sam | gzip > ${mapReduceDir}/sam/${jobName}__side1.sam.gz
	cat ${mapReduceDir}/sam/${side2ShortFileName}.c*.sam | gzip > ${mapReduceDir}/sam/${jobName}__side2.sam.gz
fi

# purge the sam files
if [ ${debugModeFlag} = 0 ] && [ ${keepSAM} = 1 ]
then
	rm ${mapReduceDir}/sam/*.c*.sam; 
fi

# now merge sort all valid pair files
cmd="sort -m -k6,6n -k12,12n -k3,3n -k9,9n "
for input in ${mapReduceDir}/validPairs/${jobName}.c*.validPair.txt.gz; 
do
    cmd="$cmd <(gunzip -c '$input')"
done
eval "$cmd" | gzip > ${mapReduceDir}/validPairs/${jobName}.validPair.txt.gz

# purge the valid pair files
if [ $debugModeFlag = 0 ]; then rm ${mapReduceDir}/validPairs/${jobName}.c*.validPair.txt.gz; fi

# collapse valid pairs
nFinalValidPairs=`perl ${collapseHiCValidPair} -jn ${mapReduceDir}/validPairs/${jobName} -i ${mapReduceDir}/validPairs/${jobName}.validPair.txt.gz`
appendToConfigFile ${configFile} "nFinalValidPairs" ${nFinalValidPairs}
cp ${mapReduceDir}/validPairs/${jobName}.pcrDupe.log ${mapReduceDir}/log/.

# combine all chunk stats + calculate dangling-end molecule length histogram
perl ${plotMoleculeSizeHistogram} -i ${mapReduceDir}/moleculeSize-log/ -jn ${mapReduceDir}/moleculeSize-log/${jobName}
Rscript ${moleculeSizeHistogram} ${mapReduceDir}/moleculeSize-log/ ${jobName}.moleculeSize.plot.log &> /dev/null
cp ${mapReduceDir}/moleculeSize-log/${jobName}.moleculeSize.log ${mapReduceDir}/log/.
cp ${mapReduceDir}/moleculeSize-log/${jobName}.moleculeSize.plot.log.png ${mapReduceDir}/plots/.

# combine all chunk stats + calculate strand bias at nearby distance + cutoff to use
distanceCutoff=`perl ${plotStrandBias} -i ${mapReduceDir}/strandBias-log/ -jn ${mapReduceDir}/strandBias-log/${jobName}`
Rscript ${strandBiasPlot} ${mapReduceDir}/strandBias-log/ ${jobName}.strandBias.log ${distanceCutoff} &> /dev/null
cp ${mapReduceDir}/strandBias-log/${jobName}.strandBias.log ${mapReduceDir}/log/.
cp ${mapReduceDir}/strandBias-log/${jobName}.strandBias.log.png ${mapReduceDir}/plots/.

# combine all chunk stats + create mapping log pie charts
perl ${aggregrateInteractionLog} -i ${mapReduceDir}/interaction-log/ -jn ${mapReduceDir}/interaction-log/${jobName}
Rscript ${plotInteractionLog} ${mapReduceDir}/interaction-log/ ${jobName}.interaction.log &> /dev/null
cp ${mapReduceDir}/interaction-log/${jobName}.interaction.log ${mapReduceDir}/log/.
cp ${mapReduceDir}/interaction-log/${jobName}.interaction.log.png ${mapReduceDir}/plots/.

# combine all chunk stats + calculate mapping log info (U + NM + MM()
perl ${aggregrateMappingLog} -i ${mapReduceDir}/mapping-log/ -jn ${mapReduceDir}/mapping-log/${jobName}
cp ${mapReduceDir}/mapping-log/${jobName}.mapping.log ${mapReduceDir}/log/.

# draw plots for iterative mapping
cat ${mapReduceDir}/stats/iterativeMapping.header ${mapReduceDir}/stats/iterativeMapping.*.i*.stats > ${mapReduceDir}/stats/${jobName}.iterativeMapping.log
Rscript ${drawIterativePlots} ${mapReduceDir}/stats/ ${jobName}.iterativeMapping.log ${iterativeMappingStart} ${iterativeMappingEnd} ${iterativeMappingStep} &> /dev/null
cp ${mapReduceDir}/stats/${jobName}.iterativeMapping.log ${mapReduceDir}/log/.
cp ${mapReduceDir}/stats/${jobName}.iterativeMapping.log.png ${mapReduceDir}/plots/.

# combine all allele-log info
perl ${aggregrateAlleleLog} -i ${mapReduceDir}/allele-log/ -jn ${mapReduceDir}/allele-log/${jobName}
cp ${mapReduceDir}/allele-log/${jobName}.alleleOverride.log ${mapReduceDir}/log/.

# summarize all data and produce main email/log file
perl ${createHiCMappingLog} -jn ${mapReduceDir}/log/${jobName} -ilf ${mapReduceDir}/log/${jobName}.interaction.log -mlf ${mapReduceDir}/log/${jobName}.mapping.log -pdlf ${mapReduceDir}/log/${jobName}.pcrDupe.log -alf ${mapReduceDir}/log/${jobName}.alleleOverride.log -cf ${configFile}

# zip of /plots and /log
cp ${configFile} ${mapReduceDir}/plots/. 
tar -czvf ${jobDir}/${jobName}.tar.gz -C ${mapReduceDir}/ log/ plots/ > /dev/null

# copy tarball back to /plots to email out
cp ${jobDir}/${jobName}.tar.gz ${mapReduceDir}/plots/. 

# send out email
perl ${completionEmail} -j ${jobID} -jn ${mapReduceDir}/log/${jobName} -lf ${mapReduceDir}/log/${jobName}.end.mappingLog.txt -cf ${configFile} -pf ${mapReduceDir}/plots/

# copy results back to main project dir
cp ${mapReduceDir}/validPairs/${jobName}.validPair.txt.gz ${jobDir}/.
cp ${mapReduceDir}/validPairs/${jobName}.validPair.itx.gz ${jobDir}/.

if [ $keepSAM = 1 ]
then
	cp ${mapReduceDir}/sam/${jobName}__side1.sam.gz ${jobDir}/.
	cp ${mapReduceDir}/sam/${jobName}__side2.sam.gz ${jobDir}/.
fi

# do clean up
if [ $debugModeFlag = 0 ]; then rm -rf ${mapReduceDir}/; fi
