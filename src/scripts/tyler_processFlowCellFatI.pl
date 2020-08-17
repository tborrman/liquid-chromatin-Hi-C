use 5.006;
use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Carp qw(carp cluck croak confess);
use POSIX qw(ceil floor strftime);
use List::Util qw[min max];

use Cwd 'abs_path';
use Cwd;

my $tool=(split(/\//,abs_path($0)))[-1];
my $version = "1.0.6";

sub check_options {
    my $opts = shift;
    
    my $ret={};
    
    my ($flowCellDirectory,$scratchDirectory,$outputDirectory,$genomeDirectory,$logDirectory,$userEmail,$verbose,$genomeName,$hicModeFlag,$fiveCModeFlag,$keepSAM,$assumeCisAllele,$enzyme,$splitSize,$shortMode,$snpModeFlag,$debugModeFlag);
    
    if( defined($opts->{ flowCellDirectory }) ) {
        $flowCellDirectory = $opts->{ flowCellDirectory };
        $flowCellDirectory =~ s/\/$//;
        croak "flowCellDirectory [".$flowCellDirectory."] does not exist" if(!(-d $flowCellDirectory));
    } else {
        print STDERR "\nERROR: Option inputFlowCellDirectory|i is required.\n";
        help();
    }
    
    if( defined($opts->{ scratchDirectory }) ) {
        $scratchDirectory = $opts->{ scratchDirectory };
        $scratchDirectory =~ s/\/$//;
        croak "scratchDirectory [".$scratchDirectory."] does not exist" if(!(-d $scratchDirectory));
    } else {
        print STDERR "\nERROR: Option scratchDirectory|s is required.\n";
        help();
    }
    
    if( defined($opts->{ outputDirectory }) ) {
        $outputDirectory = $opts->{ outputDirectory };
        $outputDirectory =~ s/\/$//;
        croak "outputDirectory [".$outputDirectory."] does not exist" if(!(-d $outputDirectory));
    } else {
        print STDERR "\nERROR: Option outputDirectory|o is required.\n";
        help();
    }
    
    if( defined($opts->{ genomeDirectory }) ) {
        $genomeDirectory = $opts->{ genomeDirectory };
        $genomeDirectory =~ s/\/$//;
        croak "genomeDirectory [".$genomeDirectory."] does not exist" if(!(-d $genomeDirectory));
    } else {
        print STDERR "\nERROR: Option genomeDirectory|gdir is required.\n";
        help();
    }
    
    if( defined($opts->{ logDirectory }) ) {
        $logDirectory = $opts->{ logDirectory };
        $logDirectory =~ s/\/$//;
        croak "logDirectory [".$logDirectory."] does not exist" if(!(-d $logDirectory));
    } else {
        $logDirectory=$outputDirectory."/cWorld-logs";
    }
    
    if( defined($opts->{ userEmail }) ) {
        $userEmail = $opts->{ userEmail };
    } else {
        $userEmail=getUserEmail();
    }
    
    if( exists($opts->{ verbose }) ) {
        $verbose = 1;
    } else {
        $verbose = 0;
    }
    
    if( defined($opts->{ genomeName }) ) {
        $genomeName = $opts->{ genomeName };
    } else {
        $genomeName="";
    }
    
    if( defined($opts->{ hicModeFlag }) ) {
        $hicModeFlag = 1;
    } else {
        $hicModeFlag=0;
    }
    
    if( defined($opts->{ fiveCModeFlag }) ) {
        $fiveCModeFlag = 1;
    } else {
        $fiveCModeFlag=0;
    }
    
    if( defined($opts->{ keepSAM }) ) {
        $keepSAM = 1;
    } else {
        $keepSAM=0;
    }
    
    if( defined($opts->{ assumeCisAllele }) ) {
        $assumeCisAllele = 1;
    } else {
        $assumeCisAllele = 0;
    }
    
    if( defined($opts->{ enzyme }) ) {
        $enzyme = $opts->{ enzyme };
    } else {
        $enzyme = "HindIII";
    }
    
    if( defined($opts->{ splitSize }) ) {
        $splitSize = $opts->{ splitSize };
        confess "split size is too small - exiting" if($splitSize < 500000);
    } else {
        $splitSize=4000000;
    }
    
    if( defined($opts->{ shortMode }) ) {
        $shortMode = 1;
    } else {
        $shortMode=0;
    }
    
    if( defined($opts->{ snpModeFlag }) ) {
        $snpModeFlag = 1;
    } else {
        $snpModeFlag=0;
    }
    
    if( $opts->{ debugModeFlag } ) {
        $debugModeFlag=1;
    } else {
        $debugModeFlag=0;
    }
    
    croak "\nERROR: must choose either -f (5C) or -h (HiC) mode option\n" if( (($fiveCModeFlag+$hicModeFlag) <= 0) or (($fiveCModeFlag+$hicModeFlag) >= 2) );
    croak "\nERROR: cannot use SNP mode with -f option\n" if(($snpModeFlag+$fiveCModeFlag >= 2) );
    
    $ret->{ flowCellDirectory }=$flowCellDirectory;
    $ret->{ scratchDirectory }=$scratchDirectory;
    $ret->{ outputDirectory }=$outputDirectory;
    $ret->{ genomeDirectory }=$genomeDirectory;
    $ret->{ logDirectory }=$logDirectory;
    $ret->{ userEmail }=$userEmail;
    $ret->{ verbose }=$verbose;
    $ret->{ genomeName }=$genomeName;
    $ret->{ hicModeFlag }=$hicModeFlag;
    $ret->{ fiveCModeFlag }=$fiveCModeFlag;
    $ret->{ keepSAM }=$keepSAM;
    $ret->{ assumeCisAllele }=$assumeCisAllele;
    $ret->{ enzyme }=$enzyme;
    $ret->{ splitSize }=$splitSize;
    $ret->{ shortMode }=$shortMode;
    $ret->{ snpModeFlag }=$snpModeFlag;
    $ret->{ debugModeFlag }=$debugModeFlag;
    
    return($flowCellDirectory,$scratchDirectory,$outputDirectory,$genomeDirectory,$logDirectory,$userEmail,$verbose,$genomeName,$hicModeFlag,$fiveCModeFlag,$keepSAM,$assumeCisAllele,$enzyme,$splitSize,$shortMode,$snpModeFlag,$debugModeFlag);
}

sub getRestrictionEnzymeSequences() {
    my %restrictionEnzymeSequences=();
    
    $restrictionEnzymeSequences{ HindIII } = "AAGCTT";
    $restrictionEnzymeSequences{ EcoRI } = "GAATTC";
    $restrictionEnzymeSequences{ NcoI } = "CCATGG";
    $restrictionEnzymeSequences{ DpnII } = "GATC";
    $restrictionEnzymeSequences{ MNase } = "MNase";
    $restrictionEnzymeSequences{ FatI } = "CATG";
    
    return(\%restrictionEnzymeSequences);
}

sub getUserEmail() {
    
    # hb67w:x:10839:1081:Houda Belaghzal [Houda.belaghzal@umassmed.edu]:/home/hb67w:/bin/bash
    my $user_info=`grep \$USER /etc/passwd`;
    chomp($user_info);
    
    my @tmp=split(/:/,$user_info);
    my $user_email=$tmp[4];
    $user_email=(split(/\[/,$user_email))[1];
    $user_email =~ s/\]//;
    
    $user_email = "" if($user_email !~ /\@/);
    
    return($user_email);
}

sub getUserHomeDirectory() {
    my $userHomeDirectory = `echo \$HOME`;
    chomp($userHomeDirectory);
    return($userHomeDirectory);
}

sub getUniqueString() {
    my $UUID = `uuidgen`;
    chomp($UUID);
    return($UUID);
}

sub getSmallUniqueString() {
    my $UUID=`uuidgen | rev | cut -d '-' -f 1`;
    chomp($UUID);
    return($UUID);
}

sub getComputeResource() {
    my $hostname = `hostname`;
    chomp($hostname);
    return($hostname);
}

sub translateFlag($) {
    my $flag=shift;
    
    my $response="no";
    $response="yes" if($flag);    
    return($response);
}

sub check_dependency($;$) {
    # required
    my $command=shift;
    # optional
    my $weblink=shift;
    
    my $repo=(split(/\//,$command))[-3];
    
    confess "missing dependency [$repo] - $command.\n\tPlease install\n\t$weblink\n\n" if(!-e($command));
    
    return($command);
    
}

sub which($;$) {
    # required
    my $command=shift;
    # optional
    my $die=1;
    $die=shift if @_;
    
    my $path="";
    $path=`which $command 2>&1`;
    chomp($path);
    
    confess "no path for $command" if(($path =~ /which: no/) and ($die == 1));
        
    return($path);
}

sub getDate() {
    my $time = strftime '%I:%M:%S %P, %m/%d/%Y', localtime;
    
    return($time);
}

sub commify {
   (my $num = shift) =~ s/\G(\d{1,3})(?=(?:\d\d\d)+(?:\.|$))/$1,/g; 
   return $num; 
}

#

sub getScriptOpts($$) {
    # required
    my $ret=shift;
    my $tool=shift;
    
    my $commentLine="";
    
    $commentLine = "## Tool:\t".$tool;
    $commentLine .= "\n## ";
    
    foreach my $opt ( sort keys %{$ret} ) {
        my $value=$ret->{$opt};
        $value="" if(!defined($value));
        if (ref $value eq 'ARRAY') {
            my @value=@$value;
            my $value_str=join(',',@value);
            $value=$value_str;
        }
        $commentLine .= "\n## ".$opt." = '".$value."'";
    }
    
    return $commentLine;
}
   
   
sub getAlignmentSoftware() {
    
    my $userHomeDirectory = getUserHomeDirectory();
    
    my %alignmentSoftware=();
    
    my $ret="";
    
    $ret=`which bowtie2 2> /dev/null`;
    chomp($ret);
    $alignmentSoftware{ bowtie2 }=$ret;
    
    $ret=`which novoalign 2> /dev/null`;
    chomp($ret);
    $alignmentSoftware{ novoalign }=$ret;
    
    return(\%alignmentSoftware);
}

sub getGenomePath($$$$) {
    my $genomeDirectory=shift;
    my $aligner=shift;
    my $genomeName=shift;
    my $restrictionSite=shift;
        
        
    my $fastaDirectory=$genomeDirectory."/fasta/".$genomeName;
    my $indexDirectory=$genomeDirectory."/".$aligner."/".$genomeName;
    my $restrictionFragmentFile=$genomeDirectory."/restrictionFragments/".$genomeName."/".$genomeName."__".$restrictionSite.".txt";
    
    confess "invalid genome directory ($indexDirectory)\n" if(!(-d($indexDirectory)));
    confess "invalid fasta directory ($fastaDirectory)\n" if(!(-d($fastaDirectory)));
    confess "invalid restriction fragment file ($restrictionFragmentFile)\n" if(!(-e($restrictionFragmentFile)));
    
    return($indexDirectory,$fastaDirectory,$restrictionFragmentFile);
}

sub getDefaultOutputFolder($$$) {
    my $flowCell=shift;
    my $laneName=shift;
    my $outputFolder=shift;
    
    my $userHomeDirectory = getUserHomeDirectory();
    
    $outputFolder=$userHomeDirectory."/scratch/cData" if($outputFolder eq "");
    croak "scratch dir [".$outputFolder."] does not exist" if(!(-d $outputFolder));
    
    $outputFolder .= "/$flowCell/$laneName";
    
    return($outputFolder);
}

sub readFlowCellDirectory($) {
    my $dataDirectory=shift;
    $dataDirectory =~ s/\/$//;

    opendir(BIN, $dataDirectory) or die "Can't open [$dataDirectory]: $!";
    my @fileNames = readdir BIN ;
    close(BIN);
    my $nFiles = @fileNames;
    
    my @lanes=();
    for(my $i=0; $i<$nFiles; $i++) {
        my $dataFileName=$fileNames[$i];        
        
        next if ($dataFileName =~ /^\.\.?$/);
        
        push(@lanes,$dataFileName) if(-d $dataDirectory."/".$dataFileName);
    }
    
    return(\@lanes);
}

sub searchForFASTQ($$$$$$$) {
    my $fastqFiles=shift;
    my $dataDirectory=shift;
    $dataDirectory =~ s/\/$//;
    my $laneNum=shift;
    my $side1FastqFile=shift;
    my $side2FastqFile=shift;
    my $readLength=shift;
    my $zippedFlag=shift;
    
    opendir(BIN, $dataDirectory) or die "Can't open $dataDirectory: $!";
    my @fileNames = readdir BIN ;
    close(BIN);
    my $nFiles = @fileNames;
    
    my $index="";
    my $last_side1_dataFileName="";
    my $last_side2_dataFileName="";
    
    for(my $i=0; $i<$nFiles; $i++) {
        my ($side1File,$side2File);
        
        my $dataFileName=$fileNames[$i];        
        
        next if ($dataFileName =~ /^\.\.?$/);
        next if ($dataFileName =~ /^\./);
        
        if(-d $dataDirectory."/".$dataFileName) {
            
            &searchForFASTQ($fastqFiles,$dataDirectory."/".$dataFileName,$laneNum,$side1FastqFile,$side2FastqFile,$readLength,$zippedFlag) if(($dataFileName =~ /bustard/i) or ($dataFileName =~ /gerald/i));
            
        } else {
        
            next if(($dataFileName !~ /.fq$/) and ($dataFileName !~ /_sequence.txt/) and ($dataFileName !~ /.fastq$/) and ($dataFileName !~ /.fastq.gz$/));
            
            if($dataFileName =~ /\.gz$/) {
                
                my $fastqFile=$dataDirectory."/".$dataFileName;
                my $sampleLine = `zcat $fastqFile | head -n 2 | tail -n 1`;
                chomp($sampleLine);
                my $readLength = length($sampleLine);
                my @tmp=split(/_/,$dataFileName);
                
                $index=$tmp[1] if(@tmp == 5);
                die("\nfile name error ($dataFileName | @tmp)\n") if(@tmp != 5);
                
                foreach(@tmp) { $laneNum = $_ if($_ =~ /L[0-9]{3}/); }
                $laneNum =~ s/L00//;
                
                my $tmp_dataFileName=$dataFileName;
                $tmp_dataFileName =~ s/_[0-9]{3}\.fastq\.gz/_\*\.fastq\.gz/;
                $tmp_dataFileName =~ s/\_$index\_/\_\*\_/;
                
                my $side="NA";
                $side=1 if($tmp_dataFileName =~ /_R1_/);
                $side=2 if($tmp_dataFileName =~ /_R2_/);
                
                die("Multiple file names detected!\n\t$last_side1_dataFileName vs $tmp_dataFileName") if(($side == 1) and (($last_side1_dataFileName ne "") and ($last_side1_dataFileName ne $tmp_dataFileName)));
                die("Multiple file names detected!\n\t$last_side2_dataFileName vs $tmp_dataFileName") if(($side == 2) and (($last_side2_dataFileName ne "") and ($last_side2_dataFileName ne $tmp_dataFileName)));
                
                $fastqFiles->{$side}->{"path"}=$dataDirectory."/".$tmp_dataFileName if(($side == 1) or ($side == 2));
                $fastqFiles->{$side}->{"readLength"}=$readLength if(($side == 1) or ($side == 2));
                $fastqFiles->{"laneNum"}=$laneNum;
                
                $zippedFlag=1;
                $last_side1_dataFileName=$tmp_dataFileName if($side == 1);
                $last_side2_dataFileName=$tmp_dataFileName if($side == 2);
                
            } else {
                
                my $fastqFile=$dataDirectory."/".$dataFileName;
                my $sampleLine = `head -n 2 $fastqFile | tail -n 1`;
                chomp($sampleLine);
                my $readLength = length($sampleLine);
                $laneNum=(split(/_/,$dataFileName))[1];
                
                my $side="NA";
                $side=1 if($dataFileName =~ /s_.+_1_sequence.txt/);
                $side=2 if($dataFileName =~ /s_.+_2_sequence.txt/);
                
                $fastqFiles->{$side}->{"path"}=$dataDirectory."/".$dataFileName if(($side == 1) or ($side == 2));
                $fastqFiles->{$side}->{"readLength"}=$readLength if(($side == 1) or ($side == 2));
                $fastqFiles->{"laneNum"}=$laneNum;
                
            }
            
        }
    }
    
    print STDERR "Warning - encountered single end reads! - skipping\n" if( (!exists($fastqFiles->{1})) or (!exists($fastqFiles->{2})) );
    next if( (!exists($fastqFiles->{1})) or (!exists($fastqFiles->{2})) );
    
    $side1FastqFile=$fastqFiles->{1}->{"path"};
    my $side1ReadLength=$fastqFiles->{1}->{"readLength"};
    $side2FastqFile=$fastqFiles->{2}->{"path"};
    my $side2ReadLength=$fastqFiles->{2}->{"readLength"};
    $laneNum=$fastqFiles->{"laneNum"};
    
    print "\n\t\tWARNING: read lengths are not equal!\n\t\t\t(side1 - $side1ReadLength | side2 - $side2ReadLength)!\n\n" if(abs($side1ReadLength-$side2ReadLength) > 1);
    $readLength=min($side1ReadLength,$side2ReadLength);
    
    return($fastqFiles,$dataDirectory,$laneNum,$side1FastqFile,$side2FastqFile,$readLength,$zippedFlag);

}

sub logConfigVariable($$$) {
    my $configFileVariables=shift;
    my $configVariableName=shift;
    my $configVariableValue=shift;
    
    $configFileVariables->{$configVariableName}=$configVariableValue;
    
    return($configFileVariables);
}

sub printConfigFile($$$) {
    my $configFileVariables=shift;
    my $tmpConfigFileVariables=shift;
    my $configFileName=shift;
    
    open(OUT,">".$configFileName);
    
    my $time=getDate();
    my $userHomeDirectory=getUserHomeDirectory();
    
    print OUT "# cWorld processFlowCell\n";
    print OUT "# my5C.umassmed.edu\n";
    print OUT "# $time\n";
    print OUT "# $userHomeDirectory\n";
    print OUT "# ".$configFileVariables->{ computeResource}."\n";
    print OUT "# initial variables\n";
    
    for my $variableName ( sort {$a cmp $b} keys %{$configFileVariables}) {
        my $variableValue=$configFileVariables->{$variableName};
        print OUT $variableName."="."\"$variableValue\"\n";
    }
    
    for my $variableName ( sort {$a cmp $b} keys %{$tmpConfigFileVariables}) {
        my $variableValue=$tmpConfigFileVariables->{$variableName};
        print OUT $variableName."="."\"$variableValue\"\n";
    }
    
    print OUT "# dynamic variables\n";
    
    close(OUT);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$version."\n";
    print STDERR "Summary:\tcMapping pipeline - stage 1\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl processFlowCell.pl [OPTIONS] -i <inputFlowCellDirectory> -o <outputDirectory> --gdir <genomeDirectory>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "flow cell directory (path)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-s", "[]", "scratch directory (path)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "output directory (path)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--gdir", "[]", "genome directory (fasta,index,restrictionSite)");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--log", "[]", "log directory");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--email", "[]", "user email address");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--split", "[]", "splitSize, # reads per chunk");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-g", "[]", "genomeName, genome to align");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-e", "[]", "enzyme name (DpnII, HindIII etc.)");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-h", "[]", "FLAG, hic flag ");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-f", "[]", "FLAG, 5C flag");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-d", "[]", "FLAg, debugMode - keep all files for debug purposes");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--ks", "[]", "FLAG, keep sam files");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--short", "[]", "FLAG, use the short queue");
    printf STDERR ("\t%-10s %-10s %-10s\n", "--sm", "[]", "FLAG, snpMode - allelic Hi-C");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    Stage 1 of the cMapping pipeline, for processing 5C/Hi-C data [UMMS specific].\n";
    
    print STDERR "\n";
    
    print STDERR "Contact:
    Bryan R. Lajoie
    Dekker Lab 2016
    https://github.com/blajoie/cMapping
    https://github.com/blajoie/cWorld-dekker
    http://my5C.umassmed.edu";
    
    print STDERR "\n";
    print STDERR "\n";
    
    exit;
}

my %options;
my $results = GetOptions( \%options,'flowCellDirectory|i=s','scratchDirectory|s=s','outputDirectory|o=s','genomeDirectory|gdir=s','logDirectory|log=s','userEmail|email=s','hicModeFlag|h','verbose|v','genomeName|g=s','fiveCModeFlag|f','keepSAM|ks','assumeCisAllele|aca','enzyme|e=s','splitSize|split=i','shortMode|short','snpModeFlag|sm','debugModeFlag|d') or croak help();
my ($flowCellDirectory,$scratchDirectory,$outputDirectory,$genomeDirectory,$logDirectory,$userEmail,$verbose,$genomeName,$hicModeFlag,$fiveCModeFlag,$keepSAM,$assumeCisAllele,$enzyme,$splitSize,$shortMode,$snpModeFlag,$debugModeFlag)=check_options( \%options );

intro();

my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
my @scriptDir=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@scriptDir);
my @gitDir=@fullScriptPathArr[0..@fullScriptPathArr-5];
my $gitPath=join("/",@gitDir);

my $configFileVariables={};
my $userHomeDirectory = getUserHomeDirectory();
my $cMapping = $scriptPath;

$configFileVariables=logConfigVariable($configFileVariables,"cMapping",$cMapping);
$configFileVariables=logConfigVariable($configFileVariables,"gitDir",$gitPath);
$configFileVariables=logConfigVariable($configFileVariables,"keepSAM",$keepSAM);
$configFileVariables=logConfigVariable($configFileVariables,"splitSize",$splitSize);
$configFileVariables=logConfigVariable($configFileVariables,"hicModeFlag",$hicModeFlag);
$configFileVariables=logConfigVariable($configFileVariables,"snpModeFlag",$snpModeFlag);
$configFileVariables=logConfigVariable($configFileVariables,"fiveCModeFlag",$fiveCModeFlag);
$configFileVariables=logConfigVariable($configFileVariables,"debugModeFlag",$debugModeFlag);

# setup scratch space
my $reduceScratchDir=$scratchDirectory;
my $mapScratchDir="/tmp";
$mapScratchDir=$scratchDirectory if($debugModeFlag == 1);

# setup queue/timelimit for LSF
my $reduceQueue="long";
$reduceQueue="short" if(($debugModeFlag == 1) or ($shortMode == 1));
my $reduceTimeNeeded="120:00";
$reduceTimeNeeded="04:00" if(($debugModeFlag == 1) or ($shortMode ==1));
$configFileVariables=logConfigVariable($configFileVariables,"reduceQueue",$reduceQueue);
$configFileVariables=logConfigVariable($configFileVariables,"reduceTimeNeeded",$reduceTimeNeeded);

my $computeResource = getComputeResource();
$configFileVariables=logConfigVariable($configFileVariables,"computeResource",$computeResource);

$flowCellDirectory =~ s/\/$//;

my $lanesRef=readFlowCellDirectory($flowCellDirectory);
my $nLanes=@$lanesRef;

for(my $i=0;$i<$nLanes;$i++) {

    my $tmpConfigFileVariables={};
    
    my $laneName=$lanesRef->[$i];
    
    print "\n\t$laneName | process? (y,n) [n]: ";
    my $processFlag = <STDIN>;
    chomp($processFlag);
    $processFlag = "n" if($processFlag eq "");
    next if($processFlag ne "y");
    
    my @tmp=split(/\//,$flowCellDirectory);
    my $flowCellName=$tmp[@tmp-1];
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"flowCellName",$flowCellName);
    
    my $emptyHashRef={};
    my ($fastqFiles,$dataDirectory,$laneNum,$side1File,$side2File,$readLength,$zippedFlag)=searchForFASTQ($emptyHashRef,$flowCellDirectory."/".$laneName,0,"NA","NA",0,0);
    
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"flowCellName",$flowCellName);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"laneName",$laneName);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"laneNum",$laneNum);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"side1File",$side1File);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"side2File",$side2File);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"readLength",$readLength);

    print "\t\t$laneName [$laneNum]\n";
    print "\t\t(1)\t$side1File\t$readLength\n";
    print "\t\t(2)\t$side2File\t$readLength\n";    
    
    my $workDirectory = $flowCellDirectory;
    $workDirectory =~ s/$flowCellName//;
    $workDirectory =~ s/\/$//;
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"workDirectory",$workDirectory);
    
    print "\n";
    
    print "\t\treduceScratchDir [$reduceScratchDir]\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"reduceScratchDir",$reduceScratchDir);
    
    print "\t\tmapScratchDir [$mapScratchDir]:\t";
    my $userScratchDir = <STDIN>;
    chomp($userScratchDir);
    $mapScratchDir = $userScratchDir if($userScratchDir ne "");
    $mapScratchDir =~ s/\/$// if($mapScratchDir =~ /\/$/); # remove trailing / 
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"mapScratchDir",$mapScratchDir);
    print "\t\t\t$mapScratchDir\n";
    
    # assume 1 byte per ASCII.
    # 40 chars per header line
    # readLength chars per SEQ/QV line
    my $mapScratchSize = (((80*2)+($readLength*2))*($splitSize/4));
    $mapScratchSize = ceil((($mapScratchSize)/1024)/1024);
    $mapScratchSize = ($mapScratchSize * 30); #assume 10 fold input data of max /tmp usage
    print "\t\t\t".$mapScratchSize."M\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"mapScratchSize",$mapScratchSize);
    
    my $reduceScratchSize = 10000;
    print "\t\t\t".$reduceScratchSize."M\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"reduceScratchSize",$reduceScratchSize);
    
    # process the output directory
    my ($outputFolder)=getDefaultOutputFolder($flowCellName,$laneName,$outputDirectory);
    print "\t\toutputFolder [$outputFolder]:\t";
    my $userOutputFolder = <STDIN>;
    chomp($userOutputFolder);
    $outputFolder = $userOutputFolder if($userOutputFolder ne "");
    $outputFolder =~ s/\/$// if($outputFolder =~ /\/$/); # remove trailing / 
    $outputFolder = $userHomeDirectory."/".$outputFolder if($outputFolder !~ /^\//);
    
    system("mkdir -p $outputFolder") if(!(-d $outputFolder));
    confess "warning - cannot use specified outputFolder ($outputFolder)\n" if(!(-d $outputFolder));
    print "\t\t\t$outputFolder\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"outputFolder",$outputFolder);
    
    print "\t\tzipModeFlag: $zippedFlag\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"zippedFlag",$zippedFlag);
    
    my $alignmentSoftware=getAlignmentSoftware();
    
    # alignment software choice
    my $aligner="bowtie2";    
    if(($hicModeFlag == 1) and ($snpModeFlag == 0) ) {
        print "\t\taligner (bowtie2,novoalign) [$aligner]: ";
        my $userAligner = <STDIN>;
        chomp($userAligner);
        confess "invalid aligner ($userAligner)!" if(($userAligner ne "") and (($userAligner ne "bowtie2") and ($userAligner ne "novoalign")));
        $aligner = $userAligner if($userAligner ne "");
    } else {  # SNP mode or 5C mode
        $aligner="novoalign";
    }
    print "\t\t\t$aligner\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"aligner",$aligner);
    
    my $alignmentSoftwarePath=which($aligner,0);
    # alignment software path choice
    $alignmentSoftwarePath=$alignmentSoftware->{ $aligner } if(exists($alignmentSoftware->{ $aligner }));
    print "\t\talignerPath [$alignmentSoftwarePath]: ";
    my $userAlignmentSoftwarePath = <STDIN>;
    chomp($userAlignmentSoftwarePath);
    $alignmentSoftwarePath=$userAlignmentSoftwarePath if($userAlignmentSoftwarePath ne "");
    
    confess "invalid path for $aligner" if(!(-e($alignmentSoftwarePath)));
    print "\t\t\t$alignmentSoftwarePath\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"alignmentSoftwarePath",$alignmentSoftwarePath);
        
    # alignment options choice
    my $alignmentOptions="";
    $alignmentOptions="--very-sensitive --no-hd --no-sq --mm --qc-filter" if($aligner eq "bowtie2");
    $alignmentOptions=" -r all 5 -R 30 -q 2" if(($aligner eq "novoalign") and ($snpModeFlag == 1));
    print "\t\talignmentOptions [$alignmentOptions]: ";
    my $userAlignmentOptions = <STDIN>;
    chomp($userAlignmentOptions);
    $alignmentOptions=$userAlignmentOptions if($userAlignmentOptions ne "");
    print "\t\t\t$alignmentOptions\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"alignmentOptions",$alignmentOptions);
    
    # alignment options choice
    my $optionalSide1AlignmentOptions="";
    print "\t\toptional side1 alignmentOptions []: ";
    my $userOptionalSide1AlignmentOptions = <STDIN>;
    chomp($userOptionalSide1AlignmentOptions);
    $optionalSide1AlignmentOptions=$userOptionalSide1AlignmentOptions if($userOptionalSide1AlignmentOptions ne "");
    print "\t\t\t$optionalSide1AlignmentOptions\n" if($optionalSide1AlignmentOptions ne "");
    print "\t\t\t[none]\n" if($optionalSide1AlignmentOptions eq "");
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"optionalSide1AlignmentOptions",$optionalSide1AlignmentOptions);
    
    # alignment options choice
    my $optionalSide2AlignmentOptions="";
    print "\t\toptional side2 alignmentOptions []: ";
    my $userOptionalSide2AlignmentOptions = <STDIN>;
    chomp($userOptionalSide2AlignmentOptions);
    $optionalSide2AlignmentOptions=$userOptionalSide2AlignmentOptions if($userOptionalSide2AlignmentOptions ne "");
    print "\t\t\t$optionalSide2AlignmentOptions\n" if($optionalSide2AlignmentOptions ne "");
    print "\t\t\t[none]\n" if($optionalSide2AlignmentOptions eq "");
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"optionalSide2AlignmentOptions",$optionalSide2AlignmentOptions);
    
    my $minimumReadDistance=5;
    if(($aligner eq "novoalign") and ($snpModeFlag == 1)) {
        print "\t\tminimumReadDistance [5]: ";
        my $userMinimumReadDistance = <STDIN>;
        chomp($userMinimumReadDistance);
        $minimumReadDistance = $userMinimumReadDistance if($userMinimumReadDistance ne "");
        print "\t\t\t$minimumReadDistance\n";
        $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"minimumReadDistance",$minimumReadDistance);
    }
    
    if($snpModeFlag == 1) {
        print "\t\tassumeCisAllele [".translateFlag($assumeCisAllele)."]: ";
        my $userAssumeCisAllele = <STDIN>;
        chomp($userAssumeCisAllele);
        $assumeCisAllele=0 if(($userAssumeCisAllele ne "on") and ($userAssumeCisAllele ne "") and ($userAssumeCisAllele != 1));
        print "\t\t\t".translateFlag($assumeCisAllele)."\n";
    }
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"assumeCisAllele",$assumeCisAllele);
    
    # enzyme choice 
    my $restrictionEnzymeSequences=getRestrictionEnzymeSequences();
    my $enzymeString=join(',', (keys %{$restrictionEnzymeSequences}));

    my $restrictionSite="NA";
    if($hicModeFlag == 1) {
        print "\t\tenzyme (".$enzymeString.") [".$enzyme."]: ";
        my $userEnzyme = <STDIN>;
        chomp($userEnzyme);
        $enzyme=$userEnzyme if($userEnzyme ne "");
        confess "invalid restriction enzyme! ($enzyme)\n" if(!(exists($restrictionEnzymeSequences->{ $enzyme })));
        $restrictionSite=$restrictionEnzymeSequences->{ $enzyme };
        print "\t\t\t$enzyme / $restrictionSite\n";
    } else {
        $enzyme="NA"; 
    }
    
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"enzyme",$enzyme);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"restrictionSite",$restrictionSite);
    
    my $iterativeMappingFlag=0;
    my $iterativeMappingStart=$readLength;
    my $iterativeMappingEnd=$readLength;
    my $iterativeMappingStep=5;
    my $iterativeMappingIterations=0;
    
    if(($hicModeFlag == 1) and ($fiveCModeFlag == 0) and ($snpModeFlag == 0)) { #HiC data
        
        $iterativeMappingFlag=1;
        print "\t\titerative mapping mode? [yes] (yes|no): ";
        my $userIterativeMappingFlag = <STDIN>;
        chomp($userIterativeMappingFlag);        
        $iterativeMappingFlag = 0 if($userIterativeMappingFlag eq "no");
        print "\t\t\t".translateFlag($iterativeMappingFlag)."\n";
        
        # default iterative mapping options - if off (use full length read)
    
        if($iterativeMappingFlag == 1) {
            print "\t\t\titerativeMappingStart [25]: ";
            $iterativeMappingStart = <STDIN>;
            chomp($iterativeMappingStart);
            $iterativeMappingStart = 25 if($iterativeMappingStart eq "");
        
            print "\t\t\titerativeMappingEnd [$readLength]: ";
            $iterativeMappingEnd = <STDIN>;
            chomp($iterativeMappingEnd);
            $iterativeMappingEnd = $readLength if(($iterativeMappingEnd eq "") or ($iterativeMappingEnd < $iterativeMappingStart) or ($iterativeMappingEnd > $readLength));
        
            print "\t\t\titerativeMappingStep [5]: ";
            $iterativeMappingStep = <STDIN>;
            chomp($iterativeMappingStep);
            $iterativeMappingStep = 5 if(($iterativeMappingStep eq "") or ($iterativeMappingStep < 2));
            
            $iterativeMappingIterations = (($iterativeMappingEnd-$iterativeMappingStart)/$iterativeMappingStep);
            print "\t\t\t\titerativeMapping $iterativeMappingStart - $iterativeMappingEnd [$iterativeMappingStep] [$iterativeMappingIterations]\n";
        }
    }
    
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"iterativeMappingFlag",$iterativeMappingFlag);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"iterativeMappingStart",$iterativeMappingStart);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"iterativeMappingEnd",$iterativeMappingEnd);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"iterativeMappingStep",$iterativeMappingStep);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"iterativeMappingIterations",$iterativeMappingIterations);
    
    # genomeName choice
    print "\t\tgenome [".$genomeName."]: ";
    my $genome = <STDIN>;
    chomp($genome);
    
    $genomeName = $genome if($genome ne "");
    print "\n\tMust select a genome! ($genome | $genomeName) - skipping lane...\n\n" if($genomeName eq "");
    print "\t\t\t$genomeName\n";
    
    my $indexPath="NA";
    my $fastaPath="NA";
    my $restrictionFragmentPath="NA";
    
    if($hicModeFlag == 1) {
        ($indexPath,$fastaPath,$restrictionFragmentPath)=getGenomePath($genomeDirectory,$aligner,$genomeName,$restrictionSite);
        $restrictionFragmentPath="" if($snpModeFlag == 1);
        if($snpModeFlag == 1) {
            print "\t\trestrictionFragmentPath [$restrictionFragmentPath]: ";
            my $userRestrictionFragmentPath = <STDIN>;
            chomp($userRestrictionFragmentPath);
            $userRestrictionFragmentPath = "" if(!(-e($userRestrictionFragmentPath)));
            print "\t\t\t$userRestrictionFragmentPath\n";
            $restrictionFragmentPath = $userRestrictionFragmentPath if($userRestrictionFragmentPath ne "");
        }
        
        confess "invalid restriction fragment file path!\n" if(!(-e($restrictionFragmentPath)));
    }
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"restrictionFragmentPath",$restrictionFragmentPath);
    
    my $genomePath=$genomeDirectory."/".$aligner."/".$genomeName;
    my $genomeDir=$genomeDirectory."/".$aligner."/".$genomeName;
    if(($hicModeFlag == 1) and ($fiveCModeFlag == 0)) {
        $genomePath .= "/".$genomeName;
        confess "invalid genome path! (".$genomePath."*)\n" if( (!(glob($genomePath))) and (!(glob($genomePath."*"))) );
    }
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"genomeName",$genomeName);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"genomePath",$genomePath);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"genomeDir",$genomeDir);
    
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"userHomeDirectory",$userHomeDirectory);
    
    print "\t\tuserEmail (email address) [$userEmail]:\t";
    my $userEmailChoice = <STDIN>;
    chomp($userEmailChoice);
    $userEmail = $userEmailChoice if(($userEmailChoice ne "") and ($userEmailChoice =~ /@/) and ($userEmailChoice =~ /\s+/));
    print "\t\t\t$userEmail\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"userEmail",$userEmail);
    
    print "\t\temailTo (email address) [none]:\t";
    my $emailTo = <STDIN>;
    chomp($emailTo);
    $emailTo = "none" if(($emailTo eq "") or ($emailTo !~ /@/) or ($emailTo =~ /\s+/));
    print "\t\t\t$emailTo\n";
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"emailTo",$emailTo);
    
    my $UUID=getUniqueString();
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"UUID",$UUID);
    my $jobName=$flowCellName."__".$laneName."__".$genomeName;
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"jobName",$jobName);
    
    print "\t\tconfigFile\n";
    my $configFilePath=$logDirectory."/".$UUID.".cWorld-stage1.cfg";
    print "\t\t\t$configFilePath\n";
    
    my $reduceID=getSmallUniqueString();
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"reduceID",$reduceID);
    
    my $cType="unknown";
    $cType="Hi-C" if($hicModeFlag == 1);
    $cType="5C" if($fiveCModeFlag == 1);
    confess "invalid cType! ($cType)\n" if($cType eq "null");
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"cType",$cType);
    
    # calculate compute resources
    my $indexSize=`du -b $genomeDir`;
    chomp($indexSize);
    $indexSize=(split(/\t/,$indexSize))[0];
    # add 8GB to each, to account for working memory per tile size (8)
    my $indexSizeMegabyte = 8198+(ceil(($indexSize*1.25) / 1000000)); # scale index size by 1.25 fold
    my $splitSizeMegabyte = 8192+(ceil(((500*($splitSize/4))/1000)/1000)); # assume 500 bytes per line of side1+side2 SAM
    
    my $intervalSizeMegabyte=0;
    if($hicModeFlag == 1) {
        my $intervalSize=`du -b $restrictionFragmentPath`;
        chomp($intervalSize);
        $intervalSize=(split(/\t/,$intervalSize))[0];
        $intervalSizeMegabyte = (ceil(($intervalSize*1.25) / 1000000)); # scale interval size by 1.25 fold
    }
    my $mapMemoryNeededMegabyte=max($indexSizeMegabyte,$splitSizeMegabyte,$intervalSizeMegabyte);
    my $reduceMemoryNeededMegabyte=max(($splitSizeMegabyte*2),$intervalSizeMegabyte);
    
    print "\t\tcompute resources:\n";
    print "\t\t\tindexSizeMegabyte (".$indexSizeMegabyte."M) memory...\n";
    print "\t\t\tsplitSizeMegabyte (".$splitSizeMegabyte."M) memory...\n";
    print "\t\t\tintervalSizeMegabyte (".$intervalSizeMegabyte."M) memory...\n";
    print "\t\t\treduceMemoryNeededMegabyte (".$reduceMemoryNeededMegabyte."M) memory...\n";
    print "\t\t\tmapMemoryNeededMegabyte (".$mapMemoryNeededMegabyte."M) memory...\n";
    
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"indexSizeMegabyte",$indexSizeMegabyte);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"splitSizeMegabyte",$splitSizeMegabyte);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"intervalSizeMegabyte",$intervalSizeMegabyte);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"reduceMemoryNeededMegabyte",$reduceMemoryNeededMegabyte);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"mapMemoryNeededMegabyte",$mapMemoryNeededMegabyte);
    
    # calculate map time needed for LSF - assume split size is # lines not # reads (4 lines per read)
    my $mapTimeNeeded=((0.00004*$splitSize)-9.345)+240;
    $mapTimeNeeded=(((0.0003*$splitSize)+57)+720)*2 if($snpModeFlag == 1);
    # this linear approximation is done using mm9 - bowtie (excel)
    
    $shortMode = 1 if($debugModeFlag);
    my $genomeSizeFactor=max(1,($indexSizeMegabyte/3000));
    $genomeSizeFactor = log($genomeSizeFactor) if($genomeSizeFactor > 1);
    $mapTimeNeeded *= $genomeSizeFactor;
    $mapTimeNeeded *= 2*($iterativeMappingIterations/10) if($iterativeMappingFlag == 1);
    $mapTimeNeeded = 120 if($shortMode == 1);
    
    my $mapTimeNeededHour = max(2,floor($mapTimeNeeded/60));
    $mapTimeNeededHour = sprintf("%02d", $mapTimeNeededHour);
    my $mapTimeNeededMinute = ($mapTimeNeeded%60);
    $mapTimeNeededMinute = sprintf("%02d", $mapTimeNeededMinute);
    my $mapQueue="long";
    $mapQueue="short" if($shortMode == 1);
    $mapTimeNeeded=$mapTimeNeededHour.":".$mapTimeNeededMinute;
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"mapTimeNeeded",$mapTimeNeeded);
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"mapQueue",$mapQueue);
    print "\t\t\treduceResources\t$reduceQueue\t$reduceTimeNeeded\n";
    print "\t\t\tmapResources\t$mapQueue\t$mapTimeNeeded [".2*($iterativeMappingIterations/10)."]\n";
    
    # get log directory
    print "\t\tlogDirectory [$logDirectory]: ";
    my $userLogDirectory = <STDIN>;
    chomp($userLogDirectory);
    $userLogDirectory =~ s/\/$//;
    $logDirectory=$userLogDirectory if(-d($userLogDirectory));
    print "\t\t\t$logDirectory\n";
    system("mkdir -p $logDirectory") if(!(-d $logDirectory));
    croak "invalid log directory [$logDirectory]\n" if(!(-d($logDirectory)));
    $tmpConfigFileVariables=logConfigVariable($tmpConfigFileVariables,"logDirectory",$logDirectory);
    
    printConfigFile($configFileVariables,$tmpConfigFileVariables,$configFilePath);
        
    if(($hicModeFlag == 1) and ($fiveCModeFlag == 0)) { #HiC data
        print "\n";
        print "\t\tsubmitting HiC (reduceMem=$reduceMemoryNeededMegabyte:mapMem=$mapMemoryNeededMegabyte:tmp=$mapScratchSize)...\n";
        my $return=`bsub -n 2 -q $reduceQueue -R span[hosts=1] -R rusage[mem=$reduceMemoryNeededMegabyte:tmp=$reduceScratchSize] -W $reduceTimeNeeded -N -u $userEmail -J submitHiC -o $userHomeDirectory/lsf_jobs/LSB_%J.log -e $userHomeDirectory/lsf_jobs/LSB_%J.err $cMapping/utilities/tyler_submitHiC.sh $configFilePath`;
        chomp($return);
        print "\t\t$return\n";
        print "\n";
    } elsif(($hicModeFlag == 0) and ($fiveCModeFlag == 1)) { #5C data
        print "\n";
        print "\t\tsubmitting 5C (reduceMem=$reduceMemoryNeededMegabyte:mapMem=$mapMemoryNeededMegabyte:tmp=$mapScratchSize)...\n";
        my $return=`bsub -n 2 -q $reduceQueue -R span[hosts=1] -R rusage[mem=$reduceMemoryNeededMegabyte:tmp=$reduceScratchSize] -W $reduceTimeNeeded -N -u $userEmail -J submit5C -o $userHomeDirectory/lsf_jobs/LSB_%J.log -e $userHomeDirectory/lsf_jobs/LSB_%J.err $cMapping/utilities/submit5C.sh $configFilePath`;
        chomp($return);
        print "\t\t$return\n";
        print "\n";
    }
    
}
