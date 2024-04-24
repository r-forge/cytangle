#! perl -w

use strict;
use warnings;
use XML::LibXML;
use Cwd;

my $debug = 0;
my $home = getcwd;
my $pathwayDir = "D:/Projects/WayFindR/HS-WP"; # WikiPathways for Homo Sapiens

my $arrowFile = "allArrows.tsv";
open(OUTPUT, ">$arrowFile") or die "Unable to create '$arrowFile'!\n";


my $typeFile = "arrowTypes.txt";
my @types = ();
if (-e $typeFile) {
    open(TYP, "<$typeFile") or die "Unable to open '$typeFile'!\n";
    while (my $line = <TYP>) {
	chomp $line;
	print STDERR "Recording type '$line'!\n" if $debug;
	push @types, $line;
    }
    close(TYP);
    print STDERR "Found ", scalar(@types), " arrow types.\n";
} else {
    print STDERR "It seems that '$typeFile' does not exist in", , "!\n";
}

if ($#types) {
    print OUTPUT "PathID\t", join("\t", @types), "\n"; # header line
}

chdir $pathwayDir or die "Cannot change directories to '$pathwayDir'!\n";
my %arrowTypes = ();
my @allFiles = glob("*.gpml");
foreach my $filename (@allFiles) {
    ## get the short version of the pathwasy identifier out of the file name
    my $wpi = index($filename, "_WP");
    my $wp = substr($filename, 1 + $wpi);
    $wp =~ s/\.gpml//;
    my ($pid, $date) = split(/_/, $wp);
    print STDERR "$wp\t$pid\t$date\n";
    
    ## Set up the XML parser to read the file
    my $domain = XML::LibXML->load_xml(location => $filename);
    my $xpc = XML::LibXML::XPathContext->new($domain);
    $xpc->registerNs(sm => "http://pathvisio.org/GPML/2013a");
    my $search = '/sm:Pathway/sm:Interaction/sm:Graphics/sm:Point/@ArrowHead';
    my $counter = $xpc->findnodes($search, $domain)->size;
    print STDERR "Found $counter points with arrowheads.\n";
    
    my %arrows = ();
    foreach my $arrow ($xpc->findnodes($search, $domain)) {
	my $name = $arrow->localname;
	my $value = $arrow->to_literal();
	$arrows{$value}++;
	$arrowTypes{$value}++;
    }

    if ($#types) {
	print OUTPUT "$pid\t";
	my @counts = map{ exists($arrows{$_}) ? $arrows{$_} : 0} @types;
	print OUTPUT join("\t", @counts), "\n";
    } else {
	foreach my $head (keys %arrows) {
	    print STDOUT "Saw '$head' in $pid $arrows{$head} times.\n";
	}
    }
}
close(OUTPUT);
chdir $home or die "Thomas Hardy was right\n";
print STDOUT "\n\nTotal:\n";
unless (-e $typeFile) {
    open(TYP,">$typeFile") or die "Unable to creaste '$typeFile'!\n";
    my @sorted = sort keys(%arrowTypes);
    foreach my $head (@sorted) {
	print STDOUT "Saw '$head' in total $arrowTypes{$head} times.\n";
	print TYP "$head\t$arrowTypes{$head}\n";
    }
    close(TYP);
}


exit;
__END__
