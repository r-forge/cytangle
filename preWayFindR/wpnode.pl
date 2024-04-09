#! perl -w

use strict;
use warnings;
use XML::LibXML;
use Cwd;

my $debug = 0;
my $home = getcwd;
my $pathwayDir = "HS-WP"; # WikiPathways for Homo Sapiens
my $typeFile = "nodeTypes.txt";

chdir $pathwayDir or die "Cannot change directories to '$pathwayDir'!\n";
my %nodeTypes = ();
my @allFiles = glob("*.gpml");
foreach my $filename (@allFiles) {
    ## get the short version of the pathwasy identifier out of the file name
    my $wpi = index($filename, "_WP");
    my $wp = substr($filename, 1 + $wpi);
    $wp =~ s/\.gpml//;
    my ($pid, $date) = split(/_/, $wp);

    ## Set up the XML parser to read the file
    my $domain = XML::LibXML->load_xml(location => $filename);
    my $xpc = XML::LibXML::XPathContext->new($domain);
    $xpc->registerNs(sm => "http://pathvisio.org/GPML/2013a");
    my $search = '/sm:Pathway/sm:DataNode/@Type';
    my $counter = $xpc->findnodes($search, $domain)->size;
    print STDERR "Found $counter items with node types in '$wp'.\n";

    ## Add the type of each node to the hash
    foreach my $node ($xpc->findnodes($search, $domain)) {
	my $value = $node->to_literal();
	$nodeTypes{$value}++;
    }

}

chdir $home or die "Thomas Hardy was right\n";
print STDOUT "\n\nTotal:\n";
unless (-e $typeFile) {
    open(TYP,">$typeFile") or die "Unable to creaste '$typeFile'!\n";
#    my @sorted = sort {$nodeTypes{$a} <=> $nodeTypes{$b}} keys(%nodeTypes);
    my @sorted = sort keys(%nodeTypes);
    foreach my $head (@sorted) {
	print STDOUT "Saw '$head' in total $nodeTypes{$head} times.\n";
	print TYP "$head\t$nodeTypes{$head}\n";
    }
    close(TYP);
}


exit;
__END__
