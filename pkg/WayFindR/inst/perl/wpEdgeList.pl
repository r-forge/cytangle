#! perl -w

use strict;
use warnings;
use XML::LibXML;
use Cwd;

my $debug = 0;
my $home = getcwd;

my $usage = "wpEdgeList.pl WP-GPML-FILE\n";
my $wpfile = shift or die $usage;


## Set up the XML parser to read the file
my $domain = XML::LibXML->load_xml(location => $wpfile);
my $xpc = XML::LibXML::XPathContext->new($domain);
$xpc->registerNs(sm => "http://pathvisio.org/GPML/2013a");
my $search = '/sm:Pathway/sm:Interaction';
my @edges = $xpc->findnodes($search, $domain);
print STDERR "Found ", scalar(@edges), " data edges in `$wpfile`.\n";
foreach my $edge (@edges) {
    my $eid = $edge->getAttribute("GraphId");
    my @pts = $xpc->findnodes("./sm:Graphics/sm:Point", $edge);
    if (scalar(@pts) != 2) {
	print STDERR "Got wrong number of points in an interaction in '$wpfile'!\n";
    }
    my $counter = 0;
    foreach my $point (@pts) {
	++$counter;
	my $node = $point->getAttribute("GraphRef");
	my $arrow = $point->getAttribute("ArrowHead");
	if ($counter == 1) {
	    if($arrow) {
		print STDERR "Source node has an arrow head!\n";
	    } else {
		print STDOUT "$node\t";
	    }
	} else {
	    $arrow = "none" unless $arrow;
	    print STDOUT "$node\t$arrow\n";
	}
    }
}
#chdir $home or die "Thomas Hardy was right\n";

exit;
__END__
