#! perl -w

use strict;
use warnings;
use XML::LibXML;
use Cwd;

my $debug = 0;
my $home = getcwd;

my $usage = "wpNodeList.pl WP-GPML-FILE\n";
my $wpfile = shift or die $usage;


## Set up the XML parser to read the file
my $domain = XML::LibXML->load_xml(location => $wpfile);
my $xpc = XML::LibXML::XPathContext->new($domain);
$xpc->registerNs(sm => "http://pathvisio.org/GPML/2013a");
my $search = '/sm:Pathway/sm:DataNode';
my @nodes = $xpc->findnodes($search, $domain);
print STDERR ref @nodes, "\n";
print STDERR "Found ", scalar(@nodes), " data nodes.\n";
foreach my $node (@nodes) {
#    print STDERR "\t", ref($node), "\n";
    my $nid = $node->getAttribute("GraphId");
    my $label = $node->getAttribute("TextLabel");
    $label =~ s/\R//;
    my $type = $node->getAttribute("Type");
    print STDERR "$nid\t$label\t$type\n";
}
#chdir $home or die "Thomas Hardy was right\n";

exit;
__END__
