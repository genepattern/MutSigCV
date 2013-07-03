use strict;
use warnings;

while (my $line = <>) {
    $line =~ s/\[[\d]{2}:[\d]{2}:[\d]{2}\]/\[timestamp\]/g;
    $line =~ s/\/[A-Za-z\/0-9._-]+\//\[system_path\]\//g;
    
    print $line;
}
