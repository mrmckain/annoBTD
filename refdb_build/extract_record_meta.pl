#!/usr/bin/perl
# Provenance for every cached GenBank record: who submitted it and what it is a
# copy of. Used to count annotation conventions once per independent research
# group rather than once per record.
#
# USAGE: extract_record_meta.pl <gb_dir> > record_meta.tsv
# OUT: accession, organism, species, source (INSDC accession a RefSeq NC_ record is
#      identical to, else "-"), sub_date, first_author, last_author, affiliation
# The submitter reference is the one whose JOURNAL is "Submitted (...)" from a
# non-NCBI address; RefSeq copies carry NCBI as submitter, so for those the
# authors come from the first reference and the affiliation is "-" (fill it from
# the source record downstream).
use strict; use warnings;
my ($dir) = @ARGV; die "usage: $0 <gb_dir>\n" unless $dir;
opendir(my $dh, $dir) or die; my @files = sort grep { /\.gb$/ } readdir $dh; closedir $dh;
print join("\t", qw(accession organism species source sub_date first_author last_author affiliation)), "\n";
for my $f (@files) {
    (my $acc = $f) =~ s/\.gb$//;
    open my $h, "<", "$dir/$f" or next;
    my ($org, $source, @refs, $cur, $in_comment, $comment) = ('', '-');
    while (<$h>) {
        chomp; last if /^FEATURES/;
        if (/^  ORGANISM\s+(.*)/) { $org = $1; next }
        if (/^COMMENT/) { $in_comment = 1; $comment = $'; next }
        if ($in_comment) { if (/^\S/) { $in_comment = 0 } else { $comment .= " $_"; next } }
        if (/^REFERENCE/) { $cur = { authors => '', journal => '', title => '' }; push @refs, $cur; $cur->{last} = ''; next }
        next unless $cur;
        if (/^  AUTHORS\s+(.*)/) { $cur->{authors} = $1; $cur->{last} = 'authors' }
        elsif (/^  TITLE\s+(.*)/) { $cur->{title} = $1; $cur->{last} = 'title' }
        elsif (/^  JOURNAL\s+(.*)/) { $cur->{journal} = $1; $cur->{last} = 'journal' }
        elsif (/^\s{12}(\S.*)/ && $cur->{last}) { $cur->{ $cur->{last} } .= " $1" }
    }
    close $h;
    if (defined $comment && $comment =~ /(?:identical to|derived from)\s+([A-Z]{1,2}_?\d+(?:\.\d+)?)/) { $source = $1 }
    # species = genus + epithet; drop infraspecific ranks and hybrids' extra tokens
    my @w = split /\s+/, $org; my $species = @w >= 2 ? "$w[0] $w[1]" : $org;
    $species = "$w[0] $w[2]" if @w >= 3 && $w[1] eq 'x';
    my ($sub) = grep { $_->{journal} =~ /^Submitted/ && $_->{journal} !~ /National Center for Biotechnology/ } @refs;
    my ($date, $aff) = ('-', '-');
    if ($sub) { ($date) = $sub->{journal} =~ /Submitted \((\S+)\)/; ($aff) = $sub->{journal} =~ /Submitted \([^)]*\)\s*(.*)/; }
    my $auth = $sub ? $sub->{authors} : (grep { $_->{authors} } @refs)[0] ? (grep { $_->{authors} } @refs)[0]{authors} : '';
    unless ($sub) { my ($ncbi) = grep { $_->{journal} =~ /^Submitted/ } @refs; ($date) = $ncbi->{journal} =~ /Submitted \((\S+)\)/ if $ncbi }
    my @a = grep { /\S/ } map { s/^\s+|\s+$//gr } split /,\s*(?=[A-Z][^,]*,)|\s+and\s+/, $auth;
    # GenBank authors are "Surname,I.N." separated by ", " and " and "
    my @names = $auth =~ /([A-Z][\w'\-]+(?:\s[\w'\-]+)*,[A-Z][\w.\-]*)/g;
    my ($first, $last) = (@names ? $names[0] : '-', @names ? $names[-1] : '-');
    $aff =~ s/\s+/ /g; $aff = '-' if $aff eq '';
    print join("\t", $acc, $org, $species, $source, $date // '-', $first, $last, $aff), "\n";
}
