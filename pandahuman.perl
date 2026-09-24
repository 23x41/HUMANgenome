#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw(time);
use Bio::SeqIO;

$| = 1;    # print immediately

# ---------------------------------------------------------------------------
# Paths — edit if your files moved
# ---------------------------------------------------------------------------
my $HUMAN = "$ENV{HOME}/Videos/ncbi_dataset/data/GCA_000001405.29/GCA_000001405.29_GRCh38.p14_genomic.fna";
my $PANDA = "$ENV{HOME}/Videos/ncbi_dataset/data/GCA_002007445.3/GCA_002007445.3_ASM200744v3_genomic.fna";

my $GAP_MIN = 10;          # NCBI-style contig split: runs of >=10 Ns
my $UNPLACED_EVERY = 2000; # progress ping for tiny scaffolds
my $LOG = "human_vs_panda_comparison.log";

open my $log, ">", $LOG or die "Cannot write $LOG: $!\n";
$log->autoflush(1);

sub say {
    my $msg = join("", @_);
    print $msg, "\n";
    print $log $msg, "\n";
}

sub comma {
    my $n = shift;
    $n = 0 unless defined $n;
    return $n if $n =~ /[^0-9.-]/;
    my $neg = $n < 0 ? "-" : "";
    $n = abs(int($n));
    $n =~ s/(?<=\d)(?=(\d{3})+$)/,/g;
    return "$neg$n";
}

sub pct {
    my ($num, $den) = @_;
    return "NA" unless $den;
    return sprintf("%.3f", 100 * $num / $den);
}

sub nx {
    my ($lens_ref, $total, $frac) = @_;
    return 0 unless $total;
    my $need = $total * $frac;
    my $acc = 0;
    my $L50 = 0;
    for my $L (@$lens_ref) {
        $acc += $L;
        $L50++;
        return ($L, $L50) if $acc >= $need;
    }
    return ($lens_ref->[-1] // 0, $L50);
}

sub classify {
    my ($id, $desc) = @_;
    my $h = lc("$id $desc");
    return "mitochondrion" if $h =~ /mitochondr/;
    return "alt_or_patch"  if $h =~ /\balt\b|alternate locus|fix patch|novel patch|\bpatch\b/;
    return "unlocalized"   if $h =~ /unlocaliz/;
    return "unplaced"      if $h =~ /unplaced|unknown/;
    return "chromosome"    if $h =~ /chromosome/;
    return "other";
}

sub analyze_file {
    my ($label, $path) = @_;
    die "Missing file for $label:\n  $path\n" unless -s $path;

    my $t0 = time;
    my $bytes = -s $path;
    say("=" x 78);
    say("SCANNING $label");
    say("file     : $path");
    say(sprintf("size     : %s bytes (%.2f GiB)", comma($bytes), $bytes/1024/1024/1024));
    say("=" x 78);

    my $in = Bio::SeqIO->new(-file => $path, -format => "fasta");

    my %tot = map { $_ => 0 } qw(
        seqs bases ACGT N GC AT lower
        gaps gap_bp contig_pieces
        AA TT GG CC AT_din TA CG GC
    );
    my %role;
    my @scaffold_len;
    my @contig_len;
    my @chroms;          # [id, len, gc, n, role]
    my $n_unplaced_seen = 0;

    while (my $seqobj = $in->next_seq) {
        my $id   = $seqobj->display_id // "unknown";
        my $desc = $seqobj->desc // "";
        my $seq  = $seqobj->seq // "";
        my $len  = length($seq);
        next unless $len;

        my $role = classify($id, $desc);

        my $N     = ($seq =~ tr/Nn//);
        my $lower = ($seq =~ tr/acgt//);
        my $A     = ($seq =~ tr/Aa//);
        my $C     = ($seq =~ tr/Cc//);
        my $G     = ($seq =~ tr/Gg//);
        my $T     = ($seq =~ tr/Tt//);
        my $acgt  = $A + $C + $G + $T;
        my $gc    = $G + $C;

        # dinucleotides on uppercase, skip N
        my $up = uc($seq);
        my %din;
        my $prev = "";
        for my $i (0 .. $len - 1) {
            my $b = substr($up, $i, 1);
            next unless $b =~ /[ACGT]/;
            if ($prev && $prev =~ /[ACGT]/) {
                $din{$prev.$b}++;
            }
            $prev = $b;
        }

        # contig split on N-runs
        my $gap_n = 0;
        my $gap_bp = 0;
        while ($up =~ /(N{$GAP_MIN,})/g) {
            $gap_n++;
            $gap_bp += length($1);
        }
        my @contigs = grep { length($_) } split(/N{$GAP_MIN,}/, $up);

        $tot{seqs}++;
        $tot{bases} += $len;
        $tot{ACGT}  += $acgt;
        $tot{N}     += $N;
        $tot{GC}    += $gc;
        $tot{AT}    += $A + $T;
        $tot{lower} += $lower;
        $tot{gaps}  += $gap_n;
        $tot{gap_bp}+= $gap_bp;
        $tot{contig_pieces} += scalar(@contigs);
        $tot{AA} += $din{AA} // 0;
        $tot{TT} += $din{TT} // 0;
        $tot{GG} += $din{GG} // 0;
        $tot{CC} += $din{CC} // 0;
        $tot{AT_din} += $din{AT} // 0;
        $tot{TA} += $din{TA} // 0;
        $tot{CG} += $din{CG} // 0;
        $tot{GC} = $tot{GC};          # base GC already stored
        $tot{CpG} += $din{CG} // 0;
        $tot{GpC} += $din{GC} // 0;

        $role{$role}{n}++;
        $role{$role}{bp} += $len;

        push @scaffold_len, $len;
        push @contig_len, map { length($_) } @contigs;

        my $rec = {
            id => $id, len => $len, gc => $gc, acgt => $acgt, N => $N,
            role => $role, lower => $lower, gaps => $gap_n, desc => $desc
        };

        if ($role eq "chromosome" || $role eq "mitochondrion") {
            push @chroms, $rec;
            say(sprintf(
                "  [%s] %-22s  %11s bp  GC %6s%%  N %6s%%  softmask %6s%%  gaps %d",
                $label, $id, comma($len),
                pct($gc, $acgt), pct($N, $len), pct($lower, $len), $gap_n
            ));
        }
        elsif ($role eq "alt_or_patch") {
            # GRCh38 has many alts; print only long ones
            if ($len >= 1_000_000) {
                say(sprintf("  [%s ALT] %-22s  %11s bp", $label, $id, comma($len)));
            }
        }
        else {
            $n_unplaced_seen++;
            if ($n_unplaced_seen % $UNPLACED_EVERY == 0) {
                my $elapsed = time - $t0;
                say(sprintf(
                    "  [%s] unplaced/other processed: %s   running total %s bp   %.1f min",
                    $label, comma($n_unplaced_seen), comma($tot{bases}), $elapsed/60
                ));
            }
        }
    }

    @scaffold_len = sort { $b <=> $a } @scaffold_len;
    @contig_len   = sort { $b <=> $a } @contig_len;
    my $scaf_total = 0; $scaf_total += $_ for @scaffold_len;
    my $ctg_total  = 0; $ctg_total  += $_ for @contig_len;

    my ($sN50, $sL50) = nx(\@scaffold_len, $scaf_total, 0.50);
    my ($sN90, $sL90) = nx(\@scaffold_len, $scaf_total, 0.90);
    my ($cN50, $cL50) = nx(\@contig_len,   $ctg_total,  0.50);

    my $cpg_oe = "NA";
    if ($tot{ACGT} > 0) {
        my $g_frac = ($tot{GC} / 2) / $tot{ACGT};   # approx G and C each ~ GC/2
        # better: we don't have separate G/C stored globally besides GC sum
        # Use dinucleotide-based O/E = CpG / (C*G) * n
        # Fallback: CpG O/E ≈ CpG_count / ((GC/2)^2 / ACGT) roughly weak
    }

    my $elapsed = time - $t0;
    say("");
    say("SUMMARY $label  (elapsed " . sprintf("%.1f", $elapsed/60) . " min)");
    say(sprintf("  sequences          %s", comma($tot{seqs})));
    say(sprintf("  total length       %s bp", comma($tot{bases})));
    say(sprintf("  ACGT (ungapped)    %s bp", comma($tot{ACGT})));
    say(sprintf("  N / gap bases      %s bp  (%s%% of assembly)", comma($tot{N}), pct($tot{N}, $tot{bases})));
    say(sprintf("  spanned gaps>=%dN  %s gaps, %s bp", $GAP_MIN, comma($tot{gaps}), comma($tot{gap_bp})));
    say(sprintf("  GC of ACGT         %s%%", pct($tot{GC}, $tot{ACGT})));
    say(sprintf("  AT of ACGT         %s%%", pct($tot{AT}, $tot{ACGT})));
    say(sprintf("  softmasked         %s%% (lowercase = repeats in NCBI eukaryotic FNA)", pct($tot{lower}, $tot{bases})));
    say(sprintf("  scaffold N50/L50   %s / %s", comma($sN50), comma($sL50)));
    say(sprintf("  scaffold N90/L90   %s / %s", comma($sN90), comma($sL90)));
    say(sprintf("  longest scaffold   %s bp", comma($scaffold_len[0] // 0)));
    say(sprintf("  contig N50/L50     %s / %s   (split on >=%d Ns)", comma($cN50), comma($cL50), $GAP_MIN));
    say(sprintf("  contig pieces      %s", comma($tot{contig_pieces})));
    say(sprintf("  CpG dinucleotides  %s", comma($tot{CpG} // 0)));
    say(sprintf("  GpC dinucleotides  %s", comma($tot{GpC} // 0)));
    if (($tot{CpG}//0) && ($tot{GpC}//0)) {
        say(sprintf("  CpG/GpC ratio      %.4f   (<1 is typical mammalian CpG suppression)",
            $tot{CpG} / $tot{GpC}));
    }

    say("  by sequence role:");
    for my $r (sort keys %role) {
        say(sprintf("    %-16s %7s seqs   %s bp  (%s%%)",
            $r, comma($role{$r}{n}), comma($role{$r}{bp}), pct($role{$r}{bp}, $tot{bases})));
    }

    @chroms = sort { $b->{len} <=> $a->{len} } @chroms;
    if (@chroms) {
        say("  placed chromosomes / organelles:");
        for my $c (@chroms) {
            say(sprintf("    %-24s %12s bp  GC %6s%%  N %6s%%",
                $c->{id}, comma($c->{len}), pct($c->{gc}, $c->{acgt}), pct($c->{N}, $c->{len})));
        }
    }

    return {
        label => $label,
        tot   => \%tot,
        role  => \%role,
        scaf_len => \@scaffold_len,
        sN50 => $sN50, sL50 => $sL50, sN90 => $sN90,
        cN50 => $cN50, cL50 => $cL50,
        chroms => \@chroms,
        seconds => $elapsed,
    };
}

# ---------------------------------------------------------------------------
say("Human GRCh38.p14  vs  giant panda ASM200744v3");
say("Live comparative scan. Log also written to $LOG");
say("Note: GRCh38.p14 genomic.fna includes alt loci + patches, so total bp");
say("is larger than the ~3.1 Gb haploid primary assembly.");
say("Panda file includes 21 chromosomes plus ~73,000 unplaced scaffolds.");
say("");

my $h = analyze_file("HUMAN", $HUMAN);
my $p = analyze_file("PANDA", $PANDA);

# ---------------------------------------------------------------------------
say("");
say("=" x 78);
say("SIDE-BY-SIDE COMPARISON");
say("=" x 78);

my @rows = (
    ["metric",                          "HUMAN", "PANDA"],
    ["sequences",                       comma($h->{tot}{seqs}), comma($p->{tot}{seqs})],
    ["total bp",                        comma($h->{tot}{bases}), comma($p->{tot}{bases})],
    ["ungapped ACGT bp",                comma($h->{tot}{ACGT}), comma($p->{tot}{ACGT})],
    ["N bases",                         comma($h->{tot}{N}), comma($p->{tot}{N})],
    ["N percent",                       pct($h->{tot}{N}, $h->{tot}{bases})."%", pct($p->{tot}{N}, $p->{tot}{bases})."%"],
    ["GC percent (of ACGT)",            pct($h->{tot}{GC}, $h->{tot}{ACGT})."%", pct($p->{tot}{GC}, $p->{tot}{ACGT})."%"],
    ["softmasked percent",              pct($h->{tot}{lower}, $h->{tot}{bases})."%", pct($p->{tot}{lower}, $p->{tot}{bases})."%"],
    ["scaffold N50",                    comma($h->{sN50}), comma($p->{sN50})],
    ["scaffold L50",                    comma($h->{sL50}), comma($p->{sL50})],
    ["scaffold N90",                    comma($h->{sN90}), comma($p->{sN90})],
    ["contig N50 (>=10 N split)",       comma($h->{cN50}), comma($p->{cN50})],
    ["CpG count",                       comma($h->{tot}{CpG}//0), comma($p->{tot}{CpG}//0)],
    ["GpC count",                       comma($h->{tot}{GpC}//0), comma($p->{tot}{GpC}//0)],
    ["chromosomes+MT listed",           comma(scalar @{$h->{chroms}}), comma(scalar @{$p->{chroms}})],
);

for my $r (@rows) {
    say(sprintf("%-28s %20s %20s", @$r));
}

if ($h->{tot}{bases} && $p->{tot}{bases}) {
    my $size_ratio = $p->{tot}{bases} / $h->{tot}{bases};
    my $gc_delta   = 100*($p->{tot}{GC}/($p->{tot}{ACGT}||1) - $h->{tot}{GC}/($h->{tot}{ACGT}||1));
    say("");
    say(sprintf("Panda assembly is %.1f%% the size of this human FNA (alts included).", 100*$size_ratio));
    say(sprintf("GC difference (panda - human): %+.3f percentage points.", $gc_delta));
    say("Expect ~2.4 Gb panda vs ~3.1 Gb human primary; extra human bp is mostly alts/patches.");
}

say("");
say("What this did NOT do (on purpose):");
say("  - whole-genome pairwise alignment (use minimap2 / LASTZ, not BioPerl)");
say("  - gene orthology (need protein/CDS FASTA + BLASTP/OrthoFinder)");
say("  - Mash / k-mer distance (run if installed, see below)");
say("");

# Optional fast global distance if mash is on PATH
if (qx(command -v mash 2>/dev/null)) {
    say("Found mash — computing k-mer distance (this is the useful whole-genome number)...");
    system("mash", "sketch", "-s", "10000", "-o", "human.msh", $HUMAN);
    system("mash", "sketch", "-s", "10000", "-o", "panda.msh", $PANDA);
    say("mash dist:");
    system("mash", "dist", "human.msh", "panda.msh");
} else {
    say("Install mash for a 1-line genome-wide distance:");
    say("  conda install -c bioconda mash");
    say("  mash sketch -o human $HUMAN");
    say("  mash sketch -o panda $PANDA");
    say("  mash dist human.msh panda.msh");
}

say("");
say("Optional chromosome mapping after this script finishes:");
say("  minimap2 -x asm20 -t 8 \\");
say("    $HUMAN \\");
say("    $PANDA > panda_on_human.paf");
say("");
say("Done. Full text is in $LOG");
