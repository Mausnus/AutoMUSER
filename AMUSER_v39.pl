#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use File::Temp qw(tempfile);

# GenBank annotation-transfer stores (added 2026-10-29).
our %GBfrag;      # uc(fragment sequence) => arrayref of { name, features=>[...] }
our %GBplasmid;   # { name, topology, features=>[...] } from a -dpgb GenBank vector
our $dp_cassette_warning = "";   # set when -dp cassette not found; rendered as a big red banner atop the report

my $version = "1.0";

###################################################################
#      	                 AMUSER version 1                         #
#																  #
#                   main program developed by:					  #
#				 Lars Ronn Olsen (lro@binf.ku.dk) 				  #
#																  #
#			 subroutines for input handling developed by:		  #
#			   Jakob Berg Jespersen (jbeje@bio.dtu.dk)			  #
#																  #
#				    html interface developed by:				  #
#			 Frederik Otzen Bagger (frederik@binf.ku.dk)		  #
#																  #
#						web server set up by:					  #
#			    Rasmus Wernersson (raz@intomics.com)			  #
#																  #
#                first version written: 15.12.2011                #
#                     latest update: 17.12.2013                   #
###################################################################


###################################################################
#####                       MAIN PROGRAM                      #####
###################################################################


#First, input/output files and other parameters are defined:

# -i = input file name (if eq "+" then file will be piped)
# -r = report file name (if eq "+" then file will be piped)
# -o = fasta output file name (if eq "+" then file will be piped)
# -h = html output file name
# -oname = custom base name for the three outputs, without a date/time stamp
#          (<name>_report.html, <name>_insert.gb|_plasmid.gb, <name>_primers.csv)
# -workdir = subdirectory path for saving outputfiles (accepts only [DIRNAME], dir is created if non-existent)
# -baseurl = URL for the individual batch report files
# -cstandard = choice of standard cassette (name of cassette)
# -cmaninfor = cassette forward strand from manual input
# -cmaninrev = cassette reverse strand from manual input
# -cnick = choice of nicking site (name of enzyme)
# -crestr = choice of restriction site (name of enzyme)
# -db1 = directional bases 1
# -db2 = directional bases 2
# -cname = cassette name
# -tm = Tm optimization (numeric value 30-90 or "eq")
# -poly = polymerase driving the Tm model: phusion (default) or taq. Sets the
#         salt and primer concentration as well as the salt correction, because
#         all three are arguments to the same equation. There is no buffer
#         sub-choice: NEB's model gives Phusion HF and GC buffers an identical
#         Tm. phusion_hf and phusion_gc are accepted as aliases of phusion.
# -sc = Na+ concentration in mM (overrides the -poly preset; 222 for Phusion,
#       55 for Taq)
# -pc = primer concentration in the reaction tube, in uM, as pipetted
#       (overrides the -poly preset; 0.5 for Phusion, 0.2 for Taq). NEB
#       divides this by 4 for the Phusion family internally and so does this
#       script, so enter what goes in the tube, not a pre-divided figure.
# -co = circular output if "-co y"
# -verbose = verbose messages for debugging if "-verbose y"



my @fastaoutarray;		#array of fastaoutputs
my @fastabatcharray;	#array of fastaoutputs as printed for batch assembly
my @outreportarray;		#array of report output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$year = $year + 1900;
$mon++;
$mon = sprintf ("%02d", $mon);
$mday = sprintf ("%02d", $mday);
$min = sprintf ("%02d", $min);
$hour = sprintf ("%02d", $hour);
$sec = sprintf ("%02d", $sec);



###################################################################
#assigning cmd line options


###################################################################
# Primer-evaluation constants. One definition each, shared by the
# candidate ranking (sub optbranalysis) and by the report's verdicts,
# so the two can no longer drift apart.
# Binding-region length window. Every candidate length from MIN to MAX
# inclusive is generated, scored and ranked; the list is ordered shortest-first.
my $BR_LEN_MIN = 18;     my $BR_LEN_MAX = 28;    # nt, inclusive
my $TM_OPT_MIN = 55;     my $TM_OPT_MAX = 72;    # "optimal range" shown in the report
my $TM_TOL_MIN = 50;     my $TM_TOL_MAX = 77;    # wider band: outside it is penalised twice

# Tm model (2026-11-14). The thermodynamics are NEB's: SantaLucia 1998 unified
# nearest-neighbour parameters with terminal and symmetry corrections, plus the
# salt correction the NEB Tm Calculator applies for the chosen polymerase.
# Validated against NEB's own calculator over 509 primers: MAE 0.0025 C,
# max 0.005 C (the recording precision of the reference), for both products.
# Replaces the Sugimoto 1996 parameters with a 7.21*ln[Na+] salt term, which
# read 7.0 C low against NEB Phusion HF and 3.5 C low against NEB Taq, with a
# GC-dependent error spanning 10 C. The 55-72 / 50-77 bands above are physical
# statements about PCR and are deliberately unchanged.
#
# Each entry: [monovalent salt mM, primer conc nM IN THE REACTION, salt
# correction]. The concentrations are arguments to the Tm equation, not
# cosmetic defaults. The primer figure is the one you pipette -- 2.5 uL of a
# 10 uM working stock in 50 uL is 500 nM -- matching what NEB's own interface
# asks for. NEB then divides it by 4 for the Phusion family before it enters
# the equation, and so does this script (see the $pc block), so neither the
# table nor -pc ever carries a pre-divided number.
# -sc / -pc still override either value.
# One entry per polymerase, no buffer sub-choice: NEB's own model gives the
# Phusion GC buffer the same 222 mM, the same /4 and the same correction as HF,
# so a "which buffer" control would imply a Tm difference that does not exist.
# phusion_gc is kept as an accepted alias of phusion_hf, not as a separate row.
my %TM_POLY = (
	phusion    => [222, 500, 'schildkraut'],   # Phusion, HF or GC buffer (default)
	taq        => [ 55, 200, 'owczarzy'   ],   # Taq, Standard Taq Buffer
);
my %TM_POLY_ALIAS = (phusion_hf => 'phusion', phusion_gc => 'phusion');
my $TM_POLY = 'phusion';                       # set from -poly during option parsing

# SantaLucia 1998 unified nearest-neighbour parameters (kcal/mol, cal/mol/K).
my %NN_H = (AA=>-7.9,TT=>-7.9,AT=>-7.2,TA=>-7.2,CA=>-8.5,TG=>-8.5,GT=>-8.4,AC=>-8.4,
            CT=>-7.8,AG=>-7.8,GA=>-8.2,TC=>-8.2,CG=>-10.6,GC=>-9.8,GG=>-8.0,CC=>-8.0);
my %NN_S = (AA=>-22.2,TT=>-22.2,AT=>-20.4,TA=>-21.3,CA=>-22.7,TG=>-22.7,GT=>-22.4,AC=>-22.4,
            CT=>-21.0,AG=>-21.0,GA=>-22.2,TC=>-22.2,CG=>-27.2,GC=>-24.4,GG=>-19.9,CC=>-19.9);
my %NN_TERM_H = (A=>2.3, T=>2.3, G=>0.1,  C=>0.1 );
my %NN_TERM_S = (A=>4.1, T=>4.1, G=>-2.8, C=>-2.8);
my $NN_R = 1.987;

# IUPAC degenerate bases. A binding region may legitimately contain them: the
# GenBank reader accepts RYKMBVDHNWS (see the -igb validation) and the
# reverse-complement at the top of the pipeline maps them. The Tm of a
# degenerate primer is not one number -- the pool contains every concrete
# variant -- so the model expands it, DISPLAYS the min-max range, and uses the
# MINIMUM for every downstream comparison. Both are what NEB does: its report
# shows "min - max", and getAnnealTemp is handed the minimum of each primer's
# expansion, not the maximum and not the mean.
my %IUPAC = (
	A=>'A', C=>'C', G=>'G', T=>'T', U=>'T',
	R=>'AG', Y=>'CT', S=>'CG', W=>'AT', K=>'GT', M=>'AC',
	B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT',
);
# NEB's own limit, from its validateInput: nonCanonicalCount(seq) < 4, i.e. at
# most 3 ambiguous bases, counted over "wsrymkbdhvn". A primer past it is
# flagged invalidseq and NEB shows "---" for both Tm and Ta rather than
# guessing. Here the length is simply marked unusable, so tm_lengths drops it
# and the no-usable-region error reports the cause.
my $IUPAC_MAX_AMBIG = 3;
# Lengths sub Tm had to drop because the binding region could not be expanded
# within that cap, keyed "fragment\tstrand". Read only by the no-usable-region
# error, so it can name degeneracy as the cause instead of blaming the overhang.
my %NN_TOO_DEGENERATE;
my $DIMER_SCAN_CAP   = 12;   # longest 3' suffix considered when looking for a
                             # dimer; a run longer than this is already far over
                             # any threshold, and capping bounds the work.
my $DIMER_MAX_END3   = 6;    # primer dimer: the heaviest 3'-anchored duplex a
                             # primer pair may form before it is flagged, scored
                             # G:C = 2 and A:T = 1, so 6 is three consecutive
                             # G/C pairs or four mixed. A duplex that pairs a
                             # 3'-terminal base is the one a polymerase can
                             # extend, which is what turns an interaction into
                             # primer-dimer product; an internal duplex of the
                             # same strength mostly does not.
                             # Measured over 565 chosen pairs: >4 flags 60.4%,
                             # >5 36.3%, >6 18.9%, >7 12.7%, >8 6.4%.
my $TM_PAIR_MAX_DIFF = 2;                        # FW/RV Tm difference still called a match
# -tm/-tm eq: lengths this close to the target form the pool the ranking
# chooses from (0 = take only the closest, the behaviour before 2026-11-13).
# The two modes get different tolerances because their targets mean different
# things.
# Each is a LADDER, tried in order per fragment: the first step that yields at
# least one FW/RV pair is the candidate pool the ranking then scores. A single
# value behaves exactly like the old single tolerance. If no step yields a
# pair, the fallback below applies (closest length per strand, no ranking).
#
# Why a ladder rather than one number: the two ends of the trade are wanted at
# once. A tight window is what hits the requested Tm, but it leaves most
# fragments with an empty pool, and an empty pool means the clamp and every
# other criterion is never consulted for them. Escalating per fragment keeps
# the tight window for fragments that can meet it and gives the rest a ranked
# pool near the target instead of a single unranked pick.
my @TM_TARGET_TOL_EQ_STEPS    = (0.5, 1.0, 1.5, 2.0);
                                 # -tm eq: the target minimises the realised Tm
                                 # spread (see optbrTm), so this window is only
                                 # how far the OTHER criteria may pull primers
                                 # off it, and it trades against $W_CLAMP.
                                 # Swept 0-3 C as a flat window over 20 fusions,
                                 # median within-run spread / runs spanning
                                 # >3 C / 3'-terminal G/C rate:
                                 #   0-0.5  0.70 C / 3 / 52.8%
                                 #   1.0    1.20 C / 3 / 65.1%
                                 #   1.5    1.70 C / 3 / 74.5%
                                 #   2.0    1.75 C / 6 / 80.2%
                                 # Below ~0.5 nothing changes - adjacent lengths
                                 # differ by a median of 1.30 C, so a window
                                 # under half a step admits one length.
                                 #
                                 # The later rungs do NOT come free here, unlike
                                 # under -tm <N>. A fragment that misses 0.5
                                 # used to take the fallback, which picks the
                                 # length CLOSEST to the target - and closest to
                                 # target is exactly what minimises spread. Give
                                 # it a ranked pool instead and the ranking will
                                 # trade Tm distance for a clamp, because
                                 # distance to target is not one of the eight
                                 # criteria. Measured against a lone 0.5 rung:
                                 # 4 primers change, 3 clamps gained and 0 lost,
                                 # but 3 of 29 runs widen (up to +0.80 C) and
                                 # none tighten; mean spread 1.63 -> 1.69 C,
                                 # median unchanged at 0.70. Monotone-good on
                                 # the clamp, monotone-bad on the spread.
                                 #
                                 # Escalating from 1.0 instead is a measured
                                 # no-op: only 8.8% of fragments fall back at a
                                 # flat 1.0 (against 75% under -tm <N> at 0.5),
                                 # because this target is chosen to sit near
                                 # every fragment's achievable range. A (1.0,
                                 # 1.5) ladder changes zero primers.

my @TM_TARGET_TOL_FIXED_STEPS = (0.5, 1.0, 1.5, 2.0);
                                 # -tm <N>: the target is the number the user
                                 # typed, so hitting it comes first; 0.5 is the
                                 # accuracy floor (nothing below it improves on
                                 # 46.5% of primers within 0.5 C at -tm 55,
                                 # since at 0 a length would have to match the
                                 # target to floating point). The later steps
                                 # exist only to rescue fragments that cannot
                                 # meet 0.5: at a flat 0.5 the pool was empty
                                 # for 75% of fragments and held a single pair
                                 # for another 15%, so the clamp was consulted
                                 # for 10% of them.
                                 #
                                 # Tightening does NOT tighten the spread
                                 # between a run's primers - measured flat at
                                 # every fixed tolerance (median 2.60 C at
                                 # -tm 55, 2.90 C at -tm 60). A run's spread is
                                 # set by its hottest and coldest fragments, and
                                 # those are the ones that cannot reach the
                                 # target at any length and fall back whatever
                                 # this window says.
                                 #
                                 # -tm 65 is largely unreachable in an 18-28 nt
                                 # window: median |Tm - 65| is 5.10 C and 74-87%
                                 # of fragments fall back at every tolerance, so
                                 # that target degrades into "as hot as
                                 # possible". A $BR_LEN_MIN/$BR_LEN_MAX limit,
                                 # not a tolerance one.

# Ranking weights (higher = worse). Each is applied per strand unless noted.
my $W_TM_OPT   = 8;      # Tm outside the optimal range. Raised from 5 on
                         # 2026-11-15 after the first weight sweep on the NEB
                         # Tm scale (296 configurations x 400 primers, both
                         # polymerases). Under the old model nothing ever sat
                         # above 72 C, so the upper Tm penalty never fired and
                         # this weight acted only at the cold end, where
                         # $W_TM_TOL was doing the work anyway. On the NEB
                         # scale the ceiling is live and the term is two-sided.
                         # Phusion, $W_CLAMP held at 8, clamp / Tm band /
                         # above 72 C / below 50 C:
                         #    0   91.5% / 64.8% / 4.5% / 1.8%
                         #    2   91.0% / 75.5% / 4.0% / 1.8%
                         #    5   89.2% / 81.5% / 2.5% / 1.8%   <- previous
                         #    8   89.0% / 82.2% / 1.8% / 1.8%   <- applied
                         #   12   85.8% / 85.8% / 1.2% / 1.8%
                         #   16   85.8% / 85.8% / 1.2% / 1.8%
                         # 5 -> 8 is free: band +0.7, above-72 -0.7, clamp
                         # -0.2, cold tail and dimer rate unchanged, 5 primers
                         # of 400 changed. Taq agrees independently (band
                         # 62.0% -> 63.2%, above-72 2.2% -> 2.0%, clamp
                         # unchanged at 86.2%, 8 primers changed). The curve
                         # saturates at 12, where 85.8% in band costs 3.4
                         # points of clamp - measured and declined, since the
                         # clamp was chosen over Tm fit deliberately.
my $W_TM_TOL   = 9;      # ...and outside the wider band as well. Raised from 4
                         # together with $W_CLAMP 5 -> 8, and for that reason.
                         # In AT-rich sequence every candidate is already
                         # outside the 55-72 C band, so $W_TM_OPT is a constant
                         # there and carries no selection information; the only
                         # Tm term still discriminating is this one. At 4 it was
                         # cheaper to leave the 50-77 C band on both strands
                         # (2 x 4 = 8) than to forgo two clamps (2 x 8 = 16), so
                         # the selector shortened AT-rich primers by 6-7 nt and
                         # dropped them ~10 C to buy a terminal G: on fus3_at,
                         # four primers went from 50-52 C to 39.5-44.8 C and the
                         # run's Tm spread went 2.2 -> 11.9 C. At 9 that trade is
                         # unprofitable and does not happen. Measured over 46
                         # inputs x 3 Tm modes: default-mode clamp 87.5% (93.3%
                         # at the old weights, 98.3% with $W_CLAMP 8 alone),
                         # -tm eq clamp 82.5% (70.0% before), Tm band membership
                         # unchanged at 63.3%, worst within-run Tm spread
                         # 14.70 -> 12.20 C, and primers below 50 C 14 -> 5.
                         # The curve is flat from 9 to 16 - 9 is where the trade
                         # flips, not a tuned value.
my $W_GC       = 1;      # GC content outside 40-60%
my $W_DIMER    = 10;     # primer dimer: a 3'-anchored duplex heavier than
                         # $DIMER_MAX_END3, over the cross-dimer and both
                         # self-dimers (sub primerdimer). Swept 0-20 against
                         # $W_CLAMP over 36 combinations, default-Tm mode,
                         # flagged pairs / clamp / Tm band / >3 G/C in last 5:
                         #    0   13.3% / 88.9% / 65.5% / 3.1%
                         #    3    8.0% / 89.4% / 65.5% / 3.1%
                         #    6    7.1% / 89.4% / 65.5% / 3.5%
                         #   10    0.0% / 88.9% / 64.6% / 5.3%   <- applied
                         #   15    0.0% / 88.9% / 64.6% / 5.3%
                         #   20    0.0% / 88.9% / 64.6% / 5.3%
                         # All-or-nothing: flat from 3 to 6, then every flagged
                         # pair disappears at 10 and nothing moves above it. So
                         # a dimer-free candidate almost always exists and
                         # merely costs more than 3 points; taking all of them
                         # costs 0.5 points of clamp, 0.9 of band, 2.2 of the
                         # (weak) last-5 criterion and 2 more primers below
                         # 50 C, over 27 changed primers.
                         #
                         # Acts in default-Tm mode only. Under -tm eq and
                         # -tm <N> the tolerance ladder leaves candidates inside
                         # the first 0.5 C rung too alike for any weight to
                         # separate them: -tm 55 is flat at 16.8% from weight 3
                         # to 20, -tm 60 moves 23.0% -> 22.1%, -tm eq not at all.
my $W_HOMOLOGY = 0;      # intra-primer homology. DEAD CRITERION - never fires;
                         # see the note below. Was 6 until 2026-11-15, when the
                         # NEB-scale sweep confirmed it structurally inert on a
                         # second, harsher corpus: driving it to 0 or to 12
                         # changes zero of 400 primers under either polymerase.
                         # Zeroed so the table no longer reads as though the
                         # scorer has a working hairpin check. This changes no
                         # output whatsoever - it is a readability fix, not a
                         # tuning change. Restoring a live homology criterion
                         # means designing a reachable threshold, not setting
                         # this number back.
my $W_CLAMP    = 8;      # no 3'-terminal G/C clamp. Swept 0-100 over 46 inputs
                         # x 3 Tm modes (360 primers): 8 gives 98.3% of
                         # default-mode primers a 3'-terminal G/C against 93.3%
                         # at 5, with Tm band membership unchanged at 63.3%;
                         # the only cost is "more than 3 G/C in the last 5"
                         # rising 2.5% -> 4.2%. The curve saturates at 10
                         # (99.2%) and nothing moves out to a weight of 100,
                         # because the remaining primers have no candidate
                         # length ending in G or C at all. Two consequences: a
                         # hard constraint with fallback would buy nothing that
                         # a weight of 10 does not, and promoting the clamp to
                         # a lexicographic tie-break would buy nothing at all,
                         # since the minimum-penalty pair already takes every
                         # clamp that costs zero. Going past 8 costs 1.6 points
                         # of band membership for 0.9 points of clamp.
my $W_GCCOUNT  = 7;      # more than 3 G/C in the last 5 bases. Swept 0-14
                         # against $W_CLAMP: any non-zero value holds the
                         # violation rate near 2-4% (12.5% at 0), and the exact
                         # value barely matters. 10 would restore 1.7% at
                         # $W_CLAMP = 8, for 1.6 points of Tm band.
my $W_POLYN    = 1;      # polyN stretch

# $W_HOMOLOGY does nothing whatsoever, and is now 0 for that reason. Measured
# twice: over 1612 candidate pairs from 60 fragments, and again on 2026-11-15
# over 400 primers x both polymerases, where driving it to 0 or to 12 leaves
# every primer byte-identical. The cause is structural, not a property of
# either corpus. sub homology enumerates sub-fragments of a primer,
# reverse-complements them, and flags the candidate only when such a fragment
# occurs in the primer AND the fragment's own Tm is above (primer Tm - 10 C).
# A fragment that hot has to be nearly as long as the primer itself, so the
# test is asking for a near-full-length perfect palindrome. Ordinary coding
# sequence never satisfies it.
#
# So the scorer has NO working hairpin check. Making one live means choosing a
# reachable threshold (a >= 6 nt self-complement above ~35 C, say), which is a
# new criterion to design and sweep, not a weight to set back to 6.
#
# NOTE: this paragraph used to name $W_DIMER alongside $W_HOMOLOGY. That
# stopped being true on 2026-11-14, when sub primerdimer was rewritten to the
# 3'-anchored duplex test against $DIMER_MAX_END3. $W_DIMER is live: on the
# 400-primer corpus the flagged-pair rate runs 17.0% at weight 0, 8.0% at 3,
# 3.5% at 10 and 2.5% at 20.


my %options;
GetOptions("i=s"=>\$options{i},
	   "r=s"=>\$options{r},
	   "o=s"=>\$options{o},
	   "h=s"=>\$options{h},
	   "oname=s"=>\$options{oname},
	   "workdir=s"=>\$options{workdir},
	   "baseurl=s"=>\$options{baseurl},
	   "cstandard:s"=>\$options{cstandard},
	   "cmaninfor:s"=>\$options{cmaninfor},
	   "cmaninrev:s"=>\$options{cmaninrev},
	   "cnick:s"=>\$options{cnick},
	   "crestr:s"=>\$options{crestr},
	   "cdb1:s"=>\$options{cdb1},
	   "cdb2:s"=>\$options{cdb2},
	   "cname:s"=>\$options{cname},
	   "tm:s"=>\$options{tm},
	   "poly:s"=>\$options{poly},
	   "sc:s"=>\$options{sc},
	   "pc:s"=>\$options{pc},
	   "co:s"=>\$options{co},
	   "verbose:s"=>\$options{verbose},
	   "fo:s"=>\$options{fo},
	   "dp:s"=>\$options{dp},
	   "dpcsv:s"=>\$options{dpcsv},
	   "igb:s"=>\$options{igb},
	   "dpgb:s"=>\$options{dpgb},
	   "ccustfull:s"=>\$options{ccustfull},
	   "ccust5:s"=>\$options{ccust5},
	   "ccust3:s"=>\$options{ccust3},
	   "autofrag:s"=>\$options{autofrag},
	   "fraglen:s"=>\$options{fraglen},
	   "tag=s@"=>\$options{tag});
	   	   


###################################################################	   
#if in debugging mode, give message


my $verbose;

if ($options{verbose} eq "y") {
	$verbose = "TRUE";
	print "\n## DEBUG MODE MESSAGES ##\n";
}




###################################################################
#Full-plasmid reconstruction input (USER cassette cloning only).
#
# -dp = destination plasmid sequence (the vector that the USER
#       cassette sits in, i.e. the sequence being cloned into).
#       When supplied, the assembled construct is spliced into this
#       plasmid in place of the cassette, and the complete plasmid
#       sequence (sans annotations) is written to the report. The
#       plasmid's own coordinates/orientation are preserved. Does
#       nothing for linear/circular (no-cassette) assembly.


###################################################################
#GenBank input handling (added 2026-10-29). If a GenBank file is
#supplied it supersedes the corresponding plain-text/FASTA input:
# -igb  = GenBank file of the fragments (overrides -i). Converted to
#         a temporary FASTA that feeds the normal pipeline; features
#         are stashed in %GBfrag (keyed by sequence) for transfer.
# -dpgb = GenBank file of the destination plasmid (overrides -dp).
#         Its sequence becomes -dp; its features go in %GBplasmid.

if (defined $options{igb} and $options{igb} ne "") {
	my @recs = &gb_parse_file($options{igb});
	die "ERROR: no sequence records found in GenBank file: $options{igb}\n" unless @recs;
	my ($tfh, $tfn) = tempfile("amuser_igbXXXXXX", SUFFIX => ".fa", TMPDIR => 1, UNLINK => 1);
	foreach my $r (@recs) {
		(my $nm = $r->{name}) =~ s/\s+/_/g;
		my $sq = $r->{seq};
		$sq =~ s/[^A-Za-z]//g;
		die "ERROR: GenBank record '$r->{name}' contains no DNA sequence.\n" if $sq eq "";
		die "ERROR: GenBank record '$r->{name}' contains non-DNA characters.\n"
			if ($sq =~ /[^ACGTURYKMBVDHNWSacgturykmbvdhnws]/);
		print $tfh ">$nm\n$sq\n";
		push @{ $GBfrag{ uc $sq } }, { name => $nm, features => $r->{features} };
	}
	close $tfh;
	$options{i} = $tfn;		#feed the normal FASTA pipeline from the GenBank-derived FASTA
}

if (defined $options{dpgb} and $options{dpgb} ne "") {
	my @recs = &gb_parse_file($options{dpgb});
	die "ERROR: no sequence records found in GenBank plasmid file: $options{dpgb}\n" unless @recs;
	my $r = $recs[0];
	warn "WARNING: -dpgb GenBank has multiple records; using the first ('$r->{name}').\n" if @recs > 1;
	$options{dp} = $r->{seq};
	%GBplasmid = ( name => $r->{name}, topology => $r->{topology}, features => $r->{features} );
}

my $destplasmid = "";
if (defined $options{dp}) {
	$destplasmid = uc $options{dp};
	$destplasmid =~ s/[^A-Za-z]//g;		#strip whitespace, numbers, FASTA headers etc.
}

#Per-cassette ADDITIONAL accepted full-cassette variants (uppercase,
#no U). The canonical $cassette[0] is ALWAYS matched automatically;
#only list here extra variants that are excised to give identical
#overhangs and scars (e.g. multi-restriction-site linker variants of
#a cassette). Keyed by cassette name ($cassette[4]).

my %cassettevariants = (
	"PacI/Nt.BbvCI" => ["GCTGAGGCTTAATTAAGGATCCTTAATTAAACCTCAGC"],
);

#Rows accumulated for the optional full-plasmid CSV (-dpcsv). Each row:
#[construct label, length (bp), status/orientation, full plasmid sequence].
my @fullplasmidrows;




###################################################################	   
#the filename inputs are handled 
	   

my $inputfasta;		#Single-fasta input filename
if ($options{i} eq "+") {
	$inputfasta = <STDIN>;
}
else {
	$inputfasta = $options{i};
}

$options{r} = "report"  if $options{r} eq "";   # default internal filename
$options{o} = "primers" if $options{o} eq "";   # default internal filename

my $outreport;		#Report output filename
my $outreportshort;		#Report output filename without path if given by wrapper
if ($options{r} eq "+") {
	$outreport = <STDIN>;
}
else {
	$outreport = $options{r};
	$outreport =~ m/.*\/(\w+\.html)/;
	#$outreportshort = $1;
	$outreportshort = "report.txt";
}

my $fastaout;		#Fasta output filename
my $fastaoutshort;		#Fasta output filename without path if given by wrapper
if ($options{o} eq "+") {
	$fastaout = <STDIN>;
}
else {
	$fastaout = $options{o};
	$fastaout =~ m/.*\/(\w+\.html)/;
	#$fastaoutshort = $1;
	$fastaoutshort = "primers.fsa";
}

my $htmlreport;		#html file of output
my $htmlreportshort;		#html file of output without path if given by wrapper
# When no -h filename is given, defer naming until the fragments are parsed so
# the report can be named after the construct: "<fragment>_YYYY-MM-DD_HH-MM-SS"
# for a single fragment, or "fusion_YYYY-MM-DD_HH-MM-SS" for a multi-fragment
# fusion (see the $auto_h block after input parsing). The .gb outputs inherit
# this stem.
#
# -oname <name> overrides that auto-naming with a user-supplied base name and
# NO date/time stamp: the run writes "<name>_report.html", "<name>_insert.gb"
# (or "<name>_plasmid.gb") and "<name>_primers.csv". The name is sanitised to
# filesystem-safe characters; a trailing ".html"/"_report" the user may have
# typed is stripped so the suffix is not doubled. An explicit -h still wins.
# Note: re-running with the same -oname overwrites the previous run's files.
if ($options{h} eq "" and defined $options{oname} and $options{oname} ne "") {
	my $on = $options{oname};
	$on =~ s/^\s+//; $on =~ s/\s+$//;
	$on =~ s/\.html?$//i;		#user may have typed an extension
	$on =~ s/_report$//i;		#...or the report suffix itself
	$on =~ s/\s+/_/g;
	$on =~ s/[^A-Za-z0-9._-]/_/g;
	$on =~ s/^[._-]+//;		#no leading dot/dash: keep it a plain visible file
	$on = substr($on, 0, 60) if length($on) > 60;
	die "ERROR: -oname contains no usable characters for a filename.\n" if $on eq "";
	$options{h} = $on . "_report.html";
}
my $auto_h = ($options{h} eq "") ? 1 : 0;
# Write all per-run outputs (HTML report, primer CSV, GenBank) into a single
# "Output files" folder located next to this script, regardless of the current
# working directory the script is run from. (Variable name kept for brevity.)
my $htmlreportdir = dirname(abs_path(__FILE__)) . "/Output files";
mkdir $htmlreportdir unless -d $htmlreportdir;
$htmlreport = $auto_h ? "" : "$htmlreportdir/" . basename($options{h});
$htmlreportshort = "report.html";


###################################################################
#The cassette input options are then handled. The array @cassette 
#will contain all information. $cassette[0] = forward strand, 
#$cassette[1] = reverse strand, $cassette[2] = forward strand 
#fusion tail, $cassette[3] = reverse strand fusion tail, 
#$cassette[4] = cassette name, $cassette[5] = pre-contruct cassette 
#sequence, $cassette[6] = post-contruct cassette sequence.

#If a cassette is selected from the standard list:


my @cassette;		#Array containing cassette data
my $cassetteprint;	#Variable which will contain cassette details for print
unless ($options{cstandard} eq "") {
	if ($options{cstandard} eq "PacI/Nt.BbvCI") {
		$cassette[0] = "GCTGAGGCTTAATTAAACCTCAGC";
		$cassette[1] = "CGACTCCGAATTAATTTGGAGTCG";
		$cassette[2] = "GGCTTAAU";
		$cassette[3] = "GGTTTAAU";
		$cassette[4] = "PacI/Nt.BbvCI";
		$cassette[5] = "GCTGA";
		$cassette[6] = "TCAGC";
	}
	elsif ($options{cstandard} eq "XbaI2/Nt.BbvCI") {
		$cassette[0] = "GCTGAGGGAAAGTCTAGAGGATCCTCTAGATGTCTCCTCAGC";
		$cassette[1] = "CGACTCCCTTTCAGATCTCCTAGGAGATCTACAGAGGAGTCG";
		$cassette[2] = "GGGAAAGU";
		$cassette[3] = "GGAGACAU";
		$cassette[4] = "XbaI(2)/Nt.BbvCI";
		$cassette[5] = "GCTGA";
		$cassette[6] = "TCAGC";
	}
	elsif ($options{cstandard} eq "PmeI/Nb.BbvCI") {
		$cassette[0] = "CCTCAGCCGTTTAAACAGCTGAGG";
		$cassette[1] = "GGAGTCGGCAAATTTGTCGACTCC";
		$cassette[2] = "GCCGTTU";
		$cassette[3] = "GCTGTTU";
		$cassette[4] = "PmeI/Nb.BbvCI";
		$cassette[5] = "CCTCA";
		$cassette[6] = "TGAGG";
	}
	elsif ($options{cstandard} eq "AsiSI/Nb.BsmI") {
		$cassette[0] = "GAATGCGTGCGATCGCGTGCATTC";
		$cassette[1] = "CTTACGCACGCTAGCGCACGTAAG";
		$cassette[2] = "CGTGCGAU";
		$cassette[3] = "CACGCGAU";
		$cassette[4] = "AsiSI/Nb.BSMI";
		$cassette[5] = "GAATG";
		$cassette[6] = "CATTC";
	}
	elsif ($options{cstandard} eq "AsiSI/Nb.BtsI") {
		$cassette[0] = "GCAGTGAGAGCGATCGCAGACACTGC";
		$cassette[1] = "CGTCACTCTCGCTAGCGTCTGTGACG";
		$cassette[2] = "AGAGCGAU";
		$cassette[3] = "TCTGCGAU";
		$cassette[4] = "AsiSI/Nb.BtsI";
		$cassette[5] = "GCAGT";
		$cassette[6] = "TGACG";
	}
	else {
		die "ERROR: unknown standard cassette input option: \"$options{cstandard}\"\n";
	}
	$cassetteprint = "standard cassette selected: $cassette[4]\n* cassette forward strand: $cassette[0]\n* cassette reverse strand: $cassette[1]";
}



###################################################################
#Determining whether input is single or multi sequence input


my $n = 0;

open IN, "<", "$options{i}" or die "ERROR: you have entered invalid filename(s): $!";
my $line = <IN>;
chomp $line;

while (defined $line) {
	if ($line =~ m/>/) {
		$n++;
	}
	$line = <IN>;
	chomp $line;
}
close IN;



###################################################################
#Setting output path for batch outputs


my $outputpath;
my $outputurl;
if ($options{workdir} ne "") {
	$outputpath = $options{workdir};
	# system("mkdir $outputpath");
}
elsif ($options{workdir} eq "") {
	$outputpath = ".";
}
if ($options{baseurl} ne "") {
	$outputurl = $options{baseurl};
}
elsif ($options{baseurl} eq "") {
	$outputurl = ".";
}

if ($verbose eq "TRUE") {
	print "output URL: $outputurl\n";
}


###################################################################
#If a cassette is designed from the list of nicking/restriction 
#sites: The input will be the names of the two selected sites, so 
#first a table with the correponding sequences is made:


my %restrictionsites = ("AsiSI"=>"GCGAT.CGC", 
						"PacI"=>"TTAAT.TAA",
						"PmeI"=>"GTTT.AAAC",
						"SwaI"=>"ATTT.AAAT",
						"XbaI"=>"T.CTAGA",
						"MssI"=>"GTTT.AAAC",
						"RgaI"=>"GCGAT.CGC",
						"SgfI"=>"GCGAT.CGC",
						"SmiI"=>"ATTT.AAAT",
						"SrfI"=>"GCCC.GGGC");

my %nickingsites = ("Nb.BbvCI"=>"CCTCA.GC",
		    		"Nb.BsmI"=>"GAATG.CT",
		    		"Nb.BsrDI"=>"GCAATG.GG",
		    		"Nb.BtsI"=>"GCAGTG.CT",
		    		"Nt.AlwI"=>"C.TCTCGATCC",
		    		"Nt.BbvCI"=>"GCTGA.GG",
		    		"Nt.BsmAI"=>"C.AGAGAC",
		    		"Nt.BspQI"=>".AGAAGAGC",
		    		"Nt.BstNBI"=>"T.TAGTCTGAT");



###################################################################
#Then the for/rev strand and fusion tails are deduced from the 
#input:

my $db1;
my $db2;
unless ($options{crestr} eq "" and $options{cnick} eq "") {
	$db1 = uc $options{cdb1};
	$db2 = uc $options{cdb2};
	my $cassette_ref = &cassettenickandrestr(\%options, \%restrictionsites, \%nickingsites, $db1, $db2);
	@cassette = @$cassette_ref;

	$cassetteprint = "custom cassette: $cassette[4]\ncassette forward strand: $cassette[0]\ncassette reverse strand: $cassette[1]";
}



###################################################################
#If a cassette is specified from manual input:


unless ($options{cmaninfor} eq "" and $options{cmaninrev} eq "") {
	my $cassette_ref = &cassettemanin(\%options, \@cassette);
	@cassette = @$cassette_ref;
	$cassetteprint = "custom cassette (manual input)\ncassette forward strand: $cassette[0]\ncassette reverse strand: $cassette[1]";
}


###################################################################
# Custom cassette from three plain fields (added 2026-10-29):
#   -ccustfull = the full cassette forward strand (located in -dp)
#   -ccust5    = 5' fusion tail (reverse-primer/fusion-tail orientation,
#                ending in U; a trailing T is accepted and converted)
#   -ccust3    = 3' fusion tail, same convention
# The nick-site flanks [5]/[6] are derived by anchoring each overhang's
# DNA form inside the full cassette, so the three fields fully define the
# cassette without the dot notation used by -cmaninfor/-cmaninrev.

if (defined $options{ccustfull} and $options{ccustfull} ne "") {
	my $full = uc $options{ccustfull}; $full =~ s/\s//g;
	my $o5   = uc(defined $options{ccust5} ? $options{ccust5} : ""); $o5 =~ s/\s//g;
	my $o3   = uc(defined $options{ccust3} ? $options{ccust3} : ""); $o3 =~ s/\s//g;

	die "ERROR: -ccustfull (full cassette sequence) must be DNA (A/C/G/T): \"$full\"\n"
		if ($full eq "" or $full =~ /[^ACGT]/);
	die "ERROR: both -ccust5 and -ccust3 (the two overhangs) are required with -ccustfull\n"
		if ($o5 eq "" or $o3 eq "");

	# Normalise each overhang to end in a single U (accept a trailing T).
	for my $ov (\$o5, \$o3) {
		$$ov =~ s/T$/U/ unless $$ov =~ /U$/;
		die "ERROR: overhang \"$$ov\" must be DNA ending in a single U (the uracil), e.g. GGGTTTAAU\n"
			unless $$ov =~ /^[ACGT]+U$/;
	}

	(my $o5dna = $o5) =~ s/U/T/;          # DNA form (forward strand) of the 5' overhang
	(my $o3dna = $o3) =~ s/U/T/;          # DNA form of the 3' fusion tail (reverse-primer orientation)
	my $o3fwd = &revcomp($o3dna);         # 3' overhang as it reads on the forward strand of the cassette

	my $i5 = index($full, $o5dna);
	die "ERROR: the 5' overhang ($o5dna) is not found in the full cassette ($full)\n" if $i5 < 0;
	my $i3 = rindex($full, $o3fwd);
	die "ERROR: the 3' overhang ($o3fwd, forward form of $o3) is not found in the full cassette ($full)\n" if $i3 < 0;

	my @c;
	$c[0] = $full;
	$c[1] = &revcomp($full);
	$c[2] = $o5;
	$c[3] = $o3;
	$c[4] = ($options{cname} ne "") ? $options{cname} : "custom cassette";
	$c[5] = substr($full, 0, $i5);                       # pre-construct flank (nick-site remnant)
	$c[6] = substr($full, $i3 + length($o3fwd));         # post-construct flank
	@cassette = @c;
	$cassetteprint = "custom cassette ($cassette[4])\n* cassette forward strand: $cassette[0]\n* cassette reverse strand: $cassette[1]";
}



###################################################################
# Default cassette: if no cassette was specified by any method and
# circular assembly is not requested, fall back to PacI/Nt.BbvCI.
# Must come before the $nocassette check so $nocassette stays FALSE.

if (!@cassette and $options{co} ne "y") {
	$cassette[0] = "GCTGAGGCTTAATTAAACCTCAGC";
	$cassette[1] = "CGACTCCGAATTAATTTGGAGTCG";
	$cassette[2] = "GGCTTAAU";
	$cassette[3] = "GGTTTAAU";
	$cassette[4] = "PacI/Nt.BbvCI";
	$cassette[5] = "GCTGA";
	$cassette[6] = "TCAGC";
	$cassetteprint = "standard cassette selected: $cassette[4]\n* cassette forward strand: $cassette[0]\n* cassette reverse strand: $cassette[1]";
}


###################################################################
#If no cassette is chosen (linear cloning).


my $nocassette = "FALSE";
if ($cassetteprint eq "") {
	$cassetteprint = "linear cloning (no cassette)";
	$nocassette = "TRUE";
}



###################################################################
#If Tm optimimization is selected, the desired Tm is stored. 
#$Tmtarget is either a numerical value or "eq".


# -poly selects the Tm model's polymerase/buffer preset (see %TM_POLY).
if (defined $options{poly} and $options{poly} ne "") {
	my $p = lc $options{poly};
	$p =~ s/[\s-]+/_/g;
	$p = $TM_POLY_ALIAS{$p} if exists $TM_POLY_ALIAS{$p};
	unless (exists $TM_POLY{$p}) {
		die "ERROR: unknown -poly '$options{poly}'. Known: "
		  . join(", ", sort keys %TM_POLY) . "\n";
	}
	$TM_POLY = $p;
}

my $Tmtarget = "";
if ($options{tm} ne "") {
	$Tmtarget = $options{tm};
	unless ($Tmtarget eq "eq") {
		if ($Tmtarget < 40 or $Tmtarget > 72) {
			die "ERROR: tm target must be between 40 and 72\n";
		}
	}	
}



###################################################################
#If circular assembly is selected, a flag is raised


my $circularoutput = "";
if ($options{co} eq "y") {
	$circularoutput = "TRUE";
	$cassetteprint = "circular assembly (no cassette)";
	if ($options{cstandard} ne "" or $options{cmaninfor} ne "" or $options{cmaninrev} ne "" or $options{cnick} ne "" or $options{crestr} ne "" or (defined $options{ccustfull} and $options{ccustfull} ne "")) {
		die "ERROR: a cassette cannot be added when circular assembly is selected\n";
	}
}





###################################################################
#Parse forced-overhang specification (-fo).
#Format: -fo <junc>:<A|B>[,<junc>:<A|B>...]
#Junction numbers are 1-based (junction 1 = between fragment 1 and 2).
#"A" forces the overhang to come from the last 24 nt of the upstream fragment.
#"B" forces the overhang to come from the first 24 nt of the downstream fragment.
#Cassette junctions are unaffected.


my %forced_overhang;
$options{fo} = "" unless defined $options{fo};
if ($options{fo} ne "") {
	foreach my $spec (split /,/, $options{fo}) {
		if ($spec =~ /^(\d+):(A|B)$/i) {
			my ($junc, $dir) = ($1 - 1, uc($2));
			if ($junc < 0) {
				die "ERROR: -fo junction numbers must be >= 1 (got \"$spec\")\n";
			}
			$forced_overhang{$junc} = $dir;
		} else {
			die "ERROR: invalid -fo argument \"$spec\". Expected format: <junction_number>:<A|B>, e.g. -fo 1:A,2:B\n";
		}
	}
}


###################################################################
#The input fasta file is opened in a filehandle and processed. The
#array @input[$i][1] will contain the names of the sequences, extracted
#from the first word in the fasta header. The array @brickseq will
#contain the corresponding sequence. Name and sequence will be
#saved so that position in the two arrays correspond, i.e. name [1]
#match sequence [1]. Also, if inserts are entered in the fasta
#input, they are saved in the array @inserts.


my ($input_ref, $batchinput, $inputcount_ref, $backbone) = &inputfile($inputfasta, $circularoutput);
my @input = @$input_ref;
my @inputcount = @$inputcount_ref;
my $numberofinputs = 1;
my $numberofbatches = 99;
my @inputs;
my @batchoutput;

# if($batchinput eq "TRUE") {
# 	print "backbone value $backbone\n";
# }
if ($batchinput eq "TRUE") {
	my ($numberofbatches_ref, $inputs_ref, $inputcount_ref) = &batchinputconversion(\@input);
	$numberofbatches = $numberofbatches_ref;
	@inputs = @$inputs_ref;
	@inputcount = @$inputcount_ref;
	$numberofinputs = scalar @inputs;
	# print "number of inputs $numberofinputs\n";
}
else {
	$inputs[0] = [@input];
}


###################################################################
# Auto-fragmentation (backend split of large sequences).
# A gene flagged via -autofrag is cut into ~$fraglen-bp sub-fragments at native
# USER overhang spans (see overhang_cuts) and fed to the normal fusion pipeline,
# which regenerates the original sequence and enforces overhang uniqueness
# (rmclashft). Done here, before primerextension, so the sub-fragments flow
# through the entire pipeline and @primerextension stays aligned.
# @fragmap[$x][old_fragment_index] = [new_start_index, sub_fragment_count] so
# per-fragment tags (-tag) and forced overhangs (-fo), which the caller keys by
# the ORIGINAL fragment/junction numbering, can be remapped onto the expanded list.
our @fragmap;
my $fraglen = (defined $options{fraglen} && $options{fraglen} =~ /^\d+$/ && $options{fraglen} > 0)
              ? $options{fraglen} + 0 : 1000;
my %autofrag_set;
if (defined $options{autofrag} && $options{autofrag} ne "") {
	for my $tok (split /,/, $options{autofrag}) {
		$autofrag_set{$tok - 1} = 1 if ($tok =~ /^\d+$/ && $tok >= 1);
	}
}
if (%autofrag_set) {
	for (my $x = 0; $x < scalar @inputs; $x++) {
		my @old = @{$inputs[$x]};
		my (@new, @map);
		for (my $f = 0; $f < scalar @old; $f++) {
			my $el    = $old[$f];
			my $seq   = $el->[2];
			my $start = scalar @new;
			if ($autofrag_set{$f} && length($seq) > $fraglen) {
				my @cuts  = &overhang_cuts($seq, $fraglen);
				my $total = scalar(@cuts) + 1;
				my $prev  = 0; my $part = 1;
				for my $c (@cuts, length $seq) {
					my @cel = @$el;                       # clone element (keeps [0] ordinal)
					$cel[1] = $el->[1] . "_p" . $part;    # gene -> gene_p1, gene_p2, ...
					$cel[2] = substr($seq, $prev, $c - $prev);
					push @new, [ @cel ];
					$prev = $c; $part++;
				}
				$map[$f] = [ $start, $total ];
			}
			else {
				push @new, $el;
				$map[$f] = [ $start, 1 ];
			}
		}
		@{$inputs[$x]} = @new;
		$fragmap[$x]   = \@map;
	}
	# Remap -fo junctions onto the expanded list (CLI = single construct $x=0).
	# Old junction j sits between old fragments j and j+1; the equivalent new
	# junction is after the last sub-fragment of old fragment j. New internal
	# split junctions are left unforced (auto-scored).
	if (%forced_overhang && defined $fragmap[0]) {
		my %remap;
		for my $oldj (keys %forced_overhang) {
			if (defined $fragmap[0][$oldj]) {
				my $newj = $fragmap[0][$oldj][0] + $fragmap[0][$oldj][1] - 1;
				$remap{$newj} = $forced_overhang{$oldj};
			}
		}
		%forced_overhang = %remap;
	}
}


# for (my $xx = 0; $xx < $numberofinputs; $xx++) {
# 	# for (my $y = 0; $y < scalar $n; $y++) {
# 		print "input: $xx\n";
# 		for (my $z = 0; $z < scalar @{$inputs[$xx]}; $z++) {
# 			print "seq$z\n$inputs[$xx][$z][2]\n";
# 		}
# 		my $z = scalar @{$inputs[$xx]};
# 		print "\n";
# 	# }
# }
# print "$n\n";



###################################################################
#Short fragments (> 80 bp) arising from linkers, mutations, or 
#insertions are split up and added to the fragment preceeding and
#following, thus extending the primers for these fragments. If the
#fragment is first or last in the total construct, the length 
#restrictions is halved, such that the length of a fragment in the
#beginning or end of a construct has to be < 40 or > 80. In 
#addition to appending the @inputs array, the array 
#@primerextension holds the length extension for each primer.


my ($inputs_ref, $primerextension_ref, $numberofinserts) = &primerextension(\@inputs, $circularoutput);
@inputs = @$inputs_ref;
my @primerextension = @$primerextension_ref;

# for (my $y = 0; $y <= scalar $n; $y++) {
# 	if ($primerextension[0][$y][0] eq "TRUE") {
# 		print "insert $y\n$primerextension[0][$y][1]\n";
# 	}
# }

###################################################################
#AMUSER is now executed once for each batch running


for (my $x = 0; $x < $numberofinputs; $x++) {
	my $n = scalar @{$inputs[$x]};


	###################################################################
	#If circular output is selected, the @input array is adjusted to
	#contain copies of gene 0 at the end and gene n in the beginning of
	#the array. If single input fragment is submitted, the fragment is 
	#copied only once.
	
	
	if ($circularoutput eq "TRUE") {
		my $inputs_ref = &circular(\@{$inputs[$x]}, $n);
		@{$inputs[$x]} = @$inputs_ref;
		$n = scalar @{$inputs[$x]};
	}
	# print "input: $x\n";
	# for (my $xx = 0; $xx < $numberofinputs; $xx++) {
	# 	# for (my $y = 0; $y < scalar $n; $y++) {
	# 	# 	print "$y\n";
	# 		for (my $z = 0; $z < scalar @{$inputs[$x]}; $z++) {
	# 			if ($primerextension[$x][$z][0] eq "TRUE") {
	# 				print "ext$z\n$primerextension[$x][$z][1]\n";
	# 			}
	# 			print "seq$z\n$inputs[$x][$z][2]\n";
	# 		}
	# 		my $z = scalar @{$inputs[$x]};
	# 		if ($primerextension[$x][$z][0] eq "TRUE") {
	# 			print "ext$z\n$primerextension[$x][$z][1]\n";
	# 		}
	# 		print "\n";
	# 	# }
	# }
	# print "$n\n";



	###################################################################
	#The @linker array is then expanded to contain all information 
	#about the linkers. $linker[$i][0] contains the linker. 
	#$linker[$i][1] contains the front half, $linker[$i][2] the length 
	#of that half, $linker[$i][3] contains the back half, and 
	#$linker[$i][4] the length of that half. 
	
	
	# if ($linkerflag == 1) {
	# 	%linker = &linkerexpand(\%linker);
	# }
	
	
	
	###################################################################
	#Two sub fragments from each of sequences; the first and last 40 
	#bases are extracted and saved in an array of two arrays - one for 
	#the first 40 and one for the last 40 for easy access purposes. If 
	#linkers are added they will be accounted for in the sub fragments
	

	my $brickseqfrags_ref = &brickseqfrags(\@{$inputs[$x]}, $n, \@{$primerextension[$x]}, $circularoutput);
	my @brickseqfrags = @$brickseqfrags_ref;

	# for (my $y = 0; $y < scalar @{$inputs[$x]}; $y++) {
	# 	print "brickseqfrags $y front\t$brickseqfrags[$y][0]\nbrickseqfrags $y end\t$brickseqfrags[$y][1]\n";
	# }
	
	
	###################################################################
	#The last seqfrag is linked together with the first from the next 
	#brick into an 80 base fragment. This is then saved in the array 
	#@joined
	
	
	my $joined_ref = &joined(\@brickseqfrags, $n);
	my @joined = @$joined_ref;
	# print "batch $x\n";
	# for (my $y = 0; $y < scalar @{$inputs[$x]}; $y++) {
	# 	print "$inputs[$x][$y][1] $y $inputs[$x][$y][2]\n\n";
	# 	print "joined $y: $joined[$y]\n\n";
	# }
	
	
	###################################################################
	#Two array tables are created. The array @A contains a list of all 
	#positions of A's in the sequence, and @T contains the positions of 
	#T's. The data in each of these two tables are saved as 
	#$A/T[$i][$Apos], where $i links it to the sequence it was found 
	#in, and [$A/Tpos] gives the position of the A/T.
	
	
	$n = scalar @{$inputs[$x]};
	my ($A_ref, $T_ref) = &ATpos(\@joined, $n, $batchinput, \%forced_overhang);
	my @A = @$A_ref;
	my @T = @$T_ref;

	
	
	###################################################################
	#Then, each possible A-T segment is scored. Since the relative 
	#penalty for an A-T segment that is too long/short, and an A-T 
	#segment not centered is the same, a simple scoring system is 
	#defined: The penalty is simply +1 for each amino a segment is 
	#longer than 9 or shorter than 7, and +1 for each amino acid the 
	#center of the segment deviates from the center of the fragment we 
	#called "$joined". The reason that the penalty is the same in each 
	#case, is that the consequence in each case is the same: an 
	#increase in length of the primer by 1 aa. The score and segment 
	#range is then saved in an array table: @scoredfts. This enables 
	#the program to quickly retrieve the next best scoring A-T segment, 
	#in case the primer from the first are no good.
	

	my $scoredfts_ref = &ATscore(\@joined, \@A, \@T, $n, \%forced_overhang);
	my @scoredfts = @$scoredfts_ref;
   
	my $scores_ref = &scores(\@scoredfts, $n);
	my @scores = @$scores_ref;

	# for (my $y = 0; $y < (scalar @{$inputs[$x]}) -1; $y++) {
	# 	for (my $z = 0; $z < scalar @{$scoredfts[$y]}; $z++) {
	# 		print "$y\t$z\t$scoredfts[$y][$z]\n";
	# 	}
	# }
	
	
	###################################################################
	#Arrays are now made containing all the allowed lengths of binding 
	#region, so that these can be varied to have the same Tm on the 
	#same sequence for PCR
	
	
	my ($Apos_ref, $Tpos_ref, $brfor_ref, $brrev_ref, $fts_ref, $brforpos_ref, $brrevpos_ref) = &bindingregion(\@scoredfts, \@joined, $n);
	
	my @Apos = @$Apos_ref;		#Array containing the best scoring non-clashing fusion tail A position in $joined for fusion tail ranked $ii
	my @Tpos = @$Tpos_ref;		#Array containing the best scoring non-clashing fusion tail T position in $joined for fusion tail ranked $ii
	my @brfor = @$brfor_ref;	#Array containing binding regions for forward primer on biobrick $i for fusion tail ranked $ii
	my @brrev = @$brrev_ref;	#Array containing binding regions for reverse primer on biobrick $i for fusion tail ranked $ii
	my @fts = @$fts_ref;		#Array containing actual best scoring non-clashing fusion tail sequence for fusion tail ranked $ii
	my @brforpos = @$brforpos_ref;
	my @brrevpos = @$brrevpos_ref;


	#print "$brrev[0][0][20]\n";

	# for (my $y = 0; $y < scalar @{$inputs[$x]}; $y++) {
	# 	for (my $z = 0; $z < $n-1; $z++) {
	# 		print "forward position for $z: $brforpos[$z][0] T position: $Tpos[$z][0]\n";
	# 		print "reverse position for $z: $brrevpos[$z][0] A position: $Apos[$z][0]\n";
	# 		# for (my $k = 0; $k < scalar @{$brfor[$y][$z]}; $k++) {
	# 		# 	print "$y\t$brfor[$y][$z][$k]\n";
	# 		# }
	# 	}
	# }

	
	
	###################################################################
	#Similarly this is done for the bases on the first and last 
	#sequence. Since we already know where these fragments will start 
	#(next to the cassette), we don't have to worry about fusion tails
	#and simply start the binding region from base 1 on sequence 1 and 
	#on base length-$i on sequence $n
	
	
	for (my $i = $BR_LEN_MIN; $i <= $BR_LEN_MAX; $i++) {
		$brfor[0][0][$i] = substr ($brickseqfrags[0][0], 0, $i);
		$brrev[$n-1][0][$i] = substr ($brickseqfrags[$n-1][1], 50 - $i, $i);
		$brrev[$n-1][0][$i] =~ tr/ATGCRYKMBVDHSW/TACGYRMKVBHDWS/;
	}
	#print "$brfor[0][0][20]\n";


	# for (my $y = 0; $y < $n; $y++) {
	# 	print "$joined[$y]\n";
	# 	for (my $z = 18; $z < 25; $z++) {
	# 		print "$y $z for $brfor[$y][0][$z]\n";
	# 		print "$y $z rev $brrev[$y][0][$z]\n";
	# 	}
	# }

	
	###################################################################
	#The following algorithm is designed to removed ft clashes. It does 
	#so by comparing all fts and removing the one with best scoring 
	#next alternative if there is a clash. In case the two alternatives
	#have the same score, the one with the most following alternatives 
	#will be removed.
	

	my $ft_ref = &rmclashft(\@fts, \@scores, \@brrev, \@brfor, \@Apos, \@Tpos, \@scoredfts, $n);
	my @ft = @$ft_ref;
	

	
	###################################################################
	#Now to replace thymine with uracil at the end of the ft and the 
	#complement ft, the former saved in the already existing array @ft, 
	#and the latter saved in the array @ftcomp
	
	
	my ($ft_ref, $ftcomp_ref) = &uracilft(\@ft, $n);
	@ft = @$ft_ref;
	my @ftcomp = @$ftcomp_ref;
	
	
		
	###################################################################
	#If a short (< 80) fragment is inserted between fragments, the
	#primers are extended to encompass the short fragment between
	#fusion tails and binding regions.
	

	my $primerextension_ref = &extensions(\@primerextension, \@joined, \@brforpos, \@brrevpos, \@ft, $n, \@Apos, \@Tpos, \@brrev, \@brfor);
	@primerextension = @$primerextension_ref;

	# for (my $y = 0; $y < scalar @inputs; $y++) {
	# 	for (my $z = 0; $z < $n; $z++) {
	# 		print "$z\tforward: $primerextension[$y][$z][4]\n";
	# 		print "$z\treverse: $primerextension[$y][$z][5]\n";
	# 	}
	# }


	###################################################################
	#Now to make sure that primers on the same DNA fragment has similar 
	#Tm, as well as falling into the optimal range of 55 to 72, Tm for 
	#all 7 binding regions (both forward and reverse), is calculated 
	#and compared to the appropriate counterpart. If Tm is out of the 
	#optimum, or not similar to counterpart binding region, the length 
	#of one or both (if needed) is adjusted within the specified 
	#lengths of 18 to 28.
	#First, the counting of A, T, G, and C in each of the 14 binding 
	#regions from each DNA fragment:
	

	
	my ($Afor_ref, $Tfor_ref, $Gfor_ref, $Cfor_ref, $Arev_ref, $Trev_ref, $Grev_ref, $Crev_ref) = &ATGCcount(\@brfor, \@brrev, $n);
	
	my @Afor = @$Afor_ref;	#Count of A's in forward primer
	my @Tfor = @$Tfor_ref;	#Count of T's in forward primer	
	my @Gfor = @$Gfor_ref;	#Count of G's in forward primer
	my @Cfor = @$Cfor_ref;	#Count of C's in forward primer
	my @Arev = @$Arev_ref;	#Count of A's in reverse primer
	my @Trev = @$Trev_ref;	#Count of T's in reverse primer
	my @Grev = @$Grev_ref;	#Count of G's in reverse primer
	my @Crev = @$Crev_ref;	#Count of C's in reverse primer
	
	
	
	###################################################################
	#Calling salt concentration and primer concentration parameters if
	#these have been selected, otherwise the variables are set to their
	#default values
	
	
	# Defaults come from the selected polymerase (%TM_POLY), because salt and
	# primer concentration are arguments to the Tm equation, not presentation.
	# -sc / -pc still override them.
	my $sc = $TM_POLY{$TM_POLY}[0] / 1000;          # mM -> M
	# -pc is the primer concentration IN THE REACTION TUBE, in uM, as pipetted.
	# $pc_bench keeps that value for the report; $pc is what the Tm equation
	# receives, which for the Phusion family is a quarter of it, exactly as
	# NEB's calculator does with the 500 nM its own interface asks for.
	my $pc_bench = $TM_POLY{$TM_POLY}[1] / 1000000000;   # nM -> M
	if (exists $options{sc} and $options{sc} ne "") {
		$sc = $options{sc} / 1000;
		if ($sc < 10/1000 or $sc > 250/1000) {
			die "ERROR: salt concentration must be between 10 and 250 mM\n";
		}
	}
	if (exists $options{pc} and $options{pc} ne "") {
		$pc_bench = $options{pc} / 1000000;
		if ($pc_bench < 0.05/1000000 or $pc_bench > 5/1000000) {
			die "ERROR: primer concentration must be between 0.05 and 5 uM\n";
		}
	}
	my $pc = $pc_bench;
	$pc /= 4 if $TM_POLY{$TM_POLY}[2] eq 'schildkraut';
	
	
	
	###################################################################
	# Finalise the auto report name now that fragments are parsed.
	# >1 input record -> "fusion"; otherwise the single fragment's name.
	if ($auto_h) {
		my $cname = ($n > 1) ? "fusion" : $inputs[0][0][1];
		$cname = "construct" unless defined $cname and $cname ne "";
		$cname =~ s/\s+/_/g;
		$cname =~ s/[^A-Za-z0-9._-]/_/g;
		$cname = substr($cname, 0, 30) if length($cname) > 30;
		my @tt = localtime(time);
		$options{h} = sprintf("%04d-%02d-%02d_%02d-%02d-%02d_%s_report.html",
			$tt[5] + 1900, $tt[4] + 1, $tt[3], $tt[2], $tt[1], $tt[0], $cname);
		$htmlreport = "$htmlreportdir/" . basename($options{h});
	}

	###################################################################
	# Per-fragment tags (from -tag "frag|end|name|seq"). %tag5/%tag3 are
	# keyed by 0-based fragment index (= assembly order). seq is sense DNA;
	# for 3' tags the reverse complement is used in the RV primer, while the
	# sense form is spliced into the construct. Tags are pure flap/insert
	# additions: primer design (binding, Tm) is unchanged.
	###################################################################
	my (%tag5, %tag3);
	if (ref $options{tag} eq 'ARRAY') {
		for my $spec (@{$options{tag}}) {
			# New 7-field form: frag|end|name|oligo|linkoligo|linker|seq  (seq last so it
			# may hold anything; linker is DNA, no pipe). linkoligo=1 includes the linker
			# in the oligo (so the USER overhang is derived from the linker end) instead
			# of putting it on the flanking gene primer.
			# Back-compatible 5-field (frag|end|name|oligo|seq) and 4-field forms.
			my @f = split /\|/, $spec, 7;
			my ($tfrag, $tend, $tname, $toligo, $tlinkoligo, $tlinker, $tseq);
			if    (@f == 7) { ($tfrag, $tend, $tname, $toligo, $tlinkoligo, $tlinker, $tseq) = @f; }
			elsif (@f == 5) { ($tfrag, $tend, $tname, $toligo, $tseq) = @f; $tlinkoligo = 0; $tlinker = ""; }
			else            { ($tfrag, $tend, $tname, $tseq) = split /\|/, $spec, 4; $toligo = 0; $tlinkoligo = 0; $tlinker = ""; }
			next unless defined $tseq;
			$tseq = uc $tseq; $tseq =~ s/\s//g;
			next if $tseq eq "";
			$tlinker = defined $tlinker ? uc $tlinker : "";
			$tlinker =~ s/[^ACGT]//g;   # linker is DNA only
			next unless defined $tfrag and $tfrag =~ /^\d+$/;
			my $of = $tfrag - 1;   # original fragment index (caller/GUI numbering)
			my $tidx = $of;
			# If this gene was auto-fragmented, a 5' tag belongs on its first
			# sub-fragment and a 3' tag on its last; internal split junctions stay
			# tag-free.
			if (defined $fragmap[$x] && defined $fragmap[$x][$of]) {
				$tidx = ($tend eq '3') ? ($fragmap[$x][$of][0] + $fragmap[$x][$of][1] - 1)
				                       : $fragmap[$x][$of][0];
			}
			next if $tidx < 0 or $tidx >= $n;
			$tname = "tag" unless defined $tname and $tname ne "";
			my $oligo = (defined $toligo and $toligo eq '1') ? 1 : 0;
			my $linkoligo = (defined $tlinkoligo and $tlinkoligo eq '1') ? 1 : 0;
			if    ($tend eq '5') { $tag5{$tidx} = { seq => $tseq, name => $tname, oligo => $oligo, linker => $tlinker, linkoligo => $linkoligo }; }
			elsif ($tend eq '3') { $tag3{$tidx} = { seq => $tseq, name => $tname, oligo => $oligo, linker => $tlinker, linkoligo => $linkoligo }; }
		}
	}

	###################################################################
	#####                        SINGLESEQ                		  #####
	###################################################################

	
	if ($n == 1) {
		# Tag flaps for the single fragment: 5' tag -> FW (sense); 3' tag -> RV
		# (reverse complement). Empty when no tag. Spliced after the bridge,
		# immediately before the binding region, in every primer rendering.
		# Optional linker between the tag and the gene: 5' -> tag.linker (linker
		# gene-facing); 3' -> linker.tag (linker gene-facing). Carried in the flap
		# unless included in the oligo (see below).
		my $link5_0 = defined $tag5{0} ? ($tag5{0}{linker} // "") : "";
		my $link3_0 = defined $tag3{0} ? ($tag3{0}{linker} // "") : "";
		my $fwtag0 = defined $tag5{0} ? ($tag5{0}{seq} . $link5_0) : "";
		my $rvtag0 = defined $tag3{0} ? &revcomp($link3_0 . $tag3{0}{seq}) : "";

		# Oligo mode (single fragment): both ends abut the cassette, so a 5' oligo
		# tag's gene-facing side is the tag's 3' end and its outer side is the
		# cassette 5' scar; a 3' oligo tag's gene-facing side is the tag's 5' end
		# and its outer side is the cassette 3' scar. When active, the gene primer
		# uses the oligo overhang (no flap) and the tag is emitted as an oligo.
		my ($fwlead0, $rvlead0, @oligo_recs);
		# Stashed overhang data for the oligo fusion diagram (set below when a tag
		# is emitted as an oligo; read at the graphic-overview push sites).
		my ($tag5o5, $tag5o3, $tag5br3, $tag3o5, $tag3o3, $tag3br5) = ("","","","","","");
		my ($tag5core, $tag3core) = ("","");
		$fwlead0 = "$cassette[2]$primerextension[$x][0][1]$fwtag0";
		$rvlead0 = "$cassette[3]$primerextension[$x][1][1]$rvtag0";
		if (defined $tag5{0} and $tag5{0}{oligo}) {
			my $inoligo = ($link5_0 ne "" && $tag5{0}{linkoligo} && &linker_in_oligo($tag5{0}{seq}, $link5_0, 5));
			warn "WARNING: the linker on 5' tag '$tag5{0}{name}' cannot form the oligo overhang; placed on the gene primer instead.\n" if ($link5_0 ne "" && $tag5{0}{linkoligo} && !$inoligo);
			my $core = $inoligo ? ($tag5{0}{seq} . $link5_0) : $tag5{0}{seq};
			if (&oligo_viable($core)) {
				my $r = &compute_oligo($core, &cassette_scar5($cassette[2]), 'tag', 8);
				($tag5o5, $tag5o3, $tag5br3, $tag5core) = ($r->{o5}, $r->{o3}, $r->{br3}, $r->{core});
				$fwlead0 = $r->{tail3} . $r->{br3};
				$fwlead0 .= $link5_0 unless $inoligo;   # linker on the gene primer when not in the oligo
				push @oligo_recs, { name=>$tag5{0}{name}.($inoligo?"+linker":""), frag=>$inputs[$x][0][1], end=>5,
				                    fw=>$r->{oligo_fw}, rv=>$r->{oligo_rv} };
			} else { $tag5{0}{oligo}=0; warn "WARNING: tag '$tag5{0}{name}'".($inoligo?"+linker":"")." cannot form a viable oligo; using a primer flap instead.\n"; }
		}
		if (defined $tag3{0} and $tag3{0}{oligo}) {
			my $inoligo = ($link3_0 ne "" && $tag3{0}{linkoligo} && &linker_in_oligo($tag3{0}{seq}, $link3_0, 3));
			warn "WARNING: the linker on 3' tag '$tag3{0}{name}' cannot form the oligo overhang; placed on the gene primer instead.\n" if ($link3_0 ne "" && $tag3{0}{linkoligo} && !$inoligo);
			my $core = $inoligo ? ($link3_0 . $tag3{0}{seq}) : $tag3{0}{seq};
			if (&oligo_viable($core)) {
				my $r = &compute_oligo($core, 'tag', &cassette_scar3($cassette[3]), 8);
				($tag3o5, $tag3o3, $tag3br5, $tag3core) = ($r->{o5}, $r->{o3}, $r->{br5}, $r->{core});
				$rvlead0 = $r->{tail5} . &revcomp($r->{br5});
				$rvlead0 .= &revcomp($link3_0) unless $inoligo;   # linker on the gene primer when not in the oligo
				push @oligo_recs, { name=>$tag3{0}{name}.($inoligo?"+linker":""), frag=>$inputs[$x][0][1], end=>3,
				                    fw=>$r->{oligo_fw}, rv=>$r->{oligo_rv} };
			} else { $tag3{0}{oligo}=0; warn "WARNING: tag '$tag3{0}{name}'".($inoligo?"+linker":"")." cannot form a viable oligo; using a primer flap instead.\n"; }
		}


		###################################################################
		#Then calculation of GC ratio for each element. The results are 
		#saved as the fourth dimension of @TmGCrev and @TmGCfor (under 
		#array element [1])
		
		
		my @TmGCrev;
		my @TmGCfor;
		my ($TmGCrev_ref, $TmGCfor_ref) = &GCratiosingle(\@Afor, \@Tfor, \@Gfor, \@Cfor, \@Arev, \@Trev, \@Grev, \@Crev, \@brrev, \@brfor);
		@TmGCrev = @$TmGCrev_ref;
		@TmGCfor = @$TmGCfor_ref;
		
		
		
		###################################################################
		#Then calculation of Tm using the nearest-neighbor method, The 
		#results are saved as the fourth dimension of @TmGCrev and @TmGCfor 
		#(under array element [0])
		
		
		my ($TmGCrev_ref, $TmGCfor_ref) = &Tm(\@brfor, \@brrev, \@TmGCfor, \@TmGCrev, $sc, $pc, $n);
		@TmGCrev = @$TmGCrev_ref;
		@TmGCfor = @$TmGCfor_ref;
		
		
		
		###################################################################
		#Then for the comparision and optimization of Tm of reverse and 
		#forward binding regions on the same DNA fragment
		#The optimal binding regions in terms of Tm. For a given DNA 
		#fragment, $i, the optimal choice for reverse and forward binding 
		#region (length) is saved in the array $optbr[$i][0][$ii] for 
		#reverse, and $optbr[$i][1][$ii] for forward (where $i is the gene 
		#entry, and $ii are the Tm matches).
		
		
		my $optbr_ref; 
		my $primerTmwarning_ref;
		my @optbr;
		my %primerTmwarning;
		
		if ($Tmtarget eq "") {
			($optbr_ref, $primerTmwarning_ref) = &optbr(\@TmGCrev, \@TmGCfor, $n);
			@optbr = @$optbr_ref;
			%primerTmwarning = %$primerTmwarning_ref;
		}
		
		elsif ($Tmtarget ne "") {
			($optbr_ref, $primerTmwarning_ref, $Tmtarget) = &optbrTm(\@TmGCrev, \@TmGCfor, $Tmtarget, $n);
			@optbr = @$optbr_ref;
			%primerTmwarning = %$primerTmwarning_ref;
		}
		
		
		
		###################################################################
		#Analyzing for risk of primer dimer formation. This is simply done 
		#by aligning the bases of each primer in a set from length 5 (the 
		#smallest dimer that can cause problems) to full alignment. If a 
		#primer dimer risk is detected, a flag is saved in the hash 
		#%primerdimer{$i}{$ii} where $i is the gene, and $ii is the number 
		#of the Tm match.
		

		my $primerdimer_ref = &primerdimer(\@brrev, \@brfor, \@optbr, \@TmGCfor, \@TmGCrev, $n, $sc, $pc);
		my %primerdimer = %$primerdimer_ref;
		
		
		
		###################################################################
		#Testing for intra-primer homology is done by aligning segments of 
		#4 bases with the rest of the primer sequence. 3 bases in each end 
		#of the 4 base segment are replace by "---" for the alignment, 
		#since "snap-backs" only occur if a small loop can be formed. If an 
		#intra-primer homology risk is detected, a flag is saved in the 
		#hash %homology{$i}{$ii} where $i is the gene, and $ii is the 
		#number of the Tm match.
		
		
		my ($homologyrev_ref, $homologyfor_ref) = &homology(\@brrev, \@brfor, \@optbr, \@TmGCfor, \@TmGCrev, $n, $sc, $pc);
		my %homologyrev = %$homologyrev_ref;
		my %homologyfor = %$homologyfor_ref;
		
		
		
		###################################################################
		#Testing for presence of GC clamp is done simply by checking for 
		#the presence of G's or C's at the 3' end of a primer 
					
		
		my ($GCclamprev_ref, $GCclampfor_ref) = &GCclamp(\@brrev, \@brfor, \@optbr, $n);
		my %GCclamprev = %$GCclamprev_ref;
		my %GCclampfor = %$GCclampfor_ref;
		
		
		
		###################################################################
		#Testing for risk of too strong binding, i.e. more than 3 G/C's in 
		#last 5 bases of 3' end
		
		
		my ($GCcountclamprev_ref, $GCcountclampfor_ref) = &GCcountclamp(\@brrev, \@brfor, \@optbr, $n);
		my %GCcountclamprev = %$GCcountclamprev_ref;
		my %GCcountclampfor = %$GCcountclampfor_ref;
		
		
		
		###################################################################
		#Testing for presence of poly N stretches. Poly is in this case 
		#defined as more than 4 dentical nucleotides in a row
		
		
		my ($polyNrev_ref, $polyNfor_ref) = &polyN(\@brrev, \@brfor, \@optbr, $n);
		my %polyNrev = %$polyNrev_ref;
		my %polyNfor = %$polyNfor_ref;
		
		
		
		###################################################################
		#Now that all putative primers are analysed for the above mentioned
		#parameters, it's time to select the most optimal primers. This is
		#done either by ranking by all parameters. If desired Tm is defined
		#by the user, this step is skipped, since there is no variability
		#in the primers
		
		my $ranking_ref;
		my %ranking;
		
		# Both paths rank now. Under -tm the candidates are the tolerance pool
		# optbrTm built, and the band-based Tm terms are skipped: the pool already
		# carries the user's Tm intent, and a +-$TM_TARGET_TOL pool is wide enough to
		# straddle the 55 C edge and drag a deliberately low target back into the band.
		$ranking_ref = &optbranalysis(\@optbr, \@TmGCrev, \@TmGCfor, \%primerdimer, \%homologyrev, \%homologyfor, \%GCclamprev, \%GCclampfor, \%GCcountclamprev, \%GCcountclampfor, \%polyNrev, \%polyNfor, $n, ($Tmtarget ne "") ? 1 : 0);
		%ranking = %$ranking_ref;
		
		
		###################################################################
		#####              PRINTING SINGLE SEQ RESULTS                #####
		###################################################################
		
		
		
		#For the printing, all reverse primers must be complemented
		
		
		
		
		
		#First a quick overview of the needed primers is printed to the 
		#screen output report
		
		
		push(@{$outreportarray[$x]}, "--------------------------------------\n");
		push(@{$outreportarray[$x]}, "AMUSER $version\n");
		push(@{$outreportarray[$x]}, "output generated: $mday-$mon-$year $hour:$min:$sec\n");
		push(@{$outreportarray[$x]}, "--------------------------------------");
		push(@{$outreportarray[$x]}, "\n\n");
		push(@{$outreportarray[$x]}, "input parameters:\n");
		push(@{$outreportarray[$x]}, "-----------------\n");
		push(@{$outreportarray[$x]}, "number of fragments: $inputcount[$x]\n");
		unless ($numberofinserts == 0) {
			push(@{$outreportarray[$x]}, "number of inserts: $numberofinserts\n");
		}
		{
			my $tempsc = $sc * 1000;
			my $temppc = $pc_bench * 1000000;
			push(@{$outreportarray[$x]}, "Tm model: NEB / $TM_POLY\n");
			push(@{$outreportarray[$x]}, "salt concentration selected: $tempsc mM\n");
			push(@{$outreportarray[$x]}, "primer concentration selected: $temppc uM\n");
		}
		if($Tmtarget ne "") {
			my $templine = "tm optimized to: $Tmtarget ";
			if ($options{tm} eq "eq") {
				$templine .= "C (all primer Tm approximated to each other)";
			}
			else {
				$templine .= "C (user defined)";
			}
			push(@{$outreportarray[$x]}, "$templine\n");
		}
		push(@{$outreportarray[$x]}, "$cassetteprint\n\n\n"); 
		push(@{$outreportarray[$x]}, "overview of the needed primers (5'-3'):\n");
		push(@{$outreportarray[$x]}, "---------------------------------------\n\n");
		push(@{$outreportarray[$x]}, "fragment\tstrand\ttm\tprimer\n");
		push(@{$outreportarray[$x]}, "$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t${fwlead0}$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
		$brrev[0][0][$optbr[0][0][$ranking{0}{0}]] = reverse $brrev[0][0][$optbr[0][0][$ranking{0}{0}]];
		$primerextension[$x][1][1] = reverse $primerextension[$x][1][1];
		$primerextension[$x][1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		push(@{$outreportarray[$x]}, "$inputs[$x][0][1]\treverse\t$TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t${rvlead0}$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");
		$brrev[0][0][$optbr[0][0][$ranking{0}{0}]] = reverse $brrev[0][0][$optbr[0][0][$ranking{0}{0}]];
		$primerextension[$x][1][1] = reverse $primerextension[$x][1][1];
		$primerextension[$x][1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		# Oligo strands (ordered as annealed FW+RV); listed here so they also flow into the primer CSV.
		for my $orec (@oligo_recs) {
			push(@{$outreportarray[$x]}, "$orec->{name} oligo\tforward\t-\t$orec->{fw}\n");
			push(@{$outreportarray[$x]}, "$orec->{name} oligo\treverse\t-\t$orec->{rv}\n");
		}
		if (@oligo_recs) {
			push(@{$outreportarray[$x]}, "\nTags above marked \"oligo\" are ordered as the annealed FW+RV strands listed and inserted into the USER fusion; the flanking gene primer already carries the matching USER overhang (no tag flap). The graphic fusion diagram below shows each oligo as a staggered sticky-end duplex (its two USER overhangs; core elided).\n");
		}


		#If there is any warnings in the primer details section, the user is given a heads up:
		
		
		my $warningflag = 0;
		
		unless (55 < $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0]) {
			$warningflag = 1;
		}
		
		unless (39.9 < $TmGCfor[0][0][$optbr[0][0][$ranking{0}{0}]][1] and 60.1 > $TmGCfor[0][0][$optbr[0][0][$ranking{0}{0}]][1]) {
			$warningflag = 1;
		}
		
		unless (exists $GCclampfor{0}{0}) {
			$warningflag = 1;
		}
		
		if (exists $GCcountclampfor{0}{0}) {
			$warningflag = 1;
		}
		
		if (exists $primerdimer{0}{0}) {
			$warningflag = 1;
		}
		
		if (exists $homologyfor{0}{0}) {
			$warningflag = 1;
		}	
		
		if (exists $polyNfor{0}{0}) {
			$warningflag = 1;
		} 
		
		
		$brrev[0][0][$optbr[0][0][$ranking{0}{0}]] = reverse $brrev[0][0][$optbr[0][0][$ranking{0}{0}]];
		unless (55 < $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][0]) {
			$warningflag = 1;
		}
		
		unless (39.9 < $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][1] and 60.1 > $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][1]) {
			$warningflag = 1;
		}
		
		unless (exists $GCclamprev{$n-1}{0}) {
			$warningflag = 1;
		}
		
		if (exists $GCcountclamprev{$n-1}{0}) {
			$warningflag = 1;
		}
		
		if (exists $primerdimer{$n-1}{0}) {
			$warningflag = 1;
		}
		
		if (exists $homologyrev{$n-1}{0}) {
			$warningflag = 1;
		}
		
		if (exists $polyNfor{$n-1}{0}) {
			$warningflag = 1;
		} 
		
		if (exists $primerTmwarning{$n-1}) {
			$warningflag = 1;
		}
		
		if ($warningflag == 1) {
			push(@{$outreportarray[$x]},  "attention: one or more of the designed primers may have undesirable properties\n");
			push(@{$outreportarray[$x]},  "please see \"primer details\" section\n\n\n");
		}
		
		
		
		
		#Then a graphic overview of the non-fusion tail, the fusion tail, and the binding region 
		#is printed to the screenout file
		
		
		
		my $cassetteftnoUpre = $cassette[2];
		$cassetteftnoUpre =~ tr/U/T/;
		my $cassetteftnoUpost = $cassette[3];
		$cassetteftnoUpost =~ tr/U/T/;
		$cassetteftnoUpost = reverse $cassetteftnoUpost;
		$cassetteftnoUpost =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		
		push(@{$outreportarray[$x]}, "graphic overview of DNA fragments and primers:\n");
		push(@{$outreportarray[$x]}, "----------------------------------------------\n\n");
		
		if ($nocassette eq "TRUE") {
			push(@{$outreportarray[$x]}, "forward primer for $inputs[$x][0][1]:\n\n");
			push(@{$outreportarray[$x]}, "5'-$primerextension[$x][0][1]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]-3'\n");		
		}
		else {
			if (defined $tag5{0} && $tag5{0}{oligo}) {
				# N-terminal oligo tag at the 5' cassette junction: cassette - tag - gene.
				# Outer (5') overhang = cassette scar; gene-facing (3') overhang = tag-derived.
				my $cctx = (length($cassette[5]) > 30) ? substr($cassette[5], -30) : $cassette[5];
				my $ltr  = (length($cassette[5]) > 30) ? 1 : 0;
				my $gseq = defined $brickseqfrags[0][0] ? $brickseqfrags[0][0] : "";
				my $gctx = (length($gseq) > 30) ? substr($gseq, 0, 30) : $gseq;
				my $rtr  = (length($gseq) > 30) ? 1 : 0;
				push(@{$outreportarray[$x]}, &oligodiag_line(
					5, $cassette[4], $inputs[$x][0][1], $tag5{0}{name}, $cassette[4], $inputs[$x][0][1],
					$ltr, $cctx, "", $tag5o5, "cassette", $tag5o3, "tag", $tag5br3, $gctx, $rtr,
					"cassette", $inputs[$x][0][1], $tag5core));
			} else {
				push(@{$outreportarray[$x]}, "fusion region and forward primer for joining of $cassette[4] and $inputs[$x][0][1]:\n\n");
			}
			my $space = "";
			for (my $ii = 0; $ii < (length $cassette[5]); $ii++) {
				$space .= " ";
			}
			push(@{$outreportarray[$x]}, "$space");
			push(@{$outreportarray[$x]}, "5'-$cassette[2]\n");
			my $space = "";
			for (my $ii = 0; $ii < 3+(length $cassette[2])+(length $cassette[5]); $ii++) {
				$space .= " ";
			}	
			push(@{$outreportarray[$x]}, "$space$primerextension[$x][0][1]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]-3'\n");
		}
		push(@{$outreportarray[$x]}, "5'-$cassette[5]$cassetteftnoUpre$primerextension[$x][0][1]$brickseqfrags[0][0]\[...]-3'\n\n");
		push(@{$outreportarray[$x]}, "OLIGOSKIP_END\n") if (defined $tag5{0} && $tag5{0}{oligo});
		
		
		if (defined $tag3{0} && $tag3{0}{oligo}) {
			# C-terminal oligo tag at the 3' cassette junction: gene - tag - cassette.
			# Gene-facing (5') overhang = tag-derived; outer (3') overhang = cassette scar.
			my $gseq = defined $brickseqfrags[$n-1][1] ? $brickseqfrags[$n-1][1] : "";
			my $gctx = (length($gseq) > 30) ? substr($gseq, -30) : $gseq;
			my $ltr  = (length($gseq) > 30) ? 1 : 0;
			my $cctx = (length($cassette[6]) > 30) ? substr($cassette[6], 0, 30) : $cassette[6];
			my $rtr  = (length($cassette[6]) > 30) ? 1 : 0;
			push(@{$outreportarray[$x]}, &oligodiag_line(
				3, $inputs[$x][$n-1][1], $cassette[4], $tag3{0}{name}, $cassette[4], $inputs[$x][$n-1][1],
				$ltr, $gctx, $tag3br5, $tag3o5, "tag", $tag3o3, "cassette", "", $cctx, $rtr,
				$inputs[$x][$n-1][1], "cassette", $tag3core));
		} else {
			push(@{$outreportarray[$x]},  "fusion region and reverse primer for joining of $inputs[$x][$n-1][1] and $cassette[4]:\n\n");
		}
		push(@{$outreportarray[$x]},  "5'-[...]$brickseqfrags[$n-1][1]$primerextension[$x][1][1]$cassetteftnoUpost$cassette[6]-3'\n");
		
		my $space = "";
		for (my $ii = 0; $ii < (length $brickseqfrags[$n-1][1])+5-(length $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]); $ii++) {
			$space .= " ";
		}
		
		$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
		$primerextension[$x][1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		my $templine = $space;
		if ($nocassette eq "TRUE") {
			$templine .= "3'-$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]$primerextension[$x][1][1]-5'";
			push(@{$outreportarray[$x]},  $templine);
		}
		else {
			$templine .= "3'-$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]$primerextension[$x][1][1]";
			push(@{$outreportarray[$x]},  "$templine\n");
			my $space = "";
			for (my $ii = 0; $ii < (length $brickseqfrags[$n-1][1])+8 + length $primerextension[$x][1][1]; $ii++) {
				$space .= " ";
			}
			$cassette[3] = reverse $cassette[3];
			push(@{$outreportarray[$x]},  "$space$cassette[3]-5'\n");
			$cassette[3] = reverse $cassette[3];
		}
		$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
		$primerextension[$x][1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		push(@{$outreportarray[$x]},  "\n\n\n");
		push(@{$outreportarray[$x]}, "OLIGOSKIP_END\n") if (defined $tag3{0} && $tag3{0}{oligo});
		
		
		
		#Lastly, the final construct is printed to the screenout for good meassures
		
		
		
		my $final = "";
		# 5' tag (sense) before the fragment and 3' tag (sense) after it, each as
		# its own alternating case-run so each is a separate coloured construct
		# segment (matched by a tag entry in the name row below). Case is display-
		# only; the .gb annotates tags via the segment map.
		my $tag5seq = defined $tag5{0} ? $tag5{0}{seq} : "";
		my $tag3seq = defined $tag3{0} ? $tag3{0}{seq} : "";
		# linker is a separate construct segment, gene-facing of its tag.
		my $flink5 = $link5_0;
		my $flink3 = $link3_0;
		my $cf = ($cassette[5] ne "") ? 0 : 1;
		if ($cassette[5] ne "") {
			$final .= uc $cassette[5];
			$final .= uc $cassetteftnoUpre;
		}
		if ($primerextension[$x][0][0] eq "TRUE") {
			if ($cf == 0) { $final .= lc $primerextension[$x][0][1]; $cf = 1; } else { $final .= uc $primerextension[$x][0][1]; $cf = 0; }
		}
		if ($tag5seq ne "") {
			if ($cf == 0) { $final .= lc $tag5seq; $cf = 1; } else { $final .= uc $tag5seq; $cf = 0; }
		}
		if ($flink5 ne "") {
			if ($cf == 0) { $final .= lc $flink5; $cf = 1; } else { $final .= uc $flink5; $cf = 0; }
		}
		if ($cf == 0) { $final .= lc $inputs[$x][0][2]; $cf = 1; } else { $final .= uc $inputs[$x][0][2]; $cf = 0; }
		if ($flink3 ne "") {
			if ($cf == 0) { $final .= lc $flink3; $cf = 1; } else { $final .= uc $flink3; $cf = 0; }
		}
		if ($tag3seq ne "") {
			if ($cf == 0) { $final .= lc $tag3seq; $cf = 1; } else { $final .= uc $tag3seq; $cf = 0; }
		}
		if ($primerextension[$x][1][0] eq "TRUE") {
			if ($cf == 0) { $final .= lc $primerextension[$x][1][1]; $cf = 1; } else { $final .= uc $primerextension[$x][1][1]; $cf = 0; }
		}
		if ($cassette[5] ne "") {
			if ($cf == 0) { $final .= lc $cassetteftnoUpost; $final .= lc $cassette[6]; $cf = 1; } else { $final .= uc $cassetteftnoUpost; $final .= uc $cassette[6]; $cf = 0; }
		}

		my $templength = length($final);
		if ($circularoutput eq "TRUE") {
			push(@{$outreportarray[$x]},  "details of final construct after cloning (circular, length: $templength bases):\n");
			push(@{$outreportarray[$x]},  "------------------------------------------------------------------------\n\n");
		}
		else {
			push(@{$outreportarray[$x]},  "details of final construct after cloning (length: $templength bases):\n");
			push(@{$outreportarray[$x]},  "--------------------------------------------------------------\n\n");
		}

		my $string = "";
		my $insertcounter = 1;
		if ($cassette[4] ne "") {
			$string .= "$cassette[4],,";
		}
		if ($primerextension[$x][0][0] eq "TRUE") {
			$string .= "insert$insertcounter,,";
			$insertcounter++;
		}
		$string .= "$tag5{0}{name},," if defined $tag5{0};
		$string .= "linker,," if $flink5 ne "";
		$string .= "$inputs[$x][0][1],,";
		$string .= "linker,," if $flink3 ne "";
		$string .= "$tag3{0}{name},," if defined $tag3{0};
		if ($primerextension[$x][1][0] eq "TRUE") {
			$string .= "insert$insertcounter,,";
			$insertcounter++;
		}
		if ($cassette[4] ne "") {
			$string .= "$cassette[4]";
		}
		else {
			chop $string;
			chop $string;
		}
		$string .= "\n\n\n";
		push(@{$outreportarray[$x]}, "$string");
		
		# my $string = "";
		# if ($cassette[4] ne "") {
		# 	$string .= "$cassette[4],,";
		# }
		# $string .= "$inputs[$x][0][1]";
		# if ($cassette[4] ne "") {
		# 	$string .= ",,$cassette[4]";
		# }
		
		# push(@{$outreportarray[$x]}, "$string\n\n");
		
		for (my $pos = 0; $pos < length($final); $pos = $pos + 60) {
			my $string = $pos+1;
			for (my $i=length $pos; $i < 10; $i++) {
				$string .= " ";
			}
			$string .= substr($final, $pos, 60);
			push(@{$outreportarray[$x]}, "$string\n");
		}
			push(@{$outreportarray[$x]}, "");   # terminate the htmlprint construct-sequence walker

		&fullplasmidreportblock($outreportarray[$x], $destplasmid, \@cassette, $final, $inputs[$x][0][1]);

		#Detailed descriptions of the primers are then printed to the screenout file
		
		
		push(@{$outreportarray[$x]},  "primer details:\n");
		push(@{$outreportarray[$x]},  "---------------\n\n");
		push(@{$outreportarray[$x]},  "primer details for $inputs[$x][0][1]\n\n");
		push(@{$outreportarray[$x]},  "forward primer: ${fwlead0}$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
		
		if ($TM_OPT_MIN <= $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0]) {
			push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
		}
		
		if (39.9 < $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][1] and 60.1 > $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][1]) {
			push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
		}
		
		if (exists $GCclampfor{0}{0}) {
			push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
		}
		
		if (exists $GCcountclampfor{0}{0}) {
			push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
		}
		
		if (exists $primerdimer{0}{0}) {
			push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
		}
		
		if (exists $homologyfor{0}{0}) {
			push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
		}	
		
		if (exists $polyNfor{0}{0}) {
			push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
		} 
		push(@{$outreportarray[$x]},  "\n");
		
		
		
		
		
		
		
		#$brrev[$n-1][0][$optbr[$n-1][0][0]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][0]];
		$primerextension[$x][1][1] = reverse $primerextension[$x][1][1];
		$primerextension[$x][1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		push(@{$outreportarray[$x]},  "reverse primer: ${rvlead0}$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
		$primerextension[$x][1][1] = reverse $primerextension[$x][1][1];
		$primerextension[$x][1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		if ($TM_OPT_MIN <= $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0]) {
			push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
		}
		
		if (39.9 < $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1] and 60.1 > $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1]) {
			push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
		}
		
		if (exists $GCclampfor{$n-1}{0}) {
			push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
		}
		
		if (exists $GCcountclampfor{$n-1}{0}) {
			push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
		}
		
		if (exists $primerdimer{$n-1}{0}) {
			push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
		}
		
		if (exists $homologyfor{$n-1}{0}) {
			push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
		}	
		
		if (exists $polyNfor{$n-1}{0}) {
			push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
		}
		else {
			push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
		}
		
		push(@{$outreportarray[$x]},  "\n");
			
		if (exists $primerTmwarning{$n-1}) {
			push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...NO\n");
		}
		else {
			push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...YES\n");
		}
		push(@{$outreportarray[$x]},  "\n\n");
		&genbank_outputs($x, \@inputs, \@primerextension, \@cassette, $final, $destplasmid, $outreportarray[$x], $htmlreportdir, basename($options{h}), \%tag5, \%tag3)
			unless ($circularoutput eq "TRUE");

		#close OUT;
		
		
		#printing output to fasta array
		
		
		@fastaoutarray;
		push(@{$fastaoutarray[$x]},  ">FW_$inputs[$x][0][1]_[Tm:$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0]]\n${fwlead0}$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
		$primerextension[$x][1][1] = reverse $primerextension[$x][1][1];
		$primerextension[$x][1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		push(@{$fastaoutarray[$x]},  ">RV_$inputs[$x][$n-1][1]_[Tm:$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0]]\n${rvlead0}$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
		$primerextension[$x][1][1] = reverse $primerextension[$x][1][1];
		$primerextension[$x][1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		#close OUT;
	
	
	}
			
	
	
	###################################################################
	#####                MULTIPLE SEQUENCE INPUT                  #####
	###################################################################
	
	
	if ($n > 1) {

		
	
		###################################################################
		#Then calculation of GC ratio for each element. The results are 
		#saved as the fourth dimension of @TmGCrev and @TmGCfor (under 
		#array element [1])
		

		my @TmGCrev;
		my @TmGCfor;
		my ($TmGCrev_ref, $TmGCfor_ref) = &GCratiomulti(\@Afor, \@Tfor, \@Gfor, \@Cfor, \@Arev, \@Trev, \@Grev, \@Crev, \@brrev, \@brfor, $n);
		@TmGCrev = @$TmGCrev_ref;
		@TmGCfor = @$TmGCfor_ref;

		
		
		###################################################################
		#Since Tm is calculated for the binding region sequence + the part
		#of the fusion tail which also binds to the given sequence, two
		#arrays - @brforft and @brrevft - of sequences for Tm calculation 
		#are generated. They are simply copies of @brfor and @brrev with 
		#the addition of the ft part. For example, $brrevft[1][0][18] 
		#contains the reverse binding region of length 18+ft part for 
		#sequence 1.
		
		
		my @brrevft;
		my @brforft;
		for (my $i = 0; $i < $n; $i++) {
			for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
				$brforft[$i][0][$iii] = $brfor[$i][0][$iii];
			}
			for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
				$brrevft[$i][0][$iii] = $brrev[$i][0][$iii];
			}
		}
		
		my @revcompft;
		my @forft;
		for (my $i = 0; $i < $n-1; $i++) {
			$revcompft[$i] = $ft[$i];
			$revcompft[$i] =~ s/U/T/;
			$revcompft[$i] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
			for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
				$brrevft[$i][0][$iii] .= $revcompft[$i];
				# print "$i $iii rev $brrev[$i][0][$iii]\n";
			}
			$forft[$i] = $ft[$i];
			$forft[$i] =~ s/U/T/;
			for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
				my $tempbr = $forft[$i];
				$tempbr .= $brforft[$i+1][0][$iii];
				$brforft[$i+1][0][$iii] = $tempbr;
				# print "$i $iii for $brfor[$i+1][0][$iii]\n";
			}
		}
		
		# @brrevft;
		# @brforft;
		# for (my $i = 0; $i < $n; $i++) {
		# 	for (my $iii = 18+$linker{$i}{2}; $iii < 25+$linker{$i}{2}; $iii++) {
		# 		$brforft[$i][0][$iii] = $brfor[$i][0][$iii];
		# 	}
		# 	for (my $iii = 18+$linker{$i+1}{2}; $iii < 25+$linker{$i+1}{2}; $iii++) {
		# 		$brrevft[$i][0][$iii] = $brrev[$i][0][$iii];
		# 	}
		# }
		
		# @revcompft;
		# @forft;
		# for (my $i = 0; $i < $n-1; $i++) {
		# 	$revcompft[$i] = $ft[$i];
		# 	$revcompft[$i] =~ s/U/T/;
		# 	$revcompft[$i] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		# 	for (my $iii = 18+$linker{$i+1}{2}; $iii < 25+$linker{$i+1}{2}; $iii++) {
		# 		$brrevft[$i][0][$iii] .= $revcompft[$i];
		# 	}
		# 	$forft[$i] = $ft[$i];
		# 	$forft[$i] =~ s/U/T/;
		# 	for (my $iii = 18+$linker{$i+1}{2}; $iii < 25+$linker{$i+1}{2}; $iii++) {
		# 		my $tempbr = $forft[$i];
		# 		$tempbr .= $brforft[$i+1][0][$iii];
		# 		$brforft[$i+1][0][$iii] = $tempbr;
		# 	}
		# }

		
		
		###################################################################
		#Then calculation of Tm using the nearest-neighbor method, The 
		#results are saved in the fourth dimension of @TmGCrev and @TmGCfor 
		#(under array element [0])
		

		my ($TmGCrev_ref, $TmGCfor_ref) = &Tm(\@brfor, \@brrev, \@TmGCfor, \@TmGCrev, $sc, $pc, $n);
		@TmGCrev = @$TmGCrev_ref;
		@TmGCfor = @$TmGCfor_ref;

		
		
		###################################################################
		#Then for the comparision and optimization of Tm of reverse and 
		#forward binding regions on the same DNA fragment
		#The optimal binding regions in terms of Tm. For a given DNA 
		#fragment, $i, the optimal choice for reverse and forward binding 
		#region (length) is saved in the array $optbr[$i][0][$ii] for 
		#reverse, and $optbr[$i][1][$ii] for forward (where $i is the gene 
		#entry, and $ii are the Tm matches).
		
		
		my $optbr_ref; 
		my $primerTmwarning_ref;
		my @optbr;
		my %primerTmwarning;
		
		if ($Tmtarget eq "") {
			($optbr_ref, $primerTmwarning_ref) = &optbr(\@TmGCrev, \@TmGCfor, $n);
			@optbr = @$optbr_ref;
			%primerTmwarning = %$primerTmwarning_ref;
		}
		
		elsif ($Tmtarget ne "") {
			($optbr_ref, $primerTmwarning_ref, $Tmtarget) = &optbrTm(\@TmGCrev, \@TmGCfor, $Tmtarget, $n);
			@optbr = @$optbr_ref;
			%primerTmwarning = %$primerTmwarning_ref;
		}
		
		
		
		###################################################################
		#Analyzing for risk of primer dimer formation. This is simply done 
		#by aligning the bases of each primer in a set from length 5 (the 
		#smallest dimer that can cause problems) to full alignment. If a 
		#primer dimer risk is detected, a flag is saved in the hash 
		#%primerdimer{$i}{$ii} where $i is the gene, and $ii is the number 
		#of the Tm match.
		

		my $primerdimer_ref = &primerdimer(\@brrev, \@brfor, \@optbr, \@TmGCfor, \@TmGCrev, $n, $sc, $pc);
		my %primerdimer = %$primerdimer_ref;
		
		
		
		###################################################################
		#Testing for intra-primer homology is done by aligning segments of 
		#4 bases with the rest of the primer sequence. 3 bases in each end 
		#of the 4 base segment are replace by "---" for the alignment, 
		#since "snap-backs" only occur if a small loop can be formed. If an 
		#intra-primer homology risk is detected, a flag is saved in the 
		#hash %homology{$i}{$ii} where $i is the gene, and $ii is the 
		#number of the Tm match.
		
		#print "hey\n";
		my ($homologyrev_ref, $homologyfor_ref) = &homology(\@brrev, \@brfor, \@optbr, \@TmGCfor, \@TmGCrev, $n, $sc, $pc);
		my %homologyrev = %$homologyrev_ref;
		my %homologyfor = %$homologyfor_ref;
		#print "hey\n";
		
		
		###################################################################
		#Testing for presence of GC clamp is done simply by checking for 
		#the presence of G's or C's at the 3' end of a primer 
					
		
		my ($GCclamprev_ref, $GCclampfor_ref) = &GCclamp(\@brrev, \@brfor, \@optbr, $n);
		my %GCclamprev = %$GCclamprev_ref;
		my %GCclampfor = %$GCclampfor_ref;
		
		
		
		###################################################################
		#Testing for risk of too strong binding, i.e. more than 3 G/C's in 
		#last 5 bases of 3' end
		
		
		my ($GCcountclamprev_ref, $GCcountclampfor_ref) = &GCcountclamp(\@brrev, \@brfor, \@optbr, $n);
		my %GCcountclamprev = %$GCcountclamprev_ref;
		my %GCcountclampfor = %$GCcountclampfor_ref;
		
		
		
		###################################################################
		#Testing for presence of poly N stretches. Poly is in this case 
		#defined as more than 4 dentical nucleotides in a row
		
		
		my ($polyNrev_ref, $polyNfor_ref) = &polyN(\@brrev, \@brfor, \@optbr, $n);
		my %polyNrev = %$polyNrev_ref;
		my %polyNfor = %$polyNfor_ref;
		
		
		
		###################################################################
		#Now that all putative primers are analysed for the above mentioned
		#parameters, it's time to select the most optimal primers. This is
		#done either by ranking by all parameters. If desired Tm is defined
		#by the user, this step is skipped, since there is no variability
		#in the primers
		
		my $ranking_ref;
		my %ranking;
		
		# Both paths rank now. Under -tm the candidates are the tolerance pool
		# optbrTm built, and the band-based Tm terms are skipped: the pool already
		# carries the user's Tm intent, and a +-$TM_TARGET_TOL pool is wide enough to
		# straddle the 55 C edge and drag a deliberately low target back into the band.
		$ranking_ref = &optbranalysis(\@optbr, \@TmGCrev, \@TmGCfor, \%primerdimer, \%homologyrev, \%homologyfor, \%GCclamprev, \%GCclampfor, \%GCcountclamprev, \%GCcountclampfor, \%polyNrev, \%polyNfor, $n, ($Tmtarget ne "") ? 1 : 0);
		%ranking = %$ranking_ref;
		
		
		
		###################################################################
		#####               PRINTING MULTISEQ RESULTS                 #####
		###################################################################
		
		


####CIRCULAR OUPUT		



		if ($circularoutput eq "TRUE") {
			
			
			
			#First a quick overview of the needed primers is printed to the 
			#screen output report
			
			
			push(@{$outreportarray[$x]}, "--------------------------------------\n");
			push(@{$outreportarray[$x]}, "AMUSER $version\n");
			push(@{$outreportarray[$x]}, "output generated: $mday-$mon-$year $hour:$min:$sec\n");
			push(@{$outreportarray[$x]}, "--------------------------------------");
			push(@{$outreportarray[$x]}, "\n\n");
			push(@{$outreportarray[$x]}, "input parameters:\n");
			push(@{$outreportarray[$x]}, "-----------------\n");
			my $tempn = $n-1;
			push(@{$outreportarray[$x]}, "number of fragments: $inputcount[$x]\n");
			unless ($numberofinserts == 0) {
				push(@{$outreportarray[$x]}, "number of inserts: $numberofinserts\n");
			}
			{
				my $tempsc = $sc * 1000;
				my $temppc = $pc_bench * 1000000;
				push(@{$outreportarray[$x]}, "Tm model: NEB / $TM_POLY\n");
				push(@{$outreportarray[$x]}, "salt concentration selected: $tempsc mM\n");
				push(@{$outreportarray[$x]}, "primer concentration selected: $temppc uM\n");
			}
			if($Tmtarget ne "") {
				my $templine = "tm optimized to: $Tmtarget ";
				if ($options{tm} eq "eq") {
					$templine .= "C (all primer Tm approximated to each other)";
				}
				else {
					$templine .= "C (user defined)";
				}
				push(@{$outreportarray[$x]}, "$templine\n");
			}
			push(@{$outreportarray[$x]}, "$cassetteprint\n\n\n");
			push(@{$outreportarray[$x]}, "overview of the needed primers (5'-3'):\n");
			push(@{$outreportarray[$x]}, "---------------------------------------\n\n");
			push(@{$outreportarray[$x]}, "fragment\tstrand\ttm\tprimer\n");
			# print "$optbr[$n-1][1][0]\n$brfor[$n-1][0][$optbr[$n-1][1][0]]\n";
			# print "$n\n";
			push(@{$outreportarray[$x]}, "$inputs[$x][0][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			for (my $i = 0; $i < $n-2; $i++) {
				$ftcomp[$i] = reverse $ftcomp[$i];
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
				#print "$i $brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i}{0}]]\n";
				push(@{$outreportarray[$x]}, "$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
				push(@{$outreportarray[$x]}, "$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n");
				$ftcomp[$i] = reverse $ftcomp[$i];
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
			}
			$ftcomp[$n-2] = reverse $ftcomp[$n-2];
			$brrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]] = reverse $brrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]];
			$primerextension[$x][$n-2][5] = reverse $primerextension[$x][$n-2][5];
			push(@{$outreportarray[$x]}, "$inputs[$x][$n-2][1]\treverse\t$TmGCrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]][3]\C\t$ftcomp[$n-2]$primerextension[$x][$n-2][5]$brrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]]\n\n");
			$ftcomp[$n-2] = reverse $ftcomp[$n-2];
			$brrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]] = reverse $brrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]];
			$primerextension[$x][$n-2][5] = reverse $primerextension[$x][$n-2][5];
			
			#If there is any warnings in the primer details section, the user is given a heads up:
			
			
			my $warningflag = 0;
			
			unless (55 < $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0]) {
				$warningflag = 1;
			}
			
			unless (39.9 < $TmGCfor[0][0][$optbr[0][0][$ranking{0}{0}]][1] and 60.1 > $TmGCfor[0][0][$optbr[0][0][$ranking{0}{0}]][1]) {
				$warningflag = 1;
			}
			
			unless (exists $GCclampfor{0}{0}) {
				$warningflag = 1;
			}
			
			if (exists $GCcountclampfor{0}{0}) {
				$warningflag = 1;
			}
			
			if (exists $primerdimer{0}{0}) {
				$warningflag = 1;
			}
			
			if (exists $homologyfor{0}{0}) {
				$warningflag = 1;
			}	
			
			if (exists $polyNfor{0}{0}) {
				$warningflag = 1;
			} 
			
			for (my $i = 0; $i < $n-1; $i++) {
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				$ftcomp[$i] = reverse $ftcomp[$i];
				
				unless (55 < $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][0]) {
					$warningflag = 1;
				}
				
				unless (39.9 < $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][1] and 60.1 > $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][1]) {
					$warningflag = 1;
				}
			
				unless (exists $GCclamprev{$i}{0}) {
					$warningflag = 1;
				}
				
				if (exists $GCcountclamprev{$i}{0}) {
					$warningflag = 1;
				}
				
				if (exists $primerdimer{$i}{0}) {
					$warningflag = 1;
				}
				
				if (exists $homologyrev{$i}{0}) {
					$warningflag = 1;
				}
				
				if (exists $polyNrev{$i}{0}) {
					$warningflag = 1;
				} 
				
				if (exists $primerTmwarning{$i}) {
					$warningflag = 1;
				}
				
				unless (55 < $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0]) {
					$warningflag = 1;
				}
				
				unless (39.9 < $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1] and 60.1 > $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1]) {
					$warningflag = 1;
				}
			
				unless (exists $GCclampfor{$i+1}{0}) {
					$warningflag = 1;
				}
				
				if (exists $GCcountclampfor{$i+1}{0}) {
					$warningflag = 1;
				}
				
				if (exists $primerdimer{$i+1}{0}) {
					$warningflag = 1;
				}
				
				if (exists $homologyfor{$i+1}{0}) {
					$warningflag = 1;
				}	
				
				if (exists $polyNfor{$i+1}{0}) {
					$warningflag = 1;
				}		
			}
			
			$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
			unless (55 < $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0]) {
				$warningflag = 1;
			}
			
			unless (39.9 < $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1] and 60.1 > $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1]) {
				$warningflag = 1;
			}
			
			unless (exists $GCclamprev{$n-1}{0}) {
				$warningflag = 1;
			}
			
			if (exists $GCcountclamprev{$n-1}{0}) {
				$warningflag = 1;
			}
			
			if (exists $primerdimer{$n-1}{0}) {
				$warningflag = 1;
			}
			
			if (exists $homologyrev{$n-1}{0}) {
				$warningflag = 1;
			}
			
			if (exists $polyNfor{$n-1}{0}) {
				$warningflag = 1;
			} 
			
			if (exists $primerTmwarning{$n-1}) {
				$warningflag = 1;
			}
			
			if ($warningflag == 1) {
				push(@{$outreportarray[$x]},  "attention: one or more of the designed primers may have undesirable properties\n");
				push(@{$outreportarray[$x]},  "please see \"primer details\" section\n\n\n");
			}
			
			
			
			#Then a graphic overview of the non-fusion tail, the fusion tail, and the binding region 
			#is printed to the screenout file
			
			
			
			my $cassetteftnoUpre = $cassette[2];
			$cassetteftnoUpre =~ tr/U/T/;
			my $cassetteftnoUpost = $cassette[3];
			$cassetteftnoUpost =~ tr/U/T/;
			$cassetteftnoUpost = reverse $cassetteftnoUpost;
			$cassetteftnoUpost =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
			
			push(@{$outreportarray[$x]}, "graphic overview of DNA fragments and primers:\n");
			push(@{$outreportarray[$x]}, "----------------------------------------------\n\n");

			for (my $i = 0; $i < $n-1; $i++) {
				if ($n-1 == 1 and $i == 0) {
					push(@{$outreportarray[$x]},  "fusion region and related primers for joining of $inputs[$x][0][1] (3' end) and $inputs[$x][1][1] (5' end):\n\n");
				}
				else {
					push(@{$outreportarray[$x]},  "fusion region and related primers for joining of $inputs[$x][$i][1] and $inputs[$x][$i+1][1]:\n\n");
				}

				my $space = "";
				for (my $ii = 0; $ii < $Apos[$i][0]+5; $ii++) {
					$space .= " ";
				}
				my $temp = "$space";
				$temp .= "5'-$ft[$i]\n";
				push(@{$outreportarray[$x]},  "$temp");

				my $space = "";
				for (my $ii = 0; $ii < $Tpos[$i][0]+8; $ii++) {
					$space .= " ";
				}
				push(@{$outreportarray[$x]}, "$space$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]-3'\n");

				push(@{$outreportarray[$x]},  "5'-[...]$joined[$i]\[...]-3'\n");
				my $space = "";
###NB NB NB NB NB: This should be refined to check for existence of $primerextension throughout printing
#print length $primerextension[$x][$i][5], "\n";
				if ($Apos[$i][0] > (length $joined[$i])/2) {	
					for (my $ii = 0; $ii < ($Apos[$i][0]-$optbr[$i][0][$ranking{$i}{0}]-($Apos[$i][0]-(length $joined[$i])/2))+5; $ii++) {
						$space .= " ";
					}
				}
				else {
					for (my $ii = 0; $ii < $Apos[$i][0]-$optbr[$i][0][$ranking{$i}{0}]+5-length $primerextension[$x][$i][5]; $ii++) {
						$space .= " ";
					}
				}
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				my $temp = "$space";
				$temp .= "3'-$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n";
				push(@{$outreportarray[$x]},  "$temp");
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				my $space = "";
				for (my $ii = 0; $ii < $Apos[$i][0]+8-(length $primerextension[$x][$i][5]); $ii++) {
					$space .= " ";
				}
				$ftcomp[$i] = reverse $ftcomp[$i];
				push(@{$outreportarray[$x]},  "$space$primerextension[$x][$i][5]$ftcomp[$i]-5'\n");
				push(@{$outreportarray[$x]},  "\n");
				$ftcomp[$i] = reverse $ftcomp[$i];
			}
			
			push(@{$outreportarray[$x]},  "\n");	
			
			
			
			#Lastly, the final construct is printed to the screenout for good meassures
			
			
			
			my $final .= $cassette[5];
			$final .= $cassetteftnoUpre;
			
			my $caseflag = 0;
			if ($cassette[5] eq "") {
				$caseflag = 1;
			}
	
			for (my $i = 0; $i < $n-1; $i++) {
				if ($primerextension[$x][$i][0] eq "TRUE") {
					if ($caseflag == 0) {
						$final .= lc $primerextension[$x][$i][1];
						$caseflag = 1;
					}
					elsif ($caseflag == 1) {
						$final .= uc $primerextension[$x][$i][1];
						$caseflag = 0;
					}
				}
				if ($caseflag == 0) {
					$final .= lc $inputs[$x][$i][2];
					$caseflag = 1;
				}
				elsif ($caseflag == 1) {
					$final .= uc $inputs[$x][$i][2];
					$caseflag = 0;
				}
			}
	
			if ($primerextension[$x][$n][0] eq "TRUE") {
				if ($caseflag == 0) {
					$final .= lc $primerextension[$x][$n][1];
					$caseflag = 1;
				}
				elsif ($caseflag == 1) {
					$final .= uc $primerextension[$x][$n][1];
					$caseflag = 0;
				}
			}

			if ($circularoutput eq "TRUE") {
				if ($primerextension[$x][$n-1][0] eq "TRUE") {
					if ($caseflag == 0) {
						$final .= lc $primerextension[$x][$n-1][1];
					}
					if ($caseflag == 1) {
						$final .= uc $primerextension[$x][$n-1][1];
					}
				}
			}
			
			if ($caseflag == 0) {
				$final .= lc $cassetteftnoUpost;
				$final .= lc $cassette[6];
			}
			elsif ($caseflag == 1) {
				$final .= uc $cassetteftnoUpost;
				$final .= uc $cassette[6];
			}
			
			my $templength = length($final);
	
			push(@{$outreportarray[$x]},  "details of final construct after cloning (circular, length: $templength bases):\n");
			push(@{$outreportarray[$x]},  "------------------------------------------------------------------------\n\n");
	

			my $string = "";
			my $insertcounter = 1;
			$string .= "$inputs[$x][0][1],,";
			for (my $i = 1; $i < $n-1; $i++) {
				if ($primerextension[$x][$i][0] eq "TRUE") {
					$string .= "insert$insertcounter,,";
					$insertcounter++;
				}
				$string .= "$inputs[$x][$i][1],,";
			}
			if ($primerextension[$x][$n-1][0] eq "TRUE") {
				$string .= "insert$insertcounter,,";
			}
			chop $string;
			chop $string;

			$string .= "\n\n\n";
			push(@{$outreportarray[$x]},  $string);
			
			
			for (my $pos = 0; $pos < length($final); $pos = $pos + 60) {
				my $string = $pos+1;
				for (my $i=length $pos; $i < 10; $i++) {
					$string .= " ";
				}
				$string .= substr($final, $pos, 60);
				push(@{$outreportarray[$x]}, "$string\n");
			push(@{$outreportarray[$x]}, "");   # terminate the htmlprint construct-sequence walker
			}

			&fullplasmidreportblock($outreportarray[$x], $destplasmid, \@cassette, $final, $inputs[$x][0][1]);
			#Detailed descriptions of the primers are then printed to the screenout file
			
			
			push(@{$outreportarray[$x]},  "primer details:\n");
			push(@{$outreportarray[$x]},  "---------------\n\n");
			push(@{$outreportarray[$x]},  "primer details for $inputs[$x][0][1]\n\n");
			push(@{$outreportarray[$x]},  "forward primer: $ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			
			if ($TM_OPT_MIN <= $TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][0]) {
				push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
			}
			
			if (39.9 < $TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][1] and 60.1 > $TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][1]) {
				push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
			}
			
			if (exists $GCclampfor{$n-1}{0}) {
				push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
			}
			
			if (exists $GCcountclampfor{$n-1}{0}) {
				push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
			}
			
			if (exists $primerdimer{$n-1}{0}) {
				push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
			}
			
			if (exists $homologyfor{$n-1}{0}) {
				push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
			}	
			
			if (exists $polyNfor{$n-1}{0}) {
				push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
			} 
			push(@{$outreportarray[$x]},  "\n");

			$primerextension[$x][0][5] = reverse $primerextension[$x][0][5];
			push(@{$outreportarray[$x]},  "reverse primer: $ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");
			$primerextension[$x][0][5] = reverse $primerextension[$x][0][5];
			if ($TM_OPT_MIN <= $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][0]) {
				push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
			}
			
			if (39.9 < $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][1] and 60.1 > $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][1]) {
				push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
			}
			
			if (exists $GCclampfor{0}{0}) {
				push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
			}
			
			if (exists $GCcountclampfor{0}{0}) {
				push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
			}
			
			if (exists $primerdimer{0}{0}) {
				push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
			}
			
			if (exists $homologyfor{0}{0}) {
				push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
			}	
			
			if (exists $polyNfor{0}{0}) {
				push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
			}
			else {
				push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
			} 
			push(@{$outreportarray[$x]},  "\n");
			
			if (exists $primerTmwarning{0}) {
				push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...NO\n");
			}
			else {
				push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...YES\n");
			}
			push(@{$outreportarray[$x]},  "\n");
			

			for (my $i = 0; $i < $n-2; $i++) {
				#$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				#$ftcomp[$i] = reverse $ftcomp[$i];
				
				push(@{$outreportarray[$x]},  "\nprimer details for $inputs[$x][$i+1][1]\n\n");
				push(@{$outreportarray[$x]},  "forward primer: $ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n");
				
				if ($TM_OPT_MIN <= $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0]) {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
				}
				
				if (39.9 < $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1] and 60.1 > $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1]) {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
				}
			
				if (exists $GCclampfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
				}
				
				if (exists $GCcountclampfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
				}
				
				if (exists $primerdimer{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
				}
				
				if (exists $homologyfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
				}	
				
				if (exists $polyNfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
				}
				push(@{$outreportarray[$x]},  "\n");

				$primerextension[$x][$i+1][5] = reverse $primerextension[$x][$i+1][5];
				push(@{$outreportarray[$x]},  "reverse primer: $ftcomp[$i+1]$primerextension[$x][$i+1][5]$brrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]]\n");
				$primerextension[$x][$i+1][5] = reverse $primerextension[$x][$i+1][5];
	
				if ($TM_OPT_MIN <= $TmGCrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]][0]) {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
				}
				
				if (39.9 < $TmGCrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]][1] and 60.1 > $TmGCrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]][1]) {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[$i+1][0][$optbr[$i+1][0][$ranking{$i+1}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
				}
				
				if (exists $GCclampfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
				}
				
				if (exists $GCcountclampfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
				}
				
				if (exists $primerdimer{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
				}
				
				if (exists $homologyfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
				}	
				
				if (exists $polyNfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
				} 
				push(@{$outreportarray[$x]},  "\n");
				
				if (exists $primerTmwarning{$i+1}) {
					push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...NO\n");
				}
				else {
					push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...YES\n");
				}
			
				push(@{$outreportarray[$x]},  "\n");
				
			}
			
			push(@{$outreportarray[$x]},  "\n\n");


			#printing output to fasta array
			
			if ($batchinput eq "TRUE") {

				# if ($x == 0) {
				# 	# push(@{$fastabatcharray[$x]}, "forward and reverse primers for backbone\n");
				# 	$primerextension[0][0][5] = reverse $primerextension[0][0][5];
				# 	# push(@{$fastabatcharray[$x]}, "$inputs[$x][0][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
				# 	$primerextension[0][0][5] = reverse $primerextension[0][0][5];						
				# }

				for (my $i = 0; $i < $n-2; $i++) {
					#$ftcomp[$i] = reverse $ftcomp[$i];
					#$brrev[$i][0][$optbr[$i][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][$ranking{$i}{0}]];
					# push(@{$fastabatcharray[$x]}, "$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n\n");	
					push(@{$fastaoutarray[$x]},  "$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n");
					
				}
				$primerextension[$x][$n-2][5] = reverse $primerextension[$x][$n-2][5];
				push(@{$fastaoutarray[$x]},  "$inputs[$x][$n-2][1]\treverse\t$TmGCrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]][3]\C\t$ftcomp[$n-2]$primerextension[$x][$n-2][5]$brrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]]\n");	
				$primerextension[$x][$n-2][5] = reverse $primerextension[$x][$n-2][5];
			}
			else {
				#push(@{$fastaoutarray[$x]},  ">FW_$inputs[$x][0][1]_[Tm:$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0]]\n$cassette[2]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
				push(@{$fastaoutarray[$x]},  ">FW_$inputs[$x][-1][1]_[Tm:$TmGCfor[-1][0][$optbr[-1][1][$ranking{-1}{0}]][0]]\n$ft[-1]$primerextension[$x][-1][4]$brfor[-1][0][$optbr[-1][1][$ranking{-1}{0}]]\n");
				for (my $i = 0; $i < $n-2; $i++) {
					#$ftcomp[$i] = reverse $ftcomp[$i];
					#$brrev[$i][0][$optbr[$i][0]] = reverse $brrev[$i][0][$optbr[$i][0]];
					$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
					push(@{$fastaoutarray[$x]},  ">RV_$inputs[$x][$i][1]_[Tm:$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][0]]\n$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");	
					$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
					push(@{$fastaoutarray[$x]},  ">FW_$inputs[$x][$i+1][1]_[Tm:$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0]]\n$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n");
				}
				#$brrev[$n-1][0][$optbr[$n-1][0]] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
				#$cassette[3] = reverse $cassette[3];
				$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
				$primerextension[$x][$n][1] = reverse $primerextension[$x][$n][1];
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				push(@{$fastaoutarray[$x]},  ">RV_$inputs[$x][$n-2][1]_[Tm:$TmGCrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]][0]]\n$ftcomp[$n-2]$primerextension[$x][$n-1][1]$brrev[$n-2][0][$optbr[$n-2][0][$ranking{$n-2}{0}]]\n");
				$primerextension[$x][$n][1] = reverse $primerextension[$x][$n][1];
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
			}
		}





















	## NON-CIRCULAR MULTI OUTPUT PRINTING
	###################################################################



		else {
		#For the printing, all reverse primers must be complemented


			
			
			
			#First a quick overview of the needed primers is printed to the 
			#screen output report
			
			push(@{$outreportarray[$x]}, "--------------------------------------\n");
			push(@{$outreportarray[$x]}, "AMUSER $version\n");
			push(@{$outreportarray[$x]}, "output generated: $mday-$mon-$year $hour:$min:$sec\n");
			push(@{$outreportarray[$x]}, "--------------------------------------");
			push(@{$outreportarray[$x]}, "\n\n");
			push(@{$outreportarray[$x]}, "input parameters:\n");
			push(@{$outreportarray[$x]}, "-----------------\n");
			push(@{$outreportarray[$x]}, "number of fragments: $inputcount[$x]\n");
			unless ($numberofinserts == 0) {
				push(@{$outreportarray[$x]}, "number of inserts: $numberofinserts\n");
			}
			{
				my $tempsc = $sc * 1000;
				my $temppc = $pc_bench * 1000000;
				push(@{$outreportarray[$x]}, "Tm model: NEB / $TM_POLY\n");
				push(@{$outreportarray[$x]}, "salt concentration selected: $tempsc mM\n");
				push(@{$outreportarray[$x]}, "primer concentration selected: $temppc uM\n");
			}
			if($Tmtarget ne "") {
				my $templine = "tm optimized to: $Tmtarget ";
				if ($options{tm} eq "eq") {
					$templine .= "C (all primer Tm approximated to each other)";
				}
				else {
					$templine .= "C (user defined)";
				}
				push(@{$outreportarray[$x]}, "$templine\n");
			}
			push(@{$outreportarray[$x]}, "$cassetteprint\n\n\n");
			push(@{$outreportarray[$x]}, "overview of the needed primers (5'-3'):\n");
			push(@{$outreportarray[$x]}, "---------------------------------------\n\n");
			push(@{$outreportarray[$x]}, "fragment\tstrand\ttm\tprimer\n");
			# per-fragment tag flaps: FW gets the 5' tag (sense); RV gets the 3'
			# tag reverse-complemented. Spliced after the bridge, before binding.
			my (@fwtag, @rvtag, @link5, @link3);
			for (my $f = 0; $f < $n; $f++) {
				$link5[$f] = defined $tag5{$f} ? ($tag5{$f}{linker} // "") : "";
				$link3[$f] = defined $tag3{$f} ? ($tag3{$f}{linker} // "") : "";
				$fwtag[$f] = defined $tag5{$f} ? ($tag5{$f}{seq} . $link5[$f]) : "";
				$rvtag[$f] = defined $tag3{$f} ? &revcomp($link3[$f] . $tag3{$f}{seq}) : "";
			}
			# Oligo mode: for each tag flagged "clone as oligo", the gene primer uses an
			# oligo overhang (no flap) and the adjacent neighbour primer is adjusted to
			# match; the tag is emitted as a standalone annealed oligo. $fwlead[f]/$rvlead[f]
			# override the normal [USER tail][bridge][flap] lead when defined.
			my (@fwlead, @rvlead, @oligo_recs);
			# Stashed overhang data for the oligo fusion diagrams, keyed by fragment
			# index (see render_oligo_diag). _5 = from a 5' tag's oligo, _3 = from a 3' tag's.
			my (%mo5_5, %mo3_5, %mbr5_5, %mbr3_5, %mo5_3, %mo3_3, %mbr5_3, %mbr3_3, %mcore_5, %mcore_3);
			my $csc5 = &cassette_scar5($cassette[2]);
			my $csc3 = &cassette_scar3($cassette[3]);
			for (my $f = 0; $f < $n; $f++) {
				if (defined $tag5{$f} and $tag5{$f}{oligo}) {
					my $inoligo = ($link5[$f] ne "" && $tag5{$f}{linkoligo} && &linker_in_oligo($tag5{$f}{seq}, $link5[$f], 5));
					warn "WARNING: the linker on 5' tag '$tag5{$f}{name}' cannot form the oligo overhang; placed on the gene primer instead.\n" if ($link5[$f] ne "" && $tag5{$f}{linkoligo} && !$inoligo);
					my $core = $inoligo ? ($tag5{$f}{seq} . $link5[$f]) : $tag5{$f}{seq};
					if (&oligo_viable($core)) {
						my $five = ($f > 0) ? 'tag' : $csc5;
						my $r = &compute_oligo($core, $five, 'tag', 8);
						($mo5_5{$f},$mo3_5{$f},$mbr5_5{$f},$mbr3_5{$f},$mcore_5{$f}) = ($r->{o5},$r->{o3},$r->{br5},$r->{br3},$r->{core});
						$fwlead[$f] = $r->{tail3} . $r->{br3};
						$fwlead[$f] .= $link5[$f] unless $inoligo;   # linker on the gene primer when not in the oligo
						$rvlead[$f-1] = $r->{tail5} . &revcomp($r->{br5}) if $f > 0;
						push @oligo_recs, { name=>$tag5{$f}{name}.($inoligo?"+linker":""), frag=>$inputs[$x][$f][1], end=>5,
						                    fw=>$r->{oligo_fw}, rv=>$r->{oligo_rv} };
					} else { $tag5{$f}{oligo}=0; warn "WARNING: tag '$tag5{$f}{name}'".($inoligo?"+linker":"")." cannot form a viable oligo; using a primer flap instead.\n"; }
				}
				if (defined $tag3{$f} and $tag3{$f}{oligo}) {
					my $inoligo = ($link3[$f] ne "" && $tag3{$f}{linkoligo} && &linker_in_oligo($tag3{$f}{seq}, $link3[$f], 3));
					warn "WARNING: the linker on 3' tag '$tag3{$f}{name}' cannot form the oligo overhang; placed on the gene primer instead.\n" if ($link3[$f] ne "" && $tag3{$f}{linkoligo} && !$inoligo);
					my $core = $inoligo ? ($link3[$f] . $tag3{$f}{seq}) : $tag3{$f}{seq};
					if (&oligo_viable($core)) {
						my $three = ($f < $n-1) ? 'tag' : $csc3;
						my $r = &compute_oligo($core, 'tag', $three, 8);
						($mo5_3{$f},$mo3_3{$f},$mbr5_3{$f},$mbr3_3{$f},$mcore_3{$f}) = ($r->{o5},$r->{o3},$r->{br5},$r->{br3},$r->{core});
						$rvlead[$f] = $r->{tail5} . &revcomp($r->{br5});
						$rvlead[$f] .= &revcomp($link3[$f]) unless $inoligo;   # linker on the gene primer when not in the oligo
						$fwlead[$f+1] = $r->{tail3} . $r->{br3} if $f < $n-1;
						push @oligo_recs, { name=>$tag3{$f}{name}.($inoligo?"+linker":""), frag=>$inputs[$x][$f][1], end=>3,
						                    fw=>$r->{oligo_fw}, rv=>$r->{oligo_rv} };
					} else { $tag3{$f}{oligo}=0; warn "WARNING: tag '$tag3{$f}{name}'".($inoligo?"+linker":"")." cannot form a viable oligo; using a primer flap instead.\n"; }
				}
			}
				push(@{$outreportarray[$x]}, "$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t@{[defined $fwlead[0] ? $fwlead[0] : qq($cassette[2]$primerextension[$x][0][4]$fwtag[0])]}$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
			for (my $i = 0; $i < $n-1; $i++) {
				$ftcomp[$i] = reverse $ftcomp[$i];
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
				#print "$i $brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i}{0}]]\n";
				push(@{$outreportarray[$x]},  "$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t@{[defined $rvlead[$i] ? $rvlead[$i] : qq($ftcomp[$i]$primerextension[$x][$i][5]$rvtag[$i])]}$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
				push(@{$outreportarray[$x]},  "$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t@{[defined $fwlead[$i+1] ? $fwlead[$i+1] : qq($ft[$i]$primerextension[$x][$i+1][4]$fwtag[$i+1])]}$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n");
				$ftcomp[$i] = reverse $ftcomp[$i];
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
			}
				$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
				$primerextension[$x][$n][1] = reverse $primerextension[$x][$n][1];
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				push(@{$outreportarray[$x]},  "$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t@{[defined $rvlead[$n-1] ? $rvlead[$n-1] : qq($cassette[3]$primerextension[$x][$n][1]$rvtag[$n-1])]}$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n\n");
				$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
				$primerextension[$x][$n][1] = reverse $primerextension[$x][$n][1];
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				# Oligo strands (ordered as annealed FW+RV); listed here so they also flow into the primer CSV.
				for my $orec (@oligo_recs) {
					push(@{$outreportarray[$x]}, "$orec->{frag} $orec->{name} oligo\tforward\t-\t$orec->{fw}\n");
					push(@{$outreportarray[$x]}, "$orec->{frag} $orec->{name} oligo\treverse\t-\t$orec->{rv}\n");
				}
				if (@oligo_recs) {
					push(@{$outreportarray[$x]}, "\nTags above marked \"oligo\" are ordered as the annealed FW+RV strands listed and inserted into the USER fusion; the flanking gene primers already carry the matching USER overhangs (no tag flap). The graphic fusion diagram below shows the direct gene junctions and does not depict the oligo.\n");
				}


			#If there is any warnings in the primer details section, the user is given a heads up:
			
			
			my $warningflag = 0;
			
			unless (55 < $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0]) {
				$warningflag = 1;
			}
			
			unless (39.9 < $TmGCfor[0][0][$optbr[0][0][$ranking{0}{0}]][1] and 60.1 > $TmGCfor[0][0][$optbr[0][0][$ranking{0}{0}]][1]) {
				$warningflag = 1;
			}
			
			unless (exists $GCclampfor{0}{0}) {
				$warningflag = 1;
			}
			
			if (exists $GCcountclampfor{0}{0}) {
				$warningflag = 1;
			}
			
			if (exists $primerdimer{0}{0}) {
				$warningflag = 1;
			}
			
			if (exists $homologyfor{0}{0}) {
				$warningflag = 1;
			}	
			
			if (exists $polyNfor{0}{0}) {
				$warningflag = 1;
			} 
			
			for (my $i = 0; $i < $n-1; $i++) {
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				$ftcomp[$i] = reverse $ftcomp[$i];
				
				unless (55 < $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][0]) {
					$warningflag = 1;
				}
				
				unless (39.9 < $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][1] and 60.1 > $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][1]) {
					$warningflag = 1;
				}
			
				unless (exists $GCclamprev{$i}{0}) {
					$warningflag = 1;
				}
				
				if (exists $GCcountclamprev{$i}{0}) {
					$warningflag = 1;
				}
				
				if (exists $primerdimer{$i}{0}) {
					$warningflag = 1;
				}
				
				if (exists $homologyrev{$i}{0}) {
					$warningflag = 1;
				}
				
				if (exists $polyNrev{$i}{0}) {
					$warningflag = 1;
				} 
				
				if (exists $primerTmwarning{$i}) {
					$warningflag = 1;
				}
				
				unless (55 < $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0]) {
					$warningflag = 1;
				}
				
				unless (39.9 < $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1] and 60.1 > $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1]) {
					$warningflag = 1;
				}
			
				unless (exists $GCclampfor{$i+1}{0}) {
					$warningflag = 1;
				}
				
				if (exists $GCcountclampfor{$i+1}{0}) {
					$warningflag = 1;
				}
				
				if (exists $primerdimer{$i+1}{0}) {
					$warningflag = 1;
				}
				
				if (exists $homologyfor{$i+1}{0}) {
					$warningflag = 1;
				}	
				
				if (exists $polyNfor{$i+1}{0}) {
					$warningflag = 1;
				}		
			}
			
			$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
			unless (55 < $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0]) {
				$warningflag = 1;
			}
			
			unless (39.9 < $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1] and 60.1 > $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1]) {
				$warningflag = 1;
			}
			
			unless (exists $GCclamprev{$n-1}{0}) {
				$warningflag = 1;
			}
			
			if (exists $GCcountclamprev{$n-1}{0}) {
				$warningflag = 1;
			}
			
			if (exists $primerdimer{$n-1}{0}) {
				$warningflag = 1;
			}
			
			if (exists $homologyrev{$n-1}{0}) {
				$warningflag = 1;
			}
			
			if (exists $polyNfor{$n-1}{0}) {
				$warningflag = 1;
			} 
			
			if (exists $primerTmwarning{$n-1}) {
				$warningflag = 1;
			}
			
			if ($warningflag == 1) {
				push(@{$outreportarray[$x]},  "attention: one or more of the designed primers may have undesirable properties\n");
				push(@{$outreportarray[$x]},  "please see \"primer details\" section\n\n\n");
			}
			
			
			#Then a graphic overview of the non-fusion tail, the fusion tail, and the binding region 
			#is printed to the screenout file
			
			
			
			my $cassetteftnoUpre = $cassette[2];
			$cassetteftnoUpre =~ tr/U/T/;
			my $cassetteftnoUpost = $cassette[3];
			$cassetteftnoUpost =~ tr/U/T/;
			$cassetteftnoUpost = reverse $cassetteftnoUpost;
			$cassetteftnoUpost =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
			
			push(@{$outreportarray[$x]}, "graphic overview of DNA fragments and primers:\n");
			push(@{$outreportarray[$x]}, "----------------------------------------------\n\n");
			
				if ($nocassette eq "TRUE") {
					push(@{$outreportarray[$x]}, "forward primer for $inputs[$x][0][1]:\n\n");
					push(@{$outreportarray[$x]}, "5'-$primerextension[$x][0][4]".((defined $tag5{0} && !$tag5{0}{oligo})?"[tag]":"")."$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]-3'\n");		
				}
			
				# else {
				# 	push(@{$outreportarray[$x]}, "fusion region and forward primer for joining of $cassette[4] and $inputs[$x][0][1]:\n\n");
				# 	my $space = "";
				# 	for (my $ii = 0; $ii < (length $cassette[2])+3+(length $cassette[5]); $ii++) {
				# 		$space .= " ";
				# 	}
				# 	push(@{$outreportarray[$x]}, "$space");
				# 	my $space = "";
				# 	push(@{$outreportarray[$x]}, "$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]-3'\n");
				# 	for (my $ii = 0; $ii < (length $cassette[5]); $ii++) {
				# 		$space .= " ";
				# 	}
				# 	push(@{$outreportarray[$x]}, "$space");
				# 	push(@{$outreportarray[$x]}, "5'-$cassette[2]\n");
				# }
	
				else {
					if (defined $tag5{0} && $tag5{0}{oligo}) {
						# N-terminal oligo tag at the 5' cassette junction (multi-fragment).
						my ($cctx,$ltr) = &oligodiag_ctx($cassette[5], 30, 1);
						my ($gctx,$rtr) = &oligodiag_ctx($brickseqfrags[0][0], 30, 0);
						push(@{$outreportarray[$x]}, &oligodiag_line(5, $cassette[4], $inputs[$x][0][1], $tag5{0}{name}, $cassette[4], $inputs[$x][0][1],
							$ltr, $cctx, "", $mo5_5{0}, "cassette", $mo3_5{0}, "tag", $mbr3_5{0}, $gctx, $rtr, "cassette", $inputs[$x][0][1], $mcore_5{0}));
					} else {
						push(@{$outreportarray[$x]}, "fusion region and forward primer for joining of $cassette[4] and $inputs[$x][0][1]:\n\n");
					}
					my $space = "";
					for (my $ii = 0; $ii < (length $cassette[5]); $ii++) {
						$space .= " ";
					}
					push(@{$outreportarray[$x]}, "$space");
					push(@{$outreportarray[$x]}, "5'-$cassette[2]$primerextension[$x][0][4]".((defined $tag5{0} && !$tag5{0}{oligo})?"[tag]":"")."\n");
					my $space = "";
					for (my $ii = 0; $ii < (length $cassette[2])+3+(length $cassette[5])+(length $primerextension[$x][0][1]); $ii++) {
						$space .= " ";
					}
					push(@{$outreportarray[$x]}, "$space$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]-3'\n");
		
				}				
				
				push(@{$outreportarray[$x]}, "5'-$cassette[5]$cassetteftnoUpre$primerextension[$x][0][4]$brickseqfrags[0][0]\[...]-3'\n\n");
			push(@{$outreportarray[$x]}, "OLIGOSKIP_END\n") if (defined $tag5{0} && $tag5{0}{oligo});

			for (my $i = 0; $i < $n-1; $i++) {
				my $oligo_int = ((defined $tag3{$i} && $tag3{$i}{oligo}) || (defined $tag5{$i+1} && $tag5{$i+1}{oligo}));
				if ($oligo_int) {
					my $t3 = (defined $tag3{$i}   && $tag3{$i}{oligo})   ? 1 : 0;
					my $t5 = (defined $tag5{$i+1} && $tag5{$i+1}{oligo}) ? 1 : 0;
					if ($t3 && $t5) {
						# two oligos abut at one internal junction (gene - tag - tag - gene):
						# keep the honest text pointer rather than draw an ambiguous diagram.
						push(@{$outreportarray[$x]}, &oligo_skip_start($inputs[$x][$i][1], $inputs[$x][$i+1][1],
							$tag3{$i}{name}, $tag5{$i+1}{name}));
					} else {
						# one oligo tag between two genes: gene[i] - tag - gene[i+1].
						my ($o5,$o3,$br5,$br3,$tname,$tcore) = $t3
							? ($mo5_3{$i},   $mo3_3{$i},   $mbr5_3{$i},   $mbr3_3{$i},   $tag3{$i}{name},   $mcore_3{$i})
							: ($mo5_5{$i+1}, $mo3_5{$i+1}, $mbr5_5{$i+1}, $mbr3_5{$i+1}, $tag5{$i+1}{name}, $mcore_5{$i+1});
						my ($lctx,$ltr) = &oligodiag_ctx($brickseqfrags[$i][1],   30, 1);
						my ($rctx,$rtr) = &oligodiag_ctx($brickseqfrags[$i+1][0], 30, 0);
						push(@{$outreportarray[$x]}, &oligodiag_line(5, $inputs[$x][$i][1], $inputs[$x][$i+1][1], $tname, $cassette[4], $inputs[$x][$i][1],
							$ltr, $lctx, $br5, $o5, "tag", $o3, "tag", $br3, $rctx, $rtr, $inputs[$x][$i][1], $inputs[$x][$i+1][1], $tcore));
					}
				} else {
					push(@{$outreportarray[$x]},  "fusion region and related primers for joining of $inputs[$x][$i][1] and $inputs[$x][$i+1][1]:\n\n");
				}

				my $space = "";
				for (my $ii = 0; $ii < $Apos[$i][0]+5; $ii++) {
					$space .= " ";
				}
				my $temp = "$space";
				# Forward fusion tail (USER overhang) PLUS any forced-overhang
				# bridge: both are 5' flap that does not anneal to this fragment,
				# so both are lifted onto this raised line.
				$temp .= "5'-$ft[$i]$primerextension[$x][$i+1][4]".((defined $tag5{$i+1} && !$tag5{$i+1}{oligo})?"[tag]":"")."\n";
				push(@{$outreportarray[$x]},  "$temp");

				my $space = "";
				# Annealing region now starts after the bridge (at the binding
				# position), so advance the offset by the bridge length.
				for (my $ii = 0; $ii < $Tpos[$i][0]+8+(length $primerextension[$x][$i+1][4]); $ii++) {
					$space .= " ";
				}
				push(@{$outreportarray[$x]}, "$space$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]-3'\n");
				push(@{$outreportarray[$x]},  "5'-[...]$joined[$i]\[...]-3'\n");
				
				
				my $space = "";
###NB NB	 NB NB NB: This should be refined to check for existence of $primerextension throughout printing
				if ($Apos[$i][0] > (length $joined[$i])/2) {	
					for (my $ii = 0; $ii < ($Apos[$i][0]-$optbr[$i][0][$ranking{$i}{0}]-length $primerextension[$x][$i][5])+5; $ii++) {
						$space .= " ";
					}
				}
				else {
					for (my $ii = 0; $ii < $Apos[$i][0]-$optbr[$i][0][$ranking{$i}{0}]+5- length $primerextension[$x][$i][5]; $ii++) {
						$space .= " ";
					}
				}
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				my $temp = "$space";
				# Reverse-primer annealing region only (bridge lifted below).
				$temp .= "3'-$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n";
				push(@{$outreportarray[$x]},  "$temp");
				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				my $space = "";
				# Raised reverse flap starts at the binding position (brrevpos =
				# Apos - bridge length), carrying the bridge then the USER overhang.
				for (my $ii = 0; $ii < $Apos[$i][0]+8-(length $primerextension[$x][$i][5]); $ii++) {
					$space .= " ";
				}
				$ftcomp[$i] = reverse $ftcomp[$i];
				my $rvtaglbl = (defined $tag3{$i} && !$tag3{$i}{oligo}) ? "[tag]" : "";
				my $rvtagsp = $space;
				$rvtagsp = substr($rvtagsp,0,length($rvtagsp)-length($rvtaglbl)) if length($rvtagsp) >= length($rvtaglbl);
				push(@{$outreportarray[$x]},  "$rvtagsp$rvtaglbl$primerextension[$x][$i][5]$ftcomp[$i]-5'\n");
				push(@{$outreportarray[$x]},  "\n");
				push(@{$outreportarray[$x]}, "OLIGOSKIP_END\n") if $oligo_int;
				$ftcomp[$i] = reverse $ftcomp[$i];
			}
			
				if ($nocassette eq "TRUE") {
					push(@{$outreportarray[$x]},  "reverse primer for $inputs[$x][$n-1][1]:\n\n");
					my $templine = "5'-[...]$brickseqfrags[$n-1][1]$primerextension[$x][$n][1]-3'";
					push(@{$outreportarray[$x]},  "$templine\n");
				}
				else {
					if (defined $tag3{$n-1} && $tag3{$n-1}{oligo}) {
						# C-terminal oligo tag at the 3' cassette junction (multi-fragment).
						my ($gctx,$ltr) = &oligodiag_ctx($brickseqfrags[$n-1][1], 30, 1);
						my ($cctx,$rtr) = &oligodiag_ctx($cassette[6], 30, 0);
						push(@{$outreportarray[$x]}, &oligodiag_line(3, $inputs[$x][$n-1][1], $cassette[4], $tag3{$n-1}{name}, $cassette[4], $inputs[$x][$n-1][1],
							$ltr, $gctx, $mbr5_3{$n-1}, $mo5_3{$n-1}, "tag", $mo3_3{$n-1}, "cassette", "", $cctx, $rtr, $inputs[$x][$n-1][1], "cassette", $mcore_3{$n-1}));
					} else {
						push(@{$outreportarray[$x]},  "fusion region and reverse primer for joining of $inputs[$x][$n-1][1] and $cassette[4]:\n\n");
					}
					push(@{$outreportarray[$x]},  "5'-[...]$brickseqfrags[$n-1][1]$primerextension[$x][$n][1]$cassetteftnoUpost$cassette[6]-3'\n");
				}
					
				
				my $space = "";
				for (my $ii = 0; $ii < ((length $brickseqfrags[$n-1][1])+5)-(length $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]); $ii++) {
					$space .= " ";
				}
				my $templine = $space;
				
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
				if ($nocassette eq "TRUE") {
					$templine .= "3'-$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]".((defined $tag3{$n-1} && !$tag3{$n-1}{oligo})?"[tag]":"")."$primerextension[$x][$n][1]-5'";
					push(@{$outreportarray[$x]},  "$templine\n");
				}
				else {
					$templine .= "3'-$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]";
					push(@{$outreportarray[$x]},  "$templine\n");
				}
				$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
	
				if ($nocassette eq "FALSE") {
					my $space = "";
					for (my $ii = 0; $ii < (length $brickseqfrags[$n-1][1])+8; $ii++) {
						$space .= " ";
					}
					$cassette[3] = reverse $cassette[3];
					$primerextension[$x][$n-1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
					my $rvtaglbl2 = (defined $tag3{$n-1} && !$tag3{$n-1}{oligo}) ? "[tag]" : "";
					my $rvtagsp2 = $space;
					$rvtagsp2 = substr($rvtagsp2,0,length($rvtagsp2)-length($rvtaglbl2)) if length($rvtagsp2) >= length($rvtaglbl2);
					push(@{$outreportarray[$x]},  "$rvtagsp2$rvtaglbl2$primerextension[$x][$n][1]$cassette[3]-5'");
					$cassette[3] = reverse $cassette[3];
					$primerextension[$x][$n-1][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				}
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;

			push(@{$outreportarray[$x]},  "\n\n\n");
			push(@{$outreportarray[$x]}, "OLIGOSKIP_END\n") if ($nocassette eq "FALSE" && defined $tag3{$n-1} && $tag3{$n-1}{oligo});
			
			

			
			#Lastly, the final construct is printed to the screenout for good meassures
			
			
			
			my $final .= $cassette[5];
			$final .= $cassetteftnoUpre;
			
			my $caseflag = 0;
			if ($cassette[5] eq "") {
				$caseflag = 1;
			}
	
			for (my $i = 0; $i < $n; $i++) {
				if ($primerextension[$x][$i][0] eq "TRUE") {
					if ($caseflag == 0) {
						$final .= lc $primerextension[$x][$i][1];
						$caseflag = 1;
					}
					elsif ($caseflag == 1) {
						$final .= uc $primerextension[$x][$i][1];
						$caseflag = 0;
					}
				}
				# Emit the 5' tag (sense), the fragment, and the 3' tag (sense)
				# each as its OWN alternating case-run, so each becomes a
				# separate coloured segment in the construct view (matched by a
				# tag entry added to the name row below).
				my $fl5 = (defined $tag5{$i} && defined $tag5{$i}{linker}) ? $tag5{$i}{linker} : "";
				my $fl3 = (defined $tag3{$i} && defined $tag3{$i}{linker}) ? $tag3{$i}{linker} : "";
				if (defined $tag5{$i}) {
					if ($caseflag == 0) { $final .= lc $tag5{$i}{seq}; $caseflag = 1; }
					else                { $final .= uc $tag5{$i}{seq}; $caseflag = 0; }
				}
				if ($fl5 ne "") {   # 5' linker: between tag and fragment
					if ($caseflag == 0) { $final .= lc $fl5; $caseflag = 1; }
					else                { $final .= uc $fl5; $caseflag = 0; }
				}
				if ($caseflag == 0) { $final .= lc $inputs[$x][$i][2]; $caseflag = 1; }
				else                { $final .= uc $inputs[$x][$i][2]; $caseflag = 0; }
				if ($fl3 ne "") {   # 3' linker: between fragment and tag
					if ($caseflag == 0) { $final .= lc $fl3; $caseflag = 1; }
					else                { $final .= uc $fl3; $caseflag = 0; }
				}
				if (defined $tag3{$i}) {
					if ($caseflag == 0) { $final .= lc $tag3{$i}{seq}; $caseflag = 1; }
					else                { $final .= uc $tag3{$i}{seq}; $caseflag = 0; }
				}
			}
	
			if ($primerextension[$x][$n][0] eq "TRUE") {
				if ($caseflag == 0) {
					$final .= lc $primerextension[$x][$n][1];
					$caseflag = 1;
				}
				elsif ($caseflag == 1) {
					$final .= uc $primerextension[$x][$n][1];
					$caseflag = 0;
				}
			}
			
			if ($caseflag == 0) {
				$final .= lc $cassetteftnoUpost;
				$final .= lc $cassette[6];
			}
			elsif ($caseflag == 1) {
				$final .= uc $cassetteftnoUpost;
				$final .= uc $cassette[6];
			}
			
			my $templength = length($final);
	
			push(@{$outreportarray[$x]},  "details of final construct after cloning (length: $templength bases):\n");
			push(@{$outreportarray[$x]},  "--------------------------------------------------------------\n\n");
		
			
			# my $linkernumber = 1;
			my $string = "";
			if ($cassette[4] ne "") {
				$string .= "$cassette[4],,";
			}
			my $insertcounter = 1;
			for (my $i = 0; $i < $n; $i++) {
				if ($primerextension[$x][$i][0] eq "TRUE") {
					$string .= "insert$insertcounter,,";
					$insertcounter++;
				}
				$string .= "$tag5{$i}{name},," if defined $tag5{$i};
				$string .= "linker,," if (defined $tag5{$i} && defined $tag5{$i}{linker} && $tag5{$i}{linker} ne "");
				$string .= "$inputs[$x][$i][1],,";
				$string .= "linker,," if (defined $tag3{$i} && defined $tag3{$i}{linker} && $tag3{$i}{linker} ne "");
				$string .= "$tag3{$i}{name},," if defined $tag3{$i};
			}
			if ($primerextension[$x][$n][0] eq "TRUE") {
				$string .= "insert$insertcounter,,";
				$insertcounter++;
			}
			if ($cassette[4] ne "") {
				$string .= "$cassette[4]";
			}
			else {
				chop $string;
				chop $string;
			}
			$string .= "\n\n\n";
			push(@{$outreportarray[$x]},  $string);
			
			for (my $pos = 0; $pos < length($final); $pos = $pos + 60) {
				my $string = $pos+1;
				for (my $i=length $pos; $i < 10; $i++) {
					$string .= " ";
				}
				$string .= substr($final, $pos, 60);
				push(@{$outreportarray[$x]}, "$string\n");
			}
			push(@{$outreportarray[$x]}, "");   # terminate the htmlprint construct-sequence walker

			&fullplasmidreportblock($outreportarray[$x], $destplasmid, \@cassette, $final, $inputs[$x][0][1]);

			#Detailed descriptions of the primers are then printed to the screenout file
			
			
			push(@{$outreportarray[$x]},  "primer details:\n");
			push(@{$outreportarray[$x]},  "---------------\n\n");
			push(@{$outreportarray[$x]},  "primer details for $inputs[$x][0][1]\n\n");
			my (@fwtagd, @rvtagd);   # tag flaps for the detail section (same as the overview)
			for (my $f = 0; $f < $n; $f++) {
				$fwtagd[$f] = defined $tag5{$f} ? $tag5{$f}{seq} : "";
				$rvtagd[$f] = defined $tag3{$f} ? &revcomp($tag3{$f}{seq}) : "";
			}
				push(@{$outreportarray[$x]},  "forward primer: @{[defined $fwlead[0] ? $fwlead[0] : qq($cassette[2]$primerextension[$x][0][4]$fwtagd[0])]}$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");

				if ($TM_OPT_MIN <= $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0]) {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
				}
				
				if (39.9 < $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][1] and 60.1 > $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][1]) {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
				}
				
				if (exists $GCclampfor{0}{0}) {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
				}
				
				if (exists $GCcountclampfor{0}{0}) {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
				}
				
				if (exists $primerdimer{0}{0}) {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
				}
				
				if (exists $homologyfor{0}{0}) {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
				}	
				
				if (exists $polyNfor{0}{0}) {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
				} 
				push(@{$outreportarray[$x]},  "\n");
			
			for (my $i = 0; $i < $n-1; $i++) {
				#$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
				#$ftcomp[$i] = reverse $ftcomp[$i];
				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
				push(@{$outreportarray[$x]},  "reverse primer: @{[defined $rvlead[$i] ? $rvlead[$i] : qq($ftcomp[$i]$primerextension[$x][$i][5]$rvtagd[$i])]}$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
	
				if ($TM_OPT_MIN <= $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][0]) {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
				}
				
				if (39.9 < $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][1] and 60.1 > $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][1]) {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
				}
				
				if (exists $GCclampfor{$i}{0}) {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
				}
				
				if (exists $GCcountclampfor{$i}{0}) {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
				}
				
				if (exists $primerdimer{$i}{0}) {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
				}
				
				if (exists $homologyfor{$i}{0}) {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
				}	
				
				if (exists $polyNfor{$i}{0}) {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
				} 
				push(@{$outreportarray[$x]},  "\n");
				
				if (exists $primerTmwarning{$i}) {
					push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...NO\n");
				}
				else {
					push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...YES\n");
				}
			
				push(@{$outreportarray[$x]},  "\n");
				
				push(@{$outreportarray[$x]},  "\nprimer details for $inputs[$x][$i+1][1]\n\n");
				push(@{$outreportarray[$x]},  "forward primer: @{[defined $fwlead[$i+1] ? $fwlead[$i+1] : qq($ft[$i]$primerextension[$x][$i+1][4]$fwtagd[$i+1])]}$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n");
				
				if ($TM_OPT_MIN <= $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0] and $TM_OPT_MAX >= $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0]) {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
				}
				
				if (39.9 < $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1] and 60.1 > $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1]) {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
				}
			
				if (exists $GCclampfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
				}
				
				if (exists $GCcountclampfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
				}
				
				if (exists $primerdimer{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
				}
				
				if (exists $homologyfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
				}	
				
				if (exists $polyNfor{$i+1}{0}) {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
				}
				push(@{$outreportarray[$x]},  "\n");
			}
			
			
			
			#$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
				$primerextension[$x][$n][1] = reverse $primerextension[$x][$n][1];
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				push(@{$outreportarray[$x]},  "reverse primer: @{[defined $rvlead[$n-1] ? $rvlead[$n-1] : qq($cassette[3]$primerextension[$x][$n][1]$rvtagd[$n-1])]}$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
				$primerextension[$x][$n][1] = reverse $primerextension[$x][$n][1];
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		
				if ($TM_OPT_MIN <= $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0] and $TM_OPT_MAX >= $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0]) {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C - in optimal range (55-72)?                   ...YES\n");
				}										
				else {
					push(@{$outreportarray[$x]},  "* Tm: $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C - in optimal range (55-72)?                   ...NO\n");
				}
				
				if (39.9 < $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1] and 60.1 > $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1]) {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1]\% - in optimal range (40-60)?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC content: $TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][1]\% - in optimal range (40-60)?          ...NO\n");
				}
				
				if (exists $GCclampfor{$n-1}{0}) {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* GC clamp present at 3' end?                             ...NO\n");
				}
				
				if (exists $GCcountclampfor{$n-1}{0}) {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* more than 3 G/C out of last 5 bases at 3' end?          ...NO\n");
				}
				
				if (exists $primerdimer{$n-1}{0}) {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of primer dimer formation in primer pair?          ...NO\n");
				}
				
				if (exists $homologyfor{$n-1}{0}) {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* risk of intra-primer homology (secondary structures)?   ...NO\n");
				}	
				
				if (exists $polyNfor{$n-1}{0}) {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...YES\n");
				}
				else {
					push(@{$outreportarray[$x]},  "* presence of polyN stretches?                            ...NO\n");
				}
				
				push(@{$outreportarray[$x]},  "\n");
					
			if (exists $primerTmwarning{$n-1}) {
				push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...NO\n");
			}
			else {
				push(@{$outreportarray[$x]},  "Tm of primers within 2C of each other?                    ...YES\n");
			}
			push(@{$outreportarray[$x]},  "\n\n");
			&genbank_outputs($x, \@inputs, \@primerextension, \@cassette, $final, $destplasmid, $outreportarray[$x], $htmlreportdir, basename($options{h}), \%tag5, \%tag3);

			close OUT;

			#printing output to fasta array

			# if ($batchinput eq "TRUE") {
			# 	if ($circularoutput eq "TRUE") {
			# 		if ($x == 0) {
			# 			push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$cassette[2]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
			# 			push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\treverse\t$TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");	
			# 		}
	
			# 		for (my $i = 1; $i < $n-1; $i++) {
			# 			my $insertno = $x+1;
			# 			$ftcomp[$i] = reverse $ftcomp[$i];
			# 			$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
			# 			$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];												  
			# 			push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\tforward\t$TmGCfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i][4]$brfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]]\n");
			# 			push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
			# 			$ftcomp[$i] = reverse $ftcomp[$i];
			# 			$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
			# 			$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];	
			# 		}
					
			# 		#double backbone
			# 		if ($backbone == 1 and $x == $numberofbatches-1) {
			# 			                                    #"$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n"
			# 			push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[0][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			# 			# push(@{$fastabatcharray[$x]}, "$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			# 			$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
			# 			# push(@{$fastaoutarray[$x]},  "$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$ftcomp[$n-1]$primerextension[$x][$n-1][5]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");	
			# 			$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
			# 		}

			# 		#single backbone, front
			# 		elsif ($backbone == 2 and $x == $numberofbatches-1) {
			# 			my $insertno = $x+1;
			# 			                                    #"$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n"
			# 			push(@fastabatcharray, "insert $insertno\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[0][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			# 			# push(@{$fastabatcharray[$x]}, "$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			# 			$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
			# 			# push(@{$fastaoutarray[$x]},  "$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$ftcomp[$n-1]$primerextension[$x][$n-1][5]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");	
			# 			$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
			# 		}

			# 		#single backbone, back
			# 	}
				
			# 	if ($circularoutput ne "TRUE") {
			# 		unless ($backbone == 3) {
			# 			if ($x == 0) {
			# 				push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$cassette[2]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
			# 				push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\treverse\t$TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");	
			# 			}
			# 		}
					
			# 		if ($backbone == 1) {
			# 			for (my $i = 1; $i < $n-1; $i++) {
			# 				my $insertno = $x+1;
			# 				$ftcomp[$i] = reverse $ftcomp[$i];
			# 				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
			# 				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];												  
			# 				push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\tforward\t$TmGCfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i][4]$brfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]]\n");
			# 				push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
			# 				$ftcomp[$i] = reverse $ftcomp[$i];
			# 				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
			# 				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];	
			# 			}
			# 		}

			# 		if ($backbone == 2) {
			# 			for (my $i = 0; $i < $n-1; $i++) {
			# 				my $insertno = $x+1;
			# 				$ftcomp[$i] = reverse $ftcomp[$i];
			# 				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
			# 				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];												  
			# 				push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\tforward\t$TmGCfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i][4]$brfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]]\n");
			# 				push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
			# 				$ftcomp[$i] = reverse $ftcomp[$i];
			# 				$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
			# 				$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];	
			# 			}
			# 		}

			# 		if ($x == $numberofbatches-1) {
			# 			if ($backbone == 1 or $backbone == 3) {
			# 				push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			# 				push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$cassette[3]$primerextension[$x][$n-1][1]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
			# 			}
			# 		}
			# 		# elsif ($backbone == 1 and $x == $numberofbatches-1) {
			# 		# 	my $insertno = $x+1;
			# 		# 	                                    #"$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n"
			# 		# 	push(@fastabatcharray, "insert $insertno\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[0][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			# 		# 	# push(@{$fastabatcharray[$x]}, "$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
			# 		# 	$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
			# 		# 	# push(@{$fastaoutarray[$x]},  "$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$ftcomp[$n-1]$primerextension[$x][$n-1][5]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");	
			# 		# 	$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
			# 		# }
			# 	}
			# }
			
			if ($batchinput ne "TRUE") {
				push(@{$fastaoutarray[$x]},  ">FW_$inputs[$x][0][1]_[Tm:$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][0]]\n$cassette[2]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
				for (my $i = 0; $i < $n-1; $i++) {
					#$ftcomp[$i] = reverse $ftcomp[$i];
					#$brrev[$i][0][$optbr[$i][0]] = reverse $brrev[$i][0][$optbr[$i][0]];
					$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
					push(@{$fastaoutarray[$x]},  ">RV_$inputs[$x][$i][1]_[Tm:$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][0]]\n$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");	
					$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];
					push(@{$fastaoutarray[$x]},  ">FW_$inputs[$x][$i+1][1]_[Tm:$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][0]]\n$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n");
				}
				#$brrev[$n-1][0][$optbr[$n-1][0]] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
				#$cassette[3] = reverse $cassette[3];
				$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]] = reverse $brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]];
				$primerextension[$x][$n][1] = reverse $primerextension[$x][$n][1];
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				push(@{$fastaoutarray[$x]},  ">RV_$inputs[$x][$n-1][1]_[Tm:$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][0]]\n$cassette[3]$primerextension[$x][$n][1]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
				$primerextension[$x][$n][1] = reverse $primerextension[$x][$n][1];
				$primerextension[$x][$n][1] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
			
				push (@batchoutput, @fastaoutarray);
			}
		}



##### BATCH PRINTING
		if ($batchinput eq "TRUE") {
			if ($circularoutput eq "TRUE") {
				if ($backbone == 1) {
					# backbone on front and back of combinatorial input Z+XY+Z

					die "sorry, when circular assembly is selected, only input in the following format is currently available in AMUSER: \"back bone\"+\"frament 1\" \"fragment n\"+. More options to be added soon\n";

					if ($x == 0) {
						my $insertno = 1;
						push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$ft[1]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
						push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\treverse\t$TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");
						
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\tforward\t$TmGCfor[1][0][$optbr[1][1][$ranking{1}{0}]][3]\C\t$ft[0]$primerextension[$x][1][4]$brfor[1][0][$optbr[1][1][$ranking{1}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\treverse\t$TmGCrev[1][0][$optbr[1][0][$ranking{1}{0}]][3]\C\t$ftcomp[1]$primerextension[$x][1][5]$brrev[1][0][$optbr[1][0][$ranking{1}{0}]]\n");
					}
					if ($x > 0 and $x < $numberofbatches) {
						my $insertno = $x+1;
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\tforward\t$TmGCfor[1][0][$optbr[1][1][$ranking{1}{0}]][3]\C\t$ft[0]$primerextension[$x][1][4]$brfor[1][0][$optbr[1][1][$ranking{1}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\treverse\t$TmGCrev[1][0][$optbr[1][0][$ranking{1}{0}]][3]\C\t$ftcomp[1]$primerextension[$x][1][5]$brrev[1][0][$optbr[1][0][$ranking{1}{0}]]\n");
					}
					if ($x == $numberofbatches-1) {
						push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
	 					push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$cassette[3]$primerextension[$x][$n-1][1]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
					}
				}
				elsif ($backbone == 2) {
					# only front backbone is present Z+XY+

					if ($x == 0) {
						my $insertno = 1;
						push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$ft[1]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
						push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\treverse\t$TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");
						
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\tforward\t$TmGCfor[1][0][$optbr[1][1][$ranking{1}{0}]][3]\C\t$ft[0]$primerextension[$x][1][4]$brfor[1][0][$optbr[1][1][$ranking{1}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\treverse\t$TmGCrev[1][0][$optbr[1][0][$ranking{1}{0}]][3]\C\t$ftcomp[1]$primerextension[$x][1][5]$brrev[1][0][$optbr[1][0][$ranking{1}{0}]]\n");
					}
					if ($x > 0 and $x < $numberofbatches) {
						my $insertno = $x+1;
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\tforward\t$TmGCfor[1][0][$optbr[1][1][$ranking{1}{0}]][3]\C\t$ft[0]$primerextension[$x][1][4]$brfor[1][0][$optbr[1][1][$ranking{1}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\treverse\t$TmGCrev[1][0][$optbr[1][0][$ranking{1}{0}]][3]\C\t$ftcomp[1]$primerextension[$x][1][5]$brrev[1][0][$optbr[1][0][$ranking{1}{0}]]\n");
					}									
				}		
				elsif ($backbone == 3) {

					die "sorry, when circular assembly is selected, only input in the following format is currently available in AMUSER: \"back bone\"+\"frament 1\" \"fragment n\"+. More options to be added soon\n";

					# only back backbone is present +XY+Z
					if ($x >= 0 and $x < $numberofbatches) {
						my $insertno = $x+1;
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][0][1]\tforward\t$TmGCfor[1][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$ft[1]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][0][1]\treverse\t$TmGCrev[1][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");
					}
					if ($x == $numberofbatches-1) {
						push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
	 					push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$cassette[3]$primerextension[$x][$n-1][1]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
					}
				}
			}


			elsif ($circularoutput ne "TRUE") {
				if ($backbone == 1) {
					# backbone on front and back of combinatorial input Z+XY+Z
					if ($x == 0) {
						my $insertno = 1;
						push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$cassette[2]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
						push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\treverse\t$TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");
						
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\tforward\t$TmGCfor[1][0][$optbr[1][1][$ranking{1}{0}]][3]\C\t$ft[0]$primerextension[$x][1][4]$brfor[1][0][$optbr[1][1][$ranking{1}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\treverse\t$TmGCrev[1][0][$optbr[1][0][$ranking{1}{0}]][3]\C\t$ftcomp[1]$primerextension[$x][1][5]$brrev[1][0][$optbr[1][0][$ranking{1}{0}]]\n");
					}
					if ($x > 0 and $x < $numberofbatches) {
						my $insertno = $x+1;
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\tforward\t$TmGCfor[1][0][$optbr[1][1][$ranking{1}{0}]][3]\C\t$ft[0]$primerextension[$x][1][4]$brfor[1][0][$optbr[1][1][$ranking{1}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\treverse\t$TmGCrev[1][0][$optbr[1][0][$ranking{1}{0}]][3]\C\t$ftcomp[1]$primerextension[$x][1][5]$brrev[1][0][$optbr[1][0][$ranking{1}{0}]]\n");
					}
					if ($x == $numberofbatches-1) {
						push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
	 					push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$cassette[3]$primerextension[$x][$n-1][1]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
					}
				}
				elsif ($backbone == 2) {
					# only front backbone is present Z+XY+
					if ($x == 0) {
						my $insertno = 1;
						push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$cassette[2]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
						push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\treverse\t$TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");
						
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\tforward\t$TmGCfor[1][0][$optbr[1][1][$ranking{1}{0}]][3]\C\t$ft[0]$primerextension[$x][1][4]$brfor[1][0][$optbr[1][1][$ranking{1}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\treverse\t$TmGCrev[1][0][$optbr[1][0][$ranking{1}{0}]][3]\C\t$ftcomp[1]$primerextension[$x][1][5]$brrev[1][0][$optbr[1][0][$ranking{1}{0}]]\n");
					}
					if ($x > 0 and $x < $numberofbatches) {
						my $insertno = $x+1;
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\tforward\t$TmGCfor[1][0][$optbr[1][1][$ranking{1}{0}]][3]\C\t$ft[0]$primerextension[$x][1][4]$brfor[1][0][$optbr[1][1][$ranking{1}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][1][1]\treverse\t$TmGCrev[1][0][$optbr[1][0][$ranking{1}{0}]][3]\C\t$ftcomp[1]$primerextension[$x][1][5]$brrev[1][0][$optbr[1][0][$ranking{1}{0}]]\n");
					}									
				}		
				elsif ($backbone == 3) {
					# only back backbone is present +XY+Z
					if ($x >= 0 and $x < $numberofbatches) {
						my $insertno = $x+1;
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][0][1]\tforward\t$TmGCfor[1][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$ft[1]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
						push(@fastabatcharray, "insert $insertno\t$inputs[$x][0][1]\treverse\t$TmGCrev[1][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");
					}
					if ($x == $numberofbatches-1) {
						push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
	 					push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$cassette[3]$primerextension[$x][$n-1][1]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
					}
				}
			}
		}
	}
}



		# for (my $i = 1; $i < $n-1; $i++) {
		# 	my $insertno = $x+1;
		# 	$ftcomp[$i] = reverse $ftcomp[$i];
		# 	$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
		# 	$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];												  
		# 	push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\tforward\t$TmGCfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i][4]$brfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]]\n");
		# 	push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
		# 	$ftcomp[$i] = reverse $ftcomp[$i];
		# 	$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
		# 	$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];	
		# }
		
		# #double backbone
		# if ($backbone == 1 and $x == $numberofbatches-1) {
		# 	                                    #"$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n"
		# 	push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[0][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
		# 	# push(@{$fastabatcharray[$x]}, "$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
		# 	$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
		# 	# push(@{$fastaoutarray[$x]},  "$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$ftcomp[$n-1]$primerextension[$x][$n-1][5]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");	
		# 	$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
	
		# #single backbone, front
		# elsif ($backbone == 2 and $x == $numberofbatches-1) {
		# 	my $insertno = $x+1;
		# 	                                    #"$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n"
		# 	push(@fastabatcharray, "insert $insertno\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[0][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
		# 	# push(@{$fastabatcharray[$x]}, "$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
		# 	$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
		# 	# push(@{$fastaoutarray[$x]},  "$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$ftcomp[$n-1]$primerextension[$x][$n-1][5]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");	
		# 	$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
	
		# #single backbone, back
	# }
	
	# if ($circularoutput ne "TRUE") {
	# 	unless ($backbone == 3) {
	# 		if ($x == 0) {
	# 			push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\tforward\t$TmGCfor[0][0][$optbr[0][1][$ranking{0}{0}]][3]\C\t$cassette[2]$primerextension[$x][0][4]$brfor[0][0][$optbr[0][1][$ranking{0}{0}]]\n");
	# 			push(@fastabatcharray, "5' backbone\t$inputs[$x][0][1]\treverse\t$TmGCrev[0][0][$optbr[0][0][$ranking{0}{0}]][3]\C\t$ftcomp[0]$primerextension[$x][0][5]$brrev[0][0][$optbr[0][0][$ranking{0}{0}]]\n");	
	# 		}
	# 	}
		
	# 	if ($backbone == 1) {
	# 		for (my $i = 1; $i < $n-1; $i++) {
	# 			my $insertno = $x+1;
	# 			$ftcomp[$i] = reverse $ftcomp[$i];
	# 			$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
	# 			$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];												  
	# 			push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\tforward\t$TmGCfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i][4]$brfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]]\n");
	# 			push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
	# 			$ftcomp[$i] = reverse $ftcomp[$i];
	# 			$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
	# 			$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];	
	# 		}
	
	# 	if ($backbone == 2) {
	# 		for (my $i = 0; $i < $n-1; $i++) {
	# 			my $insertno = $x+1;
	# 			$ftcomp[$i] = reverse $ftcomp[$i];
	# 			$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
	# 			$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];												  
	# 			push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\tforward\t$TmGCfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i][4]$brfor[$i][0][$optbr[$i][1][$ranking{$i}{0}]]\n");
	# 			push(@fastabatcharray, "insert $insertno\t$inputs[$x][$i][1]\treverse\t$TmGCrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]][3]\C\t$ftcomp[$i]$primerextension[$x][$i][5]$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]]\n");
	# 			$ftcomp[$i] = reverse $ftcomp[$i];
	# 			$brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]] = reverse $brrev[$i][0][$optbr[$i][0][$ranking{$i}{0}]];
	# 			$primerextension[$x][$i][5] = reverse $primerextension[$x][$i][5];	
	# 		}
	
	# 	if ($x == $numberofbatches-1) {
	# 		if ($backbone == 1 or $backbone == 3) {
	# 			push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-2]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
	# 			push(@fastabatcharray, "3' backbone\t$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$cassette[3]$primerextension[$x][$n-1][1]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");
	# 		}
	# 	}
	# 	# elsif ($backbone == 1 and $x == $numberofbatches-1) {
	# 	# 	my $insertno = $x+1;
	# 	# 	                                    #"$inputs[$x][$i+1][1]\tforward\t$TmGCfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]][3]\C\t$ft[$i]$primerextension[$x][$i+1][4]$brfor[$i+1][0][$optbr[$i+1][1][$ranking{$i+1}{0}]]\n"
	# 	# 	push(@fastabatcharray, "insert $insertno\t$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[0][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
	# 	# 	# push(@{$fastabatcharray[$x]}, "$inputs[$x][$n-1][1]\tforward\t$TmGCfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]][3]\C\t$ft[$n-1]$primerextension[$x][$n-1][4]$brfor[$n-1][0][$optbr[$n-1][1][$ranking{$n-1}{0}]]\n");
	# 	# 	$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
	# 	# 	# push(@{$fastaoutarray[$x]},  "$inputs[$x][$n-1][1]\treverse\t$TmGCrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]][3]\C\t$ftcomp[$n-1]$primerextension[$x][$n-1][5]$brrev[$n-1][0][$optbr[$n-1][0][$ranking{$n-1}{0}]]\n");	
	# 	# 	$primerextension[$x][$n-1][5] = reverse $primerextension[$x][$n-1][5];
	# 	# }
	# }
# }



###################################################################
#####           HTML FORMATTING AND FILE PRINTING             #####
###################################################################


# MAKE INPUT OPTIONS FOR PRINT OUT, PRINT STDOUT, OR BOTH


# Re-decide primer dimer / intra-primer homology / polyN on the COMPLETE primers
# (USER tail + bridge + tag/linker + binding region) rather than the binding
# region alone; see &fullprimerchecks. Done after the report builders so every
# path (single, linear fusion, circular) is covered, and before any output.
{
	# same salt / primer concentrations the design used (validated there)
	my $sc_f = (exists $options{sc} and $options{sc} ne "") ? $options{sc} / 1000      : 0.05;
	my $pc_f = (exists $options{pc} and $options{pc} ne "") ? $options{pc} / 1000000   : 0.00000005;
	for (my $xx = 0; $xx < scalar @outreportarray; $xx++) {
		&fullprimerchecks(\@{$outreportarray[$xx]}, $sc_f, $pc_f);
		&move_primer_overview(\@{$outreportarray[$xx]});
	}
}

if ($batchinput eq "FALSE") {
	if ($options{r} ne "") {
		open (OUT, '>', $outreport);
		print OUTFRED "-r txt report: $outreport\n";
		for (my $i = 0; $i < scalar @{$outreportarray[0]}; $i++) {
			print OUT "$outreportarray[0][$i]";
		}
		close OUT;
	}
	
	if ($options{o} ne "") {
		open (OUT, '>', $fastaout);
		print OUTFRED "-o fsa report: $fastaout\n";
		for (my $i = 0; $i < scalar @{$fastaoutarray[0]}; $i++) {
			print OUT "$fastaoutarray[0][$i]";
		}
		close OUT;
	}	
	
	my $htmlarray_ref = &htmlprint(\@outreportarray, \@fastaoutarray, $batchinput);
	my @htmlarray = @$htmlarray_ref;
	
	if ($options{h} ne "") {
		open (OUT, '>', $htmlreport);
		print OUTFRED "-h html report: $htmlreport\n";
		for (my $i = 0; $i < scalar @{$htmlarray[0]}; $i++) {
			print OUT "$htmlarray[0][$i]";
		}
		close OUT;
		print "AMUSER report: " . basename($htmlreport) . "\n";   # report the chosen name (used by the GUI)
		# Primer list CSV (kept output, like the GenBank files): name, Tm, sequence.
		my $csvdir = $htmlreportdir;   # same "Output files" folder
		mkdir $csvdir unless -d $csvdir;
		(my $csvstem = basename($htmlreport)) =~ s/\.html?$//i;
		$csvstem =~ s/_report$//;   # share the base stem with the .gb (report suffix is HTML-only)
		&write_primer_csv($outreportarray[0], "$csvdir/${csvstem}_primers.csv");
		unless ($ENV{AMUSER_NO_OPEN}) {
			if ($^O eq 'MSWin32') {
				system(1, "cmd", "/c", "start", "", $htmlreport);
			} else {
				system("xdg-open", $htmlreport);
			}
		}
	}
	else {
		for (my $i = 0; $i < scalar @{$htmlarray[0]}; $i++) {
			print "$htmlarray[0][$i]";
		}
	}
}

elsif ($batchinput eq "TRUE") {
	for (my $x = 0; $x < scalar @inputs; $x++) {
		if ($options{r} ne "") {
			if ($options{workdir} ne "") {
				my $filename = $outputpath;
				$filename .= "batch$x.$outreportshort";
				open (OUT, '>', "$filename");
				print OUTFRED "-workdir + batch outreport: $outputpath/batch$x.$outreportshort\n";
			}
			else {
				open (OUT, '>', "batch$x.$outreportshort");
				print OUTFRED "batch outreport: batch$x.$outreportshort\n";
			}
	
			for (my $i = 0; $i < scalar @{$outreportarray[$x]}; $i++) {
				print OUT "$outreportarray[$x][$i]";
			}
			close OUT;
		}
		
		if ($options{o} ne "") {
			if ($options{workdir} ne "") {
				my $filename = $outputpath;
				$filename .= "batch$x.$fastaoutshort";
				open (OUT, '>', "$filename");
				print OUTFRED "-workdir + batch fsa: $outputpath/batch$x.$fastaoutshort\n";
			}
			else {
				open (OUT, '>', "batch$x.$fastaout");
				print OUTFRED "batch fsa: batch$x.$fastaoutshort\n";
			}
			for (my $i = 0; $i < scalar @{$fastaoutarray[$x]}; $i++) {
				print OUT "$fastaoutarray[$x][$i]";
			}
			close OUT;
		}
	}
	# for (my $p = 0; $p < scalar @fastabatcharray; $p++) {
	# 	print "\n$fastabatcharray[$p]\n";
	# }
	my $htmlarray_ref = &htmlprintbatch(\@fastaoutarray, \@fastabatcharray);
	my @htmlarray = @$htmlarray_ref;		
	
	if ($options{h} ne "") {

		open (OUT, '>', $htmlreport);
		print OUTFRED "-h html report: $htmlreport\n";
		for (my $i = 0; $i < scalar @htmlarray; $i++) {
			print OUT "$htmlarray[$i]";
		}
		close OUT;

		for (my $x = 0; $x < $numberofinputs; $x++) {
			my $htmlarray_batch_ref = &htmlprint(\@outreportarray, \@fastaoutarray, $batchinput, $x);
			my @htmlarray_batch = @$htmlarray_batch_ref;
			my $filename = $outputpath;
			$filename .= "batch$x.$htmlreportshort";
			open (OUT, '>', "$filename");
			print OUTFRED "-h + batch html: $filename\n";
			for (my $i = 0; $i < scalar @{$htmlarray_batch[$x]}; $i++) {
				print OUT "$htmlarray_batch[$x][$i]";
			}
			close OUT;
		}
	}
	else {
		for (my $i = 0; $i < scalar @htmlarray; $i++) {
			print "$htmlarray[$i]";
		}
	}	
}



###################################################################
##one last line break in debug messages for good measures
if ($verbose eq "TRUE") {
	print "\n";
} 



###################################################################
#####                       SUBROUTINES                       #####
###################################################################
				
				
sub cassettenickandrestr {

    my %options = %{$_[0]};
    my %restrictionsites = %{$_[1]};
    my %nickingsites = %{$_[2]};
    my $db1 = $_[3];
    my $db2 = $_[4];

	unless (exists $nickingsites{$options{cnick}}) {
		die "ERROR: unknown nicking input option: \"$options{cnick}\"\n";
	}
	unless (exists $restrictionsites{$options{crestr}}) {
		die "ERROR: unknown restriction input option: \"$options{crestr}\"\n";
	}
	if ($db1 eq "" or $db2 eq "" or $db1 =~ m/[^ATGC]/ or $db2 =~ m/[^ATGC]/) {
		die "ERROR: please enter valid directional bases\n";
	}
	
	my $tempfor;
	my $temprev;
	my $nickpos = 0;
	my $restrpos = 0;
	my $temp;
	
	$tempfor = $nickingsites{$options{cnick}};
	$tempfor =~ s/\.//;
	$tempfor .= "$db1";
	$tempfor .= $restrictionsites{$options{crestr}};
	$tempfor .= "$db2";
	$temp = reverse $nickingsites{$options{cnick}};
	$temp =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
	$tempfor .= $temp;
	$cassette[0] = $tempfor;
	
	$temprev = $nickingsites{$options{cnick}};
	$temprev =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
	$temprev .= "$db2";
	$temp = $restrictionsites{$options{crestr}};
	$temp = reverse $temp;
	$temprev .= $temp;
	$temprev .= "$db1";
	$temp = reverse $nickingsites{$options{cnick}};
	$temp =~ s/\.//;
	$temprev .= $temp;
	$cassette[1] = $temprev;

	while (substr($cassette[1], $nickpos, 1) ne ".") {
		$nickpos++;
	}
	
	while (substr($cassette[0], $restrpos, 1) ne ".") {
		$restrpos++;
	}
	
	$cassette[2] = substr($cassette[0], $nickpos, $restrpos-$nickpos-1);
	$cassette[2] .= "U";

	$temprev = reverse $cassette[1];

	$cassette[3] = substr($temprev, $nickpos, $restrpos-$nickpos-1);
	$cassette[3] .= "U";
	
	$cassette[4] = "$options{crestr}/$options{cnick}";

	$cassette[5] = substr($cassette[0], 0, $nickpos);
	$cassette[6] = reverse $cassette[5];
	$cassette[6] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
	
	return \@cassette;
}


sub cassettemanin {

	my %options = %{@_[0]};
	my @cassette = @{@_[1]};

	my $nickpos = 0;
	my $restrpos = 0;
	my $temprev;
	
	$cassette[0] = uc $options{cmaninfor};
	$cassette[1] = uc $options{cmaninrev};
	
	if ($cassette[0] =~ m/[^ATGC.]/g or $cassette[0] !~ m/.*\..+\./g) {
		die "ERROR: unknown custom forward cassette input option: \"$options{cmaninfor}\"\n";
	}

	if ($cassette[1] =~ m/[^ATGC.]/g or $cassette[1] !~ m/.*\..+\./g) {
		die "ERROR: unknown custom reverse cassette input option: \"$options{cmaninrev}\"\n";
	}


	if ($options{cname} ne "") {
		$cassette[4] = $options{cname};
	}
	elsif ($cassette[4] eq "") {
		$cassette[4] = "custom cassette";
	}

	while (substr($cassette[1], $nickpos, 1) ne ".") {
		$nickpos++;
	}
	
	while (substr($cassette[0], $restrpos, 1) ne ".") {
		$restrpos++;
	}
	
	$cassette[2] = substr($cassette[0], $nickpos, $restrpos-$nickpos-1);
	$cassette[2] .= "U";

	$temprev = reverse $cassette[1];

	$cassette[3] = substr($temprev, $nickpos, $restrpos-$nickpos-1);
	$cassette[3] .= "U";

	$cassette[5] = substr($cassette[0], 0, $nickpos);
	$cassette[6] = reverse $cassette[5];
	$cassette[6] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
	
	return \@cassette;
}



sub inputfile {
	# print "Entering sub inputfile\n" if ($debug eq "TRUE");
	###
	# Input file handeling
	###

	# Patch 2014.01.08 Boston, Jakob:
	my $backbone = 0; 
	# 0 =  uninitialized, 1 = front and back, 2 = only front, 3 = only back, 4 = none


	my $inputfasta = $_[0];
	my $circular_dna = $_[1];
	my $fasta_count = 0;
	my $fasta_header_in_input = 0;
	my $no = 0;
	my $batchinput = "FALSE";
	my $flag_for_varying_fasta = "FALSE";
	my $flag_long_insert = 0;
	my $end_insert;
	my @splittype;
	## splittype [fasta count][0] identifyer of which input number we came from
	## splittype [fasta count][1] identifyer: 2=not touched; 1=parted; 0=inserted
	my @input;
	my @namecount;
	my @partcount;
	my $name_i=-1;
	#my $inputcount[0] = 0;
	my $namepatch;
	my $flag_first;
	my $prev_length = 0;	
	my $cur_do_not_cross = 0;
	my $prev_do_not_cross = 0;
	my $cur_name;
	my @array_of_names;
	
	# read into variables

	open (IN, '<', $inputfasta) or die "ERROR: you have entered invalid filename(s): $!";
	my $line = <IN>;
	
	# Patch 2014.01.08 Boston, Jakob:
	my $lastline = 'void';

	if($line=~m/\+/){
		$backbone = 3;
	}else{
		$backbone = 1;
	}

	# LOOP
	while (defined $line){
		$line =~ tr/\r/\n/; # translitterate makes carrrige return to newlines
		chomp($line);
		# &sanity_check_line($line);
		# Skip empty lines in input
		if ($line =~ m/^\s*$/){
			# Patch 2014.01.08 Boston, Jakob:
			$lastline = $line;

			$line = <IN>;

		# save name if line starts with >
		} elsif ( $line =~ m/^>(.*)/ ){
			push (@array_of_names,$1);
			$flag_first = 1;
			$name_i++;
			$inputcount[0]++;
			# "$no" is a number to determine if input is a part of combinatorial input.
			$cur_name = $1;
			if ($flag_for_varying_fasta eq "FALSE"){
				$no=0;
			} else {
				$no++;
				$cur_do_not_cross++;
			}
			$input[$fasta_count][0] = $no;
			$input[$fasta_count][1] = $1;

			# Patch 2014.01.08 Boston, Jakob:
			$lastline = $line;

			# Raise flag if input contains fasta header
			$fasta_header_in_input = 1;

			$line = <IN>;

		# change state of flag if line starts with +
		} elsif ( $line =~ m/^\+(.*)$/ ){
			$cur_do_not_cross++;
			$batchinput = "TRUE";
			if ($flag_for_varying_fasta eq "FALSE") {
				$flag_for_varying_fasta = "TRUE";
			} else {
				$flag_for_varying_fasta = "FALSE";
				$no=0;
			}
			# Patch 2014.01.08 Boston, Jakob:
			$lastline = $line;

			$line = <IN>;
		}

		# if neither line starts with > nor + DNA in fasta is assumed, 
		### Now we check if the DNA should be modified, that is the case if we have squared brackets present.
		else {
			# all DNA will be concatenated into $string
			my $string;
			while (defined $line and ($line =~ m/^[>\+]/) == 0 ){
				$line =~ tr/\r/\n/; # translitterate makes carrrige return to newlines
				chomp($line);
				$line =~ s/\s//g;			#removing all whitespace
				$string .= uc($line);		#converting DNA to uppercase

				# Patch 2014.01.08 Boston, Jakob:
				$lastline = $line;

				$line = <IN>;
			}

			if ($fasta_header_in_input == 0 and $batchinput eq "FALSE") {
				die "Error: no fasta input detected. AMUSER accepts only fasta inputs. Please check input\n";
			}
		#	my $fragment_no=0;
		#	my $insertion_no=0;
		#	&sanity_check_string($string);
	
		# check id the DNA contains an equal sign: = 
		# in this case, DNA must be split up according to rules
			if ($string =~ m/=/){
				while ($string =~m/^(.*?)\[(.*?)=(.*?)\](.*)$/){
					my ($dollar1,$dollar2,$dollar3,$dollar4) = ($1,$2,$3,$4);
					my ($length1,$length2,$length3,$length4) = (length($dollar1),length($dollar2),length($dollar3),length($dollar4));
					
					#Pseudo Algorithm

					#a + may not have been crossed
					# if varyig fasta flag is up, a > must not have been crossed
					#we are not supposed to merge strings to inserts
					if ($length1 > 0){
						if ($prev_do_not_cross == $cur_do_not_cross and $fasta_count > 0 and ($length1 + $prev_length) < 80){
						#we are okay to merge strings of short length to inserts
							$input[$fasta_count-1][2] .= $dollar1; #string
						} else {
							if ($length1 < 80){
								$splittype[$fasta_count][0] = 0; # type 0 = insert
							}else{
								$splittype[$fasta_count][0] = 1; # type 1 = part x of y
								if (defined $partcount[$name_i]){
									$partcount[$name_i] += 1;
								} else {
									$partcount[$name_i] = 1;
								}
							}
							$input[$fasta_count][0] = $no; 
							$input[$fasta_count][1] = $cur_name;
							$input[$fasta_count][2] = $dollar1; #string
							$splittype[$fasta_count][1] = $name_i;
							$fasta_count++;
						}
						$string = "[=$dollar3]$dollar4";
					}elsif($length3 > 0){
						if ($prev_do_not_cross == $cur_do_not_cross and $fasta_count > 0 and ($length3 + $prev_length) < 80){
						#we are okay to merge strings of short length to inserts
							$input[$fasta_count-1][2] .= $dollar3; #string
						} else {
							if ($length3 < 80){
								$splittype[$fasta_count][0] = 0; # type 0 = insert
							}else{
								#WARNING LONG INSERT
								$flag_long_insert = 1;
								$splittype[$fasta_count][0] = 1; # type 1 = part x of y
							}
							$input[$fasta_count][0] = $no; 
							$input[$fasta_count][1] = $cur_name;
							$input[$fasta_count][2] = $dollar3; #string
							$splittype[$fasta_count][1] = $name_i;
							$fasta_count++;
						}
						$string = $dollar4;
					}else{
						$string = $dollar4;
					}
				

					# if previous string is short (<80) and this one is as well, they can be merged as an insert

					# if length1 != 0 dollar1 is saved
					# dollar3 and 4 are put into string
					# else if length3 !=0 length3 is saved
					#dollar4 is put into string.
					
					$flag_first = 0;
					$prev_length = length($input[$fasta_count-1][2]);
					$prev_do_not_cross = $cur_do_not_cross;

				}
				if (length($string)>0){
					if ($prev_do_not_cross == $cur_do_not_cross and $fasta_count > 0 and (length($string) + $prev_length) < 80){
						#we are okay to merge strings of short length to inserts
						$input[$fasta_count-1][2] .= $string; #string
					} else {
						if (length($string) < 80){
							$splittype[$fasta_count][0] = 0; # type 0 = insert
						}else{
							$splittype[$fasta_count][0] = 1; # type 1 = part x of y
							if (defined $partcount[$name_i]){
								$partcount[$name_i] += 1;
							} else {
								$partcount[$name_i] = 1;
							}
						}
						$input[$fasta_count][0] = $no; 
						$input[$fasta_count][1] = $cur_name;
						$input[$fasta_count][2] = $string; #string
						$splittype[$fasta_count][1] = $name_i;
						$fasta_count++;
					}
					if (defined $partcount[$name_i]){
						$partcount[$name_i] += 1;
					} else {
						$partcount[$name_i] = 1;
					}

					$flag_first = 0;
					$prev_length = length($input[$fasta_count-1][2]);
					$prev_do_not_cross = $cur_do_not_cross;
				}
	
			} 

			else {
	
				if (length($string)>0){
					if ($prev_do_not_cross == $cur_do_not_cross and $fasta_count > 0 and (length($string) + $prev_length) < 80){
						#we are okay to merge strings of short length to inserts
						$input[$fasta_count-1][2] .= $string; #string
					} else {
						if (length($string) < 80){
							$splittype[$fasta_count][0] = 0; # type 0 = insert
						}else{
							$splittype[$fasta_count][0] = 2; # type 2 = not modyfied
						}
						$input[$fasta_count][0] = $no; 
						$input[$fasta_count][1] = $cur_name;
						$input[$fasta_count][2] = $string; #string
						$splittype[$fasta_count][1] = $name_i;

						$fasta_count++;
					}
					if (defined $partcount[$name_i]){
						$partcount[$name_i] += 1;
					} 

					else {
						$partcount[$name_i] = 1;
					}

					$flag_first = 0;
					$prev_length = length($input[$fasta_count-1][2]);
					$prev_do_not_cross = $cur_do_not_cross;
				}

				if (defined $namecount[$name_i]){
					$namecount[$name_i] += 1;
				} 

				else {
					$namecount[$name_i] = 1;
				}
			}
		#	$prev_length = length($input[$fasta_count-1][2]);
		#	$prev_do_not_cross = $cur_do_not_cross;
	
		}
		# if line starts with > elsif +
		# define approperiate arrays/strings	
	}
	
	# Patch 2014.01.08 Boston, Jakob:
	if($lastline=~m/\+/){
		$backbone++;
	}

	close IN; 
	
	
	#sanity final check
	
	# if ($sanity_batch_odd == 1){
	# 	die "odd number of plusses";
	# }

	### Disabled on request from Lars
	#if ($circular_dna eq "TRUE"){
	#	$fasta_count++;
	#	$name_i++;
	#	$input[$fasta_count][0] = 0;
	#	$input[$fasta_count][1] = "Virtual circular DNA (based on: $input[0][1])";
	#	$input[$fasta_count][2] = $input[0][2];
	#	$splittype[$fasta_count][0] = 3;
	#	$splittype[$fasta_count][1] =  $name_i;
	#} 
	
	
	#names
	my $j=1;
	my $prev_input=0;
	for (my $i=0;$i<$fasta_count;$i++){
		#loop over all entries
		# If they are of type fragment, add  part x of y if more than 1 part
		if ($splittype[$i][1]!=$prev_input){
			$prev_input++;
			$j=1;
		}
	
		# add the name
		if($splittype[$i][0]==1 and $partcount[$prev_input]!=1){
			my $tempcount = $partcount[$prev_input]-1;
			$input[$i][1] .= " - part $j of $tempcount";
			$j++;
		}elsif($splittype[$i][0]==0){
			$input[$i][1] .= " - insert";
		}
	
	#	#flag for long fasta insert
	#	if ($splittype[$i][0]==0 or $splittype[$i][0]==0){
	#		if (length($input[$i][2])>80){
	#			$flag_long_insert=1;
	#		}
	#	}

	}
	

	# print "flag_long_insert = $flag_long_insert\n" if ($verbose eq "TRUE");





###JAKOB option
#	 # printing to test of program is working
#	 if ($verbose eq "TRUE"){
#	 	print "test\n";
#	 	for (my $x = 0; $x <= $#input; $x++){
#	 		for (my $y = 0; $y <= $#{$input[$x]}; $y++){
#	 			print "x=$x y=$y ",$input[$x][$y],"\t";
#	 		}
#	 		print"$splittype[$x][1]\n\n";
#	 	}
#	 }
#	 print "tesr2\n";
#	 print(@array_of_names);
#	 print"\n";
#
#	 print "\$backbone_around_combinatorial =  $backbone_around_combinatorial\n";


###LARS option
	return (\@input, $batchinput, \@inputcount, $backbone);


}

# sub inputfile {
# 	# print "Entering sub inputfile\n" if ($debug eq "TRUE");
# 	###
# 	# Input file handeling
# 	###
	
# 	my $inputfasta = $_[0];
# 	my $circular_dna = $_[1];
# 	my $fasta_count = 0;
# 	my $no = 0;
# 	my $batchinput = "FALSE";
# 	my $flag_for_varying_fasta = "FALSE";
# 	my $flag_long_insert = 0;
# 	my $end_insert;
# 	my @splittype;
# 	## splittype [fasta count][0] identifyer of which input number we came from
# 	## splittype [fasta count][1] identifyer: 2=not touched; 1=parted; 0=inserted
# 	my @input;
# 	my @namecount;
# 	my @partcount;
# 	my $name_i=-1;
# 	#my $inputcount[0] = 0;
# 	my $namepatch;
# 	my $flag_first;
# 	my $prev_length = 0;	
# 	my $cur_do_not_cross = 0;
# 	my $prev_do_not_cross = 0;
# 	my $cur_name;
# 	my @array_of_names;
	
# 	# read into variables

# 	open (IN, '<', $inputfasta) or die "ERROR: you have entered invalid filename(s): $!";
# 	my $line = <IN>;
	
# 	# LOOP
# 	while (defined $line){
# 		$line =~ tr/\r/\n/; # translitterate makes carrrige return to newlines
# 		chomp($line);
# 		# &sanity_check_line($line);
# 		# Skip empty lines in input
# 		if ($line =~ m/^\s*$/){
# 			$line = <IN>;
# 		# save name if line starts with >
# 		} elsif ( $line =~ m/^>(.*)/ ){
# 			push (@array_of_names,$1);
# 			$flag_first = 1;
# 			$name_i++;
# 			$inputcount[0]++;
# 			# "$no" is a number to determine if input is a part of combinatorial input.
# 			$cur_name = $1;
# 			if ($flag_for_varying_fasta eq "FALSE"){
# 				$no=0;
# 			} else {
# 				$no++;
# 				$cur_do_not_cross++;
# 			}
# 			$input[$fasta_count][0] = $no;
# 			$input[$fasta_count][1] = $1;

# 			$line = <IN>;
# 		# change state of flag if line starts with +
# 		} elsif ( $line =~ m/^\+(.*)$/ ){
# 			$cur_do_not_cross++;
# 			$batchinput = "TRUE";
# 			if ($flag_for_varying_fasta eq "FALSE") {
# 				$flag_for_varying_fasta = "TRUE";
# 			} else {
# 				$flag_for_varying_fasta = "FALSE";
# 				$no=0;
# 			}
# 			$line = <IN>;
# 		} 
# 		# if neither line starts with > nor + DNA in fasta is assumed, 
# 		### Now we check if the DNA should be modified, that is the case if we have squared brackets present.
# 		else {
# 			# all DNA will be concatenated into $string
# 			my $string;
# 			while (defined $line and ($line =~ m/^[>\+]/) == 0 ){
# 				$line =~ tr/\r/\n/; # translitterate makes carrrige return to newlines
# 				chomp($line);
# 				$line =~ s/\s//g;			#removing all whitespace
# 				$string .= uc($line);		#converting DNA to uppercase
# 				$line = <IN>;
# 			}
# 		#	my $fragment_no=0;
# 		#	my $insertion_no=0;
# 		#	&sanity_check_string($string);
	
# 		# check id the DNA contains an equal sign: = 
# 		# in this case, DNA must be split up according to rules
# 			if ($string =~ m/=/){
# 				while ($string =~m/^(.*?)\[(.*?)=(.*?)\](.*)$/){
# 					my ($dollar1,$dollar2,$dollar3,$dollar4) = ($1,$2,$3,$4);
# 					my ($length1,$length2,$length3,$length4) = (length($dollar1),length($dollar2),length($dollar3),length($dollar4));
					
# 					#Pseudo Algorithm

# 					#a + may not have been crossed
# 					# if varyig fasta flag is up, a > must not have been crossed
# 					#we are not supposed to merge strings to inserts
# 					if ($length1 > 0){
# 						if ($prev_do_not_cross == $cur_do_not_cross and $fasta_count > 0 and ($length1 + $prev_length) < 80){
# 						#we are okay to merge strings of short length to inserts
# 							$input[$fasta_count-1][2] .= $dollar1; #string
# 						} else {
# 							if ($length1 < 80){
# 								$splittype[$fasta_count][0] = 0; # type 0 = insert
# 							}else{
# 								$splittype[$fasta_count][0] = 1; # type 1 = part x of y
# 								if (defined $partcount[$name_i]){
# 									$partcount[$name_i] += 1;
# 								} else {
# 									$partcount[$name_i] = 1;
# 								}
# 							}
# 							$input[$fasta_count][0] = $no; 
# 							$input[$fasta_count][1] = $cur_name;
# 							$input[$fasta_count][2] = $dollar1; #string
# 							$splittype[$fasta_count][1] = $name_i;
# 							$fasta_count++;
# 						}
# 						$string = "[=$dollar3]$dollar4";
# 					}elsif($length3 > 0){
# 						if ($prev_do_not_cross == $cur_do_not_cross and $fasta_count > 0 and ($length3 + $prev_length) < 80){
# 						#we are okay to merge strings of short length to inserts
# 							$input[$fasta_count-1][2] .= $dollar3; #string
# 						} else {
# 							if ($length3 < 80){
# 								$splittype[$fasta_count][0] = 0; # type 0 = insert
# 							}else{
# 								#WARNING LONG INSERT
# 								$flag_long_insert = 1;
# 								$splittype[$fasta_count][0] = 1; # type 1 = part x of y
# 							}
# 							$input[$fasta_count][0] = $no; 
# 							$input[$fasta_count][1] = $cur_name;
# 							$input[$fasta_count][2] = $dollar3; #string
# 							$splittype[$fasta_count][1] = $name_i;
# 							$fasta_count++;
# 						}
# 						$string = $dollar4;
# 					}else{
# 						$string = $dollar4;
# 					}
				

# 					# if previous string is short (<80) and this one is as well, they can be merged as an insert

# 					# if length1 != 0 dollar1 is saved
# 					# dollar3 and 4 are put into string
# 					# else if length3 !=0 length3 is saved
# 					#dollar4 is put into string.
					
# 					$flag_first = 0;
# 					$prev_length = length($input[$fasta_count-1][2]);
# 					$prev_do_not_cross = $cur_do_not_cross;

# 				}
# 				if (length($string)>0){
# 					if ($prev_do_not_cross == $cur_do_not_cross and $fasta_count > 0 and (length($string) + $prev_length) < 80){
# 						#we are okay to merge strings of short length to inserts
# 						$input[$fasta_count-1][2] .= $string; #string
# 					} else {
# 						if (length($string) < 80){
# 							$splittype[$fasta_count][0] = 0; # type 0 = insert
# 						}else{
# 							$splittype[$fasta_count][0] = 1; # type 1 = part x of y
# 							if (defined $partcount[$name_i]){
# 								$partcount[$name_i] += 1;
# 							} else {
# 								$partcount[$name_i] = 1;
# 							}
# 						}
# 						$input[$fasta_count][0] = $no; 
# 						$input[$fasta_count][1] = $cur_name;
# 						$input[$fasta_count][2] = $string; #string
# 						$splittype[$fasta_count][1] = $name_i;
# 						$fasta_count++;
# 					}
# 					if (defined $partcount[$name_i]){
# 						$partcount[$name_i] += 1;
# 					} else {
# 						$partcount[$name_i] = 1;
# 					}

# 					$flag_first = 0;
# 					$prev_length = length($input[$fasta_count-1][2]);
# 					$prev_do_not_cross = $cur_do_not_cross;
# 				}
	
# 			} 

# 			else {
	
# 				if (length($string)>0){
# 					if ($prev_do_not_cross == $cur_do_not_cross and $fasta_count > 0 and (length($string) + $prev_length) < 80){
# 						#we are okay to merge strings of short length to inserts
# 						$input[$fasta_count-1][2] .= $string; #string
# 					} else {
# 						if (length($string) < 80){
# 							$splittype[$fasta_count][0] = 0; # type 0 = insert
# 						}else{
# 							$splittype[$fasta_count][0] = 2; # type 2 = not modyfied
# 						}
# 						$input[$fasta_count][0] = $no; 
# 						$input[$fasta_count][1] = $cur_name;
# 						$input[$fasta_count][2] = $string; #string
# 						$splittype[$fasta_count][1] = $name_i;

# 						$fasta_count++;
# 					}
# 					if (defined $partcount[$name_i]){
# 						$partcount[$name_i] += 1;
# 					} 

# 					else {
# 						$partcount[$name_i] = 1;
# 					}

# 					$flag_first = 0;
# 					$prev_length = length($input[$fasta_count-1][2]);
# 					$prev_do_not_cross = $cur_do_not_cross;
# 				}

# 				if (defined $namecount[$name_i]){
# 					$namecount[$name_i] += 1;
# 				} 

# 				else {
# 					$namecount[$name_i] = 1;
# 				}
# 			}
		
# 	#	$prev_length = length($input[$fasta_count-1][2]);
# 	#	$prev_do_not_cross = $cur_do_not_cross;
	
# 		}
# 		# if line starts with > elsif +
# 		# define approperiate arrays/strings	

# 	}
# 	close IN; 
	
	
# 	#sanity final check
	
# 	# if ($sanity_batch_odd == 1){
# 	# 	die "odd number of plusses";
# 	# }

# 	### Disabled on request from Lars
# 	#if ($circular_dna eq "TRUE"){
# 	#	$fasta_count++;
# 	#	$name_i++;
# 	#	$input[$fasta_count][0] = 0;
# 	#	$input[$fasta_count][1] = "Virtual circular DNA (based on: $input[0][1])";
# 	#	$input[$fasta_count][2] = $input[0][2];
# 	#	$splittype[$fasta_count][0] = 3;
# 	#	$splittype[$fasta_count][1] =  $name_i;
# 	#} 
	
	
# 	#names
# 	my $j=1;
# 	my $prev_input=0;
# 	for (my $i=0;$i<$fasta_count;$i++){
# 		#loop over all entries
# 		# If they are of type fragment, add  part x of y if more than 1 part
# 		if ($splittype[$i][1]!=$prev_input){
# 			$prev_input++;
# 			$j=1;
# 		}
	
# 		# add the name
# 		if($splittype[$i][0]==1 and $partcount[$prev_input]!=1){
# 			my $tempcount = $partcount[$prev_input]-1;
# 			$input[$i][1] .= " - part $j of $tempcount";
# 			$j++;
# 		}elsif($splittype[$i][0]==0){
# 			$input[$i][1] .= " - insert";
# 		}
	
# 	#	#flag for long fasta insert
# 	#	if ($splittype[$i][0]==0 or $splittype[$i][0]==0){
# 	#		if (length($input[$i][2])>80){
# 	#			$flag_long_insert=1;
# 	#		}
# 	#	}

# 	}
	
	
# 	# print "flag_long_insert = $flag_long_insert\n" if ($verbose eq "TRUE");





# ###JAKOB option
# 	# # printing to test of program is working
# 	# if ($verbose eq "TRUE"){
# 	# 	print "test\n";
# 	# 	for (my $x = 0; $x <= $#input; $x++){
# 	# 		for (my $y = 0; $y <= $#{$input[$x]}; $y++){
# 	# 			print "x=$x y=$y ",$input[$x][$y],"\t";
# 	# 		}
# 	# 		print"$splittype[$x][1]\n\n";
# 	# 	}
# 	# }

# 	# print "tesr2\n";
# 	# print(@array_of_names);
# 	# print"\n";


# ###LARS option
# 	return (\@input, $batchinput, \@inputcount);


# }
			
		
				
sub batchinputconversion {
	my @input = @{$_[0]};
	my @inputcount;
	# print "number of input sequences ", scalar @input, "\n";
	my @inputs;
	my $numberofbatches = 0;
	for (my $i = 0; $i < scalar @input; $i++) {
		if ($input[$i][0] > $numberofbatches) {
			$numberofbatches = $input[$i][0];
		}
	}
	for (my $i = 1; $i <= $numberofbatches; $i++) {
		$inputcount[$i-1] = 0;
		my @temparray;
		for (my $ii = 0; $ii < scalar @input; $ii++) {
			if ($input[$ii][0] == 0 or $input[$ii][0] == $i) {
				$temparray[$inputcount[$i-1]][0] = $input[$ii][0];
				$temparray[$inputcount[$i-1]][1] = $input[$ii][1];
				$temparray[$inputcount[$i-1]][2] = $input[$ii][2];
				$inputcount[$i-1]++;
			}
		}
		$inputs[$i-1] = [@temparray];
	}
	# for (my $i = 0; $i <= $numberofbatches; $i++) {
	# 	#for (my $ii = 0; $ii < scalar @{$input[$i]}; $ii++) {
	# 		print "seq $i $input[$i][0]\n$input[$i][1]\n\n";
	# 	#}
	# }
	# print "number of batches: $numberofbatches\n";
	return ($numberofbatches, \@inputs, \@inputcount);
}



sub primerextension {
	my @inputs = @{$_[0]};
	my $circularoutput = $_[1];
	my @primerextension;
	my $numberofinserts = 0;
	# print scalar @{$inputs[0]}, "\n";
	for (my $i = 0; $i < scalar @inputs; $i++) {
		for (my $ii = 0; $ii < scalar @{$inputs[$i]}; $ii++) {
			if (length $inputs[$i][$ii][2] < 81) {
				# print "$ii $inputs[$i][$ii][2]\n";
				if ($ii == (scalar @{$inputs[$i]})-1) {
					my $shortseq = $inputs[$i][$ii][2];
					$primerextension[$i][$ii][0] = "TRUE";
					$primerextension[$i][$ii][1] = $shortseq;
					$numberofinserts++;
					splice (@{$inputs[$i]}, $ii, 1);
				}
				else {
					my $shortseq = $inputs[$i][$ii][2];
					my $startlength;
					my $endlength;
					if ($shortseq % 2 == 0) {
						$startlength = length ($shortseq)/2;
						$endlength = length ($shortseq)/2;
					} 
					else {
						$startlength = (length ($shortseq)/2)+0.5;
						$endlength = (length ($shortseq)/2)-0.5;
					}
					$primerextension[$i][$ii][0] = "TRUE";
					$primerextension[$i][$ii][1] = $shortseq;
					$primerextension[$i][$ii][2] = $startlength;
					$primerextension[$i][$ii][3] = $endlength;
					# print "$primerextension[$i][$ii][1]\n$primerextension[$i][$ii][2]\n$primerextension[$i][$ii][3]\n";
					$numberofinserts++;
					splice (@{$inputs[$i]}, $ii, 1);
				}
			}
		}
		my $n = scalar @{$inputs[$i]};
		if ($circularoutput eq "TRUE") {
			if ($primerextension[$i][0][0] eq "TRUE" and $primerextension[$i][$n][0] eq "TRUE") {
				my $tempseq = $primerextension[$i][$n][1];
				$tempseq .= $primerextension[$i][0][1];
				$primerextension[$i][0][1] = $tempseq;
				if ((length $primerextension[$i][0][1]) > 80) {
					die "ERROR: two insertions are attempted to be introduced next to each other in circular construct. The combined length of these exceed 80 bp, which may cause primer instability. Please check input\n";
				}

				if ($primerextension[$i][0][1] % 2 == 0) {
					$primerextension[$i][0][2] = length ($primerextension[$i][0][1])/2;
					$primerextension[$i][0][3] = length ($primerextension[$i][0][1])/2;
				} 
				else {
					$primerextension[$i][0][2] = (length ($primerextension[$i][0][1])/2)+0.5;
					$primerextension[$i][0][3] = (length ($primerextension[$i][0][1])/2)-0.5;
				}

				$primerextension[$i][$n][0] = "";
				$primerextension[$i][$n][1] = "";
				$primerextension[$i][$n][2] = "";
				$primerextension[$i][$n][3] = "";
			}
			if ($primerextension[$i][0][0] eq "TRUE") {
				$primerextension[$i][$n][0] = $primerextension[$i][0][0];
				$primerextension[$i][$n][1] = $primerextension[$i][0][1];
				$primerextension[$i][$n][2] = $primerextension[$i][0][2];
				$primerextension[$i][$n][3] = $primerextension[$i][0][3];
				$primerextension[$i][0][0] = "";
				$primerextension[$i][0][1] = "";
				$primerextension[$i][0][2] = "";
				$primerextension[$i][0][3] = "";
			}	
		}
	}

	#print scalar @{$inputs[0]}, "\n";
	return (\@inputs, \@primerextension, $numberofinserts);

	# for (my $i = 0; $i < scalar @inputs; $i++) {
	# 	for (my $ii = 0; $ii < scalar @{$inputs[$i]}; $ii++) {
	# 		if ($primerextension[$i][$ii][1] ne "TRUE") {
	# 			$tempinputs[$i][$ii] = [$inputs[$i][$ii]];
	# 		}
	# 	}
	# }

	# for (my $i = 0; $i < scalar @inputs; $i++) {
	# 	for (my $ii = 0; $ii < scalar @{$inputs[$i]}; $ii++) {
	# 		if ($ii == 0 or $ii == (scalar @{$inputs[$i]})-1) {
	# 			#oh, what to do if fragment is first or last...
	# 		}
	# 		else {
	# 			if (length $inputs[$i][$ii][2] < 81) {
	# 				my $shortseq = $inputs[$i][$ii][2];
	# 				splice (@{$inputs[$i]}, $ii, 1);
	# 				my $frontlength;
	# 				my $backlength;
	# 				if (length $shortseq % 2) {
	# 					$frontlength = length ($shortseq)/2;
	# 					$backlength = length ($shortseq)/2;
	# 				} 
	# 				else {
	# 					$frontlength = (length ($shortseq)/2)+0.5;
	# 					$backlength = (length ($shortseq)/2)-0.5;
	# 				}
	# 				my $front = substr($shortseq, 0, $frontlength);
	# 				my $back = substr($shortseq, $frontlength, $backlength);
	# 				$inputs[$i][$ii-1][2] .= $front;
	# 				$back .= $inputs[$i][$ii][2];
	# 				$inputs[$i][$ii][2] = $back; 
	# 				$primerextension[$i][$ii-1][1] = $frontlength;
	# 				$primerextension[$i][$ii][0] = $backlength;
	# 			}
	# 		}
	# 	}
	# }
}

# sub linkerexpand {
# 	for (my $i = 0; $i < $n; $i++) {
# 		unless ($linker{$i}{0} eq "") {
# 			my $templength = length $linker{$i}{0};
# 			if ($templength%2 == 1) {
# 				$linker{$i}{1} = substr ($linker{$i}{0}, 0, ($templength / 2)-0.5);
# 				$linker{$i}{2} = ($templength / 2)-0.5;
# 				$linker{$i}{3} = substr ($linker{$i}{0}, ($templength / 2), ($templength / 2)+0.5);
# 				$linker{$i}{4} = ($templength / 2)+0.5;
# 			}
# 			else {
# 				$linker{$i}{1} = substr ($linker{$i}{0}, 0, ($templength / 2));
# 				$linker{$i}{2} = ($templength / 2);
# 				$linker{$i}{3} = substr ($linker{$i}{0}, ($templength / 2), ($templength / 2));
# 				$linker{$i}{4} = ($templength / 2);
# 			}
# 		}
# 	}
# 	return %linker;
# }



sub circular {
	my @inputs = @{$_[0]};
	my $n = $_[1];
	push (@inputs, $inputs[0]);
	return (\@inputs);
}



sub brickseqfrags {
	my @input = @{$_[0]};
	my $n = $_[1];
	my @primerextension = @{$_[2]};
	my $circularoutput = $_[3];
	my @brickseqfrags;

	if ($primerextension[0][0] eq "TRUE") {
		my $frag1 = $primerextension[0][1];
		$frag1 .= substr($input[0][2], 0, 50);
		my $frag2 = substr($input[0][2], (length $input[0][2])-50, 50);
		@brickseqfrags[0] = ([$frag1,$frag2]);
	}
	else {
		my $frag1 = substr($input[0][2], 0, 50);
		my $frag2 = substr($input[0][2], (length $input[0][2])-50, 50);
		@brickseqfrags[0] = ([$frag1,$frag2]);
	}

	for (my $i = 1; $i < $n; $i++) {
		if ($primerextension[$i][0] eq "TRUE") {
			my $frag1 = $primerextension[$i][1];
			$frag1 .= substr($input[$i][2], 0, 50);
			my $frag2 = substr($input[$i][2], (length $input[$i][2])-50, 50);
			@brickseqfrags[$i] = ([$frag1,$frag2]);
		}
		else {
			my $frag1 = substr($input[$i][2], 0, 50);
			my $frag2 = substr($input[$i][2], (length $input[$i][2])-50, 50);
			@brickseqfrags[$i] = ([$frag1,$frag2]);
		}
	}

	# if ($primerextension[$n-1][0] eq "TRUE") {
	# 	my $tempseq = $primerextension[$n-1][1];
	# 	$tempseq .= $brickseqfrags[$n-1][1];
	# 	print "$tempseq\n";
	# 	$brickseqfrags[$n-1][1] = $tempseq;
	# }

	# if ($circularoutput eq "TRUE") {

	# 	if ($primerextension[0][0] eq "TRUE") {
	# 		my $tempseq = $primerextension[0][1];
	# 		$tempseq .= $brickseqfrags[$n-1][1];
	# 		$brickseqfrags[$n-1][0] = $tempseq;
	# 	}	
	# }

	# for (my $i = 0; $i < $n; $i++) {
	# 	# print "$i $input[$i][2]\n";
	# 	# print "$i $primerextension[$i][1]\n";
	# 	# if ($i == 0) {
	# 	# 	if ($primerextension[$i][0] eq "TRUE") {
	# 	# 		# Oh, what to do if insert is first...
	# 	# 	}
	# 	# 	else {
	# 	# 		my $frag1 = substr($input[$i][2], 0, 40);
	# 	# 		my $frag2 = substr($input[$i][2], (length $input[$i][2])-40, 40);
	# 	# 		@brickseqfrags[$i] = ([$frag1,$frag2]);
	# 	# 	}
	# 	# }
	# 	# elsif ($i > 0 and $i < $n) {
	# 	if ($i > 0) {
	# 		if ($primerextension[$i][0] eq "TRUE") {
	# 			my $frag0 = substr($primerextension[$i][1], 0, $primerextension[$i][2]);
	# 			$brickseqfrags[$i-1][1] .= $frag0;
	# 			# print "$brickseqfrags[$i-1][1]\n";
	# 			my $frag1 = substr($primerextension[$i][1], $primerextension[$i][2], $primerextension[$i][3]);
	# 			$frag1 .= substr($input[$i][2], 0, 40);
	# 			# print "$frag2\n";
	# 			my $frag2 = substr($input[$i][2], (length $input[$i][2])-40, 40);
	# 			@brickseqfrags[$i] = ([$frag1,$frag2]);
	# 			print "$i\n$frag1\n$frag2\n";
	# 		}
	# 		else {
	# 			my $frag1 = substr($input[$i][2], 0, 40);
	# 			my $frag2 = substr($input[$i][2], (length $input[$i][2])-40, 40);
	# 			@brickseqfrags[$i] = ([$frag1,$frag2]);
	# 			#print "$i\n$frag1\n$frag2\n";
	# 		}
	# 	}
	# 	else {
	# 		my $frag1 = substr($input[$i][2], 0, 40);
	# 		my $frag2 = substr($input[$i][2], (length $input[$i][2])-40, 40);
	# 		@brickseqfrags[$i] = ([$frag1,$frag2]);
	# 		#print "$i\n$frag1\n$frag2\n";
	# 	}

	# 	# }
	# 	# elsif ($i == $n-1) {
	# 	# 	if ($primerextension[$i][0] eq "TRUE") {
	# 	# 	# Oh, what to do if insert is last...
	# 	# 	}
	# 	# 	else {
	# 	# 		my $frag1 = substr($input[$i][2], 0, 40);
	# 	# 		my $frag2 = substr($input[$i][2], (length $input[$i][2])-40, 40);
	# 	# 		@brickseqfrags[$i] = ([$frag1,$frag2]);
	# 	# 	}
	# 	# }
	# }
	return \@brickseqfrags;
}



sub joined {
	my @brickseqfrags = @{$_[0]};
	my $n = $_[1];
    my @joined;

	for (my $i = 0; $i < $n-1; $i++) {	
		my $join = "";	
		$join .= $brickseqfrags[$i][1];
		$join .= $brickseqfrags[$i+1][0];
		$joined[$i] = $join;
	}
	return \@joined;
}



sub ATpos {

	my @joined = @{$_[0]};
	my $n = $_[1];
	my $batchinput = $_[2];
	my %forced_overhang = %{$_[3]} if defined $_[3];
	my @A;
	my @T;


	###################################################################
	#All positions of A's and T's for each joining region are stored
	#in the arrays @A and @T.

	# print "$n\n";
	for (my $i = 0; $i < $n-1; $i++) {
		my $Apos = 0;
		for (my $ii = 0; $ii < (length $joined[$i]); $ii++) {
			if (substr($joined[$i], $ii, 1) eq "A") {
				$A[$i][$Apos] = $ii;
				$Apos++;
				# print "$ii\t";
			}
		}
		# print "\n";
		# print "A: $Apos\n";
		my $Tpos = 0;
		for (my $ii = 0; $ii < (length $joined[$i]); $ii++) {
			if (substr($joined[$i], $ii, 1) eq "T") {
				$T[$i][$Tpos] = $ii+1;
				$Tpos++;
			}
		}
		# print "T: $Tpos\n";
	}


	###################################################################
	#For forced-overhang junctions, restrict A/T positions to the
	#relevant 24-nt window on the forced side of the boundary (pos 50).
	#"A": span must lie entirely in positions 26-49 (last 24 nt of A half).
	#"B": span must lie entirely in positions 50-73 (first 24 nt of B half).

	for (my $i = 0; $i < $n-1; $i++) {
		if (exists $forced_overhang{$i}) {
			if ($forced_overhang{$i} eq "A") {
				@{$A[$i]} = grep { $_ >= 26           } @{$A[$i]};
				@{$T[$i]} = grep { $_ <= 50           } @{$T[$i]};
			}
			elsif ($forced_overhang{$i} eq "B") {
				@{$A[$i]} = grep { $_ >= 50           } @{$A[$i]};
				@{$T[$i]} = grep { $_ <= 74           } @{$T[$i]};
			}
		}
	}


	###################################################################
	#For batch inputs,

	# print "A's before\n";
	# for (my $a = 0; $a < scalar @{$A[0]}; $a++) {
	# 	print "$A[0][$a]\t";
	# }
	# print "\n\n";

	# print "T's before\n";
	# for (my $t = 0; $t < scalar @{$T[0]}; $t++) {
	# 	print "$T[0][$t]\t";
	# }
	# print "\n\n";

	if ($batchinput eq "TRUE") {
		if ($n < 2) {
			for (my $ii = 0; $ii < scalar @{$T[0]}; $ii++) {
					if ($T[0][$ii] > 50) {
						splice (@{$T[0]}, $ii, 1);
						$ii--;
					}
			}
			for (my $ii = 0; $ii < scalar @{$A[-1]}; $ii++) {
				if ($A[-1][$ii] < 50) {
					splice (@{$A[-1]}, $ii, 1);
					$ii--;
				}
			}
		}
		for (my $ii = 0; $ii < scalar @{$A[0]}; $ii++) {
				if ($A[0][$ii] > 50) {
					splice (@{$A[0]}, $ii, 1);
					$ii--;
				}
		}
		for (my $ii = 0; $ii < scalar @{$T[-1]}; $ii++) {
				if ($T[-1][$ii] < 50) {
					splice (@{$T[-1]}, $ii, 1);
					$ii--;
				}
		}
		# print "A: ", scalar @{$A[0]}, "\n";
		# print "T: ", scalar @{$T[0]}, "\n";
	}

	# print "A's after\n";
	# for (my $a = 0; $a < scalar @{$A[0]}; $a++) {
	# 	print "$A[0][$a]\t";
	# }
	# print "\n\n";

	# print "T's after\n";
	# for (my $t = 0; $t < scalar @{$T[0]}; $t++) {
	# 	print "$T[0][$t]\t";
	# }
	# print "\n\n";

	# if ($A[0][-1] < 24) {
	# 	die "ERROR: A-to-T segment in joining region of $inputs[0][0][1] and $inputs[0][0+1][1] is $A[0][-1] bp from joining region, USER fusion might not be possible\n";
	# }
	# print "HERE $T[-1][0]\n";
	# if ($T[-1][0] > 56) {
	# 	die "ERROR: A-to-T segment in joining region of $inputs[0][-2][1] and $inputs[0][-1][1] is $T[0][0] bp from joining region, USER fusion might not be possible\n";	
	# }

	for (my $i = 0; $i < $n-1; $i++) {
		my $nicelengthflag = 0;
		# print "$i $joined[$i]\n";
		for (my $a = 0; $a < scalar @{$A[$i]}; $a++) {
			for (my $t = 0; $t < scalar @{$T[$i]}; $t++) {
				if ($A[$i][$a] < $T[$i][$t]) {
					my $tempAT = substr($joined[$i], $A[$i][$a], $T[$i][$t]-$A[$i][$a]);
					# print "$A[$i][$a] $T[$i][$t] $tempAT\n";
					if (length $tempAT > 7 and length $tempAT < 13) {
						$nicelengthflag = 1;
						# print "$A[$i][$a] $T[$i][$t] $tempAT\n";
					}
				}
			}
		}
		if ($nicelengthflag == 0) {
			if (exists $forced_overhang{$i}) {
				my $side = $forced_overhang{$i};
				my $fragname = ($side eq "A") ? $inputs[0][$i][1] : $inputs[0][$i+1][1];
				my $window = ($side eq "A") ? "last 24 nt of $inputs[0][$i][1]" : "first 24 nt of $inputs[0][$i+1][1]";
				die "ERROR: no A-to-T segment of 8-12 bp found in the $window (forced overhang junction " . ($i+1) . "). USER fusion is not possible with this forced-overhang setting.\n";
			} else {
				die "ERROR: no A-to-T segment in joining region of $inputs[0][$i][1] and $inputs[0][$i+1][1] is between 7 and 13 bp in length, USER fusion might not be possible\n";
			}
		}
	}

	return (\@A, \@T);
}



sub ATscore {

	my @joined = @{$_[0]};
    my @A = @{$_[1]};
    my @T = @{$_[2]};
    my $n = $_[3];
    my %forced_overhang = %{$_[4]} if defined $_[4];
    my @scoredfts;
	my $count;

	for (my $i = 0; $i < $n-1; $i++) { 
		$count = 0;
		my $penlen;	#penalty for length
		my $penpos;	#penalty for position
		my $score;	#score of fusion tail
		my $ft;		#fusion tail
		for (my $ia = 0; $ia < scalar @{$A[$i]}; $ia++) {
			for (my $it = 0; $it < scalar @{$T[$i]}; $it++) {
				if ($A[$i][$ia] < $T[$i][$it] and $A[$i][$ia] ne "" and $T[$i][$it] ne "") {
					$ft = substr($joined[$i], $A[$i][$ia], $T[$i][$it] - $A[$i][$ia]);	
					#if ($i == 0) {
						#print "$ft\n";
					#}
					if (length $ft < 7) {
						$penlen = 1000;
					}
					if (length $ft > 13) {
						$penlen = 1000;
					}
					elsif (length $ft > 9) {
						$penlen = (length $ft) - 9;
					}
					elsif (length $ft == 7 or length $ft == 8 or length $ft == 9) {
						$penlen = 0;
					}
					if    ( ( (length $joined[$i])/2 ) > ( $A[$i][$ia] + ($T[$i][$it]-$A[$i][$ia])/2 ) ) { 
						$penpos = (((length $joined[$i])/2) - ($A[$i][$ia]+(($T[$i][$it]-$A[$i][$ia])/2)));
					}
					elsif ( ( (length $joined[$i])/2 ) < ( $A[$i][$ia] + ($T[$i][$it]-$A[$i][$ia])/2 ) ) {
						$penpos = (($A[$i][$ia]+(($T[$i][$it]-$A[$i][$ia])/2)) - (length $joined[$i])/2);
					}
					elsif (((length $joined[$i])/2) == ($A[$i][$ia]+(($T[$i][$it]-$A[$i][$ia])/2))) {
						$penpos = 0;
					}
					$score = 0 + $penlen + $penpos;
					if ($i == 0) {
						# print "$i ", length $ft, " $ft\n$A[$i][$ia]\t$T[$i][$it]\t$penpos\t$penlen\t$score\n";
					}
					#if ($score < 1000) {
						# if ($i == 3) {
						# 	print length $ft, " $ft\n$A[$i][$ia]\t$T[$i][$it]\t$penpos\t$penlen\t$score\n";
						# }


					###################################################################
					#If the joining regions contain degenerate nucleotides, these must
					#be in the overhangs of the primers, not the fusion tails. All
					#fusion tails containing dg nts must be removed from the pool
				
					unless (substr ($joined[$i], $A[$i][$ia], $T[$i][$it] - $A[$i][$ia]) =~ m/[RYKMBVDHN]/) {
						$scoredfts[$i][$count] = "$score $A[$i][$ia] $T[$i][$it]";
						$count++;
					}
						#print "$i\n";
						#print "$i ", length $ft, " $ft\n$A[$i][$ia]\t$T[$i][$it]\t$penpos\t$penlen\t$score\n";
					#}
				}
			}
		}
		#print "$i $count\n";
		# if ($n > 1) {
		# 	if ($count == 0) {
		# 		die "ERROR\n";
		# 	}
		# }
	}
	

	for (my $i = 0; $i < $n-1; $i++) {
		@{$scoredfts[$i]} = sort {$a <=> $b} @{$scoredfts[$i]};
		# print "$i\n";
		for (my $ii= 0; $ii < scalar @{$scoredfts[$i]}; $ii++) {
			# print "$scoredfts[$i][$ii]\n";
		}

		# Forced-overhang: warn if the best-scoring span requires a bridge > 8 nt
		if (exists $forced_overhang{$i} and scalar @{$scoredfts[$i]} > 0) {
			$scoredfts[$i][0] =~ m/\d+\.*\d*\s(\d+)\s(\d+)/;
			my ($best_Apos, $best_Tpos) = ($1, $2);
			my $bridge_len;
			if ($forced_overhang{$i} eq "A") {
				$bridge_len = 50 - $best_Tpos;   # nt from span end to A/B boundary
			} else {
				$bridge_len = $best_Apos - 50;   # nt from A/B boundary to span start
			}
			if ($bridge_len > 8) {
				my $span = substr($joined[$i], $best_Apos, $best_Tpos - $best_Apos);
				print STDERR "WARNING: junction " . ($i+1) . " forced-overhang span \"$span\" requires a bridge of $bridge_len nt (> 8 nt threshold). Primers will be longer than optimal.\n";
			}
		}
	}
	return \@scoredfts;
}



sub scores {

	my @scoredfts = @{$_[0]};
	my $n = $_[1];
	my @scores;

	for (my $i = 0; $i < $n-1; $i++) {
		for (my $ii = 0; $ii < $#{$scoredfts[$i]}; $ii++) {
			$scoredfts[$i][$ii] =~ m/(\d+\.*\d*)\s\d+\s\d+/;
			$scores[$i][$ii] = $1;
			# print "$scores[$i][$ii]\n";
		}
	}
	return \@scores;
}



# A junction window of 100 nt (50 nt of each neighbour) cannot always supply a
# binding region of every length in the window: the reverse region needs
# Apos >= len and the forward region needs Tpos + len <= 100. bindingregion
# leaves those lengths unset, so every consumer must ask whether a length exists
# rather than assume it. br_have: is this candidate sequence present and non-empty.
sub br_have { return (defined $_[0] and length $_[0]) ? 1 : 0; }

# The binding-region lengths this strand actually has, shortest first, read from
# the Tm array (which is only filled for lengths whose sequence exists).
sub tm_lengths {
	my $slot = $_[0];
	return () unless defined $slot;
	return grep { defined $slot->[$_] and defined $slot->[$_][0] } ($BR_LEN_MIN .. $BR_LEN_MAX);
}

sub bindingregion {
    my @scoredfts = @{$_[0]};
    my @joined = @{$_[1]};
    my $n = $_[2];
	my @Apos;
	my @Tpos;
	my @brforpos;
	my @brrevpos;
	my @brfor;
	my @brrev;
	my @fts;

	for (my $i = 0; $i < $n-1; $i++) {
		for (my $ii = 0; $ii < scalar @{$scoredfts[$i]}; $ii++) {
			# sprint "$i $ii\n";
			my $temp;
			my $insertlimstart = 50;
			my $insertlimend = (length $joined[$i])-50;
			# print "$i start $insertlimstart\n";
			# print "$i end $insertlimend\n";
			# print length $joined[$i], " $insertlimstart $insertlimend\n";
			if ($scoredfts[$i][$ii] =~ m/\d+\.*\d*\s(\d+)\s(\d+)/) {
				$Apos[$i][$ii] = $1;
				$Tpos[$i][$ii] = $2;
				# print "$Apos[$i][$ii]\t$Tpos[$i][$ii]\n";
				# if (length $joined[$i] > 80) {
				# 	if (((length $joined[$i])/2)
				# }
			}
			for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
				my $fwstart = ($Tpos[$i][$ii] > $insertlimend) ? $Tpos[$i][$ii] : $insertlimend;
				$brforpos[$i][$ii] = $fwstart;
				# the window holds 50 nt of each neighbour; a longer region than
				# remains to its right does not exist here, so skip that length
				next if $fwstart + $iii > length $joined[$i];
				$temp = substr ($joined[$i], $fwstart, $iii);
				if (length $temp == $iii) {
					$brfor[$i+1][$ii][$iii] = $temp;
				}
			}
			for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
				my $rvend = ($insertlimstart > $Apos[$i][$ii]) ? $Apos[$i][$ii] : $insertlimstart;
				$brrevpos[$i][$ii] = $rvend;
				# a region longer than the sequence to its left does not exist in
				# this window; a negative substr start would silently wrap
				next if $rvend - $iii < 0;
				$temp = substr ($joined[$i], $rvend-$iii, $iii);
				$temp =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				if (length $temp == $iii) {
					$brrev[$i][$ii][$iii] = $temp;
					#print "$i $ii $iii $brrev[$i][$ii][$iii]\n";
				}
			}

			$fts[$i][$ii] = substr ($joined[$i], $Apos[$i][$ii], $Tpos[$i][$ii] - $Apos[$i][$ii]);
		}	
	}
	return (\@Apos, \@Tpos, \@brfor, \@brrev, \@fts, \@brforpos, \@brrevpos);	


	# for (my $i = 0; $i < $n; $i++) {
	# 	for (my $ii = 0; $ii < $#{$scoredfts[$i]}; $ii++) {
	# 		my $temp;
	# 		if ($scoredfts[$i][$ii] =~ m/\d+\.*\d*\s(\d+)\s(\d+)/) {
	# 			$Apos[$i][$ii] = $1;
	# 			$Tpos[$i][$ii] = $2;
	# 		}
	# 		for (my $iii = 7; $iii < 25; $iii++) {
	# 			if ($Tpos[$i][$ii] < (length $joined[$i])/2) {
	# 				$brfor[$i+1][$ii][$iii] = substr ($joined[$i], $Tpos[$i][$ii], ($iii + ((length $joined[$i])/2)-$Tpos[$i][$ii]));
	# 			}
	# 			else {
	# 				$brfor[$i+1][$ii][$iii] = substr ($joined[$i], $Tpos[$i][$ii], $iii);
	# 			}
	# 		}
			
	# 		for (my $iii = 7; $iii < 25; $iii++) {
	# 			if ($Apos[$i][$ii] > (length $joined[$i])/2) {
	# 				$temp = substr ($joined[$i], $Apos[$i][$ii]-($iii + ($Apos[$i][$ii] - (length $joined[$i])/2)), $iii + ($Apos[$i][$ii] - (length $joined[$i])/2));
	# 			}
	# 			else {		
	# 				$temp = substr ($joined[$i], $Apos[$i][$ii]-$iii, $iii);
	# 			}
	# 			$temp =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
	# 			$brrev[$i][$ii][$iii] = $temp;
	# 			print "$brfor[$i+1][$ii][$iii]\n";
	# 		}
			
	# 		$fts[$i][$ii] = substr ($joined[$i], $Apos[$i][$ii], $Tpos[$i][$ii] - $Apos[$i][$ii]);
	# 	}
	# }
	# return (\@Apos, \@Tpos, \@brfor, \@brrev, \@fts);
}



sub rmclashft {
	
	my @fts = @{$_[0]};
    my @scores = @{$_[1]};
    my @brrev = @{$_[2]};
    my @brfor = @{$_[3]};
    my @Apos = @{$_[4]};
    my @Tpos = @{$_[5]};
    my @scoredfts = @{$_[6]};
    my $n = $_[7];
	my @ft;
	my $flag = 1;

	while ($flag == 1) {
		$flag = 0;
		for (my $i = 0; $i < $n-1; $i++) {
			for (my $ii = 1; $ii < $n; $ii++) {
				if ($fts[$i][0] eq $fts[$i+$ii][0]) {
					if ($scores[$i][0] > $scores[$i+$ii][0]) {	
						shift (@{$fts[$i+$ii]});
						shift (@{$scores[$i+$ii]});
						shift (@{$brrev[$i+$ii]});
						shift (@{$brfor[$i+$ii+1]});
						shift (@{$Apos[$i+$ii]});
						shift (@{$Tpos[$i+$ii]});
						if ((scalar @{$fts[$i]}) == 0) {
							die "ERROR: it is not possible to generate unique fusion tails for all fragments. This can, for example, occur when attempting to fuse several identical fragments together. Please check input\n";
						}
					}
					elsif ($scores[$i][0] < $scores[$i+$ii][0]) {
						shift (@{$fts[$i]});
						shift (@{$scores[$i]});
						shift (@{$brrev[$i]});
						shift (@{$Apos[$i]});
						shift (@{$Tpos[$i]});
						shift (@{$brfor[$i+1]});
						if ((scalar @{$fts[$i]}) == 0) {
							die "ERROR: it is not possible to generate unique fusion tails for all fragments. This can, for example, occur when attempting to fuse several identical fragments together. Please check input\n";
						}						
					}
					elsif ($scores[$i][0] == $scores[$i+$ii][0]) {
						if ($#{$scoredfts[$i]} > $#{$scoredfts[$i+$ii]}) {
							shift (@{$fts[$i]});
							shift (@{$scores[$i]});
							shift (@{$brrev[$i]});
							shift (@{$brfor[$i+1]});
							shift (@{$Apos[$i]});
							shift (@{$Tpos[$i]});
							if ((scalar @{$fts[$i]}) == 0) {
								die "ERROR: it is not possible to generate unique fusion tails for all fragments. This can, for example, occur when attempting to fuse several identical fragments together. Please check input\n";
							}
						}
						elsif ($#{$scoredfts[$i]} < $#{$scoredfts[$i+1]}) {
							shift (@{$fts[$i+$ii]});
							shift (@{$scores[$i+$ii]});
							shift (@{$brrev[$i+$ii]});
							shift (@{$brfor[$i+$ii+1]});
							shift (@{$Apos[$i+$ii]});
							shift (@{$Tpos[$i+$ii]});
							if ((scalar @{$fts[$i]}) == 0) {
								die "ERROR: it is not possible to generate unique fusion tails for all fragments. This can, for example, occur when attempting to fuse several identical fragments together. Please check input\n";
							}								
						}
						elsif ($#{$scoredfts[$i]} == $#{$scoredfts[$i+1]}) {
							shift (@{$fts[$i]});
							shift (@{$scores[$i]});
							shift (@{$brrev[$i]});
							shift (@{$brfor[$i+1]});
							shift (@{$Apos[$i]});
							shift (@{$Tpos[$i]});
							if ((scalar @{$fts[$i]}) == 0) {
								die "ERROR: it is not possible to generate unique fusion tails for all fragments. This can, for example, occur when attempting to fuse several identical fragments together. Please check input\n";
							}
						}
						else {
							die "ERROR: it is not possible to generate unique fusion tails for all fragments. This can, for example, occur when attempting to fuse several identical fragments together. Please check input\n";
						}
					}
					else {
						die "ERROR: it is not possible to generate unique fusion tails for all fragments. This can, for example, occur when attempting to fuse several identical fragments together. Please check input\n";
					}
					$flag = 1;
			 	}
				elsif ($fts[$i][0] ne $fts[$i+$ii][0]) {
					$ft[$i] = $fts[$i][0];
				}
				else {
					die "ERROR: it is not possible to generate unique fusion tails for all fragments. This can, for example, occur when attempting to fuse several identical fragments together. Please check input\n";
				}

			}
		}
	}
	return (\@ft);
}



sub uracilft {

	my @ft = @{$_[0]};
	my $n = $_[1];
	my @ftcomp;

	for (my $i = 0; $i < $n-1; $i++) {
		$ftcomp[$i] = $ft[$i];
	}
	for (my $i = 0; $i < $n-1; $i++) {
		$ftcomp[$i] =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		$ftcomp[$i] =~ s/T/U/;
	}
	for (my $i = 0; $i < $n-1; $i++) {
		chop $ft[$i];
		$ft[$i] .= "U";
	}
	return (\@ft, \@ftcomp);
}



sub extensions {
	my @primerextension = @{$_[0]};
	my @joined = @{$_[1]};
	my @brforpos = @{$_[2]};
	my @brrevpos = @{$_[3]};
	my @ft = @{$_[4]};
	my $n = $_[5];
	my @Apos = @{$_[6]};
	my @Tpos = @{$_[7]};
	my @brfor = @{$_[8]};
	my @brrev = @{$_[9]};
	my $x = 0;
	for (my $i = 0; $i < scalar @inputs; $i++) {
		my $fwextension;
		my $rvextension;
		if ($primerextension[$i][0][0] eq "TRUE") {
			$fwextension = $primerextension[$i][0][1];
		}
		$primerextension[$i][0][4] = $fwextension;

		for (my $ii = 1; $ii < $n; $ii++) {
			$fwextension = "";
			$rvextension = "";
			if ($primerextension[$i][$ii][0] eq "TRUE") {
				# print "$ii T: $Tpos[$ii-1][0] pos: $brforpos[$ii-1][0]\n";
				# print "$ii A: $Apos[$ii-1][0] pos: $brrevpos[$ii-1][0]\n";
				# print "$joined[$ii-1]\n";
				if ($Tpos[$ii-1][0] < $brforpos[$ii-1][0]) {
					$fwextension = substr ($joined[$ii-1], $Tpos[$ii-1][0], $brforpos[$ii-1][0] - $Tpos[$ii-1][0]);
				}
				
				if ($Apos[$ii-1][0] > $brrevpos[$ii-1][0]) {
					$rvextension = substr ($joined[$ii-1], $brrevpos[$ii-1][0], $Apos[$ii-1][0] - $brrevpos[$ii-1][0]);
					$rvextension =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				}
			}
			else {
				if ($Tpos[$ii-1][0] < $brforpos[$ii-1][0]) {
					$fwextension = substr ($joined[$ii-1], $Tpos[$ii-1][0], $brforpos[$ii-1][0] - $Tpos[$ii-1][0]);
				}
				if ($Apos[$ii-1][0] > $brrevpos[$ii-1][0]) {
					$rvextension = substr ($joined[$ii-1], $brrevpos[$ii-1][0], $Apos[$ii-1][0] - $brrevpos[$ii-1][0]);
					$rvextension =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				}
			}
			$primerextension[$i][$ii][4] = $fwextension;
			$primerextension[$i][$ii-1][5] = $rvextension;
		}

		$rvextension = "";
		if ($primerextension[$i][$n][0] eq "TRUE") {
			$rvextension = $primerextension[$i][$n][1];
			$rvextension =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
		}
		$primerextension[$i][$n][5] = $rvextension;
	}
	return (\@primerextension)
}



sub ATGCcount {

    my @brfor = @{$_[0]};
    my @brrev = @{$_[1]};
    my $n = $_[2];	
	my @Afor;
	my @Tfor;
	my @Gfor;
	my @Cfor;
	my @Arev;
	my @Trev;
	my @Grev;
	my @Crev;

	for (my $i = 0; $i < $n; $i++) {
		for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
			next unless &br_have($brfor[$i][0][$iii]) or &br_have($brrev[$i][0][$iii]);
			for (my $ii = 0; $ii < length $brfor[$i][0][$iii]; $ii++) {
				if (substr ($brfor[$i][0][$iii], $ii, 1) eq "A") {
					$Afor[$i][0][$iii]++;
				}
				elsif (substr ($brfor[$i][0][$iii], $ii, 1) eq "T") {
					$Tfor[$i][0][$iii]++;
				}
				elsif (substr ($brfor[$i][0][$iii], $ii, 1) eq "G") {
					$Gfor[$i][0][$iii]++;
				}
				elsif (substr ($brfor[$i][0][$iii], $ii, 1) eq "C") {
					$Cfor[$i][0][$iii]++;
				}	
			}
			for (my $ii = 0; $ii < length $brrev[$i][0][$iii]; $ii++) {
				if (substr ($brrev[$i][0][$iii], $ii, 1) eq "A") {
					$Arev[$i][0][$iii]++;
				}
				elsif (substr ($brrev[$i][0][$iii], $ii, 1) eq "T") {
					$Trev[$i][0][$iii]++;
				}
				elsif (substr ($brrev[$i][0][$iii], $ii, 1) eq "G") {
					$Grev[$i][0][$iii]++;
				}
				elsif (substr ($brrev[$i][0][$iii], $ii, 1) eq "C") {
					$Crev[$i][0][$iii]++;
				}	
			}
		}
	}
	return (\@Afor, \@Tfor, \@Gfor, \@Cfor, \@Arev, \@Trev, \@Grev, \@Crev);
}


sub GCratiosingle {

	my @Afor = @{$_[0]};
    my @Tfor = @{$_[1]};
    my @Gfor = @{$_[2]};
    my @Cfor = @{$_[3]};
    my @Arev = @{$_[4]};
    my @Trev = @{$_[5]};
    my @Grev = @{$_[6]};
    my @Crev = @{$_[7]};
    my @brrev = @{$_[8]};
    my @brfor = @{$_[9]};
	my @TmGCrev;
	my @TmGCfor;

	for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
		next unless &br_have($brrev[0][0][$iii]) and &br_have($brfor[0][0][$iii]);
		$TmGCrev[0][0][$iii][1] = ($Grev[0][0][$iii]+$Crev[0][0][$iii]) / (length $brrev[0][0][$iii]) * 100;
		$TmGCrev[0][0][$iii][1] = sprintf ("%.2f", $TmGCrev[0][0][$iii][1]);
		$TmGCfor[0][0][$iii][1] = ($Gfor[0][0][$iii]+$Cfor[0][0][$iii]) / (length $brfor[0][0][$iii]) * 100;
		$TmGCfor[0][0][$iii][1] = sprintf ("%.2f", $TmGCfor[0][0][$iii][1]);
	}
	return (\@TmGCrev, \@TmGCfor);
}


sub GCratiomulti {

	my @Afor = @{$_[0]};
    my @Tfor = @{$_[1]};
    my @Gfor = @{$_[2]};
    my @Cfor = @{$_[3]};
    my @Arev = @{$_[4]};
    my @Trev = @{$_[5]};
    my @Grev = @{$_[6]};
    my @Crev = @{$_[7]};
    my @brrev = @{$_[8]};
    my @brfor = @{$_[9]};
    my $n = $_[10];
	my @TmGCrev;
	my @TmGCfor;

	for (my $i = 0; $i < $n; $i++) {
		for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
			next unless &br_have($brrev[$i][0][$iii]);
			$TmGCrev[$i][0][$iii][1] = ($Grev[$i][0][$iii]+$Crev[$i][0][$iii]) / (length $brrev[$i][0][$iii]) * 100;
			$TmGCrev[$i][0][$iii][1] = sprintf ("%.2f", $TmGCrev[$i][0][$iii][1]);
		}
		for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
			next unless &br_have($brfor[$i][0][$iii]);
			$TmGCfor[$i][0][$iii][1] = ($Gfor[$i][0][$iii]+$Cfor[$i][0][$iii]) / (length $brfor[$i][0][$iii]) * 100;
			$TmGCfor[$i][0][$iii][1] = sprintf ("%.2f", $TmGCfor[$i][0][$iii][1]);
			#print "$i\t$iii\t$brfor[$i][0][$iii]\t$TmGCfor[$i][0][$iii][1]\n";
		}
	}
	# for (my $i = 0; $i < $n; $i++) {
	# 	for (my $iii = 18+$linker{$i+1}{2}; $iii < 25+$linker{$i+1}{2}; $iii++) {
	# 		$TmGCrev[$i][0][$iii][1] = ($Grev[$i][0][$iii]+$Crev[$i][0][$iii]) / (length $brrev[$i][0][$iii]) * 100;
	# 		$TmGCrev[$i][0][$iii][1] = sprintf ("%.2f", $TmGCrev[$i][0][$iii][1]);
	# 	}
	# 	for (my $iii = 18+$linker{$i}{2}; $iii < 25+$linker{$i}{2}; $iii++) {
	# 		$TmGCfor[$i][0][$iii][1] = ($Gfor[$i][0][$iii]+$Cfor[$i][0][$iii]) / (length $brfor[$i][0][$iii]) * 100;
	# 		$TmGCfor[$i][0][$iii][1] = sprintf ("%.2f", $TmGCfor[$i][0][$iii][1]);
	# 	}
	# }
	return (\@TmGCrev, \@TmGCfor);
}


# Tm of one binding region, 5'->3'. NEB Tm Calculator equivalent: SantaLucia
# 1998 unified NN parameters, terminal corrections, symmetry correction for
# self-complementary sequences, then the salt correction belonging to the
# selected polymerase. $sc is molar monovalent salt, $pc molar primer, as
# everywhere else in this script. Degenerate bases are handled by nn_tm, which
# expands them before calling this.
# Every concrete sequence a degenerate one stands for, or () if that would
# exceed $IUPAC_MAX_COMBOS. A sequence of plain A/C/G/T expands to itself.
sub iupac_expand {
	my $seq = uc $_[0];
	my $ambig = 0;
	foreach my $b (split //, $seq) {
		return () unless exists $IUPAC{$b};
		# U is a plain T here, as it is to NEB's nonCanonicalCount, which counts
		# only wsrymkbdhvn.
		$ambig++ if length($IUPAC{$b}) > 1;
		return () if $ambig > $IUPAC_MAX_AMBIG;
	}
	my @out = ('');
	foreach my $b (split //, $seq) {
		my @bases = split //, $IUPAC{$b};
		@out = map { my $p = $_; map { $p . $_ } @bases } @out;
	}
	return @out;
}

# How a Tm is shown: one value, or NEB's "min - max" when the binding region is
# degenerate and the variants disagree once rounded to 0.1 C.
sub tm_display {
	my ($min, $max) = @_;
	return "" unless defined $min;
	my $lo = sprintf("%.1f", $min);
	my $hi = defined $max ? sprintf("%.1f", $max) : $lo;
	return $lo eq $hi ? $lo : "$lo - $hi";
}

# Tm of a binding region. In list context returns (min, max) across the
# degenerate expansion; in scalar context the min. undef if the sequence
# cannot be expanded within the cap, or holds a character that is not DNA.
sub nn_tm {
	my ($seq, $sc, $pc) = @_;
	return undef unless defined $seq and length($seq) > 1;
	my @seqs = &iupac_expand($seq);
	return undef unless @seqs;
	my ($min, $max);
	foreach my $one (@seqs) {
		my $t = &nn_tm_one($one, $sc, $pc);
		return undef unless defined $t;
		$min = $t if !defined $min or $t < $min;
		$max = $t if !defined $max or $t > $max;
	}
	return wantarray ? ($min, $max) : $min;
}

# The model itself, on one concrete A/C/G/T sequence.
sub nn_tm_one {
	my ($seq, $sc, $pc) = @_;
	return undef unless defined $seq and length($seq) > 1;
	$seq = uc $seq;

	my ($dH, $dS) = (0, 0);
	for (my $i = 0; $i < (length($seq) - 1); $i++) {
		my $d = substr($seq, $i, 2);
		return undef unless exists $NN_H{$d};
		$dH += $NN_H{$d};
		$dS += $NN_S{$d};
	}
	foreach my $e (substr($seq, 0, 1), substr($seq, -1)) {
		return undef unless exists $NN_TERM_H{$e};
		$dH += $NN_TERM_H{$e};
		$dS += $NN_TERM_S{$e};
	}
	# Self-complementary duplex: SantaLucia symmetry correction. Omitting this
	# is the one place the polymerase-tm package departs from NEB (+0.7..1.0 C
	# on perfect palindromes); the concentration term must NOT also change.
	my $rc = reverse $seq;
	$rc =~ tr/ACGT/TGCA/;
	$dS += -1.4 if $seq eq $rc;

	my $denom = $dS + $NN_R * log($pc);
	return undef unless $denom != 0;
	my $tm = ($dH * 1000) / $denom - 273.15;

	if ($TM_POLY{$TM_POLY}[2] eq 'schildkraut') {
		$tm += 16.6 * (log($sc) / log(10));
	}
	else {
		my $fgc = ($seq =~ tr/GC//) / length($seq);
		my $ln  = log($sc);
		my $c   = (4.29e-5 * $fgc - 3.95e-5) * $ln + 9.40e-6 * $ln ** 2;
		$tm = 1 / (1 / ($tm + 273.15) + $c) - 273.15;
	}
	return $tm;
}


sub Tm {
	my @brfor = @{$_[0]};
    my @brrev = @{$_[1]};
    my @TmGCfor = @{$_[2]};
    my @TmGCrev = @{$_[3]};
	my $sc = $_[4];
	my $pc = $_[5];
	my $n = $_[6];

	for (my $i = 0; $i < $n; $i++) {
		for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
			next unless &br_have($brrev[$i][0][$iii]);
			# brrev is stored 3'->5' and already complemented, so reversing it
			# gives the primer's real 5'->3' sequence.
			my $tempseq = reverse $brrev[$i][0][$iii];
			my ($tmin, $tmax) = &nn_tm($tempseq, $sc, $pc);
			if (defined $tmin) {
				$TmGCrev[$i][0][$iii][0] = sprintf ("%.1f", $tmin);
				$TmGCrev[$i][0][$iii][3] = &tm_display($tmin, $tmax);
			}
			else { $NN_TOO_DEGENERATE{"$i\trev"}++ }

			# print "rev $i $iii $TmGCrev[$i][0][$iii][0]\t$brrev[$i][0][$iii]\n";
		}

		for (my $iii = $BR_LEN_MIN; $iii <= $BR_LEN_MAX; $iii++) {
			next unless &br_have($brfor[$i][0][$iii]);
			my ($tmin, $tmax) = &nn_tm($brfor[$i][0][$iii], $sc, $pc);
			if (defined $tmin) {
				$TmGCfor[$i][0][$iii][0] = sprintf ("%.1f", $tmin);
				$TmGCfor[$i][0][$iii][3] = &tm_display($tmin, $tmax);
			}
			else { $NN_TOO_DEGENERATE{"$i\tfor"}++ }

			#print "for $i $iii $TmGCfor[$i][0][$iii][0]\t$brfor[$i][0][$iii]\n";
		}
	}
	return (\@TmGCrev, \@TmGCfor);
}


sub simpleTm {
	my $tempseq = $_[0];
	my $sc = $_[1];
	my $pc = $_[2];
	my $temptm = &nn_tm($tempseq, $sc, $pc);
	return undef unless defined $temptm;
	return sprintf ("%.1f", $temptm);
}

	
sub optbr {

    my @TmGCrev = @{$_[0]};
    my @TmGCfor = @{$_[1]}; 
	my $n = $_[2];

	my @optbr;
	my %primerTmwarning;
	for (my $i = 0; $i < $n; $i++) {
		my @revlens = &tm_lengths($TmGCrev[$i][0]);
		my @forlens = &tm_lengths($TmGCfor[$i][0]);
		if (!@revlens or !@forlens) {
			my $strand = !@revlens ? "reverse" : "forward";
			my $degen  = $NN_TOO_DEGENERATE{"$i\t" . (!@revlens ? "rev" : "for")} || 0;
			die "ERROR: fragment " . ($i+1) . " has no usable binding region. Every length\n"
			  . "       between $BR_LEN_MIN and $BR_LEN_MAX nt on the $strand strand carries more than\n"
			  . "       $IUPAC_MAX_AMBIG ambiguous bases, which is NEB's limit too, so its Tm is not a\n"
			  . "       single value. Use a less degenerate sequence near that junction.\n"
				if $degen;
			die "ERROR: fragment " . ($i+1) . " has no usable binding region. The 100 nt\n"
			  . "       junction window could not supply a region of any length between\n"
			  . "       $BR_LEN_MIN and $BR_LEN_MAX nt on the $strand strand,\n"
			  . "       because the chosen USER overhang sits too far from the fragment\n"
			  . "       boundary. Try a different fusion order, or -fo to force the overhang.\n";
		}
		my $counter = 0;
		foreach my $rev (@revlens) {
			foreach my $for (@forlens) {
				if (($TmGCrev[$i][0][$rev][0] - $TmGCfor[$i][0][$for][0]) <= $TM_PAIR_MAX_DIFF and 
				    ($TmGCrev[$i][0][$rev][0] - $TmGCfor[$i][0][$for][0]) >= -$TM_PAIR_MAX_DIFF) {
				    	$optbr[$i][0][$counter] = $rev;
					$optbr[$i][1][$counter] = $for;
					$counter++
				}
			}
		}

		if ($optbr[$i][0][0] eq "" and $optbr[$i][1][0] eq "") {
		
			$primerTmwarning{$i} = 1;
			
			# the window edges are the shortest/longest lengths this junction
			# actually has, which are the constants unless the window was short
			my ($rvmin, $rvmax) = ($revlens[0], $revlens[-1]);
			my ($fwmin, $fwmax) = ($forlens[0], $forlens[-1]);

			my $mintemprev = $TmGCrev[$i][0][$rvmin][0];
			my $maxtemprev = $TmGCrev[$i][0][$rvmax][0];

			my $mintempfor = $TmGCfor[$i][0][$fwmin][0];
			my $maxtempfor = $TmGCfor[$i][0][$fwmax][0];

			if ($TmGCrev[$i][0][$rvmin][0] > $TmGCfor[$i][0][$fwmax][0]) {
				$optbr[$i][0][0] = $rvmin;
				$optbr[$i][1][0] = $fwmax;

				my $alt = 1;
				foreach my $f (grep { $_ < $fwmax } reverse @forlens) {
					last if $alt > 2;
					if (($TmGCrev[$i][0][$rvmin][0] - $TmGCfor[$i][0][$f][0]) < 5) {
						$optbr[$i][0][$alt] = $rvmin;
						$optbr[$i][1][$alt] = $f;
						$alt++;
					}
					else { last }
				}
			}
			elsif ($TmGCrev[$i][0][$rvmax][0] < $TmGCfor[$i][0][$fwmin][0]) {
				$optbr[$i][0][0] = $rvmax;
				$optbr[$i][1][0] = $fwmin;
				my $alt = 1;
				foreach my $f (grep { $_ > $fwmin } @forlens) {
					last if $alt > 2;
					if (($TmGCfor[$i][0][$f][0] - $TmGCrev[$i][0][$rvmax][0]) < 5) {
						$optbr[$i][0][$alt] = $rvmax;
						$optbr[$i][1][$alt] = $f;
						$alt++;
					}
					else { last }
				}
			}
		}
	}
	return (\@optbr, \%primerTmwarning);
}


# The available length whose Tm is nearest $target, lowest length winning a tie --
# what optbrTm returned on its own before the tolerance pool was introduced.
sub closest_to_target {
	my ($slot, $target) = @_;
	my ($best, $bestdiff);
	foreach my $L (&tm_lengths($slot)) {
		my $d = abs($slot->[$L][0] - $target);
		if (!defined $bestdiff or $d < $bestdiff) { $bestdiff = $d; $best = $L; }
	}
	return $best;
}

sub optbrTm {
    my @TmGCrev = @{$_[0]};
    my @TmGCfor = @{$_[1]};
    my $Tmtarget = $_[2];
	my $n = $_[3];

	my @optbr;
	my %primerTmwarning;

	# Which pool tolerance applies has to be decided here, before the "eq"
	# branch below overwrites $Tmtarget with the candidate-pool median: a
	# user-supplied target is owed accuracy, a median is not.
	my @tolsteps = ($Tmtarget eq "eq") ? @TM_TARGET_TOL_EQ_STEPS
	                                  : @TM_TARGET_TOL_FIXED_STEPS;



	###################################################################
	#If target Tm is "eq" (i.e. no specific target value, but Tm for 
	#all primers to be as close to each other as possible), the median
	#Tm for all possible binding regions is set as $Tmtarget


	if ($Tmtarget eq "eq") {

		###################################################################
		#"eq" asks for primers whose Tms are as close to each other as
		#possible. The target used to be the median over ALL candidate Tms
		#pooled together, which optimises nothing in particular and lets a
		#strand-fragment offering 11 usable lengths outvote one offering 3.
		#It is now the target that minimises the spread actually realised:
		#every strand-fragment takes its nearest achievable Tm, and the
		#candidate target with the smallest (max - min) wins. Ties break
		#toward the old median, so a run whose median already minimises the
		#spread keeps exactly the primers it had.
		#
		#Each strand-fragment's achievable set is kept separate (@sets);
		#@Tmvalues stays pooled only to compute the median tie-break.
		#The optimum is attained at an achievable Tm, so scanning the
		#distinct achievable values is exact.

		my @sets;
		my @Tmvalues;

		for (my $i = 0; $i < $n; $i++) {
			my @r;
			foreach my $ii (&tm_lengths($TmGCrev[$i][0])) {
				push(@r, $TmGCrev[$i][0][$ii][0]);
				push(@Tmvalues, $TmGCrev[$i][0][$ii][0]);
			}
			push(@sets, \@r) if scalar @r;

			my @f;
			foreach my $ii (&tm_lengths($TmGCfor[$i][0])) {
				push(@f, $TmGCfor[$i][0][$ii][0]);
				push(@Tmvalues, $TmGCfor[$i][0][$ii][0]);
			}
			push(@sets, \@f) if scalar @f;
		}

		my $Tmcount = scalar @Tmvalues;
		my @sortedTmvalues = sort { $a <=> $b } @Tmvalues;

		my $Tmmedian;
		if ($Tmcount % 2) {
			$Tmmedian = $sortedTmvalues[int($Tmcount/2)];
		} 
		else {
			$Tmmedian = ($sortedTmvalues[$Tmcount/2] + $sortedTmvalues[$Tmcount/2 - 1]) / 2;
		}

		$Tmtarget = $Tmmedian;

		if (scalar @sets > 1) {

			my @cands;
			my $prev;
			foreach my $t (@sortedTmvalues) {
				push(@cands, $t) unless defined $prev and abs($t - $prev) < 1e-9;
				$prev = $t;
			}

			my $bestspread;
			my $bestdist;
			my $besttarget;

			foreach my $cand (@cands) {
				my $lo;
				my $hi;
				foreach my $set (@sets) {
					my $near;
					foreach my $t (@$set) {
						$near = $t if !defined $near or abs($t - $cand) < abs($near - $cand);
					}
					$lo = $near if !defined $lo or $near < $lo;
					$hi = $near if !defined $hi or $near > $hi;
				}
				my $spread = $hi - $lo;
				my $dist = abs($cand - $Tmmedian);
				if (!defined $bestspread
				    or $spread < $bestspread - 1e-9
				    or (abs($spread - $bestspread) < 1e-9 and $dist < $bestdist)) {
					$bestspread  = $spread;
					$bestdist    = $dist;
					$besttarget  = $cand;
				}
			}

			$Tmtarget = $besttarget if defined $besttarget;
		}
	}


	###################################################################
	#The variable $Tmtarget is now set (either to median, or user
	#defined), and the optimal binding regions adapted to this Tm can
	#be stored in the array @optbr	
	
	
	for (my $i = 0; $i < $n; $i++) {
		unless (&tm_lengths($TmGCrev[$i][0]) and &tm_lengths($TmGCfor[$i][0])) {
			die "ERROR: fragment " . ($i+1) . " has no usable binding region. The 100 nt\n"
			  . "       junction window could not supply a region of any length between\n"
			  . "       $BR_LEN_MIN and $BR_LEN_MAX nt, because the chosen USER overhang\n"
			  . "       sits too far from the fragment boundary. Try a different fusion\n"
			  . "       order, or -fo to force the overhang.\n";
		}
		# the length closest to the target on each strand: the answer this routine
		# used to return on its own, and still the fallback when nothing else fits
		my $r0 = &closest_to_target($TmGCrev[$i][0], $Tmtarget);
		my $f0 = &closest_to_target($TmGCfor[$i][0], $Tmtarget);

		###################################################################
		#Every length within the current tolerance step of the target is
		#kept, and the surviving FW/RV combinations become the candidate pool
		#that optbranalysis ranks on the remaining criteria. Taking only the
		#single closest length threw that choice away: at +-1 C a quarter of
		#-tm <N> fragments and over half of -tm eq ones have more than one
		#length to choose between.
		#
		#The steps are tried in order and the FIRST one that yields a pair
		#wins, per fragment: a fragment that can sit within the tightest step
		#of the target is held there, and only a fragment that cannot is given
		#a wider window - rather than being dropped to the unranked fallback.
		#$tol_used is kept for the report/diagnostics.

		my @pairs;
		my $tol_used;

		foreach my $tol (@tolsteps) {

			my @revpool = grep { abs($TmGCrev[$i][0][$_][0] - $Tmtarget) <= $tol }
			              &tm_lengths($TmGCrev[$i][0]);
			my @forpool = grep { abs($TmGCfor[$i][0][$_][0] - $Tmtarget) <= $tol }
			              &tm_lengths($TmGCfor[$i][0]);

			# pair only combinations the FW/RV filter accepts, so primerdimer and
			# the pair-Tm row are never handed a mismatched pair
			@pairs = ();
			foreach my $rev (@revpool) {
				foreach my $for (@forpool) {
					next if abs($TmGCrev[$i][0][$rev][0] - $TmGCfor[$i][0][$for][0]) > $TM_PAIR_MAX_DIFF;
					push @pairs, [$rev, $for];
				}
			}

			if (scalar @pairs) { $tol_used = $tol; last }
		}

		# closest to the target first, so a tie on every other criterion still
		# returns the primer this routine would have returned before
		@pairs = sort {
			   ( abs($TmGCrev[$i][0][$a->[0]][0] - $Tmtarget) + abs($TmGCfor[$i][0][$a->[1]][0] - $Tmtarget) )
			<=> ( abs($TmGCrev[$i][0][$b->[0]][0] - $Tmtarget) + abs($TmGCfor[$i][0][$b->[1]][0] - $Tmtarget) )
			|| $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]
		} @pairs;

		# nothing in tolerance, or nothing survived the pair filter: the closest
		# pair, exactly as before
		@pairs = ([$r0, $f0]) unless @pairs;

		for (my $c = 0; $c < scalar @pairs; $c++) {
			$optbr[$i][0][$c] = $pairs[$c][0];
			$optbr[$i][1][$c] = $pairs[$c][1];
		}

	###################################################################
	#A warning flag is raised if the Tm for the primer pair is not
	#within two degrees C of each other. Computed from the closest pair,
	#as before; the printed row is recomputed from the primers the report
	#actually shows, by fullprimerchecks.

		my $pairTmdiff = $TmGCrev[$i][0][$r0][0] - $TmGCfor[$i][0][$f0][0];
		$pairTmdiff =~ s/-//;

		if ($pairTmdiff > $TM_PAIR_MAX_DIFF) {
			$primerTmwarning{$i} = 1;
		}
	}
	# print "$Tmtarget\n";
	return (\@optbr, \%primerTmwarning, $Tmtarget);
}


# Longest CONTIGUOUS perfectly-paired run that includes a 3'-terminal base of
# either strand, weighted G:C = 2 / A:T = 1. Such a run of k bases exists
# exactly when the reverse complement of one strand's k-base 3' suffix occurs
# somewhere in the other strand, which makes this a suffix search rather than a
# register scan - verified identical to a full O(n^2) scan over 8000 random
# sequence pairs, and fast enough to run on every candidate.
#
# Only the 3'-anchored run is scored. The obvious companion measure - the best
# total complementarity in any register, Primer3's SELF_ANY - was measured and
# dropped: over 565 chosen pairs it ran from 10 to 30 with a median of 18, so
# any useful threshold flagged 60-85% of all primers. Two natural 18-28mers
# always share a dozen scattered complementary bases; that is base composition,
# not a dimer, and a criterion that fires on four primers in five carries no
# selection information (the same failure the old "G or C anywhere in the last
# 5" clamp had).
sub dimer_rc {
	my $t = reverse uc($_[0] // "");
	$t =~ tr/AUTGC/TAACG/;
	return $t;
}

sub dimer_end3 {
	my ($a, $b) = @_;
	$a = uc($a // ""); $b = uc($b // "");
	$a =~ tr/U/T/;    $b =~ tr/U/T/;
	return 0 if $a eq "" or $b eq "";
	my $best = 0;
	foreach my $ord ([$a, $b], [$b, $a]) {
		my ($x, $y) = @$ord;
		my $max = (length($x) < length($y)) ? length($x) : length($y);
		$max = $DIMER_SCAN_CAP if $max > $DIMER_SCAN_CAP;
		for (my $k = $max; $k >= 1; $k--) {
			my $suf = substr($x, -$k, $k);
			next if index($y, &dimer_rc($suf)) < 0;
			my $w = 0;
			foreach my $c (split //, $suf) { $w += ($c eq 'G' or $c eq 'C') ? 2 : 1 }
			$best = $w if $w > $best;
			last;
		}
	}
	return $best;
}

# Worst 3'-anchored duplex this pair can form, over the cross-dimer and BOTH
# self-dimers. Self-dimers were not examined at all before; measured over 565
# chosen pairs they are the worst interaction for 34% of them (FW x FW 19%,
# RV x RV 16%), so a cross-only check misses a third of the real cases.
sub dimer_worst {
	my ($fw, $rv) = @_;
	my $w = &dimer_end3($fw, $rv);
	my $s = &dimer_end3($fw, $fw);  $w = $s if $s > $w;
	   $s = &dimer_end3($rv, $rv);  $w = $s if $s > $w;
	return $w;
}

sub primerdimer {

	# NOTE: this never flags anything in practice, so $W_DIMER = 3 is dead
	# weight, by the same mechanism as sub homology above - the (primer Tm
	# - 10 C) threshold asks for a near-full-length perfect cross-match
	# between the two primers. Measured: 0 hits over 1612 candidates from 60
	# fragments. See the note at the weight constants near the top.

    my @brrev = @{$_[0]};
    my @brfor = @{$_[1]};
    my @optbr = @{$_[2]};
    my @TmGCfor = @{$_[3]};
    my @TmGCrev = @{$_[4]};
    my $n = $_[5];
    my $sc = $_[6];
    my $pc = $_[7];
	my %primerdimer;

	# Cross-annealing between the two primers of a pair, evaluated per CANDIDATE
	# (candidate $ii of the reverse list pairs with candidate $ii of the forward
	# list, which is how the ranking treats them); it used to be computed for
	# candidate 0 only and stored under a fixed key, so no other candidate was
	# ever examined.
	# @brrev holds the complement in gene order, so the reverse primer as it
	# would be ordered (5'->3') is its reverse. Both strands are handed to
	# dimer_worst in that 5'->3' sense, which is what the pairing geometry
	# assumes.
	for (my $i = 0; $i < $n; $i++) {
		my $nc = (scalar @{$optbr[$i][0]} <= scalar @{$optbr[$i][1]})
		       ?  scalar @{$optbr[$i][0]} :  scalar @{$optbr[$i][1]};
		for (my $ii = 0; $ii < $nc; $ii++) {
			my $rvraw = $brrev[$i][0][$optbr[$i][0][$ii]];
			my $fw    = $brfor[$i][0][$optbr[$i][1][$ii]];
			next unless defined $rvraw and defined $fw and $rvraw ne "" and $fw ne "";
			my $rv = reverse $rvraw;
			$primerdimer{$i}{$ii} = 1 if &dimer_worst($fw, $rv) > $DIMER_MAX_END3;
		}
	}
	return \%primerdimer;
}


	
sub homology {

	# NOTE: this never flags anything in practice, so $W_HOMOLOGY = 6 is dead
	# weight. The snapback test below requires a self-complementary stretch
	# whose own Tm is within 10 C of the whole primer's, i.e. a near-full-
	# length perfect palindrome. Measured: 0 hits over 1612 candidates from 60
	# fragments; zeroing $W_HOMOLOGY leaves every primer of the test corpus
	# byte-identical. See the note at the weight constants near the top.

    my @brrev = @{$_[0]};
    my @brfor = @{$_[1]};
    my @optbr = @{$_[2]};
    my @TmGCfor = @{$_[3]};
    my @TmGCrev = @{$_[4]};
    my $n = $_[5];
    my $sc = $_[6];
    my $pc = $_[7];
	my %homologyrev;
	my %homologyfor;

	# Segments of a primer that can fold back on the primer itself ("snap-backs").
	# Evaluated per CANDIDATE binding-region length, so the ranking can tell
	# candidates apart; the fragment list is rebuilt for each candidate (it used
	# to be built once from candidate 0 and to accumulate across genes).
	my $snapback = sub {
		my ($seq, $tm) = @_;
		return 0 unless defined $seq and $seq ne "" and defined $tm;
		my @frags;
		for (my $a = 0; $a <= length $seq; $a++) {
			for (my $b = $a+5; $b <= (length $seq)-$a; $b++) {
				my $t = substr($seq, $a, $b);
				$t = reverse $t;
				$t =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
				push @frags, $t;
			}
		}
		foreach my $fr (@frags) {
			next if $fr eq "";
			if ($seq =~ m/\Q$fr\E/) {
				my $ftm = &simpleTm($fr, $sc, $pc);
				return 1 if defined $ftm and ($tm - 10) < $ftm;
			}
		}
		return 0;
	};

	for (my $i = 0; $i < $n; $i++) {
		for (my $ii = 0; $ii < scalar @{$optbr[$i][0]}; $ii++) {
			my $c = $optbr[$i][0][$ii];
			$homologyrev{$i}{$ii} = 1
				if &$snapback($brrev[$i][0][$c], $TmGCrev[$i][0][$c][0]);
		}
		for (my $ii = 0; $ii < scalar @{$optbr[$i][1]}; $ii++) {
			my $c = $optbr[$i][1][$ii];
			$homologyfor{$i}{$ii} = 1
				if &$snapback($brfor[$i][0][$c], $TmGCfor[$i][0][$c][0]);
		}
	}
	return (\%homologyrev, \%homologyfor);	
}
	

	
sub GCclamp {

    my @brrev = @{$_[0]};
    my @brfor = @{$_[1]};
    my @optbr = @{$_[2]};
    my $n = $_[3];
	my %GCclamprev;
	my %GCclampfor;

	# @brrev holds the complement in gene order, so the primer as ordered
	# (5'->3') is its reverse; gcclamp_full reads its 3'-terminal base. One
	# helper decides the clamp for both the penalty and the printed verdict.
	for (my $i = 0; $i < $n; $i++) {
		for (my $ii = 0; $ii < scalar @{$optbr[$i][0]}; $ii++) {
			my $p = reverse $brrev[$i][0][$optbr[$i][0][$ii]];
			$GCclamprev{$i}{$ii} = 1 if &gcclamp_full($p);
		}
	}

	for (my $i = 0; $i < $n; $i++) { 
		for (my $ii = 0; $ii < scalar @{$optbr[$i][1]}; $ii++) {
			$GCclampfor{$i}{$ii} = 1 if &gcclamp_full($brfor[$i][0][$optbr[$i][1][$ii]]);
		}
	}
	return (\%GCclamprev, \%GCclampfor);	
}
	

sub GCcountclamp {

    my @brrev = @{$_[0]};
    my @brfor = @{$_[1]};
    my @optbr = @{$_[2]};
    my $n = $_[3];
	my %GCcountclamprev;
	my %GCcountclampfor;
	my $count = 0;

	for (my $i = 0; $i < $n; $i++) {
		for (my $ii = 0; $ii < scalar @{$optbr[$i][0]}; $ii++) {
			$brrev[$i][0][$optbr[$i][0][$ii]] = reverse $brrev[$i][0][$optbr[$i][0][$ii]];
			my $temp = substr($brrev[$i][0][$optbr[$i][0][$ii]], length ($brrev[$i][0][$optbr[$i][0][$ii]])-5, 5);
			for (my $iii = 0; $iii < 6; $iii++) {
				if (substr ($temp, $iii, 1) =~ m/[C,G]/) {
					$count++
				}
			}
			if ($count > 3) {
				$GCcountclamprev{$i}{$ii} = 1;
			}
			$brrev[$i][0][$optbr[$i][0][$ii]] = reverse $brrev[$i][0][$optbr[$i][0][$ii]];
			$count = 0;
		}
	}

	for (my $i = 0; $i < $n; $i++) {
		for (my $ii = 0; $ii < scalar @{$optbr[$i][1]}; $ii++) {
			my $temp = substr($brfor[$i][0][$optbr[$i][1][$ii]], length ($brfor[$i][0][$optbr[$i][1][$ii]])-5, 5);
			for (my $iii = 0; $iii < length $temp; $iii++) {
				if (substr ($temp, $iii, 1) =~ m/[C,G]/) {
					$count++
				}
			}
			if ($count > 3) {
				$GCcountclampfor{$i}{$ii} = 1;
			}
			$count = 0;
		}
	}
	return (\%GCcountclamprev, \%GCcountclampfor);
}


###################################################################
# Full-primer versions of the three sequence-composition checks.
#
# The original checks (subs polyN/homology/primerdimer above) run on the BINDING
# REGION only, because that is what the binding-region-length optimiser ranks.
# But primer dimer, intra-primer homology and polyN stretches are properties of
# the molecule that is actually synthesised, which also carries the 5' USER tail
# (with its dU), any forced-overhang bridge, and any tag/linker flap. These
# helpers therefore take the COMPLETE primer string and are applied to the
# finished report by &fullprimerchecks below. The algorithms are unchanged -
# only the sequence they see is wider.
#
# Tm threshold: a self/cross structure matters when it is stable close to the
# temperature at which the primer anneals, so the threshold stays the reported
# (binding-region) Tm minus 10 C, exactly as in the original subs.

# 3'-end checks, mirroring subs GCclamp / GCcountclamp. Both look only at the
# last five bases, which are identical in the full primer and in the binding
# region (every addition is 5'), so running them on the printed primer gives
# the intended binding-region answer while automatically describing the primer
# that was actually chosen for the report.
# A GC clamp is a G or C at the 3'-TERMINAL base: that is the base the polymerase
# extends from, so 3' duplex stability is what the criterion is about. The
# previous test asked only whether any G or C appeared anywhere in the last 5
# bases, which "GAAAA" passes while ending in four A's.
# $_[0] is the primer written 5'->3'.
sub gcclamp_full {
	my $q = uc($_[0] // "");
	return 0 if length($q) < 1;
	return (substr($q, -1, 1) =~ m/[CG]/) ? 1 : 0;
}

sub gccountclamp_full {
	my $q = uc($_[0] // "");
	return 0 if length($q) < 1;
	my $t = substr($q, (length($q) >= 5 ? length($q)-5 : 0), 5);
	my $c = () = $t =~ m/[CG]/g;
	return ($c > 3) ? 1 : 0;   # "more than 3 G/C out of last 5", as the row says
}

sub polyN_full {
	my $s = uc(shift // "");
	return ($s =~ m/G{4}/ or $s =~ m/C{4}/ or $s =~ m/A{4}/ or $s =~ m/T{4}/) ? 1 : 0;
}

sub homology_full {
	my ($s, $tm, $sc, $pc) = @_;
	$s = uc($s // "");
	return 0 if $s eq "" or !defined $tm;
	my @frags;
	for (my $ii = 0; $ii <= length $s; $ii++) {
		for (my $iii = $ii+5; $iii <= (length $s)-$ii; $iii++) {
			my $t = substr($s, $ii, $iii);
			$t = reverse $t;
			$t =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
			push @frags, $t;
		}
	}
	foreach my $fr (@frags) {
		next if $fr eq "";
		if ($s =~ m/\Q$fr\E/) {
			my $ftm = &simpleTm($fr, $sc, $pc);
			return 1 if defined $ftm and ($tm - 10) < $ftm;
		}
	}
	return 0;
}

# Report-side verdict. Same measure as the ranking, and deliberately on the
# same sequence: the binding region only. The printed primer carries the USER
# tail (and any bridge or tag flap) in front of the dU, and those are shared
# boilerplate rather than a property of the design - every pair in a run would
# score alike on them - so the tail is stripped back to the last U before
# scoring. Keeping the scope identical is what lets the printed verdict and the
# penalty that drove the choice agree.
sub primerdimer_full {
	my ($fw, $rv, $tm, $sc, $pc) = @_;
	$fw = uc($fw // ""); $rv = uc($rv // "");
	return 0 if $fw eq "" or $rv eq "";
	my $f = ($fw =~ m/U([^U]*)$/) ? $1 : $fw;
	my $r = ($rv =~ m/U([^U]*)$/) ? $1 : $rv;
	return 0 if $f eq "" or $r eq "";
	return (&dimer_worst($f, $r) > $DIMER_MAX_END3) ? 1 : 0;
}

# Walk a finished construct's report array and re-decide the three full-primer
# checks from the primer sequences the report itself prints, so every report
# builder (single, linear fusion, circular) is covered in one place and the
# GenBank/CSV scraping of these same lines is untouched. A "uracil" U in a
# primer is read as the T it replaces for the purpose of base-pairing.
# Move the primer overview ("overview of the needed primers") so it sits
# directly above the primer evaluation section instead of near the top of the
# report: the list of primers and the table judging them are read together.
# Works on the report array before any output, so the HTML, the text report and
# the CSV/GenBank scrapes all see the same (re-ordered) lines.
sub move_primer_overview {
	my $arr = shift;
	my ($from, $to, $details);
	for (my $i = 0; $i < scalar @$arr; $i++) {
		$from    = $i if !defined $from    and $arr->[$i] =~ m/^overview of the needed primers/;
		$to      = $i if  defined $from and !defined $to and $arr->[$i] =~ m/^graphic overview of DNA fragments/;
		$details = $i if !defined $details and $arr->[$i] =~ m/^primer details:/;
	}
	return unless defined $from and defined $to and defined $details;
	return unless $to > $from and $details > $to;          # already in place / unexpected layout
	my @block = splice(@$arr, $from, $to - $from);          # the overview section
	$details -= scalar @block;                              # it shifted left by that much
	splice(@$arr, $details, 0, @block);
	return;
}

sub fullprimerchecks {
	my ($arr, $sc, $pc) = @_;
	my ($fw, $rv, $fwtm, $rvtm);
	my @idx;            # line indices of the three checks, per primer currently open
	my $cur = "";       # "fw" or "rv"
	my $apply = sub {
		my ($seq, $tm, $other, $iref) = @_;
		return unless defined $seq and defined $tm and @$iref;
		my $dim  = &primerdimer_full($seq, $other, $tm, $sc, $pc);
		my $hom  = &homology_full($seq, $tm, $sc, $pc);
		my $poly = &polyN_full($seq);
		my $clam = &gcclamp_full($seq);
		my $gcnt = &gccountclamp_full($seq);
		foreach my $e (@$iref) {
			my ($i, $kind) = @$e;
			my $yes = $kind eq "dimer" ? $dim : $kind eq "hom"   ? $hom
			         : $kind eq "clamp" ? $clam : $kind eq "gcnt"  ? $gcnt : $poly;
			$arr->[$i] =~ s/\.\.\.(?:YES|NO)/"..." . ($yes ? "YES" : "NO")/e;
		}
		@$iref = ();
	};
	my @pend;
	for (my $i = 0; $i < scalar @$arr; $i++) {
		my $l = $arr->[$i];
		if ($l =~ m/^primer details for/) {
			($fw, $rv, $fwtm, $rvtm, $cur) = (undef, undef, undef, undef, "");
			next;
		}
		if ($l =~ m/^forward primer:\s*(\S+)/) {
			&$apply($cur eq "rv" ? $rv : $fw, $cur eq "rv" ? $rvtm : $fwtm, $cur eq "rv" ? $fw : $rv, \@pend) if @pend;
			($fw, $fwtm, $cur) = (uc($1), undef, "fw"); $fw =~ tr/U/T/;
			next;
		}
		if ($l =~ m/^reverse primer:\s*(\S+)/) {
			&$apply($fw, $fwtm, $rv, \@pend) if @pend;
			($rv, $rvtm, $cur) = (uc($1), undef, "rv"); $rv =~ tr/U/T/;
			next;
		}
		if ($l =~ m/^\* Tm:\s*([\d.]+)C/) {
			if ($cur eq "fw") { $fwtm = $1 } elsif ($cur eq "rv") { $rvtm = $1 }
			next;
		}
		if ($l =~ m/^Tm of primers within 2C of each other\?/) {
			# the flag the builders set is computed on candidate pair 0, not on the
			# pair that was actually picked; recompute it from the two primers above.
			if (defined $fwtm and defined $rvtm) {
				my $yes = (abs($fwtm - $rvtm) <= $TM_PAIR_MAX_DIFF) ? 1 : 0;
				$arr->[$i] =~ s/\.\.\.(?:YES|NO)/"..." . ($yes ? "YES" : "NO")/e;
			}
			next;
		}
		push @pend, [$i, "dimer"] if $l =~ m/^\* risk of primer dimer formation/;
		push @pend, [$i, "hom"]   if $l =~ m/^\* risk of intra-primer homology/;
		push @pend, [$i, "poly"]  if $l =~ m/^\* presence of polyN stretches/;
		push @pend, [$i, "clamp"] if $l =~ m/^\* GC clamp present/;
		push @pend, [$i, "gcnt"]  if $l =~ m/^\* more than 3 G\/C out of last 5/;
		if ($l =~ m/^\s*$/ and @pend) {
			if ($cur eq "fw") { &$apply($fw, $fwtm, $rv, \@pend) }
			elsif ($cur eq "rv") { &$apply($rv, $rvtm, $fw, \@pend) }
		}
	}
	if (@pend) {
		if ($cur eq "fw") { &$apply($fw, $fwtm, $rv, \@pend) }
		elsif ($cur eq "rv") { &$apply($rv, $rvtm, $fw, \@pend) }
	}
	return;
}

sub polyN {

    my @brrev = @{$_[0]};
    my @brfor = @{$_[1]};
    my @optbr = @{$_[2]};
    my $n = $_[3];
	my %polyNrev;
	my %polyNfor;

	for (my $i = 0; $i < $n; $i++) {
		for (my $ii = 0; $ii < scalar @{$optbr[$i][0]}; $ii++) {
			if ($brrev[$i][0][$optbr[$i][0][$ii]] =~ m/[G]{4}/ or $brrev[$i][0][$optbr[$i][0][$ii]] =~ m/[C]{4}/ or $brrev[$i][0][$optbr[$i][0][$ii]] =~ m/[A]{4}/ or $brrev[$i][0][$optbr[$i][0][$ii]] =~ m/[T]{4}/) {
				$polyNrev{$i}{$ii} = 1;
			}
		}
	}

	for (my $i = 0; $i < $n; $i++) {
		for (my $ii = 0; $ii < scalar @{$optbr[$i][1]}; $ii++) {
			if ($brfor[$i][0][$optbr[$i][1][$ii]] =~ m/[G]{4}/ or $brfor[$i][0][$optbr[$i][1][$ii]] =~ m/[C]{4}/ or $brfor[$i][0][$optbr[$i][1][$ii]] =~ m/[A]{4}/ or $brfor[$i][0][$optbr[$i][1][$ii]] =~ m/[T]{4}/) {
				$polyNfor{$i}{$ii} = 1;
			}
		}
	}
	return (\%polyNrev, \%polyNfor);
}


sub optbranalysis {

    my @optbr = @{$_[0]};
	my @TmGCrev = @{$_[1]};
    my @TmGCfor = @{$_[2]};
    my %primerdimer = %{$_[3]};
    my %homologyrev = %{$_[4]};
    my %homologyfor = %{$_[5]};
    my %GCclamprev = %{$_[6]};
    my %GCclampfor = %{$_[7]};
    my %GCcountclamprev = %{$_[8]};
    my %GCcountclampfor = %{$_[9]};
    my %polyNrev = %{$_[10]};
    my %polyNfor = %{$_[11]};
    my $n = $_[12];   # (the %ranking argument was dropped: it was always empty here)
    my $tm_target_mode = $_[13] ? 1 : 0;   # -tm/-tm eq: the pool is already Tm-selected
    my %penalty;

	for (my $i = 0; $i < $n; $i++) {
		# both candidate lists are indexed below, so bound by the shorter one
		my $ncand = (scalar @{$optbr[$i][0]} <= scalar @{$optbr[$i][1]})
		          ?  scalar @{$optbr[$i][0]} :  scalar @{$optbr[$i][1]};
		for (my $ii = 0; $ii < $ncand; $ii++) {
			$penalty{$i}{$ii} = 0;
			unless ($tm_target_mode or ($TM_OPT_MIN <= $TmGCrev[$i][0][$optbr[$i][0][$ii]][0] and $TM_OPT_MAX >= $TmGCrev[$i][0][$optbr[$i][0][$ii]][0])) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_TM_OPT;
				#print "$i $ii Tm rev! length: $optbr[$i][0][$ii]\n";
			}
			unless ($TM_TOL_MIN <= $TmGCrev[$i][0][$optbr[$i][0][$ii]][0] and $TM_TOL_MAX >= $TmGCrev[$i][0][$optbr[$i][0][$ii]][0]) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_TM_TOL;
				#print "$i $ii Tm rev! length: $optbr[$i][0][$ii]\n";
			}
			unless ($tm_target_mode or ($TM_OPT_MIN <= $TmGCfor[$i][0][$optbr[$i][1][$ii]][0] and $TM_OPT_MAX >= $TmGCfor[$i][0][$optbr[$i][1][$ii]][0])) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_TM_OPT;   # was ++ : the only asymmetric criterion in the block
				#print "$i $ii Tm for! length: $optbr[$i][1][$ii]\n";
			}
			unless ($TM_TOL_MIN <= $TmGCfor[$i][0][$optbr[$i][1][$ii]][0] and $TM_TOL_MAX >= $TmGCfor[$i][0][$optbr[$i][1][$ii]][0]) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_TM_TOL;   # was ++
				#print "$i $ii Tm for! length: $optbr[$i][1][$ii]\n";
			}
			
			unless (39.9 < $TmGCrev[$i][0][$optbr[$i][0][$ii]][1] and 60.1 > $TmGCrev[$i][0][$optbr[$i][0][$ii]][1]) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_GC;
				#print "$i $ii GC rev! length: $optbr[$i][0][$ii]\n";
			}
			unless (39.9 < $TmGCfor[$i][0][$optbr[$i][1][$ii]][1] and 60.1 > $TmGCfor[$i][0][$optbr[$i][1][$ii]][1]) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_GC;
				#print "$i $ii GC for! length: $optbr[$i][1][$ii]\n";
			}
			
			if (exists $primerdimer{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_DIMER;
				#print "$i $ii primer dimer! length: $optbr[$i][0][$ii]\n";
			}
			if (exists $homologyrev{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_HOMOLOGY;
				#print "$i $ii hom rev! length: $optbr[$i][0][$ii]\n";
			}
			if (exists $homologyfor{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_HOMOLOGY;
				#print "$i $ii hom for! length: $optbr[$i][1][$ii]\n";
			}
			unless (exists $GCclamprev{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_CLAMP;
				#print "$i $ii GC clamp rev! length: $optbr[$i][0][$ii]\n";
			}
			unless (exists $GCclampfor{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_CLAMP;
				#print "$i $ii GC clamp for! length: $optbr[$i][1][$ii]\n";
			}
			if (exists $GCcountclamprev{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_GCCOUNT;
				#print "$i $ii GC count clamp rev! length: $optbr[$i][0][$ii]\n";
			}
			if (exists $GCcountclampfor{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_GCCOUNT;
				#print "$i $ii GC count clamp for! length: $optbr[$i][1][$ii]\n";
			}
			if (exists $polyNrev{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_POLYN;
				#print "$i $ii polyN rev! length: $optbr[$i][0][$ii]\n";
			}
			if (exists $polyNfor{$i}{$ii}) {
				$penalty{$i}{$ii} = $penalty{$i}{$ii} + $W_POLYN;
				#print "$i $ii polyN for! length: $optbr[$i][1][$ii]\n";
			}
		}
	}
	
	my %rank;
	for (my $i = 0; $i < $n; $i++) {
		my $ranking = 0;
		foreach my $ii (sort { $penalty{$i}{$a} <=> $penalty{$i}{$b} || $a <=> $b } keys %{$penalty{$i}}) {
			$rank{$i}{$ranking} = $ii;
			$ranking++;
		}
	}
	return \%rank;
}


sub htmlprint {
	
	my @outreportarray = @{$_[0]};
	my @fastaoutarray = @{$_[1]};
	my $batchinput = $_[2];
	my $x = 0;
	my $tableflag = 0;
	if ($batchinput eq "TRUE") {
		$x = $_[3];
	}
	my @htmlarray;
	my %htmlcolor;
	$htmlcolor{0} = "#99cccc";
	$htmlcolor{1} = "#cc9999";
	$htmlcolor{2} = "#99cc66";
	$htmlcolor{3} = "#6699cc";
	$htmlcolor{4} = "#ffcc66";
	$htmlcolor{5} = "#9999cc";
	$htmlcolor{6} = "#66cccc";
	$htmlcolor{7} = "#cccc99";
	my $htmlcount = 0;

	if ($dp_cassette_warning ne "") {
		my $m = $dp_cassette_warning;
		$m =~ s/&/&amp;/g; $m =~ s/</&lt;/g; $m =~ s/>/&gt;/g;
		push(@{$htmlarray[$x]}, "<div style=\"background:#ffe3e3;border:3px solid #c00000;color:#b00000;font-family:Arial,Helvetica,sans-serif;font-size:24px;font-weight:bold;padding:18px 20px;margin:0 0 18px 0;text-align:center;border-radius:6px;\">&#9888; USER cassette not found in the destination plasmid &mdash; full plasmid NOT generated.<br><span style=\"font-size:15px;font-weight:normal;\">$m. The designed primers below are unaffected; check the pasted destination plasmid sequence.</span></div>\n");
	}

	push(@{$htmlarray[$x]}, "<font size=2 face=\"courier new\">\n");

	# Canonical name -> colour map, taken from the "details of final construct"
	# segment list (the authoritative colouring), so the graphic fusion diagram
	# colours each fragment/cassette with the SAME colour as the construct view.
	# Repeated names (e.g. the cassette at both ends) reuse their first colour.
	my %namecolor;
	{
		my $ncol = scalar keys %htmlcolor;
		my $cc = 0;
		for (my $k = 0; $k < scalar @{$outreportarray[$x]}; $k++) {
			if ($outreportarray[$x][$k] =~ m/details of final construct after cloning/) {
				my $row = defined $outreportarray[$x][$k+2] ? $outreportarray[$x][$k+2] : "";
				for my $name (split(/,,/, $row)) {
					$name =~ s/^\s+//; $name =~ s/\s+$//;
					next if $name eq "" or $name eq "upstream backbone" or $name eq "downstream backbone";
					next if exists $namecolor{$name};
					$namecolor{$name} = $htmlcolor{$cc % $ncol};
					$cc++;
				}
				last;
			}
		}
	}
	my ($fusA, $fusB);   # current fusion-region fragment names (from the "joining of X and Y" headers)
	my $oligoskip = 0;   # while set, suppress the stale flap-geometry diagram of an oligo junction

	for (my $i = 0; $i < scalar @{$outreportarray[$x]}; $i++) {
		my $temp = $outreportarray[$x][$i];

		# Oligo-junction sentinels. OLIGODIAG renders a staggered sticky-end diagram
		# (two USER overhangs, core elided) coloured from the construct's name->colour
		# map; both set $oligoskip to suppress the stale flap-geometry lines up to END.
		if ($outreportarray[$x][$i] =~ m/^OLIGODIAG\t(.*)/s) {
			$oligoskip = 1;
			push(@{$htmlarray[$x]}, &render_oligo_diag($1, \%namecolor, \%htmlcolor) . "\n");
			next;
		}
		# emit an honest header+pointer and skip the
		# stale fusion diagram lines between START and END (see oligo_skip_start).
		if ($outreportarray[$x][$i] =~ m/^OLIGOSKIP_START\t(.*)/s) {
			$oligoskip = 1;
			# the payload was built without the colour map; colour-chip the bold
			# fragment/tag names it mentions, exactly as the other headers do.
			my $pl = $1;
			$pl =~ s{<b>([^<]+)</b>}{ exists $namecolor{$1} ? &name_chip($1, \%namecolor) : "<b>$1</b>" }ge;
			push(@{$htmlarray[$x]}, "$pl\n");
			next;
		}
		if ($outreportarray[$x][$i] =~ m/^OLIGOSKIP_END/) {
			$oligoskip = 0;
			next;
		}
		next if $oligoskip;
	
		if ($outreportarray[$x][$i] =~ m/(AMUSER .*)/) {
			$temp = "<b>$1</b>\n";
		}
	
		if ($outreportarray[$x][$i] =~ m/(output generated: .+)/) {
			$temp = "<b>$1</b>\n";
		}
	
		if ($outreportarray[$x][$i] =~ m/(input parameters:)/) {
			$temp = "<b>$1</b>\n";
		}
	
		if ($outreportarray[$x][$i] =~ m/(overview of the needed primers \(5'-3'\):)/) {
			$temp = "<b>$1</b>\n";
		}
	
		if ($outreportarray[$x][$i] =~ m/(primer details for .+)/) {
			$temp = "<b>$1</b>\n\n";
		}
	
		if ($outreportarray[$x][$i] =~ m/(primer details:)/) {
			push(@{$htmlarray[$x]}, "</table>") if $tableflag == 0;   # the primer overview's table, if still open
			$tableflag = 1;
			push(@{$htmlarray[$x]}, "<b>$1</b><br>\n");
			$i++;
			push(@{$htmlarray[$x]}, "$outreportarray[$x][$i]<br>\n");
			push(@{$htmlarray[$x]}, "click <a href=\"http://www.cbs.dtu.dk/services/AMUSER/instructions.php#primereval\" target=\"_blank\">here</a> for detailed descriptions of primer evaluation parameters<br><br>\n");
			# The per-primer lines that follow are rendered as one paired-column
			# table (a column pair per sequence) instead of a long list; the text
			# lines themselves are left untouched for the GenBank/CSV scraping.
			my ($ptable, $pnext) = &primer_eval_table($outreportarray[$x], $i+1, \%namecolor);
			if ($ptable ne "") {
				push(@{$htmlarray[$x]}, $ptable);
				$i = $pnext - 1;
				next;
			}
		}
	
		if ($outreportarray[$x][$i] =~ m/(\* Tm: (.+)C - in optimal range \(55-72\)\?\s+)...YES/) {
			$temp = "$1<b><font color=\"#347C17\">...YES</font></b><br>";
		}
		
		if ($outreportarray[$x][$i] =~ m/(\* Tm: (.+)C - in optimal range \(55-72\)\?\s+)...NO/) {
			$temp = "$1<b><font color=\"#C11B17\">...NO</font></b><br>";
		}
	
		if ($outreportarray[$x][$i] =~ m/(\* GC content: (.+)% - in optimal range \(40-60\)\?\s+)...YES/) {
			$temp = "$1<b><font color=\"#347C17\">...YES</font></b><br>";
		}
		if ($outreportarray[$x][$i] =~ m/(\* GC content: (.+)% - in optimal range \(40-60\)\?\s+)...NO/) {
			$temp = "$1<b><font color=\"#C11B17\">...NO</font></b><br>";
		}
		
		if ($outreportarray[$x][$i] =~ m/(\* GC clamp present at 3' end\?\s+)...YES/) {
			$temp = "$1<b><font color=\"#347C17\">...YES</font></b><br>";
		}
		if ($outreportarray[$x][$i] =~ m/(\* GC clamp present at 3' end\?\s+)...NO/) {
			$temp = "$1<b><font color=\"#C11B17\">...NO</font></b><br>";
		}
		
		if ($outreportarray[$x][$i] =~ m/(\* more than 3 G\/C out of last 5 bases at 3' end\?\s+)...NO/) {
			$temp = "$1<b><font color=\"#347C17\">...NO</font></b><br>";
		}
		if ($outreportarray[$x][$i] =~ m/(\* more than 3 G\/C out of last 5 bases at 3' end\?\s+)...YES/) {
			$temp = "$1<b><font color=\"#C11B17\">...YES</font></b><br>";
		}
		
		if ($outreportarray[$x][$i] =~ m/(\* risk of primer dimer formation in primer pair\?\s+)...NO/) {
			$temp = "$1<b><font color=\"#347C17\">...NO</font></b><br>";
		}
		if ($outreportarray[$x][$i] =~ m/(\* risk of primer dimer formation in primer pair\?\s+)...YES/) {
			$temp = "$1<b><font color=\"#C11B17\">...YES</font></b><br>";
		}
		
		if ($outreportarray[$x][$i] =~ m/(\* risk of intra-primer homology \(secondary structures\)\?\s+)...NO/) {
			$temp = "$1<b><font color=\"#347C17\">...NO</font></b><br>";
		}
		if ($outreportarray[$x][$i] =~ m/(\* risk of intra-primer homology \(secondary structures\)\?\s+)...YES/) {
			$temp = "$1<b><font color=\"#C11B17\">...YES</font></b><br>";
		}
		
		if ($outreportarray[$x][$i] =~ m/(\* presence of polyN stretches\?\s+)...NO/) {
			$temp = "$1<b><font color=\"#347C17\">...NO</font></b><br>";
		}
		if ($outreportarray[$x][$i] =~ m/(\* presence of polyN stretches\?\s+)...YES/) {
			$temp = "$1<b><font color=\"#C11B17\">...YES</font></b><br>";
		}
		
		if ($outreportarray[$x][$i] =~ m/(\Tm of primers within 2C of each other\?\s+)...YES/) {
			$temp = "$1<b><font color=\"#347C17\">...YES</font></b><br><br>";
		}
		if ($outreportarray[$x][$i] =~ m/(\Tm of primers within 2C of each other\?\s+)...NO/) {
			$temp = "$1<b><font color=\"#C11B17\">...NO</font><br><br></b>";
		}
	
		if ($outreportarray[$x][$i] =~ m/(please see "primer details" section)/) {
			$temp = "<i>$1</i>\n\n\n";
		}
	
		if ($outreportarray[$x][$i] =~ m/(overview of your final construct after cloning.*:)/) {
			if ($tableflag == 0) {
				push(@{$htmlarray[$x]}, "</table>");
			}
			$temp = "<b>$1</b><br>";
			push(@{$htmlarray[$x]}, $temp);
			$i++;
			$temp = $outreportarray[$x][$i];
			$temp =~ s/\n/<br>/g;
			push(@{$htmlarray[$x]}, $temp);
			$i++;
			my @fragmentnames = split (",,", $outreportarray[$x][$i]);
			my $html = 0;
			for (my $ii = 0; $ii < scalar @fragmentnames; $ii++) {
				unless ($fragmentnames[$ii] eq "") {
					push(@{$htmlarray[$x]}, "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$html}\">&nbsp;$fragmentnames[$ii]&nbsp;</FONT>");
					if ($html+1 < keys %htmlcolor) {
						$html++;
					}
					else {
						$html = 0;
					}
				}
			}
			push(@{$htmlarray[$x]}, "<br><br><br>\n");
			$i++;
			$temp = $outreportarray[$x][$i];
		}
	
		if ($outreportarray[$x][$i] =~ m/(details of final construct after cloning \(.*length: \d+ bases\):)/) {
			$temp = "<b>$1</b><br>";
			push(@{$htmlarray[$x]}, $temp);
			$i++;
			$temp = $outreportarray[$x][$i];
			$temp =~ s/\n/<br>/g;
			push(@{$htmlarray[$x]}, $temp);
			$i++;
			my @fragmentnames = split (",,", $outreportarray[$x][$i]);
			for (@fragmentnames) { s/^\s+//; s/\s+$//; }	#strip stray whitespace (e.g. trailing newlines on the last name)

			#Per-segment colours from the fragment names. Named backbone
			#segments ("upstream backbone"/"downstream backbone") render
			#white with black text and do NOT consume a colour-cycle slot.
			#Colours come from %namecolor (built once above from this same
			#name row), so a given name keeps ONE colour across both this
			#construct view and the graphic fusion diagram, and a repeated
			#name (e.g. the cassette at both ends) reuses its first colour.
			my $ncol = scalar keys %htmlcolor;
			my @segcolor;
			my $fallbackcc = 0;
			foreach my $nm (@fragmentnames) {
				next if ($nm eq "");
				if ($nm eq "upstream backbone" or $nm eq "downstream backbone") {
					push @segcolor, "WHITE";
				}
				elsif (exists $namecolor{$nm}) {
					push @segcolor, $namecolor{$nm};
				}
				else {
					push @segcolor, $htmlcolor{$fallbackcc % $ncol};
					$fallbackcc++;
				}
			}
			my $segc = 0;
			foreach my $nm (@fragmentnames) {
				next if ($nm eq "");
				if ($segcolor[$segc] eq "WHITE") {
					push(@{$htmlarray[$x]}, "<FONT style=\"BACKGROUND-COLOR: #ffffff; COLOR: #000000\">&nbsp;$nm&nbsp;</FONT>");
				}
				else {
					push(@{$htmlarray[$x]}, "<FONT style=\"BACKGROUND-COLOR: $segcolor[$segc]\">&nbsp;$nm&nbsp;</FONT>");
				}
				$segc++;
			}
			push(@{$htmlarray[$x]}, "<br><br><br>\n");
			$i++; $temp = $outreportarray[$x][$i];
	
			my $numbering;
			my $tempseq;
			my $totalnumbering = "";
			my $totalseq = "";
			my $segidx = 0;
			my $prevcase = "";   # case of the previous run, tracked ACROSS lines so a
			                     # segment is advanced on every real case transition
			                     # (not the old per-line "last run" heuristic, which
			                     # mis-coloured whenever a segment ended on a 60-char
			                     # line boundary).
			while ($outreportarray[$x][$i] ne "") {
				if ($outreportarray[$x][$i] =~ m/(\d*\s+)(.*)/) {
					$numbering = $1;
					$tempseq = $2;
				}
				my $temp = "";

				my @line;
				my @caseof;
				for (my $j = 0; $j < length $tempseq; $j++) {
					my $uc = "";
					while (defined $tempseq and substr($tempseq, $j, 1) =~ m/[A-Z]/) {
						$uc .= substr($tempseq, $j, 1);
						$j++;
					}
					if ($uc ne "") { push(@line, $uc); push(@caseof, 'U'); }

					my $lc = "";
					while (defined $tempseq and substr($tempseq, $j, 1) =~ m/[a-z]/) {
						$lc .= substr($tempseq, $j, 1);
						$j++;
					}
					if ($lc ne "") { push(@line, $lc); push(@caseof, 'L'); }

					$j--;
				}

				for (my $j = 0; $j < scalar @line; $j++) {
					if ($prevcase ne "" and $caseof[$j] ne $prevcase) { $segidx++; }
					$prevcase = $caseof[$j];
					my $col = defined $segcolor[$segidx] ? $segcolor[$segidx] : $htmlcolor{$segidx % $ncol};
					if ($col eq "WHITE") {
						$temp .= "<FONT style=\"BACKGROUND-COLOR: #ffffff; COLOR: #000000\">$line[$j]</FONT>";
					}
					else {
						$temp .= "<FONT style=\"BACKGROUND-COLOR: $col\">$line[$j]</FONT>";
					}
				}
				$temp = uc $temp;
				$totalseq .= "$temp<br>";
				$totalnumbering .= "$numbering<br>";

				$i++;
			}
			push(@{$htmlarray[$x]}, "<table>");
			push(@{$htmlarray[$x]}, "<tr><td><font size=2 face=\"courier new\">$totalnumbering</font></td>");
			push(@{$htmlarray[$x]}, "<td>&nbsp;&nbsp</td>");
			push(@{$htmlarray[$x]}, "<td><font size=2 face=\"courier new\">$totalseq</font></td></tr>");
			push(@{$htmlarray[$x]}, "</table>\n");
			$temp = $outreportarray[$x][$i];
		}


				
		# 		my $ustring = "";
		# 		my $lstring = "";
				
		# 		my $ii;
		# 		my $iii;

		# 		my $temptotal = $numbering;
		# 		my $case = "u";
		# 		my @templine;
		# 		my $elem = 0;
		# 		for ($ii = 0; $ii < 61; $ii++) {
		# 			if (substr($tempseq, $ii, 1) =~ m/[ATGC]/) {
		# 				while (substr($tempseq, $ii, 1) =~ m/[ATGC]/ and $ii < 61) {
		# 					@templine[$elem] .= substr($tempseq, $ii, 1);
		# 					$ii++;
		# 				}
		# 			}
		# 				for ($ii; $ii < 61; $ii++) {
		# 					$temppart .= substr($tempseq, $ii, 1);
		# 					$ii++;
		# 				}
		# 				$temptotal .= "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$html}\">$temppart</FONT>";
		# 				if ($case eq "l") {	
		# 					if ($html+1 < keys %htmlcolor) {
		# 						$html++;
		# 					}
		# 					else {
		# 						$html = 0;
		# 					}
		# 				}
		# 				$case = "u";
		# 			}



		# 			if (substr($tempseq, $ii, 1) =~ m/[atgc]/) {
		# 				my $temppart = "";
		# 				while (substr($tempseq, $ii, 1) =~ m/[atgc]/ and $ii < 61) {
		# 					$temppart .= substr($tempseq, $ii, 1);
		# 					$ii++;
		# 				}
		# 				$temptotal .= "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$html}\">$temppart</FONT>";
		# 				if ($case eq "u") {
		# 					if ($html+1 < keys %htmlcolor) {
		# 						$html++;
		# 					}
		# 					else {
		# 						$html = 0;
		# 					}
		# 				}
		# 				$case = "l";
		# 			}
		# 			$ii--;
		# 		}
		# 		$temptotal =~ s/\s\s/&nbsp;&nbsp;/g;
		# 		push(@{$htmlarray[$x]}, "$temptotal<br>");

		# 		$i++; $temp = $outreportarray[$x][$i];
		# 	}
		# }
		# 		if ($case eq "u") {
		# 			for ($ii = 0; $ii < 61; $ii++) {
		# 				if (substr($tempseq, $ii, 1) =~ m/[ATGC]/) {
		# 					$ustring .= substr($tempseq, $ii, 1);
		# 				}
		# 				elsif (substr($tempseq, $ii, 1) =~ m/[atgc]/) {
		# 					$iii = $ii;
		# 					$ii = 62;
		# 					$case = "l";
		# 				}
		# 			}
					
		# 			for ($iii = $iii; $iii < 61; $iii++) {
		# 				if (substr($tempseq, $iii, 1) =~ m/[atgc]/) {
		# 					$lstring .= substr($tempseq, $iii, 1);
		# 				}
		# 				elsif (substr($tempseq, $iii, 1) =~ m/[ATGC]/) {
		# 					$iii = 62;
		# 				}
		# 			}
					
		# 			unless ($ustring eq "") {
		# 				$temp = "$numbering<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$html}\">$ustring</FONT>";
		# 				if ($lstring ne "") {
		# 					if ($html+1 < keys %htmlcolor) {
		# 						$html++;
		# 					}
		# 					else {
		# 						$html = 0;
		# 					}
		# 				}
		# 			}
						
		# 			unless ($lstring eq "") {
		# 				$lstring = uc $lstring;
		# 				$temp .= "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$html}\">$lstring</FONT>";
		# 				$case = "l";
		# 			}
		# 		}
				
		# 		elsif ($case eq "l") {
		# 			for ($ii = 0; $ii < 61; $ii++) {
		# 				if (substr($tempseq, $ii, 1) =~ m/[atgc]/) {
		# 					$lstring .= substr($tempseq, $ii, 1);
		# 				}
		# 				elsif (substr($tempseq, $ii, 1) =~ m/[ATGC]/) {
		# 					$iii = $ii;
		# 					$ii = 62;
		# 					$case = "u";
		# 				}
		# 			}
					
		# 			for ($iii = $iii; $iii < 61; $iii++) {
		# 				if (substr($tempseq, $iii, 1) =~ m/[ATGC]/) {
		# 					$ustring .= substr($tempseq, $iii, 1);
		# 				}
		# 				elsif (substr($tempseq, $iii, 1) =~ m/[atgc]/) {
		# 					$iii = 62;
		# 				}
		# 			}
						
		# 			unless ($lstring eq "") {
		# 				$lstring = uc $lstring;
		# 				$temp = "$numbering<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$html}\">$lstring</FONT>";
		# 				if ($ustring ne "") {
		# 					if ($html+1 < keys %htmlcolor) {
		# 						$html++;
		# 					}
		# 					else {
		# 						$html = 0;
		# 					}
		# 				}
		# 			}
		# 			unless ($ustring eq "") {
		# 				$temp .= "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$html}\">$ustring</FONT>";
		# 				$case = "u";
		# 			}
		# 		}			
				
		# 		$temp =~ s/\s\s/&nbsp;&nbsp;/g;
		# 		push(@{$htmlarray[$x]}, "$temp<br>");
		# 		$i++; $temp = $outreportarray[$x][$i];
		# 	}
		# }
	
		if ($outreportarray[$x][$i] =~ m/(graphic overview of DNA fragments and primers:)/) {
			$temp = "<b>$1</b>\n";

		}

		# if ($outreportarray[$x][$i] =~ m/forward primer for (.+):/) {
		# 	$temp = "forward primer for <b>$1</b>:<br><br>";
		# }

		# if ($outreportarray[$x][$i] =~ m/reverse primer for (.+):/) {
		# 	$temp = "reverse primer for <b>$1</b>:<br><br>";
		# }		

		if ($outreportarray[$x][$i] =~ m/fragment\tstrand\ttm\tprimer/) {
			push(@{$htmlarray[$x]} ,"<table><tr><td><font size=2 face=\"courier new\"><b>fragment</b></font></td><td>&nbsp;&nbsp;&nbsp;</td><td><font size=2 face=\"courier new\"><b>strand</b></font></td><td>&nbsp;&nbsp;&nbsp;</td><td><font size=2 face=\"courier new\"><b>tm</b></font></td><td>&nbsp;&nbsp;&nbsp;</td><td><font size=2 face=\"courier new\"><b>primer</b></font></td></tr>");
			$i++;
			if ($outreportarray[$x][$i] =~ m/(.+)\t(.+)\t(.+)\t(.+)/) {
				my $fragments = "";
				my $strands = "";
				my $tms = "";
				my $primers = "";
				while ($outreportarray[$x][$i] =~ m/(.+)\t(.+)\t(.+)\t(.+)/) {
					$fragments .= "$1<br>";
					$strands .= "$2<br>";
					$tms .= "$3<br>";
					$primers .= "$4<br>";
					$i++;
				}
				push(@{$htmlarray[$x]}, "<tr><td><font size=2 face=\"courier new\">$fragments</font></td><td></td><td><font size=2 face=\"courier new\">$strands</font></td><td></td><td><font size=2 face=\"courier new\">$tms</font></td><td></td><td><font size=2 face=\"courier new\">$primers</font></td></tr>");
				# the while loop above stops ON the first non-row line, and the
				# enclosing loop's own $i++ would then skip it (it used to be a
				# blank line; since the overview was moved it is "primer details:").
				$i--;
				next;
			}
		}	

		# Oligo note: sits just after the primer table, so close the table first.
		if ($outreportarray[$x][$i] =~ m/Tags above marked "oligo"/) {
			push(@{$htmlarray[$x]}, "</table>") if $tableflag == 0;
			$tableflag = 1;
			my $note = $outreportarray[$x][$i];
			$note =~ s/^\s+//; $note =~ s/\s+$//; $note =~ s/\n/<br>/g;
			$temp = "<br><i>$note</i><br>\n";
		}

		if ($outreportarray[$x][$i] =~ m/(attention: one or more of the designed primers may have undesirable properties)/) {
			push (@{$htmlarray[$x]}, "</table>") if $tableflag == 0;
			$temp = "<br><i>$1</i>\n";
			$tableflag = 1;
		}

		# if ($outreportarray[$x][$i] =~ m/forward primer for (.+):(.+)/) {
		# 	$temp = "<tr><td><font size=2 face=\"courier new\">$1</font></td><td><font size=2 face=\"courier new\">forward</font></td><td><font size=2 face=\"courier new\">tm</font></td><td><font size=2 face=\"courier new\">$2</font></td></tr>";
		# }

		# if ($outreportarray[$x][$i] =~ m/reverse primer for (.+):(.+)/) {
		# 	$temp = "<tr><td><font size=2 face=\"courier new\">$1</font></td><td><font size=2 face=\"courier new\">reverse</font></td><td><font size=2 face=\"courier new\">tm</font></td><td><font size=2 face=\"courier new\">$2</font></td></tr>";
		# }
	
		if ($outreportarray[$x][$i] =~ m/fusion region and forward primer for joining of (.+) and (.+):/) {
			($fusA, $fusB) = ($1, $2);   # capture BEFORE name_chip, whose own regexes reset $1/$2
			$temp = "fusion region and forward primer for joining of " . &name_chip($fusA, \%namecolor) . " and " . &name_chip($fusB, \%namecolor) . ":<br><br>";
			$fusA =~ s/\s*\((?:3'|5') end\)\s*$//; $fusB =~ s/\s*\((?:3'|5') end\)\s*$//;
		}

		if ($outreportarray[$x][$i] =~ m/(\s*5'-\w+)U(.*)/) {
			$temp = $1;
			$temp .= "<b>U</b>";
			$temp .= "$2<br>";
		}

		if ($outreportarray[$x][$i] =~ m/fusion region and reverse primer for joining of (.+) and (.+):/) {
			($fusA, $fusB) = ($1, $2);   # capture BEFORE name_chip, whose own regexes reset $1/$2
			$temp = "fusion region and reverse primer for joining of " . &name_chip($fusA, \%namecolor) . " and " . &name_chip($fusB, \%namecolor) . ":<br><br>";
			$fusA =~ s/\s*\((?:3'|5') end\)\s*$//; $fusB =~ s/\s*\((?:3'|5') end\)\s*$//;
		}

		if ($outreportarray[$x][$i] =~ m/(\s+(?:\[tag\])?\w*)U(\w+-5')/) {
			$temp = $1;
			$temp .= "<b>U</b>";
			$temp .= "$2<br>";
		}
	
		if ($outreportarray[$x][$i] =~ m/5'-(\w+)\[\.\.\.\]-3'/) {
			my $tempseq = $1;
			my $tempseq1 = substr ($tempseq, 0, (length $cassette[2]) + (length $cassette[5]));
			my $tempseq2 = substr ($tempseq, (length $cassette[2]) + (length $cassette[5]), length $primerextension[$x][0][4]);
			my $tempseq3 = substr ($tempseq, (length $cassette[2]) + (length $cassette[5]) + (length $primerextension[$x][0][4]), 50);
			# Colour each fragment by its name, matching the construct view.
			# tempseq1 = upstream cassette ($fusA); tempseq3 = downstream gene ($fusB).
			my $colA = exists $namecolor{$fusA} ? $namecolor{$fusA} : $htmlcolor{0};
			my $colB = exists $namecolor{$fusB} ? $namecolor{$fusB} : $htmlcolor{0};

			$temp = "5'-";
			if ($tempseq1 ne "") {
				$temp .= "<FONT style=\"BACKGROUND-COLOR: $colA\">$tempseq1</font>";
			}
			if ($tempseq2 ne "") {
				$temp .= "<FONT style=\"BACKGROUND-COLOR: $colB\">$tempseq2</font>";
			}
			if ($tempseq3 ne "") {
				$temp .= "<FONT style=\"BACKGROUND-COLOR: $colB\">$tempseq3</font>[...]-3'<br><br>";
			}

			# unless ($nocassette eq "TRUE") {
			# 	$temp = "5'-<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$htmlcount}\">$tempseq1</font>";
			# 	if ($htmlcount+1 < keys %htmlcolor) {
			# 		$htmlcount++;
			# 	}
			# 	else {
			# 		$htmlcount = 0;
			# 	}
			# 	$temp .= "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$htmlcount}\">$tempseq3</font>[...]-3'<br><br>";
			# }
			# else {
			# 	$temp = "5'-<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$htmlcount}\">$tempseq3</font>[...]-3'<br><br>";
			# }
		}
	
		if ($outreportarray[$x][$i] =~ m/5'-\[\.\.\.\](\w+)\[\.\.\.\]-3'/) {
			my $tempseq = $1;
			# my $linker = "";
			my $tempseq1 = substr ($tempseq, 0, 50);
			# if ((length $tempseq) > 80) {
			# 	$linker = substr ($tempseq, 40, (length $tempseq) - 80);
			# }
			my $tempseq2 = "";
			if (length $tempseq > 100) {
				$tempseq2 = substr ($tempseq, 50, (length $tempseq)-100);				
			}
			my $tempseq3 = substr ($tempseq, (length $tempseq)-50, 50);
			# tempseq1 = upstream gene ($fusA); tempseq3 = downstream gene ($fusB).
			my $colA = exists $namecolor{$fusA} ? $namecolor{$fusA} : $htmlcolor{0};
			my $colB = exists $namecolor{$fusB} ? $namecolor{$fusB} : $htmlcolor{0};
			$temp = "5'-[...]<FONT style=\"BACKGROUND-COLOR: $colA\">$tempseq1</font>";
			unless ($tempseq2 eq "") {
				$temp .= "<FONT style=\"BACKGROUND-COLOR: $colA\">$tempseq2</font>";
			}
			$temp .= "<FONT style=\"BACKGROUND-COLOR: $colB\">$tempseq3</font>[...]-3'<br>";
		}
	
		if ($outreportarray[$x][$i] =~ m/5'-\[\.\.\.\](\w+)-3'/) {
			my $tempseq = $1;
			my $tempseq1 = substr ($tempseq, 0, 50);
			my $tempseq2 = substr ($tempseq, 50, length $primerextension[$x][-1][1]);
			my $tempseq3 = substr ($tempseq, 50 + (length $primerextension[$x][-1][1]), (length $cassette[3]) + (length $cassette[6]));
			# tempseq1 = upstream gene ($fusA); tempseq3 = downstream cassette ($fusB).
			my $colA = exists $namecolor{$fusA} ? $namecolor{$fusA} : $htmlcolor{0};
			my $colB = exists $namecolor{$fusB} ? $namecolor{$fusB} : $htmlcolor{0};

			$temp = "5'-[...]";
			if ($tempseq1 ne "") {
				$temp .= "<FONT style=\"BACKGROUND-COLOR: $colA\">$tempseq1</font>";
			}
			if ($tempseq2 ne "") {
				$temp .= "<FONT style=\"BACKGROUND-COLOR: $colA\">$tempseq2</font>";
			}
			if ($tempseq3 ne "") {
				$temp .= "<FONT style=\"BACKGROUND-COLOR: $colB\">$tempseq3</font>";
			}
			$temp .= "-3'<br>";


			# my $tempseq = $1;
			# my $tempseq1 = substr ($tempseq, 0, 40);
			# # if ((length $tempseq) > 80) {
			# # 	$linker = substr ($tempseq, 40, (length $tempseq) - 80);
			# # }
			# my $tempseq2 = "";
			# if (length $tempseq > 80) {
			# 	$tempseq2 = substr ($tempseq, 40, (length $tempseq)-80);				
			# }
			# my $tempseq3 = substr ($tempseq, (length $tempseq)-40, 40);
			# $temp = "5'-[...]<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$htmlcount}\">$tempseq1</font>";
			# if ($htmlcount+1 < keys %htmlcolor) {
			# 	$htmlcount++;
			# }
			# else {
			# 	$htmlcount = 0;
			# }
			# # unless ($linker eq "") {
			# # 	$temp .= "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$htmlcount}\">$linker</font>";
			# # 	if ($htmlcount+1 < keys %htmlcolor) {
			# # 		$htmlcount++;
			# # 	}
			# # 	else {
			# # 		$htmlcount = 0;
			# # 	}
			# # }
			# unless ($tempseq2 eq "") {
			# 	$temp .= "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$htmlcount}\">$tempseq2</font>";
			# 	if ($htmlcount+1 < keys %htmlcolor) {
			# 		$htmlcount++;
			# 	}
			# 	else {
			# 		$htmlcount = 0;
			# 	}
			# }
			# $temp .= "<FONT style=\"BACKGROUND-COLOR: $htmlcolor{$htmlcount}\">$tempseq3</font>[...]-3'<br>";
		}

		if ($outreportarray[$x][$i] =~ m/fusion region and related primers for joining of (.+) and (.+):/) {
			($fusA, $fusB) = ($1, $2);   # capture BEFORE name_chip, whose own regexes reset $1/$2
			$temp = "fusion region and related primers for joining of " . &name_chip($fusA, \%namecolor) . " and " . &name_chip($fusB, \%namecolor) . ":<br><br>";
			$fusA =~ s/\s*\((?:3'|5') end\)\s*$//; $fusB =~ s/\s*\((?:3'|5') end\)\s*$//;
		}

	
		$temp =~ s/\n/<br>/g;
		$temp =~ s/\s\s/&nbsp;&nbsp;/g;
		push(@{$htmlarray[$x]},  "$temp\n");
	}

	# unless ($batchinput eq "TRUE") {
	# 	push(@{$htmlarray[$x]},  "<br><br><b>primers in FASTA format (5'-3')</b><br>");
	# 	push(@{$htmlarray[$x]},  "-------------------------------<br><br>");
	# 	for (my $i = 0; $i < scalar @{$fastaoutarray[0]}; $i++) {
	# 		my $temp = $fastaoutarray[$x][$i];
	# 		$temp =~ s/\n/<br>/g;
	# 		push(@{$htmlarray[$x]},  "$temp");
	# 	}
	# }

	return \@htmlarray;
}

sub htmlprintbatch {
	my @fastaoutarray = @{$_[0]};
	my @fastabatcharray = @{$_[1]};
	my @htmlarray;
	
	push (@htmlarray, "<font size=2 face=\"courier new\">\n");
	push (@htmlarray, "--------------------------------------<br>\n");
	push (@htmlarray, "<b>AMUSER $version</b><br>\n");
	push (@htmlarray, "<b>output generated: $mday-$mon-$year $hour:$min:$sec</b><br>\n");
	push (@htmlarray, "--------------------------------------<br><br>\n");
	push (@htmlarray, "<b>input parameters:</b><br>\n");
	push (@htmlarray, "<b>-----------------</b><br>\n");
	push (@htmlarray, "number of batches: $numberofinputs<br>\n");
	$cassetteprint =~ s/\n/<br>/g;
	push (@htmlarray, "$cassetteprint<br><br><br>\n");
	push (@htmlarray, "<b>overview of the needed primers (5'-3'):</b><br>");
	push (@htmlarray, "<b>---------------------------------------</b><br><br>");

	push (@htmlarray, "<table>\n");
	push (@htmlarray, "<tr><td><b><font size=2 face=\"courier new\">description</b></td><td>&nbsp;&nbsp;&nbsp;</td><td><b><font size=2 face=\"courier new\">name</b></td><td>&nbsp;&nbsp;&nbsp;</td><td><b><font size=2 face=\"courier new\">strand</b></td><td>&nbsp;&nbsp;&nbsp;</td><td><b><font size=2 face=\"courier new\">tm</b></td><td>&nbsp;&nbsp;&nbsp;</td><td><b><font size=2 face=\"courier new\">primer</b></td><td>&nbsp;&nbsp;&nbsp;</td></td><td><b><font size=2 face=\"courier new\">details for batch</b></td></td></tr>\n");	


	# my $temp = $fastabatcharray[0][0];
	# if ($temp =~ m/(.+)\t(.+)\t(.+)\t(.+)\t(.+)/) {
	# 	my $fragment = $1;
	# 	my $strand = $2;
	# 	my $tm = $3;
	# 	my $primer = $4;
	# 	push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">backbone</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td></td></tr>\n");	
	# }
	# my $temp = $fastabatcharray[0][1];
	# if ($temp =~ m/(.+)\t(.+)\t(.+)\t(.+)\t(.+)/) {
	# 	my $fragment = $1;
	# 	my $strand = $2;
	# 	my $tm = $3;
	# 	my $primer = $4;
	# 	push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">backbone</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td></td></tr>\n");	
	# }

	# for (my $i = 2; $i < scalar @{$fastabatcharray[0]}; $i++) {
	# 	my $temp = $fastabatcharray[0][$i];
	# 	$temp =~ s/\n/<br>/g;
	# 	push (@htmlarray, "$temp\n");
	# }
	# for (my $i = 0; $i < scalar @inputs; $i++) {
	# 	my $batchno = $i+1;
	# 	# push (@htmlarray, "<b>forward and reverse primers for insert in batch $batchno ");
	# 	# push (@htmlarray, "<a href=\"batch$i.report.html\" target=\"_blank\">details</a>");
	# 	push (@htmlarray, "</b>");
		# my $i = 0;
		my $batchnumber = 0;
		my $tempcounter = 0;
		for (my $ii = 0; $ii < scalar @fastabatcharray; $ii++) {
			my $temp = $fastabatcharray[$ii];
			#$temp =~ s/\n/<br>/g;
			if ($temp =~ m/(.+)\t(.+)\t(.+)\t(.+)\t(.+)/) {
				my $batch = $1;
				my $fragment = $2;
				my $strand = $3;
				my $tm = $4;
				my $primer = $5;
				my $batchnumberprint = $batchnumber+1;
				if ($batch =~ m/backbone/) {
					push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">$batch</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td></td></tr>\n");	
				}
				else {
					push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">$batch</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td><font size=2 face=\"courier new\"><a href=\"$outputurl/batch$batchnumber.report.html\" target=\"_blank\">batch $batchnumberprint</a></td></tr>\n");
					if ($tempcounter == 1) {
						$batchnumber++;
						$tempcounter = 0;
					}
					else {
						$tempcounter++;
					}
				}



				# if ($ii < 2) {
				# 	push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">$batch</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td></td></tr>\n");	
				# }
				# elsif ($ii > 1 and $ii < (scalar @fastabatcharray)-2) {
				# 	if ($ii % 2 == 0) {
				# 		my $batchnumberprint = $batchnumber+1;
				# 		push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">$batch</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td><font size=2 face=\"courier new\"><a href=\"./batch$batchnumber.report.html\" target=\"_blank\">batch $batchnumberprint</a></td></tr>\n");
				# 		$batchnumber++;
				# 	}
				# 	else {
				# 		push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">$batch</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td></td></tr>\n");
				# 	}
				# }
				# elsif ($ii > (scalar @fastabatcharray)-3) {
				# 	if ($backbone == 2) {
				# 		if ($ii % 2 == 0) {
				# 			my $batchnumberprint = $batchnumber+1;
				# 			push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">$batch</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td><font size=2 face=\"courier new\"><a href=\"./batch$batchnumber.report.html\" target=\"_blank\">batch $batchnumberprint</a></td></tr>\n");		
				# 		}
				# 	}	
				# 	push (@htmlarray, "<tr><td><font size=2 face=\"courier new\">$batch</td><td></td><td><font size=2 face=\"courier new\">$fragment</td><td></td><td><font size=2 face=\"courier new\">$strand</td><td></td><td><font size=2 face=\"courier new\">$tm</td><td></td><td><font size=2 face=\"courier new\">$primer</td><td></td><td></td></tr>\n");	
				# }
			}
		}
		#push (@htmlarray, "<br>");
	# }
	push (@htmlarray, "</table>\n");
	return \@htmlarray;
}




###################################################################
#####          FULL-PLASMID RECONSTRUCTION SUBROUTINES        #####
###################################################################


#Reverse-complement a nucleotide string (case preserved, IUPAC aware).
sub revcomp {
	my ($s) = @_;
	$s = reverse $s;
	$s =~ tr/ACGTURYKMBVDHNacgturykmbvdhn/TGCAAYRMKVBHDNtgcaayrmkvbhdn/;
	return $s;
}

###################################################################
# Tag-as-oligo design.
# A tag flagged "clone as oligo" is realised as a short synthetic
# dsDNA (ordered as FW + RV strands that anneal) inserted into the
# USER fusion, instead of riding as a flap on a gene primer.
#
# Overhangs (default 8 nt) are taken from the tag's own ends (the
# engineered 5' A / 3' T give the deoxyuridine positions on the
# flanking gene primers), so the same oligo is reusable across genes.
# When an end abuts the fixed USER cassette instead of a gene, that
# end carries the cassette's fixed overhang scar instead.
#
# Args: $tag (sense DNA, incl. any ATG/stop codon), $five_src,
#       $three_src  (each either the string 'tag' for a gene-facing,
#       tag-derived overhang, or a cassette scar sense sequence for a
#       cassette-facing end), $ovh (overhang length, default 8).
# Returns a hashref:
#   oligo_fw / oligo_rv : the two strands to order (5'->3')
#   tail3 : FW-primer USER tail (sense, dU) for a gene on the 3' side
#   tail5 : RV-primer USER tail (antisense, dU) for a gene on the 5' side
#   br5 / br3 : tag bridge bases (sense) outside the 5' / 3' overhang
#   o5 / o3 : the 5' / 3' overhang sequences (sense)
###################################################################
sub compute_oligo {
	my ($tag, $five_src, $three_src, $ovh) = @_;
	$ovh = 8 unless defined $ovh and $ovh > 0;
	$tag = uc $tag;
	my $a = index($tag, 'A');
	my $t = rindex($tag, 'T');
	my $L5 = ($five_src  eq 'tag') ? '' : uc $five_src;
	my $L3 = ($three_src eq 'tag') ? '' : uc $three_src;
	my $scaf = $L5 . $tag . $L3;
	my $ovh5s = ($five_src  eq 'tag') ? (length($L5) + $a)     : 0;
	my $ovh3e = ($three_src eq 'tag') ? (length($L5) + $t + 1) : length($scaf);
	my $o5   = substr($scaf, $ovh5s, $ovh);
	my $o3   = substr($scaf, $ovh3e - $ovh, $ovh);
	my $core = substr($scaf, $ovh5s + $ovh, $ovh3e - $ovh - ($ovh5s + $ovh));
	my $oligo_fw = $core . $o3;                 # sense strand
	my $oligo_rv = &revcomp($o5 . $core);       # antisense strand
	my $tail3 = substr($o3, 0, $ovh - 1) . 'U';             # dU replaces the 3' overhang's terminal T (sense)
	my $tr    = &revcomp($o5);
	my $tail5 = substr($tr, 0, $ovh - 1) . 'U';             # dU on the antisense (= A on sense at the 5' end)
	my $br5 = ($five_src  eq 'tag') ? substr($tag, 0, $a)  : '';
	my $br3 = ($three_src eq 'tag') ? substr($tag, $t + 1) : '';
	return { oligo_fw=>$oligo_fw, oligo_rv=>$oligo_rv, tail3=>$tail3, tail5=>$tail5,
	         br5=>$br5, br3=>$br3, o5=>$o5, o3=>$o3, core=>$core };
}

# The two cassette overhang scars (sense), derived from the fusion tails.
# Auto-fragmentation cut chooser. Given a sequence and a target piece length,
# return a sorted list of cut positions that split it into ~target-bp pieces,
# with each cut placed at the CENTRE of a native USER fusion-tail span (an A..T
# span of 7-9 nt, ideally; up to 13 as a fallback). Because the cut sits at the
# span centre, the normal pipeline's ATpos/ATscore re-finds that overhang with a
# near-zero position penalty at the junction, and rmclashft guarantees the
# internal overhangs are unique. Splitting then fusing regenerates the original
# sequence scarlessly. Cuts are kept monotonic and at least $minpiece from each
# end and from each other.
sub overhang_cuts {
	my ($seq, $target) = @_;
	$target = 1000 unless (defined $target && $target =~ /^\d+$/ && $target > 0);
	my $L = length $seq;
	my $N = int($L / $target + 0.5); $N = 2 if $N < 2;
	return () if $L < 2;
	my @s = split //, uc $seq;
	my @spans;                       # [centre, len]
	for (my $a = 0; $a < $L; $a++) {
		next unless $s[$a] eq 'A';
		for (my $len = 7; $len <= 13; $len++) {
			my $t = $a + $len - 1; last if $t >= $L;
			next unless $s[$t] eq 'T';
			next if substr($seq, $a, $len) =~ /[^ACGTacgt]/;   # overhangs must be non-degenerate
			push @spans, [ ($a + $t) / 2, $len ];
		}
	}
	my @cuts; my $minpiece = 100;
	for (my $k = 1; $k < $N; $k++) {
		my $ideal = $k * $L / $N;
		my ($best, $bestscore);
		for my $sp (@spans) {
			my $cut = int($sp->[0] + 0.5);
			next if $cut < $minpiece or $cut > $L - $minpiece;
			next if (@cuts and $cut <= $cuts[-1] + $minpiece);
			my $lenpen = ($sp->[1] >= 7 && $sp->[1] <= 9) ? 0 : ($sp->[1] - 9) * 3;
			my $score  = abs($sp->[0] - $ideal) + $lenpen;
			if (!defined $bestscore or $score < $bestscore) { $bestscore = $score; $best = $cut; }
		}
		if (defined $best) { push @cuts, $best; }
		else {
			my $cut = int($ideal + 0.5);
			$cut = $cuts[-1] + $minpiece if (@cuts and $cut <= $cuts[-1]);
			next if $cut > $L - $minpiece;
			push @cuts, $cut;
			warn "WARNING: no clean USER overhang span found near position $cut during auto-fragmentation; cutting there anyway (the fusion primer at this junction may be suboptimal).\n";
		}
	}
	return @cuts;
}

sub cassette_scar5 { my $s = uc $_[0]; $s =~ tr/U/T/; return $s; }                 # from $cassette[2]
sub cassette_scar3 { my $s = uc $_[0]; $s =~ tr/U/T/; return &revcomp($s); }       # from $cassette[3]

# Build the "OLIGOSKIP_START" sentinel for a fusion junction that has one or
# more tags ordered as separate oligos. ($left, $right, <names...>); undef/""
# names are ignored. The stale flap-geometry fusion diagram for this junction
# would show the gene annealing directly to its neighbour/cassette (with the
# cassette overhang, not the real tag-derived one), which is wrong for an oligo.
# So htmlprint prints this honest header + pointer instead and SKIPS the stale
# diagram lines up to the matching OLIGOSKIP_END. The correct oligo primers and
# strands are in the primer section; the full sequence is in the construct view.
sub oligo_skip_start {
	my ($left, $right) = (shift, shift);
	my @names = grep { defined $_ and $_ ne "" } @_;
	my $list  = join(" + ", @names);
	my $chain = join(" - ", $left, @names, $right);
	my $noun  = (@names > 1) ? "are" : "is";
	my $blist = join(" + ", map { "<b>$_</b>" } @names);   # each name bold on its own so htmlprint can colour-chip it
	my $html  = "fusion of <b>$left</b> to <b>$right</b> through the $blist oligo:<br><br>"
	          . "<i>the $list $noun supplied as a separate pre-annealed oligo (forward/reverse strands "
	          . "listed in the primer section). The flanking primer(s) carry a USER overhang to the oligo, "
	          . "not directly to the other partner, so the actual order is $chain. This oligo fusion is not "
	          . "drawn as a diagram here - see the oligo strands and the final construct for the exact "
	          . "sequence.</i><br><br>";
	return "OLIGOSKIP_START\t$html";
}

# Build an OLIGODIAG sentinel line (parsed by render_oligo_diag). Fields are the
# junction, left->right in sense: neighbour labels, the tag/cassette/gene names
# used for colour lookup, the left context (+trunc flag) and its kind, an optional
# tag-remainder bridge on each side, the two 8-nt overhangs and their kinds, and
# the right context (+trunc flag) and its kind. See the emit sites.
sub oligodiag_line {
	my ($end,$label1,$label4,$tagname,$cassettename,$genename,
	    $ltr,$Lctx,$Lbr,$Lovh,$LovhKind,$Rovh,$RovhKind,$Rbr,$Rctx,$rtr,
	    $LctxKind,$RctxKind,$core) = @_;
	return join("\t", "OLIGODIAG", $end, $label1, $label4, $tagname, $cassettename, $genename,
		$ltr, $Lctx, $Lbr, $Lovh, $LovhKind, $Rovh, $RovhKind, $Rbr, $Rctx, $rtr,
		$LctxKind, $RctxKind, (defined $core ? $core : "")) . "\n";
}

# Trim a flanking context string to the last/first $len bases; returns (trimmed, trunc).
sub oligodiag_ctx { my ($s,$len,$tail) = @_; $s //= "";
	return (length($s) > $len) ? ($tail ? substr($s,-$len) : substr($s,0,$len), 1) : ($s, 0); }

# Render the oligo fusion diagram from an OLIGODIAG sentinel. The tag is ordered
# as a pre-annealed oligo; we draw the assembled junction as a staggered sticky-
# end duplex showing ONLY the two 8-nt USER overhangs (the double-stranded core
# is elided as "[...]" to save width). Four strands, 5'->3' sense left to right:
#   line1  left-neighbour top        5'- ctx + outer-overhang -3'
#   line2  oligo antisense (RV)      3'- comp(outer-ovh) + [...] -5'
#   line3  oligo sense (FW)          5'- [...] + gene-facing-ovh -3'
#   line4  right-neighbour bottom    3'- comp(gene-facing-ovh) + comp(ctx) -5'
# Overhang/context segments are coloured from the construct's name->colour map so
# the diagram matches the "details of final construct" view. All alignment uses
# &nbsp; (this bypasses the later whitespace pass). See the OLIGODIAG emit sites.
# Colour-chip a fragment/cassette/tag NAME exactly as the "details of final
# construct" legend does, so a name mentioned in a fusion-region header carries
# the same colour block as its sequence further down the report. A trailing
# " (3' end)"/" (5' end)" qualifier is kept outside the chip. A name that is not
# in the construct's name->colour map (so has no sequence colour) stays bold-only.
# Blend a #rrggbb colour halfway to white, for the FW/RV sub-header of a
# sequence (its own construct colour, one shade lighter).
sub lighter_colour {
	my $c = shift // "";
	return "#e8e8e8" unless $c =~ m/^#([0-9a-fA-F]{6})$/;
	my @v = map { hex } (substr($1,0,2), substr($1,2,2), substr($1,4,2));
	return sprintf("#%02x%02x%02x", map { int(($_ + 255) / 2) } @v);
}

# Render the primer evaluation section as one paired-column table: a column
# pair per sequence (FW, RV), a row per evaluation parameter. Reads the same
# "primer details" lines the text report builds (so the GenBank/CSV scraping of
# those lines is untouched) starting at $start, and returns the HTML plus the
# index of the first line after the section.
#
# Verdict colours carry the meaning the list form used: the parameter decides which
# answer is the good one (a GC clamp is wanted, a polyN stretch is not), shown
# as a light green / light red cell fill.
sub primer_eval_table {
	my ($arr, $start, $ncref) = @_;
	my %good = (                    # answer that counts as good, per parameter
		"tm"    => "YES", "gc"    => "YES", "clamp" => "YES",
		"gcnt"  => "NO",  "dimer" => "NO",  "hom"   => "NO", "poly" => "NO",
	);
	my %label = (
		"tm"    => "Tm in optimal range (55-72)?",
		"gc"    => "GC content in optimal range (40-60)?",
		"clamp" => "GC clamp present at 3' end?",
		"gcnt"  => "more than 3 G/C out of last 5 bases at 3' end?",
		"dimer" => "risk of primer dimer formation in primer pair?",
		"hom"   => "risk of intra-primer homology (secondary structures)?",
		"poly"  => "presence of polyN stretches?",
	);
	my @roworder = qw(tm gc clamp gcnt dimer hom poly);

	my (@cols, $cur, $dir);
	my $i = $start;
	for (; $i < scalar @$arr; $i++) {
		my $l = $arr->[$i];
		last if $l =~ m/^(graphic overview|details of final construct|overview of the needed primers|output generated|input parameters)/;
		if ($l =~ m/^\s*primer details for (.+?)\s*$/) {
			$cur = { name => $1, fw => {}, rv => {}, pair => "" };
			push @cols, $cur; $dir = "";
			next;
		}
		next unless defined $cur;
		if ($l =~ m/^forward primer:/) { $dir = "fw"; next; }
		if ($l =~ m/^reverse primer:/) { $dir = "rv"; next; }
		if ($l =~ m/^Tm of primers within 2C of each other\?\s*\.\.\.(YES|NO)/) { $cur->{pair} = $1; next; }
		next unless $dir;
		my $t = $cur->{$dir};
		if    ($l =~ m/^\* Tm: ([\d.]+)C - in optimal range \(55-72\)\?\s*\.\.\.(YES|NO)/)          { $t->{tm}    = [$1, $2] }
		elsif ($l =~ m/^\* GC content: ([\d.]+)% - in optimal range \(40-60\)\?\s*\.\.\.(YES|NO)/)  { $t->{gc}    = [$1, $2] }
		elsif ($l =~ m/^\* GC clamp present at 3' end\?\s*\.\.\.(YES|NO)/)                          { $t->{clamp} = [$1, $1] }
		elsif ($l =~ m/^\* more than 3 G\/C out of last 5 bases at 3' end\?\s*\.\.\.(YES|NO)/)      { $t->{gcnt}  = [$1, $1] }
		elsif ($l =~ m/^\* risk of primer dimer formation in primer pair\?\s*\.\.\.(YES|NO)/)       { $t->{dimer} = [$1, $1] }
		elsif ($l =~ m/^\* risk of intra-primer homology \(secondary structures\)\?\s*\.\.\.(YES|NO)/) { $t->{hom} = [$1, $1] }
		elsif ($l =~ m/^\* presence of polyN stretches\?\s*\.\.\.(YES|NO)/)                         { $t->{poly}  = [$1, $1] }
	}
	return ("", $start+1) unless @cols;

	my $cellbase = "padding:5px 10px;text-align:center;border:1px solid #ccc;font-weight:normal;";
	my $okfill   = "background:#d8f0d0;color:#1e5c12;";
	my $badfill  = "background:#f9d6d3;color:#9c1511;";
	my $gap      = "border-left:8px solid #ffffff;";

	my $h = "<table style=\"border-collapse:collapse;font-family:'courier new',monospace;font-size:13px;margin:4px 0 10px 0\">\n";
	# header row 1: sequence names, in their construct colour
	$h .= "<tr><td style=\"border:none;background:#ffffff\" rowspan=\"2\"></td>";
	foreach my $c (@cols) {
		my $col = (exists $ncref->{$c->{name}}) ? $ncref->{$c->{name}} : "#e8e8e8";
		$h .= "<td colspan=\"2\" style=\"$cellbase$gap" . "background:$col\">" . $c->{name} . "</td>";
	}
	$h .= "</tr>\n<tr>";
	foreach my $c (@cols) {
		my $col = (exists $ncref->{$c->{name}}) ? $ncref->{$c->{name}} : "#e8e8e8";
		my $lt  = &lighter_colour($col);
		$h .= "<td style=\"$cellbase$gap" . "background:$lt;font-weight:bold\">FW</td>";
		$h .= "<td style=\"$cellbase" . "background:$lt;font-weight:bold\">RV</td>";
	}
	$h .= "</tr>\n";
	# one row per parameter
	foreach my $k (@roworder) {
		next unless grep { exists $_->{fw}{$k} or exists $_->{rv}{$k} } @cols;
		$h .= "<tr><td style=\"$cellbase" . "background:#f2f2f2;text-align:left;white-space:nowrap\">$label{$k}</td>";
		foreach my $c (@cols) {
			foreach my $d ("fw", "rv") {
				my $e = $c->{$d}{$k};
				my $fill = "";
				my $txt  = "&nbsp;";
				if (defined $e) {
					$txt  = $e->[0];
					$fill = ($e->[1] eq $good{$k}) ? $okfill : $badfill;
				}
				$h .= "<td style=\"$cellbase$fill" . ($d eq "fw" ? $gap : "") . "\">$txt</td>";
			}
		}
		$h .= "</tr>\n";
	}
	# pair-level row: one verdict per sequence, spanning its two columns
	if (grep { $_->{pair} ne "" } @cols) {
		$h .= "<tr><td style=\"$cellbase" . "background:#f2f2f2;text-align:left;white-space:nowrap\">Tm of primers within 2C of each other?</td>";
		foreach my $c (@cols) {
			my $fill = $c->{pair} eq "" ? "" : ($c->{pair} eq "YES" ? $okfill : $badfill);
			my $txt  = $c->{pair} eq "" ? "&nbsp;" : $c->{pair};
			$h .= "<td colspan=\"2\" style=\"$cellbase$fill$gap\">$txt</td>";
		}
		$h .= "</tr>\n";
	}
	$h .= "</table>\n";
	return ($h, $i);
}

sub name_chip {
	my ($name, $ncref) = @_;
	return "" unless defined $name;
	my $suffix = "";
	if ($name =~ s/(\s*\((?:3'|5') end\))\s*$//) { $suffix = $1; }
	$name =~ s/^\s+//; $name =~ s/\s+$//;
	my $col = (defined $ncref and exists $ncref->{$name}) ? $ncref->{$name} : "";
	# inside a chip the colour block marks the name (as in the construct legend),
	# so the text is NOT bold; an uncoloured name keeps bold as its only emphasis.
	my $chip = ($col eq "")
		? "<b>$name</b>"
		: "<FONT style=\"BACKGROUND-COLOR: $col\">&nbsp;$name&nbsp;</FONT>";
	return $chip . $suffix;
}

sub render_oligo_diag {
	my ($payload, $ncref, $hcref) = @_;
	my %nc = %$ncref; my %hc = %$hcref;
	$payload =~ s/\n\z//;   # the sentinel line's trailing newline is captured by the /s handler regex
	my ($end,$label1,$label4,$tagname,$cassettename,$genename,
	    $ltrunc,$Lctx,$Lbr,$Lovh,$LovhKind,$Rovh,$RovhKind,$Rbr,$Rctx,$rtrunc,
	    $LctxKind,$RctxKind,$core) = split(/\t/, $payload, -1);

	my $casc  = exists $nc{$cassettename} ? $nc{$cassettename} : $hc{0};
	my $genec = exists $nc{$genename}     ? $nc{$genename}     : $hc{0};
	my $tagc  = exists $nc{$tagname}      ? $nc{$tagname}      : $hc{1};
	my %kc = (cassette=>$casc, gene=>$genec, tag=>$tagc);
	# Overhangs are only ever cassette- or tag-coloured. A context kind is 'cassette',
	# 'tag', or a fragment NAME to look up directly (so two different genes flanking an
	# internal junction get their own construct colours).
	my $ctxcol = sub { my $k = shift;
		return $casc if $k eq 'cassette'; return $tagc if $k eq 'tag';
		return (exists $nc{$k} ? $nc{$k} : $hc{0}); };
	my $cLctx = $ctxcol->($LctxKind);
	my $cLovh = $kc{$LovhKind} // $hc{0};
	my $cRovh = $kc{$RovhKind} // $hc{0};
	my $cRctx = $ctxcol->($RctxKind);

	my $M  = "[...]"; my $ML = length $M;
	my $sp = sub { my ($s,$c) = @_; return "" if !defined $s || $s eq ""; return "<FONT style=\"BACKGROUND-COLOR: $c\">$s</FONT>"; };
	my $nb = sub { return "&nbsp;" x $_[0]; };
	my $cp = sub { my $s = uc($_[0] // ""); $s =~ tr/ACGT/TGCA/; return $s; };   # complement (no reverse)

	# Show the two annealed oligo strands with their overhangs plus a few core bases
	# next to each overhang, so the double-stranded (annealed) core and the single-
	# stranded overhangs are both visible; the rest of the core is elided as "[...]".
	$core = uc($core // "");
	my $cL   = length $core;
	my $show = 10;                                     # core bases shown next to each overhang
	my ($leftn, $rightn, $mid);
	if ($cL > 2*$show) { ($leftn,$rightn,$mid) = ($show, $show, $M); }       # long core: 10 + [...] + 10
	else               { $leftn = int($cL/2); $rightn = $cL - $leftn; $mid = ""; }  # short core: show it all, no [...]
	my $midlen = length $mid;
	my $W = $leftn + $midlen + $rightn;                        # core-display width (columns), same on both strands
	my $coreL = substr($core, 0, $leftn);
	my $coreR = substr($core, $cL - $rightn, $rightn);
	# Both oligo strands show BOTH core stretches, vertically paired (the annealed
	# double-stranded core); the overhangs stay single-stranded on their own strand.
	my $coreFW = $sp->($coreL,$tagc) . $mid . $sp->($coreR,$tagc);                 # top strand, 5'->3'
	my $coreRV = $sp->($cp->($coreL),$tagc) . $mid . $sp->($cp->($coreR),$tagc);   # bottom strand, 3'->5'

	# Columns before the outer overhang (sense), counting the optional left-trunc
	# marker as $ML literal columns so every strand stays in register.
	my $P = ($ltrunc ? $ML : 0) + length($Lctx) + length($Lbr);

	my $l1 = "5'-" . ($ltrunc ? $M : "") . $sp->($Lctx,$cLctx) . $sp->($Lbr,$tagc) . $sp->($Lovh,$cLovh) . "-3'";
	my $l2 = $nb->($P)       . "3'-" . $sp->($cp->($Lovh),$cLovh) . $coreRV . "-5'";
	my $l3 = $nb->($P + 8)   . "5'-" . $coreFW . $sp->($Rovh,$cRovh) . "-3'";
	my $l4 = $nb->($P + 8 + $W) . "3'-" . $sp->($cp->($Rovh),$cRovh) . $sp->($cp->($Rbr),$tagc) . $sp->($cp->($Rctx),$cRctx) . ($rtrunc ? $M : "") . "-5'";

	my $hdr  = "oligo fusion: " . &name_chip($label1, \%nc) . " &gt; " . &name_chip($tagname, \%nc)
	         . " &gt; " . &name_chip($label4, \%nc) . "<br>";
	my $coretxt = ($mid eq "") ? "the core is double-stranded (shown in full)"
	                           : "the core is double-stranded (partly shown, rest elided as [...])";
	my $note = "<i>the $tagname tag is supplied as a pre-annealed oligo (two annealed strands); the two USER "
	         . "overhangs are single-stranded, $coretxt. "
	         . "full FW/RV strands are in the primer section.</i><br><br>";
	return $hdr . $note . $l1 . "<br>" . $l2 . "<br>" . $l3 . "<br>" . $l4 . "<br><br>";
}

# Can this tag be made into a viable oligo? (A near 5', T near 3', and a >=6 bp
# double-stranded core between the two 8-nt overhangs.) Mirrors the GUI gate so a
# too-short tag flagged -tag ...|1|... from the CLI falls back to a flap.
sub oligo_viable {
	my $s = uc $_[0]; my $ovh = defined $_[1] ? $_[1] : 8; my $minc = 6;
	my $a = index($s, 'A'); my $t = rindex($s, 'T');
	return 0 if $a < 0 or $t < 0;
	return (($t - $a - (2*$ovh - 1)) >= $minc) ? 1 : 0;
}

# True when a linker can genuinely be carried INSIDE the oligo, i.e. the oligo's
# gene-facing USER overhang is derived from the linker (not the tag). For a 5'
# tag the gene-facing end is the 3' end, so the terminal T (rindex) must fall in
# the linker; for a 3' tag it is the 5' end, so the leading A (index) must fall
# in the linker. Also requires the combined tag+linker to be oligo-viable. When
# false the linker cannot sit in the oligo and must ride on the gene primer.
sub linker_in_oligo {
	my ($tagseq, $linker, $end, $ovh) = @_;
	$ovh = 8 unless (defined $ovh && $ovh > 0);
	return 0 if (!defined $linker || $linker eq "");
	my $five = ($end eq '5' || $end eq 5);
	my $combined = $five ? (uc($tagseq) . uc($linker)) : (uc($linker) . uc($tagseq));
	return 0 unless &oligo_viable($combined, $ovh);
	if ($five) {
		my $t = rindex($combined, 'T');
		return ($t >= length($tagseq)) ? 1 : 0;        # terminal T must be in the linker
	}
	my $a = index($combined, 'A');
	return ($a >= 0 && $a < length($linker)) ? 1 : 0;  # leading A must be in the linker
}


#Locate the USER cassette in the destination plasmid and splice the
#assembled construct ($final) in its place.
#
#Recognition set = canonical $cassette[0] plus any additional
#variants listed in %cassettevariants for this cassette name. Each
#variant is searched on both strands (forward and reverse-complement)
#of the plasmid. Exactly one hit (across all variants, both strands)
#is required; zero or multiple hits return an empty sequence with an
#explanatory message so the caller can warn without aborting.
#
#Returns ($fullplasmid, $message).
sub fullplasmidreconstruct {
	my ($plasmid, $cassette_ref, $final) = @_;
	my @cas = @$cassette_ref;

	my $pl = uc $plasmid;
	$pl =~ s/[^A-Za-z]//g;
	if ($pl eq "") {
		return ("", "no destination plasmid sequence provided");
	}

	my $name = $cas[4];

	#Build the recognition set: canonical cassette first, then variants.
	my @variants;
	my %seen;
	my @candidates = (uc $cas[0]);
	if (exists $cassettevariants{$name}) {
		push @candidates, @{ $cassettevariants{$name} };
	}
	foreach my $v (@candidates) {
		$v = uc $v;
		$v =~ tr/U/T/;
		next if ($v eq "");
		next if ($seen{$v}++);
		push @variants, $v;
	}

	#Search every variant on both strands; collect all hits.
	my @hits;		#each hit: [pos, len, rev(0/1)]
	foreach my $v (@variants) {
		my $len = length $v;
		my $p = 0;
		while (($p = index($pl, $v, $p)) != -1) {
			push @hits, [$p, $len, 0];
			$p++;
		}
		my $rcv = &revcomp($v);
		if ($rcv ne $v) {
			$p = 0;
			while (($p = index($pl, $rcv, $p)) != -1) {
				push @hits, [$p, $len, 1];
				$p++;
			}
		}
	}

	if (@hits == 0) {
		return ("", "USER cassette ($name) not found in the destination plasmid");
	}
	if (@hits > 1) {
		return ("", "USER cassette ($name) found " . scalar(@hits) . " times in the destination plasmid (ambiguous)");
	}

	my ($pos, $len, $rev) = @{ $hits[0] };
	my $insert = $rev ? &revcomp($final) : $final;
	my $up   = substr($pl, 0, $pos);
	my $down = substr($pl, $pos + $len);
	my $full = $up . $insert . $down;
	my $orient = $rev ? "cassette found on reverse-complement strand" : "cassette found on forward strand";
	#Returns: full plasmid, status message, upstream backbone, oriented
	#construct, downstream backbone.
	return ($full, $orient, $up, $insert, $down);
}


#Turn the construct block just written to one report (array ref) into a
#FULL-PLASMID block: the same numbered/coloured "details of final construct"
#block, but extended with the destination-vector backbone up- and downstream
#of the cassette flanks (named "upstream backbone"/"downstream backbone").
#Does nothing unless -dp was supplied AND a USER cassette is in use. On
#failure to locate the cassette it leaves the construct block intact and
#appends a WARNING (never aborts). Also records a row for the -dpcsv CSV.
sub fullplasmidreportblock {
	my ($outref, $plasmid, $cassette_ref, $final, $label) = @_;
	my @cas = @$cassette_ref;
	$label = "" unless defined $label;

	return if ($plasmid eq "");		#-dp not supplied
	return if ($cas[5] eq "");		#no USER cassette (linear/circular assembly) -> ignore

	my ($full, $msg, $up, $construct, $down) = &fullplasmidreconstruct($plasmid, $cassette_ref, $final);

	if ($full eq "") {
		#Cassette not found / ambiguous: keep the construct block as written,
		#record a CSV row, and append a warning. The empty-string element
		#terminates the html final-construct sequence walker before the text.
		push(@fullplasmidrows, [$label, 0, $msg, ""]);
		push(@$outref, "");		#terminate the html final-construct sequence walker
		#Surface the problem prominently: a big red banner at the very top of the
		#report (rendered in htmlprint), instead of a small note at the bottom.
		$dp_cassette_warning = $msg;
		return;
	}

	#Record the clean (uppercase) full plasmid for the CSV.
	push(@fullplasmidrows, [$label, length($full), $msg, uc $full]);

	#Replace the construct block just written with a full-plasmid block.
	#Pop back through that block (numbered sequence, name row, dashes,
	#header), capturing the construct's fragment-name row on the way.
	my $namerow = "";
	while (@$outref) {
		my $line = pop @$outref;
		$namerow = $line if ($namerow eq "" && index($line, ",,") >= 0);
		last if ($line =~ /details of final construct after cloning/);
	}
	$namerow =~ s/\s+$//;

	#Build the backbone-cased sequence: each backbone is forced to a single
	#case opposite its neighbouring construct segment, so the html colouring
	#(which segments by case runs) treats each backbone as one segment, and
	#the construct keeps its existing per-fragment colours.
	my $seq = "";
	my @names;
	if ($up ne "") {
		my $c0 = substr($construct, 0, 1);
		$seq .= ($c0 =~ /[A-Z]/) ? lc($up) : uc($up);
		push @names, "upstream backbone";
	}
	$seq .= $construct;
	push @names, split(/,,/, $namerow);
	if ($down ne "") {
		my $cN = substr($construct, (length($construct) - 1), 1);
		$seq .= ($cN =~ /[A-Z]/) ? lc($down) : uc($down);
		push @names, "downstream backbone";
	}
	my $fullnamerow = join(",,", @names);
	my $len = length $seq;

	push(@$outref, "details of final construct after cloning (full plasmid, length: $len bases):\n");
	push(@$outref, "-------------------------------------------------------------------------------\n\n");
	push(@$outref, "$fullnamerow\n\n\n");
	for (my $p = 0; $p < length($seq); $p += 60) {
		my $num = $p + 1;
		my $line = $num;
		for (my $k = length($num); $k < 10; $k++) { $line .= " "; }
		$line .= substr($seq, $p, 60);
		push(@$outref, "$line\n");
	}
	push(@$outref, "");		#terminate the html sequence walker
	return;
}


#Minimal CSV field quoting (RFC-4180): quote if the field contains a
#comma, double-quote or newline; double any embedded quotes.
sub csvfield {
	my ($f) = @_;
	$f = "" unless defined $f;
	if ($f =~ /[",\r\n]/) {
		$f =~ s/"/""/g;
		$f = "\"$f\"";
	}
	return $f;
}

#Write the accumulated full plasmid(s) to the -dpcsv file.
#One row per construct: construct,length_bp,status,full_plasmid_sequence.
$options{dpcsv} = "plasmid" unless defined $options{dpcsv} and $options{dpcsv} ne "";
if (@fullplasmidrows) {
	if (open(my $CSV, '>', $options{dpcsv})) {
		print $CSV "construct,length_bp,status,full_plasmid_sequence\n";
		foreach my $row (@fullplasmidrows) {
			print $CSV join(",", map { &csvfield($_) } @$row) . "\n";
		}
		close $CSV;
	}
	else {
		print STDERR "WARNING: could not open -dpcsv file \"$options{dpcsv}\" for writing: $!\n";
	}
}


close OUTFRED;

###################################################################
# Silently delete internal output files.
# -h HTML is the sole user-facing output; -r/-o/-dpcsv are written
# internally for processing and cleaned up here unconditionally.

{
	my @_cleanup;
	if ($batchinput eq "FALSE") {
		push @_cleanup, $outreport if defined $outreport and $outreport ne "";
		push @_cleanup, $fastaout  if defined $fastaout  and $fastaout  ne "";
	}
	else {
		for my $_x (0 .. $#inputs) {
			if ($options{r} ne "") {
				push @_cleanup, $options{workdir} ne ""
					? "$outputpath/batch$_x.$outreportshort"
					: "batch$_x.$outreportshort";
			}
			if ($options{o} ne "") {
				push @_cleanup, $options{workdir} ne ""
					? "$outputpath/batch$_x.$fastaoutshort"
					: "batch$_x.$fastaout";
			}
		}
	}
	push @_cleanup, $options{dpcsv}
		if defined $options{dpcsv} and $options{dpcsv} ne "" and -e $options{dpcsv};
	unlink $_ for @_cleanup;
}
###################################################################
#####     GENBANK ANNOTATION TRANSFER  (added 2026-10-29)     #####
###################################################################
# Carry annotations through to the assembled insert (and, with -dp,
# the full plasmid) and write them as GenBank. Linear assembly only.
# Coordinates are transferred by offset arithmetic from the segment
# map built while mirroring the $final assembly; a substring self-
# check guards every anchor. Restriction/nick recognition sites are
# NOT annotated (per project decision 2026-10-29).

# --- GenBank flat-file parser --------------------------------------
sub gb_parse_file {
	my ($path) = @_;
	open(my $fh, "<", $path) or die "ERROR: cannot read GenBank file '$path': $!\n";
	local $/;
	my $text = <$fh>;
	close $fh;
	return gb_parse_text($text);
}

sub gb_parse_text {
	my ($text) = @_;
	$text =~ s/\r\n?/\n/g;
	my @records;
	my @blocks = split(/^\/\/[ \t]*$/m, $text);
	foreach my $blk (@blocks) {
		next unless $blk =~ /\S/;
		next unless $blk =~ /^LOCUS\s/m;
		my %rec = (name => "", topology => "linear", features => [], seq => "");
		if ($blk =~ /^LOCUS\s+(\S+)(.*)$/m) {
			$rec{name} = $1;
			$rec{topology} = "circular" if ($2 =~ /\bcircular\b/i);
		}
		if ($rec{name} eq "" and $blk =~ /^ACCESSION\s+(\S+)/m) { $rec{name} = $1; }
		$rec{name} = "sequence" if $rec{name} eq "";
		if ($blk =~ /^ORIGIN\b(.*)$/ms) {
			my $o = $1;
			$o =~ s/[^A-Za-z]//g;
			$rec{seq} = uc $o;
		}
		if ($blk =~ /^FEATURES[^\n]*\n(.*?)^ORIGIN\b/ms) {
			my $ftext = $1;
			my @feats;
			my $cur;
			foreach my $ln (split /\n/, $ftext) {
				if ($ln =~ /^\s{1,10}(\S+)\s+((?:complement\(|join\(|order\(|<|>|\d).*)$/) {
					push @feats, $cur if $cur;
					$cur = { type => $1, loc => $2, quals => {} };
				}
				elsif ($ln =~ /^\s{10,}\/(\w+)(?:=(.*))?$/ and $cur) {
					my ($k, $v) = ($1, defined $2 ? $2 : "");
					$v =~ s/^"//; $v =~ s/"[ \t]*$//;
					$cur->{quals}{$k} = $v unless exists $cur->{quals}{$k};
				}
				elsif ($ln =~ /^\s{10,}([^\/].*)$/ and $cur and !%{$cur->{quals}}) {
					# continuation of a multi-line location
					$cur->{loc} .= $1;
				}
			}
			push @feats, $cur if $cur;
			foreach my $f (@feats) {
				my ($s, $e, $strand) = gb_parse_location($f->{loc});
				next unless defined $s;
				my $label = $f->{quals}{label};
				$label = $f->{quals}{gene}    unless defined $label and $label ne "";
				$label = $f->{quals}{product} unless defined $label and $label ne "";
				$label = $f->{quals}{note}    unless defined $label and $label ne "";
				$label = $f->{type}           unless defined $label and $label ne "";
				push @{$rec{features}}, {
					type => $f->{type}, start => $s, end => $e,
					strand => $strand, label => $label,
				};
			}
		}
		push @records, \%rec;
	}
	return @records;
}

sub gb_parse_location {
	my ($loc) = @_;
	return () unless defined $loc;
	$loc =~ s/\s//g;
	my $strand = 1;
	$strand = -1 if ($loc =~ /complement\(/i);
	my @nums = ($loc =~ /(\d+)/g);
	return () unless @nums;
	my ($s, $e) = ($nums[0], $nums[-1]);
	($s, $e) = ($e, $s) if $s > $e;
	return ($s, $e, $strand);
}

# --- GenBank flat-file writer --------------------------------------
sub gb_write {
	my ($name, $seq, $circular, $feats) = @_;
	my $len = length $seq;
	$name =~ s/\s+/_/g;
	$name =~ s/-/_/g;   # underscores only in the LOCUS name (hyphens are non-standard)
	my $namew = length($name) > 16 ? length($name) : 16;   # pad short names to the
	                                                       # canonical 16-col field;
	                                                       # longer names (the full
	                                                       # output-file stem) just
	                                                       # extend it, fields stay
	                                                       # single-space separated.
	my $topo = $circular ? "circular" : "linear";
	my @mon = qw(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC);
	my @t = localtime(time);
	my $date = sprintf("%02d-%s-%04d", $t[3], $mon[$t[4]], $t[5] + 1900);
	my $out = "";
	# Canonical GenBank LOCUS columns (name 13, length ->40, bp 42, molecule
	# 48, topology 56, division 65, date 69) so strict parsers read topology.
	$out .= sprintf("LOCUS       %-${namew}s %11d bp    %-7s %-8s %-3s %s\n",
	                $name, $len, "DNA", $topo, "SYN", $date);
	$out .= "DEFINITION  AMUSER assembled construct with transferred annotations.\n";
	$out .= "ACCESSION   .\n";
	$out .= "VERSION     .\n";
	$out .= "KEYWORDS    .\n";
	# CLC Main Workbench requires a populated SOURCE/ORGANISM block before
	# FEATURES, otherwise it refuses to open the file. Emit a real organism
	# (Arabidopsis, the group's default) rather than ".".
	$out .= "SOURCE      Arabidopsis\n";
	$out .= "  ORGANISM  Arabidopsis\n";
	$out .= "FEATURES             Location/Qualifiers\n";
	# source feature spanning the whole molecule
	$out .= sprintf("     %-16s%s\n", "source", "1..$len");
	foreach my $f (sort { $a->{start} <=> $b->{start} or $a->{end} <=> $b->{end} } @$feats) {
		my $s = $f->{start}; my $e = $f->{end};
		next if ($s < 1 or $e > $len or $s > $e);
		my $locstr = ($f->{strand} && $f->{strand} == -1) ? "complement($s..$e)" : "$s..$e";
		my $type = $f->{type} || "misc_feature";
		$out .= sprintf("     %-16s%s\n", $type, $locstr);
		my $label = defined $f->{label} ? $f->{label} : "";
		$label =~ s/"/'/g;
		$out .= "                     /label=\"$label\"\n" if $label ne "";
		if (defined $f->{note} and $f->{note} ne "") {
			my $nt = $f->{note}; $nt =~ s/"/'/g;
			$out .= "                     /note=\"$nt\"\n";
		}
	}
	$out .= "ORIGIN\n";
	my $lc = lc $seq;
	for (my $i = 0; $i < $len; $i += 60) {
		$out .= sprintf("%9d", $i + 1);
		for (my $j = $i; $j < $i + 60 and $j < $len; $j += 10) {
			$out .= " " . substr($lc, $j, 10);
		}
		$out .= "\n";
	}
	$out .= "//\n";
	return $out;
}

# --- segment offset map (mirrors the $final assembly) --------------
sub gb_build_segmap {
	my ($x, $inputs, $pe, $cas, $final, $tag5, $tag3) = @_;
	$tag5 = {} unless ref $tag5 eq 'HASH';
	$tag3 = {} unless ref $tag3 eq 'HASH';
	my $pre  = (defined $cas->[5]) ? $cas->[5] : "";
	my $post = (defined $cas->[6]) ? $cas->[6] : "";
	my $ftpre  = "";
	my $ftpost = "";
	if (defined $cas->[2] and $cas->[2] ne "") {
		$ftpre = uc $cas->[2];
		$ftpre =~ tr/U/T/;
	}
	if (defined $cas->[3] and $cas->[3] ne "") {
		$ftpost = uc $cas->[3];
		$ftpost =~ tr/U/T/;
		$ftpost = reverse $ftpost;
		$ftpost =~ tr/ATGCRYKMBVDH/TACGYRMKVBHD/;
	}
	if ($pre eq "") { $ftpre = ""; $ftpost = ""; $post = ""; }   # no-cassette assembly has no flanks
	my @seg;
	my $p = 0;
	$p += length($pre);
	my $ftpre_start = $p;
	$p += length($ftpre);
	my $n = scalar @{$inputs->[$x]};
	for (my $i = 0; $i < $n; $i++) {
		if (defined $pe->[$x][$i][0] and $pe->[$x][$i][0] eq "TRUE") {
			push @seg, { kind => 'insert', name => undef, start => $p,
			             len => length($pe->[$x][$i][1]), seq => uc $pe->[$x][$i][1] };
			$p += length($pe->[$x][$i][1]);
		}
		if (defined $tag5->{$i}) {   # 5' tag (sense) sits right before the fragment
			my $ts = uc $tag5->{$i}{seq};
			push @seg, { kind => 'tag', name => $tag5->{$i}{name}, start => $p, len => length($ts), seq => $ts };
			$p += length($ts);
		}
		if (defined $tag5->{$i} and defined $tag5->{$i}{linker} and $tag5->{$i}{linker} ne "") {
			my $lk = uc $tag5->{$i}{linker};   # 5' linker: between tag and fragment
			push @seg, { kind => 'linker', name => 'linker', start => $p, len => length($lk), seq => $lk };
			$p += length($lk);
		}
		push @seg, { kind => 'fragment', name => $inputs->[$x][$i][1], start => $p,
		             len => length($inputs->[$x][$i][2]), seq => uc $inputs->[$x][$i][2] };
		$p += length($inputs->[$x][$i][2]);
		if (defined $tag3->{$i} and defined $tag3->{$i}{linker} and $tag3->{$i}{linker} ne "") {
			my $lk = uc $tag3->{$i}{linker};   # 3' linker: between fragment and tag
			push @seg, { kind => 'linker', name => 'linker', start => $p, len => length($lk), seq => $lk };
			$p += length($lk);
		}
		if (defined $tag3->{$i}) {   # 3' tag (sense) sits right after the fragment
			my $ts = uc $tag3->{$i}{seq};
			push @seg, { kind => 'tag', name => $tag3->{$i}{name}, start => $p, len => length($ts), seq => $ts };
			$p += length($ts);
		}
	}
	if (defined $pe->[$x][$n][0] and $pe->[$x][$n][0] eq "TRUE") {
		push @seg, { kind => 'insert', name => undef, start => $p,
		             len => length($pe->[$x][$n][1]), seq => uc $pe->[$x][$n][1] };
		$p += length($pe->[$x][$n][1]);
	}
	my $ftpost_start = $p;
	$p += length($ftpost);
	$p += length($post);
	return {
		seg => \@seg, total => $p,
		ftpre_start => $ftpre_start, ftpre_len => length($ftpre),
		ftpost_start => $ftpost_start, ftpost_len => length($ftpost),
	};
}

# --- primer annealing features (parsed from the report) ------------
sub gb_primer_anneal_fwd {
	my ($p, $frags) = @_;
	$p =~ tr/U/T/;
	my ($bestk, $bestseg);
	foreach my $f (@$frags) {
		my $fs = $f->{seq};
		my $max = (length($p) < length($fs)) ? length($p) : length($fs);
		for (my $k = $max; $k >= 12; $k--) {
			if (substr($p, length($p) - $k) eq substr($fs, 0, $k)) {
				if (!defined $bestk or $k > $bestk) { $bestk = $k; $bestseg = $f; }
				last;
			}
		}
	}
	return unless defined $bestk;
	return ($bestseg->{start} + 1, $bestseg->{start} + $bestk);
}

sub gb_primer_anneal_rev {
	my ($p, $frags) = @_;
	$p =~ tr/U/T/;
	my ($bestk, $bestseg);
	foreach my $f (@$frags) {
		my $fs = $f->{seq};
		my $max = (length($p) < length($fs)) ? length($p) : length($fs);
		for (my $k = $max; $k >= 12; $k--) {
			if (&revcomp(substr($p, length($p) - $k)) eq substr($fs, length($fs) - $k)) {
				if (!defined $bestk or $k > $bestk) { $bestk = $k; $bestseg = $f; }
				last;
			}
		}
	}
	return unless defined $bestk;
	return ($bestseg->{start} + $bestseg->{len} - $bestk + 1, $bestseg->{start} + $bestseg->{len});
}

sub gb_primer_features {
	my ($map, $report) = @_;
	my @feats;
	return @feats unless $report;
	my @frags = grep { $_->{kind} eq 'fragment' } @{$map->{seg}};
	return @feats unless @frags;
	my $curname = "";
	foreach my $ln (@$report) {
		if ($ln =~ /primer details for (.+?)\s*$/)       { $curname = $1; next; }
		if ($ln =~ /^forward primer:\s*([A-Za-z]+)/) {
			my ($s, $e) = gb_primer_anneal_fwd(uc $1, \@frags);
			push @feats, { type => 'primer_bind', start => $s, end => $e, strand => 1,
			               label => ($curname ne "" ? "$curname FW" : "forward primer") } if defined $s;
		}
		if ($ln =~ /^reverse primer:\s*([A-Za-z]+)/) {
			my ($s, $e) = gb_primer_anneal_rev(uc $1, \@frags);
			push @feats, { type => 'primer_bind', start => $s, end => $e, strand => -1,
			               label => ($curname ne "" ? "$curname RV" : "reverse primer") } if defined $s;
		}
	}
	return @feats;
}

# --- insert feature list (de novo + carry-through) -----------------
sub gb_insert_features {
	my ($map, $final, $report) = @_;
	my @feats;
	foreach my $s (@{$map->{seg}}) {
		my $st = $s->{start} + 1;
		my $en = $s->{start} + $s->{len};
		if ($s->{kind} eq 'fragment') {
			push @feats, { type => 'misc_feature', start => $st, end => $en, strand => 1,
			               label => $s->{name}, note => 'fragment' };
		}
		elsif ($s->{kind} eq 'tag') {
			push @feats, { type => 'misc_feature', start => $st, end => $en, strand => 1,
			               label => $s->{name}, note => 'tag' };
		}
		elsif ($s->{kind} eq 'linker') {
			push @feats, { type => 'misc_feature', start => $st, end => $en, strand => 1,
			               label => ($s->{name} // 'linker'), note => 'linker' };
		}
		else {
			push @feats, { type => 'misc_feature', start => $st, end => $en, strand => 1,
			               label => 'insert', note => 'USER insert/linker' };
		}
		# carry-through: match the segment sequence to a parsed GenBank record
		my $recs = $GBfrag{ $s->{seq} };
		if ($recs and @$recs) {
			my $rec = shift @$recs;
			foreach my $f (@{$rec->{features}}) {
				next if lc($f->{type}) eq 'source';   # whole-molecule tag; not meaningful in the construct
				my $ns = $s->{start} + $f->{start};   # 0-based start + 1-based a = 1-based mapped
				my $ne = $s->{start} + $f->{end};
				$ns = $st if $ns < $st;
				$ne = $en if $ne > $en;
				next if $ns > $ne;
				push @feats, { type => $f->{type}, start => $ns, end => $ne,
				               strand => $f->{strand}, label => $f->{label} };
			}
		}
	}
	if ($map->{ftpre_len} > 0) {
		push @feats, { type => 'misc_feature', start => $map->{ftpre_start} + 1,
		               end => $map->{ftpre_start} + $map->{ftpre_len}, strand => 1,
		               label => "USER overhang 5'", note => 'USER fusion overhang' };
	}
	if ($map->{ftpost_len} > 0) {
		push @feats, { type => 'misc_feature', start => $map->{ftpost_start} + 1,
		               end => $map->{ftpost_start} + $map->{ftpost_len}, strand => 1,
		               label => "USER overhang 3'", note => 'USER fusion overhang' };
	}
	push @feats, gb_primer_features($map, $report);
	return @feats;
}

# --- map backbone + insert features onto the full plasmid ----------
sub gb_plasmid_features {
	my ($insert_feats, $final, $destplasmid, $cas) = @_;
	my ($full, $msg, $up, $insert, $down) = &fullplasmidreconstruct($destplasmid, $cas, $final);
	return () if (!defined $full or $full eq "");
	my $rev  = ($insert ne $final) ? 1 : 0;
	my $L    = length($final);
	my $ulen = length($up);
	my $clen = length($destplasmid) - length($up) - length($down);   # excised cassette length
	my $shift = length($insert) - $clen;
	my @out;
	# insert-derived features -> plasmid coordinates
	foreach my $f (@$insert_feats) {
		my ($a, $b, $sdir) = ($f->{start}, $f->{end}, ($f->{strand} || 1));
		my ($na, $nb, $ns);
		if (!$rev) { $na = $ulen + $a;            $nb = $ulen + $b;            $ns = $sdir; }
		else       { $na = $ulen + ($L - $b + 1); $nb = $ulen + ($L - $a + 1); $ns = -$sdir; }
		push @out, { type => $f->{type}, start => $na, end => $nb, strand => $ns,
		             label => $f->{label}, note => $f->{note} };
	}
	# backbone (destination-vector) features, clipped around the cassette window
	if ($GBplasmid{features}) {
		my $pos = $ulen;              # cassette occupies 1-based [pos+1 .. pos+clen]
		foreach my $f (@{$GBplasmid{features}}) {
			next if lc($f->{type}) eq 'source';
			my ($a, $b) = ($f->{start}, $f->{end});
			# up-side surviving portion
			if ($a <= $pos) {
				my $ue = ($b < $pos) ? $b : $pos;
				push @out, { type => $f->{type}, start => $a, end => $ue, strand => $f->{strand},
				             label => $f->{label}, note => 'backbone' } if $a <= $ue;
			}
			# down-side surviving portion
			if ($b > $pos + $clen) {
				my $ds = ($a > $pos + $clen) ? $a : ($pos + $clen + 1);
				push @out, { type => $f->{type}, start => $ds + $shift, end => $b + $shift,
				             strand => $f->{strand}, label => $f->{label}, note => 'backbone' } if $ds <= $b;
			}
		}
	}
	return ($full, @out);
}

# --- top-level: build and write the GenBank output file(s) ---------
sub genbank_outputs {
	my ($x, $inputs, $pe, $cas, $final, $destplasmid, $report, $htmldir, $base, $tag5, $tag3) = @_;
	return unless defined $final and $final ne "";
	my $map = gb_build_segmap($x, $inputs, $pe, $cas, $final, $tag5, $tag3);
	if ($map->{total} != length($final)) {
		warn "WARNING: GenBank segment map length (" . $map->{total} .
		     ") != construct length (" . length($final) . "); skipping GenBank output.\n";
		return;
	}
	foreach my $s (@{$map->{seg}}) {
		my $chk = uc substr($final, $s->{start}, $s->{len});
		if ($chk ne $s->{seq}) {
			warn "WARNING: GenBank anchor mismatch for segment '" .
			     (defined $s->{name} ? $s->{name} : "insert") . "'; skipping GenBank output.\n";
			return;
		}
	}
	my @ifeats = gb_insert_features($map, $final, $report);

	my $gbdir = $htmldir;   # same "Output files" folder
	mkdir $gbdir unless -d $gbdir;
	(my $stem = $base) =~ s/\.html?$//i;
	$stem =~ s/_report$//;   # share the base stem (the report suffix is HTML-only)

	# GenBank output is exactly one file per run:
	#   - with a destination plasmid (-dp/-dpgb): the fully assembled plasmid
	#     (insert cloned into the backbone, circular) -> <stem>_plasmid.gb
	#   - without one: the linear assembled insert -> <stem>_insert.gb
	if (defined $destplasmid and $destplasmid ne "") {
		my ($full, @pfeats) = gb_plasmid_features(\@ifeats, $final, $destplasmid, $cas);
		if (defined $full and $full ne "") {
			my $pname = "${stem}_plasmid";   # match the output filename stem
			open(my $pfh, ">", "$gbdir/${stem}_plasmid.gb") or warn "WARNING: cannot write plasmid GenBank: $!\n";
			if ($pfh) { print $pfh gb_write($pname, uc $full, 1, \@pfeats); close $pfh; }
		}
	}
	else {
		my $iname = "${stem}_insert";   # match the output filename stem
		open(my $gfh, ">", "$gbdir/${stem}_insert.gb") or warn "WARNING: cannot write insert GenBank: $!\n";
		if ($gfh) { print $gfh gb_write($iname, uc $final, 0, \@ifeats); close $gfh; }
	}
	return;
}

###################################################################
# Primer-list CSV output.
# Scrapes the already-assembled text-report primer overview
# ("fragment\tstrand\ttm\tprimer" + rows) so the CSV matches the
# report 1:1, and writes three columns: primer name, Tm, sequence.
# Primer name = <fragment>_FW (forward) / <fragment>_RV (reverse).
###################################################################
sub write_primer_csv {
	my ($rows, $csvpath) = @_;
	return 0 unless ref $rows eq 'ARRAY';
	my @out;
	my $in = 0;
	foreach my $line (@$rows) {
		if ($line =~ /^fragment\tstrand\ttm\tprimer/) { $in = 1; next; }
		next unless $in;
		if ($line =~ /\t(?:forward|reverse)\t/) {
			my @f = split /\t/, $line;
			next unless scalar(@f) >= 4;
			my ($frag, $strand, $tm, $seq) = @f[0 .. 3];
			$seq  =~ s/[\r\n]+//g;
			$frag =~ s/^\s+|\s+$//g;
			$tm   =~ s/\s*C\s*$//i;    # strip the trailing "C" from e.g. "59.8C"
			$tm   =~ s/^\s+|\s+$//g;
			my $suffix = ($strand eq 'forward') ? 'FW' : 'RV';
			push @out, join(",", map { &_csv_field($_) } ("${frag}_${suffix}", $tm, $seq));
		}
		else { last if @out; }   # overview block has ended
	}
	return 0 unless @out;
	open(my $fh, '>', $csvpath) or do { warn "WARNING: cannot write primer CSV: $!\n"; return 0; };
	print $fh "primer name,Tm,primer sequence\n";
	print $fh "$_\n" for @out;
	close $fh;
	return 1;
}
sub _csv_field {
	my $v = shift;
	$v = "" unless defined $v;
	if ($v =~ /[",\r\n]/) { $v =~ s/"/""/g; return "\"$v\""; }
	return $v;
}
