#!/usr/bin/env perl
#######################################################################
# amuser_gui.pl - local web GUI launcher for AMUSER_v39.pl
#
# Starts a tiny web server on 127.0.0.1 (loopback only), serves a form,
# and runs AMUSER_v39.pl when you submit. The HTML report is shown
# straight back in the browser. No files leave your machine, nothing is
# installed (core Perl only), and it works on Windows, macOS and Linux.
#
# Usage:  perl amuser_gui.pl        (or double-click a launcher)
# Then open the printed http://127.0.0.1:<port>/ address in a browser.
#######################################################################
use strict;
use warnings;
use IO::Socket::INET;
use IPC::Open3;
use Symbol 'gensym';
use File::Temp qw(tempfile);
use File::Basename qw(dirname basename);
use Cwd qw(abs_path);
use JSON::PP;

$| = 1;   # autoflush so the launcher window shows the URL immediately
$SIG{PIPE} = 'IGNORE';   # a client that disconnects mid-response must not kill us

my $SCRIPT_DIR = dirname(abs_path(__FILE__));
my $AMUSER     = "$SCRIPT_DIR/AMUSER_v39.pl";
my @CASSETTES  = (
    "PacI/Nt.BbvCI", "XbaI2/Nt.BbvCI", "PmeI/Nb.BbvCI",
    "AsiSI/Nb.BsmI", "AsiSI/Nb.BtsI",
);

die "ERROR: cannot find AMUSER_v39.pl next to this launcher (looked in $SCRIPT_DIR)\n"
    unless -f $AMUSER;

#--- test hook: load subs without starting the server ----------------
return 1 if $ENV{AMUSER_GUI_TEST};

#--- bind a loopback socket, trying a small range of ports ------------
my $server;
my $PORT;
for my $p (8347 .. 8377) {
    $server = IO::Socket::INET->new(
        LocalAddr => '127.0.0.1', LocalPort => $p,
        Proto => 'tcp', Listen => 5, ReuseAddr => 1,
    );
    if ($server) { $PORT = $p; last; }
}
die "ERROR: could not bind any port in 8347-8377. Is a previous GUI still running?\n"
    unless $server;

my $URL = "http://127.0.0.1:$PORT/";
print "\n  AMUSER GUI is running.\n";
print "  Open this in your browser:  $URL\n";
print "  (Leave this window open while you work. Press Ctrl+C here to stop.)\n\n";

# open the GUI in the default browser (best-effort; harmless if it fails)
unless ($ENV{AMUSER_GUI_NOOPEN}) {
    my @open = $^O eq 'MSWin32' ? ('cmd', '/c', 'start', '', $URL)
             : $^O eq 'darwin'  ? ('open', $URL)
             :                     ('xdg-open', $URL);
    system(@open);
}

#--- accept loop ------------------------------------------------------
while (my $client = $server->accept) {
    $client->timeout(30);
    my $req = _read_request($client);
    if (!$req) { close $client; next; }

    my ($method, $path, $body) = @$req;

    if ($method eq 'GET' && ($path eq '/' || $path eq '/index.html')) {
        _send($client, 200, 'text/html; charset=utf-8', page_html());
    }
    elsif ($method eq 'GET' && $path eq '/favicon.ico') {
        _send($client, 204, 'text/plain', '');
    }
    elsif ($method eq 'POST' && $path eq '/run') {
        my $out = handle_run($body);
        _send($client, 200, 'application/json; charset=utf-8', encode_json($out));
    }
    elsif ($method eq 'GET' && $path =~ m{^/report/([\w.\-]+\.html)$}) {
        my $file = "$SCRIPT_DIR/Output files/$1";
        if (-f $file) {
            open my $rh, '<:raw', $file;
            local $/; my $html = <$rh>; close $rh;
            _send($client, 200, 'text/html; charset=utf-8', $html);
        }
        else {
            _send($client, 404, 'text/plain; charset=utf-8', "Report not found\n");
        }
    }
    else {
        _send($client, 404, 'text/plain; charset=utf-8', "Not found\n");
    }
    close $client;
}

#######################################################################
# Request handling
#######################################################################
sub _read_request {
    my $c = shift;
    local $/ = "\r\n";
    my $line = <$c>;
    return undef unless defined $line;
    $line =~ s/\r\n$//;
    my ($method, $path) = split /\s+/, $line;
    return undef unless defined $method && defined $path;

    my %h;
    while (my $hl = <$c>) {
        $hl =~ s/\r\n$//;
        last if $hl eq '';
        $h{lc $1} = $2 if $hl =~ /^([^:]+):\s*(.*)$/;
    }
    my $body = '';
    if (defined $h{'content-length'} && $h{'content-length'} > 0) {
        my $need = $h{'content-length'};
        while (length($body) < $need) {
            my $chunk;
            my $got = read($c, $chunk, $need - length($body));
            last if !defined $got || $got == 0;
            $body .= $chunk;
        }
    }
    return [$method, $path, $body];
}

sub _send {
    my ($c, $code, $ctype, $body) = @_;
    my %msg = (200 => 'OK', 204 => 'No Content', 404 => 'Not Found', 500 => 'Internal Server Error');
    my $status = $msg{$code} || 'OK';
    # Serialise to UTF-8 bytes so Content-Length matches the bytes actually
    # written (a character string with any non-ASCII char would otherwise make
    # Content-Length too small and the browser would truncate the response,
    # breaking the inline <script>). Byte strings (e.g. a report read raw) are
    # left untouched.
    utf8::encode($body) if utf8::is_utf8($body);
    my $len = length($body);
    print $c "HTTP/1.1 $code $status\r\n";
    print $c "Content-Type: $ctype\r\n";
    print $c "Content-Length: $len\r\n";
    print $c "Connection: close\r\n";
    print $c "Cache-Control: no-store\r\n";
    print $c "\r\n";
    print $c $body;
}

#######################################################################
# Run AMUSER with the submitted parameters
#######################################################################
sub handle_run {
    my $body = shift;
    my $in = eval { decode_json($body) };
    return { ok => 0, error => "Could not read the submitted form." } if !$in;

    # --- validate / normalise ---------------------------------------
    # GenBank fragment input (-igb) supersedes the FASTA text box when present.
    my $igb = defined $in->{igb} ? $in->{igb} : '';
    $igb =~ s/\r\n/\n/g; $igb =~ s/^\s+//; $igb =~ s/\s+$//;

    my $fasta = defined $in->{fasta} ? $in->{fasta} : '';
    $fasta =~ s/\r\n/\n/g;
    $fasta =~ s/^\s+//; $fasta =~ s/\s+$//;
    if ($igb eq '') {
        return { ok => 0, error => "Please enter at least one sequence." } if $fasta eq '';
        # single headerless sequence -> add a header
        $fasta = ">seq1\n$fasta\n" unless $fasta =~ /^>/m;
        $fasta .= "\n" unless $fasta =~ /\n$/;
    }
    else {
        return { ok => 0, error => "The uploaded GenBank file has no LOCUS record." }
            unless $igb =~ /^LOCUS\s/m;
    }

    my $cassette = $in->{cassette} // $CASSETTES[0];
    my @args;
    if ($cassette eq '__CUSTOM__') {
        my $full = uc(defined $in->{ccustfull} ? $in->{ccustfull} : ''); $full =~ s/\s//g;
        my $o5   = uc(defined $in->{ccust5}    ? $in->{ccust5}    : ''); $o5   =~ s/\s//g;
        my $o3   = uc(defined $in->{ccust3}    ? $in->{ccust3}    : ''); $o3   =~ s/\s//g;
        return { ok => 0, error => "Custom cassette: enter the full cassette and both overhangs." }
            if ($full eq '' or $o5 eq '' or $o3 eq '');
        return { ok => 0, error => "Custom cassette: the full cassette must be DNA (A/C/G/T only)." }
            if ($full =~ /[^ACGT]/);
        return { ok => 0, error => "Custom cassette: the overhangs must be DNA ending in U (or T)." }
            if ($o5 =~ /[^ACGTU]/ or $o3 =~ /[^ACGTU]/);
        push @args, '-ccustfull', $full, '-ccust5', $o5, '-ccust3', $o3;
        my $cn = defined $in->{cname} ? $in->{cname} : ''; $cn =~ s/^\s+//; $cn =~ s/\s+$//;
        push @args, '-cname', $cn if $cn ne '';
    }
    else {
        return { ok => 0, error => "Unknown cassette: $cassette" }
            unless grep { $_ eq $cassette } @CASSETTES;
        push @args, '-cstandard', $cassette;
    }

    # tm: blank = off; else number 40-72 or "eq"
    my $tm = defined $in->{tm} ? $in->{tm} : '';
    $tm =~ s/^\s+//; $tm =~ s/\s+$//;
    if ($tm ne '') {
        unless (lc($tm) eq 'eq' || ($tm =~ /^\d+(\.\d+)?$/ && $tm >= 40 && $tm <= 72)) {
            return { ok => 0, error => "Tm must be blank, 'eq', or a number 40-72 (got '$tm')." };
        }
        push @args, '-tm', (lc($tm) eq 'eq' ? 'eq' : $tm);
    }

    # poly: which polymerase the Tm model uses. No buffer sub-choice and no salt
    # field -- the salt belongs to the polymerase, exactly as at NEB, and a GUI
    # salt box would let the user produce a Tm matching no NEB product.
    my %POLY_DEFAULT_PC = (phusion => '0.5', taq => '0.2');
    my $poly = defined $in->{poly} ? lc $in->{poly} : 'phusion';
    $poly =~ s/^\s+//; $poly =~ s/\s+$//;
    return { ok => 0, error => "Polymerase must be 'phusion' or 'taq' (got '$poly')." }
        unless exists $POLY_DEFAULT_PC{$poly};
    push @args, '-poly', $poly;

    # pc: primer concentration IN THE REACTION TUBE, in uM, range 0.05-5.
    # Default follows the polymerase: 0.5 uM Phusion, 0.2 uM Taq.
    my $pc = $in->{pc};
    $pc = $POLY_DEFAULT_PC{$poly} unless defined $pc && $pc ne '';
    $pc =~ s/^\s+//; $pc =~ s/\s+$//;
    return { ok => 0, error => "-pc must be a number (got '$pc')." }
        unless $pc =~ /^\d+(\.\d+)?$/;
    return { ok => 0, error => "-pc must be between 0.05 and 5 uM (got '$pc')." }
        if $pc < 0.05 || $pc > 5;
    push @args, '-pc', $pc;

    # fo: array of {junction=>N(1-based), side=>'A'|'B'}
    my @fo;
    if (ref $in->{fo} eq 'ARRAY') {
        for my $j (@{$in->{fo}}) {
            next unless defined $j->{side} && ($j->{side} eq 'A' || $j->{side} eq 'B');
            my $num = $j->{junction};
            return { ok => 0, error => "Bad forced-overhang junction number." }
                unless defined $num && $num =~ /^\d+$/ && $num >= 1;
            push @fo, "$num:$j->{side}";
        }
    }
    push @args, '-fo', join(',', @fo) if @fo;

    # tags: array of {frag=>N(1-based), end=>5|3, name=>str, oligo=>0|1, seq=>sense DNA}.
    # Passed to AMUSER as repeatable -tag "frag|end|name|oligo|seq" (seq last so it
    # may contain anything; name has no pipe; oligo is 1 to order it as a standalone
    # annealed oligo instead of a primer flap).
    if (ref $in->{tags} eq 'ARRAY') {
        for my $t (@{$in->{tags}}) {
            my $frag = $t->{frag};
            my $end  = $t->{end};
            my $seq  = defined $t->{seq}  ? $t->{seq}  : '';
            my $name = defined $t->{name} ? $t->{name} : '';
            my $oligo = ($t->{oligo} ? 1 : 0);
            my $linkoligo = ($t->{linkoligo} ? 1 : 0);
            my $linker = defined $t->{linker} ? uc $t->{linker} : '';
            $linker =~ s/[^ACGT]//g;   # linker is DNA only; no pipe
            $seq =~ s/\s//g;
            next if $seq eq '';
            return { ok => 0, error => "Bad tag fragment number." }
                unless defined $frag && $frag =~ /^\d+$/ && $frag >= 1;
            return { ok => 0, error => "Bad tag end (must be 5 or 3)." }
                unless defined $end && ($end eq '5' || $end eq '3' || $end == 5 || $end == 3);
            $name =~ s/\|//g;   # keep the delimiter unambiguous
            push @args, '-tag', "$frag|$end|$name|$oligo|$linkoligo|$linker|$seq";
        }
    }

    # autofrag: array of 1-based fragment indices to auto-split (>2 kb genes
    # the user marked). Passed to AMUSER as -autofrag "i,j"; the backend cuts
    # each into ~1 kb sub-fragments at native USER overhangs and fuses them.
    if (ref $in->{autofrag} eq 'ARRAY' && @{$in->{autofrag}}) {
        my @af;
        for my $f (@{$in->{autofrag}}) {
            return { ok => 0, error => "Bad auto-fragment number." }
                unless defined $f && $f =~ /^\d+$/ && $f >= 1;
            push @af, $f;
        }
        push @args, '-autofrag', join(',', @af) if @af;
    }

    # Destination plasmid: a GenBank upload (-dpgb) supersedes the plain
    # text box (-dp). Written to a temp file when GenBank is supplied.
    my $dpgb = defined $in->{dpgb} ? $in->{dpgb} : '';
    $dpgb =~ s/\r\n/\n/g; $dpgb =~ s/^\s+//; $dpgb =~ s/\s+$//;
    my $dpgbfile;
    if ($dpgb ne '') {
        return { ok => 0, error => "The uploaded plasmid GenBank file has no LOCUS record." }
            unless $dpgb =~ /^LOCUS\s/m;
        my ($dfh, $dfn) = tempfile("amuser_dp_XXXXXX", SUFFIX => '.gb', TMPDIR => 1, UNLINK => 0);
        binmode $dfh; print $dfh $dpgb; close $dfh;
        $dpgbfile = $dfn;
        push @args, '-dpgb', $dfn;
    }
    else {
        # dp: plasmid sequence; accept an optional FASTA header, then strip it
        my $dp = defined $in->{dp} ? $in->{dp} : '';
        $dp =~ s/^>.*$//mg;   # drop any FASTA header line(s)
        $dp =~ s/\s+//g;      # then remove all whitespace/newlines
        if ($dp ne '') {
            return { ok => 0, error => "The full-plasmid sequence contains non-DNA characters." }
                unless $dp =~ /^[ACGTURYKMSWBDHVN]+$/i;
            push @args, '-dp', $dp;
        }
    }

    # oname: optional custom base name for the three output files. Passed
    # through as -oname; AMUSER sanitises it and builds
    # "<name>_report.html" / "_insert.gb"|"_plasmid.gb" / "_primers.csv"
    # with no date/time stamp. Blank => AMUSER's automatic timestamped name.
    my $oname = defined $in->{oname} ? $in->{oname} : '';
    $oname =~ s/^\s+//; $oname =~ s/\s+$//;
    if ($oname ne '') {
        (my $probe = $oname) =~ s/\.html?$//i;
        $probe =~ s/_report$//i;
        $probe =~ s/\s+/_/g;
        $probe =~ s/[^A-Za-z0-9._-]/_/g;
        $probe =~ s/^[._-]+//;
        return { ok => 0, error => "The output file name has no usable characters (use letters, digits, . _ or -)." }
            if $probe eq '';
        push @args, '-oname', $oname;
    }

    # --- write the fragment input to a temp file --------------------
    my $infile;
    if ($igb ne '') {
        my ($gfh, $gfn) = tempfile("amuser_igb_XXXXXX", SUFFIX => '.gb', TMPDIR => 1, UNLINK => 0);
        binmode $gfh; print $gfh $igb; close $gfh;
        $infile = $gfn;
        push @args, '-igb', $gfn;
    }
    else {
        my ($fh, $fname) = tempfile("amuser_in_XXXXXX", SUFFIX => '.fa', TMPDIR => 1, UNLINK => 0);
        binmode $fh; print $fh $fasta; close $fh;
        $infile = $fname;
        push @args, '-i', $fname;
    }

    # --- run it ------------------------------------------------------
    # No -h: AMUSER names the report after the construct
    # ("<fragment>_YYYY-MM-DD_HH-MM-SS" or "fusion_...") and reports the
    # chosen filename on stdout as "AMUSER report: <name>".
    my ($stdout, $stderr, $exit) = _run_amuser(@args);

    unlink $infile  if defined $infile;
    unlink $dpgbfile if defined $dpgbfile;

    my ($hname) = $stdout =~ /AMUSER report:\s*(\S+)/;
    $hname = '' unless defined $hname;
    if ($hname eq '') {
        my $err = $stderr || $stdout || "AMUSER did not report an output file (exit $exit).";
        $err =~ s/\s+$//;
        return { ok => 0, error => $err, command => _display_command(@args) };
    }

    my $report_path = "$SCRIPT_DIR/Output files/$hname";
    if (-f $report_path && -s $report_path) {
        open my $rh, '<:raw', $report_path or
            return { ok => 0, error => "Report was generated but could not be read." };
        local $/; my $html = <$rh>; close $rh;
        (my $stem = $hname) =~ s/\.html?$//i;
        $stem =~ s/_report$//;   # .gb/.csv share the base stem; report suffix is HTML-only
        my $gbdir = "$SCRIPT_DIR/Output files";
        my @gb;   # one of these per run: _plasmid.gb with -dp, else _insert.gb
        push @gb, "Output files/${stem}_insert.gb"  if -f "$gbdir/${stem}_insert.gb";
        push @gb, "Output files/${stem}_plasmid.gb" if -f "$gbdir/${stem}_plasmid.gb";
        my $csvrel = "Output files/${stem}_primers.csv";
        my $saved = "Output files/$hname";
        $saved .= "\nPrimer CSV: $csvrel" if -f "$SCRIPT_DIR/$csvrel";
        $saved .= "\nGenBank: " . join(", ", @gb) if @gb;
        return {
            ok      => 1,
            name    => $hname,
            command => _display_command(@args),
            warn    => $stderr,
            saved   => $saved,
        };
    }
    else {
        my $err = $stderr || $stdout || "AMUSER did not produce a report (exit $exit).";
        $err =~ s/\s+$//;
        return { ok => 0, error => $err, command => _display_command(@args) };
    }
}

sub _run_amuser {
    my @args = @_;
    local $ENV{AMUSER_NO_OPEN} = 1;   # tell AMUSER not to spawn a browser
    my $err = gensym;
    my $pid = open3(my $in, my $out, $err, $^X, $AMUSER, @args);
    close $in;
    local $/;
    my $stdout = <$out>; $stdout = '' unless defined $stdout;
    my $stderr = <$err>; $stderr = '' unless defined $stderr;
    waitpid($pid, 0);
    my $exit = $? >> 8;
    return ($stdout, $stderr, $exit);
}

sub _display_command {
    my @a = @_;
    my @show = ('perl', 'AMUSER_v39.pl');
    for my $x (@a) {
        if ($x =~ /^-/) { push @show, $x; }
        elsif ($x =~ /amuser_igb_\w+\.gb$/) { push @show, '<your GenBank>'; }
        elsif ($x =~ /amuser_dp_\w+\.gb$/)  { push @show, '<your plasmid GenBank>'; }
        elsif ($x =~ /\.fa$/) { push @show, '<your FASTA>'; }
        elsif (length($x) > 40) { push @show, '"' . substr($x, 0, 20) . '...' . '"'; }
        elsif ($x =~ /[\s\/]/) { push @show, "\"$x\""; }
        else { push @show, $x; }
    }
    return join(' ', @show);
}

#######################################################################
# The GUI page
#######################################################################
sub page_html {
    my $opts = join('', map { "<option value=\"$_\">$_</option>" } @CASSETTES);
    $opts .= "<option value=\"__CUSTOM__\">Custom\x{2026}</option>";
    my $html = <<'HTML';
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AMUSER</title>
<style>
  :root{--bg:#f6f7f9;--card:#fff;--ink:#1c2128;--mut:#6b7681;--line:#e2e6ea;
        --accent:#2f6f4f;--accent2:#39825c;--warn:#8a5a00;--warnbg:#fff6e5;--err:#b42318;--errbg:#fdecec;}
  @media (prefers-color-scheme:dark){
    :root{--bg:#15181c;--card:#1d2127;--ink:#e6e9ec;--mut:#9aa4ad;--line:#2b313a;
          --accent:#5bbd8a;--accent2:#6fcf9c;--warn:#e0b25a;--warnbg:#2a2412;--err:#ff7b6b;--errbg:#2a1715;}
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--ink);
       font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif}
  .wrap{max-width:960px;margin:0 auto;padding:24px 16px 64px}
  h1{font-size:22px;margin:0 0 2px} .sub{color:var(--mut);margin:0 0 20px;font-size:13px}
  .card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px;margin-bottom:16px}
  label{display:block;font-weight:600;margin:0 0 6px}
  .hint{color:var(--mut);font-weight:400;font-size:12.5px}
  textarea,input,select{width:100%;font:inherit;color:inherit;background:var(--bg);
       border:1px solid var(--line);border-radius:8px;padding:9px 10px}
  textarea{resize:vertical;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px}
  .row{display:flex;gap:14px;flex-wrap:wrap}
  .row>div{flex:1;min-width:150px}
  .juncrow{display:flex;align-items:center;gap:12px;margin:8px 0}
  .juncrow .jl{min-width:190px;color:var(--mut);font-size:13px}
  .juncrow select{max-width:240px}
  button{appearance:none;border:0;border-radius:8px;padding:11px 20px;font:inherit;font-weight:600;
       background:var(--accent);color:#fff;cursor:pointer}
  button:hover{background:var(--accent2)} button:disabled{opacity:.6;cursor:progress}
  .foempty{color:var(--mut);font-size:13px}
  .msg{padding:11px 13px;border-radius:8px;margin:12px 0;font-size:13.5px;white-space:pre-wrap;
       font-family:ui-monospace,Menlo,Consolas,monospace}
  .err{background:var(--errbg);color:var(--err);border:1px solid var(--err)}
  .warn{background:var(--warnbg);color:var(--warn);border:1px solid var(--warn)}
  .ok{background:color-mix(in srgb,var(--accent) 12%,transparent);color:var(--ink);border:1px solid var(--accent)}
  .ok a{color:var(--accent);font-weight:600}
  .saved{color:var(--mut);font-size:12.5px;margin-top:8px}
  #results{display:none}
  iframe{width:100%;height:70vh;border:1px solid var(--line);border-radius:10px;background:#fff}
  .topbar{display:flex;justify-content:space-between;align-items:center;gap:10px;margin-bottom:10px}
  .linkbtn{background:none;color:var(--accent);font-weight:600;padding:6px 0;cursor:pointer;border:0}
  .lblrow{display:flex;justify-content:space-between;align-items:flex-end;gap:10px;flex-wrap:wrap}
  .uploadbtn{background:var(--bg);border:1px solid var(--line);border-radius:7px;color:var(--accent);
       font:inherit;font-weight:600;font-size:13px;padding:6px 12px;cursor:pointer;white-space:nowrap}
  .uploadbtn:hover{border-color:var(--accent)}
  .tmopt label.radio{display:block;font-weight:400;margin:7px 0;cursor:pointer}
  .tmopt input[type=radio]{width:auto;margin-right:8px;vertical-align:middle}
  .tmopt input[type=number]{width:120px;display:inline-block;margin:0 6px;vertical-align:middle}
  .tmopt .radio:has(input:disabled) input[type=number]{opacity:.5}
  .note{color:var(--mut);font-style:italic;font-size:12.5px;margin:10px 0 0}
  /* ordered fragment diagram */
  .orderbox{margin-top:12px;border:1px dashed var(--line);border-radius:10px;padding:12px 12px 14px;background:var(--bg)}
  .chainlbl{color:var(--mut);font-size:12.5px;margin-bottom:10px}
  .chain{display:flex;align-items:center;flex-wrap:wrap;gap:4px}
  .chip{display:inline-flex;align-items:center;gap:3px;background:var(--card);border:1px solid var(--accent);
        border-radius:8px;padding:4px 5px;cursor:grab}
  .chip:active{cursor:grabbing}
  .chip.dragging{opacity:.35}
  .chip.over{outline:2px solid var(--accent2);outline-offset:1px}
  .chip .nm{font-weight:600;font-size:13px;padding:0 3px;
        font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
        max-width:200px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .chip .mv,.chip .rm{appearance:none;border:0;background:none;color:var(--mut);cursor:pointer;
        padding:0 2px;font-size:12px;line-height:1;font-weight:700}
  .chip .mv:hover,.chip .rm:hover{color:var(--accent)}
  .chip .mv:disabled{opacity:.3;cursor:default}
  .chip .rm{color:var(--err);font-size:15px}
  .conn{color:var(--mut);font-weight:700;padding:0 1px}
  /* collapsible advanced section */
  details.advcard{padding:0}
  details.advcard>summary{list-style:none;cursor:pointer;font-weight:600;padding:16px 18px;
        display:flex;align-items:center;gap:8px;user-select:none}
  details.advcard>summary::-webkit-details-marker{display:none}
  details.advcard>summary::before{content:"\25B8";color:var(--mut);font-size:12px;
        transition:transform .15s ease}
  details.advcard[open]>summary::before{transform:rotate(90deg)}
  details.advcard>.row{padding:0 18px 16px}
  /* per-sequence options (collapsible cards) */
  .fragopts{margin-top:12px}
  .fragopts>.cap{color:var(--mut);font-size:12.5px;margin:0 0 8px}
  details.fragopt{background:var(--bg);border:1px solid var(--line);border-radius:9px;margin:6px 0}
  details.fragopt>summary{list-style:none;cursor:pointer;display:flex;align-items:center;gap:10px;
        padding:9px 11px;user-select:none}
  details.fragopt>summary::-webkit-details-marker{display:none}
  details.fragopt>summary::before{content:"\25B8";color:var(--mut);font-size:11px;transition:transform .15s ease}
  details.fragopt[open]>summary::before{transform:rotate(90deg)}
  details.fragopt .idx{display:inline-flex;align-items:center;justify-content:center;min-width:20px;height:20px;
        border-radius:6px;background:color-mix(in srgb,var(--accent) 14%,transparent);color:var(--accent);
        font-weight:700;font-size:11.5px}
  details.fragopt .nm{font-weight:600;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px;
        overflow:hidden;text-overflow:ellipsis;white-space:nowrap;max-width:200px}
  details.fragopt .ends{color:var(--mut);font-size:11px;font-family:ui-monospace,Menlo,Consolas,monospace;margin-left:auto}
  details.fragopt .badges{color:var(--accent);font-size:11.5px;font-weight:600;font-family:ui-monospace,Menlo,Consolas,monospace}
  details.fragopt .optbody{display:flex;flex-wrap:wrap;align-items:center;gap:9px 16px;padding:2px 12px 12px 30px}
  details.fragopt .chk{display:inline-flex;align-items:center;gap:7px;font-size:13px;cursor:pointer}
  details.fragopt .chk input{width:16px;height:16px;accent-color:var(--accent)}
  details.fragopt .chk.off{color:var(--mut)} details.fragopt .chk.off input{cursor:not-allowed}
  details.fragopt .codon{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:11.5px;color:var(--mut)}
  details.fragopt .tagbody{border-top:1px dashed var(--line);padding-top:10px;gap:9px 24px}
  details.fragopt .tagwrap{display:inline-flex;flex-wrap:wrap;align-items:center;gap:7px;font-size:13px}
  details.fragopt .tagwrap .tagcodon{flex-basis:100%;margin-left:49px;font-size:12.5px}
  details.fragopt .taglbl{color:var(--mut);font-weight:600;min-width:42px}
  details.fragopt .tagwrap select{font:inherit;font-size:12.5px;border:1px solid var(--line);border-radius:6px;
        padding:4px 7px;background:var(--card);color:inherit;max-width:170px}
  details.fragopt .tagcustom{font:inherit;font-size:12px;font-family:ui-monospace,Menlo,Consolas,monospace;
        border:1px solid var(--line);border-radius:6px;padding:4px 7px;background:var(--bg);color:inherit;width:230px}
  details.fragopt .tagcustomname{font:inherit;font-size:12px;border:1px solid var(--line);border-radius:6px;
        padding:4px 7px;background:var(--bg);color:inherit;width:130px}
  details.fragopt .soon{font-size:11px;color:var(--mut);font-style:italic}
</style>
</head>
<body>
<div class="wrap">
  <h1>AMUSER</h1>
  <p class="sub">USER-cloning primer design &middot; local interface for AMUSER_v39.pl</p>

  <div id="form">
    <div class="card">
      <div class="lblrow">
        <label for="fasta" style="margin:0">Sequences (<code>-i</code>)
          <span class="hint">&mdash; FASTA, or upload GenBank (<code>.gb</code>) to carry its annotations through to the output. Order = fragment order.</span>
        </label>
        <span style="display:flex;gap:8px">
          <button type="button" class="uploadbtn" onclick="document.getElementById('fastafile').click()">Upload FASTA</button>
          <button type="button" class="uploadbtn" onclick="document.getElementById('gbfile').click()">Upload GenBank</button>
        </span>
        <input type="file" id="fastafile" accept=".fa,.fasta,.fna,.ffn,.frn,.seq,.txt" multiple style="display:none" onchange="loadFastaFiles(this)">
        <input type="file" id="gbfile" accept=".gb,.gbk,.genbank,.ape,.txt" multiple style="display:none" onchange="loadFragGenbank(this)">
      </div>
      <textarea id="fasta" rows="9" placeholder="ATGAGT...   (or paste multiple &gt;headed records, or upload several files)" oninput="onFastaInput()" style="margin-top:8px"></textarea>
      <div id="gbfragstatus" class="saved"></div>
      <div id="orderdiagram" class="orderbox" style="display:none"></div>
      <div id="fragopts" class="fragopts" style="display:none"></div>
    </div>

    <div class="card">
      <div>
        <label for="cassette">Cassette (<code>-cstandard</code>)</label>
        <select id="cassette" onchange="cassetteChanged()">__CASSETTE_OPTIONS__</select>
      </div>
    </div>

    <details class="card advcard">
      <summary>Advanced PCR parameters</summary>
      <div class="row">
        <div>
          <label for="poly">Polymerase (<code>-poly</code>)</label>
          <select id="poly" onchange="polyChanged()" title="Sets the Tm model: buffer salt and salt correction follow the polymerase, as at NEB">
            <option value="phusion" selected>Phusion</option>
            <option value="taq">Taq (Standard Taq Buffer)</option>
          </select>
        </div>
        <div>
          <label for="pc">Primer &micro;M in the reaction (<code>-pc</code>) <span class="hint">0.05&ndash;5</span></label>
          <input id="pc" value="0.5" inputmode="decimal" oninput="this.dataset.touched='1'" title="Primer concentration in the reaction tube, in uM, as pipetted. 2.5 uL of a 10 uM stock in 50 uL is 0.5 uM.">
        </div>
      </div>
    </details>

    <div class="card" id="customcard" style="display:none">
      <label>Custom cassette <span class="hint">&mdash; paste the full cassette and its two fusion-tail overhangs (reverse-primer orientation, ending in <code>U</code>; a trailing <code>T</code> is accepted).</span></label>
      <div style="margin-top:8px">
        <label for="ccustfull" style="font-weight:600">Full cassette sequence <span class="hint">&mdash; located in the destination plasmid for <code>-dp</code></span></label>
        <textarea id="ccustfull" rows="2" placeholder="e.g. GCTGAGGGTTTAATTAAGGATCCTTAATTAAGACCTCAGC" style="margin-bottom:10px"></textarea>
        <div class="row">
          <div>
            <label for="ccust5" style="font-weight:600">5&prime; overhang</label>
            <input id="ccust5" placeholder="e.g. GGGTTTAAU">
          </div>
          <div>
            <label for="ccust3" style="font-weight:600">3&prime; overhang</label>
            <input id="ccust3" placeholder="e.g. GGTCTTAAU">
          </div>
          <div>
            <label for="cname" style="font-weight:600">Name <span class="hint">(optional)</span></label>
            <input id="cname" placeholder="e.g. myPacI">
          </div>
        </div>
      </div>
    </div>

    <div class="card">
      <label>Temperature optimization (<code>-tm</code>)</label>
      <div class="tmopt">
        <label class="radio"><input type="radio" name="tmmode" value="individual" checked onchange="tmModeChanged()"> Optimize individual primer pair Tm</label>
        <label class="radio"><input type="radio" name="tmmode" value="eq" onchange="tmModeChanged()"> Optimize equal primer Tm for all fragments</label>
        <label class="radio"><input type="radio" name="tmmode" value="value" onchange="tmModeChanged()"> Optimize all primer Tm to
          <input id="tmvalue" type="number" min="40" max="72" step="0.1" disabled
                 onfocus="selectTmValue()" oninput="selectTmValue()" placeholder="40&ndash;72"> &deg;C</label>
      </div>
      <p class="note">Please note: optimizing for one temperature might result in primers with sub-optimal properties due to priority of temperature.</p>
    </div>

    <div class="card">
      <label>Forced overhangs (<code>-fo</code>)
        <span class="hint">&mdash; one choice per junction between fragments. Updates automatically with your sequences.</span>
      </label>
      <div id="fo"><div class="foempty">Add two or more sequences to choose forced overhangs.</div></div>
    </div>

    <div class="card">
      <div class="lblrow">
        <label for="dp" style="margin:0">Full plasmid for reconstruction (<code>-dp</code>)
          <span class="hint">&mdash; optional. Paste the destination plasmid sequence, or upload GenBank (<code>.gb</code>) to transfer the backbone annotations too.</span>
        </label>
        <span style="display:flex;gap:8px">
          <button type="button" class="uploadbtn" onclick="document.getElementById('dpfile').click()">Upload FASTA</button>
          <button type="button" class="uploadbtn" onclick="document.getElementById('dpgbfile').click()">Upload GenBank</button>
        </span>
        <input type="file" id="dpfile" accept=".fa,.fasta,.fna,.ffn,.frn,.seq,.txt" style="display:none" onchange="loadFile(this,'dp')">
        <input type="file" id="dpgbfile" accept=".gb,.gbk,.genbank,.ape,.txt" style="display:none" onchange="loadGenbank(this,'plasmid')">
      </div>
      <textarea id="dp" rows="4" placeholder="(optional) full circular plasmid sequence" style="margin-top:8px"></textarea>
      <div id="gbplasmidstatus" class="saved"></div>
    </div>

    <div class="card">
      <label for="oname">Output file name (<code>-oname</code>)
        <span class="hint">&mdash; optional. Base name for the three output files (<code>_report.html</code>, <code>_insert.gb</code> / <code>_plasmid.gb</code>, <code>_primers.csv</code>), with no date/time stamp. Leave blank for the automatic time-stamped name. Re-using a name overwrites that run's files.</span>
      </label>
      <input type="text" id="oname" placeholder="(optional) e.g. pMR012_GFP_fusion" style="margin-top:8px">
    </div>

    <div id="msg"></div>
    <button id="go" onclick="run()">Design primers</button>
  </div>
</div>

<script>
var foState = {};  // junctionIndex(0-based) -> 'none'|'A'|'B'
var gbRecords = [];    // ordered GenBank fragment records: {name, block}. Non-empty => GenBank mode (-igb).
var gbPlasmid = null;  // uploaded GenBank text for the plasmid (-dpgb), or null

// ------------------------------------------------------------------
// Fragment model. The ordered fragment list is the single source of
// truth for cloning order. In FASTA mode the <textarea> text IS that
// list (one record per ">" header); reordering rewrites the textarea.
// In GenBank mode gbRecords[] is the list; reordering reorders it and
// the -igb file is rebuilt in that order at submit time.
// ------------------------------------------------------------------
function inGbMode(){ return gbRecords.length > 0; }

// Parse FASTA text into [{header, name, body}]. A headerless single
// sequence becomes one record named "seq 1"; otherwise one record per ">"
// header. `header` keeps the FULL header line (everything after the ">",
// so any id/description/coordinates survive a reorder); `name` is just its
// first token, used only for the short chip/junction labels.
function getFastaRecords(t){
  t = t.replace(/\r/g,'');
  if(t.trim() === '') return [];
  if(!/^>/m.test(t)) return [{header:null, name:'seq 1', body:t.trim()}];
  var recs = [], cur = null, lines = t.split('\n');
  for(var i=0;i<lines.length;i++){
    var ln = lines[i];
    if(/^>/.test(ln)){
      if(cur) recs.push(cur);
      var hdr = ln.replace(/^>/,'').replace(/\s+$/,'');   // full header, trailing ws only
      var nm = hdr.trim().split(/\s+/)[0];
      cur = {header: hdr, name: nm || ('seq '+(recs.length+1)), body:''};
    } else if(cur){
      cur.body += (cur.body ? '\n' : '') + ln;
    }
  }
  if(cur) recs.push(cur);
  recs.forEach(function(r){ r.body = r.body.replace(/^\s+|\s+$/g,''); });
  return recs;
}
// Serialise FASTA records back to text, preserving each full header line.
function fastaText(recs){
  return recs.map(function(r){
    var h = (r.header != null && r.header !== '') ? r.header : r.name;
    return '>'+h+'\n'+r.body;
  }).join('\n') + '\n';
}

// Split GenBank text into per-record blocks [{name, block}]. Each block
// runs from a LOCUS line to its "//" terminator, so a multi-record file
// (or several files) becomes individually re-orderable records.
function splitGenbank(text){
  text = text.replace(/\r/g,'');
  var recs = [], re = /^LOCUS[\s\S]*?^\/\/[ \t]*$/gm, m;
  while((m = re.exec(text)) !== null){
    var block = m[0];
    var nmM = block.match(/^LOCUS\s+(\S+)/);
    recs.push({name: nmM ? nmM[1] : ('seq '+(recs.length+1)), block: block});
  }
  if(recs.length === 0 && /^LOCUS\s/m.test(text)){   // no // terminator: treat whole file as one record
    var nm = (text.match(/^LOCUS\s+(\S+)/) || [])[1] || 'seq 1';
    recs.push({name: nm, block: text.trim()});
  }
  return recs;
}
// Concatenate the ordered GenBank records into one multi-record file.
function assembleIgb(){
  if(!gbRecords.length) return '';
  return gbRecords.map(function(r){ return r.block.replace(/\s+$/,''); }).join('\n') + '\n';
}

// The ordered fragment names (both modes), used by the diagram and junctions.
function orderedNames(){
  if(inGbMode()) return gbRecords.map(function(r){ return r.name; });
  return getFastaRecords(document.getElementById('fasta').value).map(function(r){ return r.name; });
}

// ---- multi-file uploads ------------------------------------------
// FASTA: append every selected file's records to the textarea (a
// headerless file is given a header from its filename so it stays a
// distinct record). GenBank upload takes precedence, so ignore while
// in GenBank mode.
function loadFastaFiles(input){
  var files = input.files; if(!files || !files.length) return;
  if(inGbMode()){ input.value=''; showMsg('warn','Clear the GenBank upload first to use plain FASTA.'); return; }
  var arr = Array.prototype.slice.call(files);
  var texts = new Array(arr.length), pending = arr.length;
  arr.forEach(function(f, idx){
    var rd = new FileReader();
    rd.onload = function(){
      var txt = String(rd.result).replace(/\r/g,'').trim();
      if(txt && !/^>/m.test(txt)){ txt = '>' + f.name.replace(/\.[^.]+$/,'') + '\n' + txt; }
      texts[idx] = txt;
      if(--pending === 0) appendFasta(texts);
    };
    rd.onerror = function(){ texts[idx] = ''; if(--pending === 0) appendFasta(texts); };
    rd.readAsText(f);
  });
  input.value = '';
}
function appendFasta(texts){
  var add = texts.filter(function(x){ return x && x.trim(); }).join('\n');
  if(!add) return;
  var ta = document.getElementById('fasta');
  ta.value = ta.value.trim() ? (ta.value.replace(/\s+$/,'') + '\n' + add + '\n') : (add + '\n');
  renderAll();
}

// GenBank fragments: append all records from every selected file, in order.
function loadFragGenbank(input){
  var files = input.files; if(!files || !files.length) return;
  var arr = Array.prototype.slice.call(files);
  var parts = new Array(arr.length), pending = arr.length;
  arr.forEach(function(f, idx){
    var rd = new FileReader();
    rd.onload = function(){ parts[idx] = splitGenbank(String(rd.result)); if(--pending === 0) addGbRecords(parts); };
    rd.onerror = function(){ parts[idx] = []; if(--pending === 0) addGbRecords(parts); };
    rd.readAsText(f);
  });
  input.value = '';
}
function addGbRecords(parts){
  var firstTime = gbRecords.length === 0, added = 0;
  parts.forEach(function(recs){ recs.forEach(function(r){ gbRecords.push(r); added++; }); });
  if(added === 0){
    if(firstTime) showMsg('err','No LOCUS records found in the uploaded GenBank file(s).');
    return;
  }
  var fa = document.getElementById('fasta'); fa.disabled = true; fa.style.opacity = .5;
  updateGbStatus();
  renderAll();
}
function updateGbStatus(){
  var nms = gbRecords.map(function(r){ return escapeHtml(r.name); }).join(', ');
  document.getElementById('gbfragstatus').innerHTML =
    'GenBank mode: <b>'+gbRecords.length+'</b> record'+(gbRecords.length>1?'s':'')+' &mdash; '+nms+
    '. Annotations will be transferred; the sequence box is ignored. '+
    '<button type="button" class="linkbtn" onclick="clearFragGenbank()">clear &amp; use plain text</button>';
}
function clearFragGenbank(){
  gbRecords = [];
  var fa = document.getElementById('fasta'); fa.disabled = false; fa.style.opacity = 1;
  document.getElementById('gbfragstatus').innerHTML = '';
  renderAll();
}

// Plasmid GenBank upload (-dpgb). One plasmid, so single file.
function loadGenbank(input, which){
  var f = input.files && input.files[0];
  if(!f) return;
  var rd = new FileReader();
  rd.onload = function(){
    gbPlasmid = rd.result;
    var dp = document.getElementById('dp'); dp.disabled = true; dp.style.opacity = .5;
    document.getElementById('gbplasmidstatus').innerHTML =
      'Plasmid GenBank loaded: <b>'+escapeHtml(f.name)+'</b> &mdash; backbone annotations will be transferred; the sequence box is ignored. '
      +'<button type="button" class="linkbtn" onclick="clearGenbank(\'plasmid\')">use plain text instead</button>';
  };
  rd.readAsText(f);
  input.value = '';
}
function clearGenbank(which){
  gbPlasmid = null;
  var dp = document.getElementById('dp'); dp.disabled = false; dp.style.opacity = 1;
  document.getElementById('gbplasmidstatus').innerHTML = '';
}

function onFastaInput(){ if(inGbMode()) return; renderAll(); }
function renderAll(){ renderDiagram(); rebuildJunctions(); renderFragOpts(); }

// ---- per-sequence options (trim start/stop; tag reserved for later) ----
var fragOpts = {};        // fragment name -> {trimStart, trimStop, tag}
var openFragCards = {};   // fragment name -> bool (expanded state, preserved across re-render)
var currentFragOrder = [];// names in current order, so option handlers avoid embedding names in markup
var STOPS = { TAA:1, TAG:1, TGA:1 };

// Predefined tags (sense strand). A 5' tag goes into the FW primer as-is; a 3'
// tag goes into the RV primer reverse-complemented. "Custom..." reveals a box.
// Sequences carry an A near the 5' end and a T near the 3' end so that, when a
// tag is ordered as a standalone annealed oligo, the flanking genes' USER
// (deoxyuridine) primers can regenerate the overhangs from the tag's own ends.
var TAGS = [
  {name:"His6", seq:"CACCATCATCACCATCAC"},
  {name:"FLAG", seq:"GACTACAAGGACGACGACGATAAG"},
  {name:"Strep-tag II", seq:"TGGAGTCACCCACAGTTTGAGAAG"},
  {name:"S-tag", seq:"AAGGAGACCGCCGCCGCCAAGTTCGAGCGTCAGCACATGGACTCT"},
  {name:"AviTag", seq:"GGATTGAACGACATCTTCGAGGCCCAGAAGATCGAGTGGCATGAG"},
  {name:"CBP", seq:"AAGCGTCGTTGGAAGAAGAACTTCATCGCCGTTTCTGCCGCCAACCGTTTCAAGAAGATCTCTTCTTCTGGTGCCTTG"},
  {name:"HA", seq:"TACCCATACGACGTTCCAGACTACGCT"},
  {name:"c-Myc", seq:"GAGCAGAAGTTGATCTCTGAGGAGGACTTG"},
  {name:"V5", seq:"GGAAAGCCAATCCCAAACCCATTGTTGGGTTTGGACTCTACT"},
  {name:"T7", seq:"ATGGCCTCTATGACCGGTGGTCAGCAGATGGGT"},
  {name:"ALFA-tag", seq:"AGTCGTTTGGAGGAGGAGTTGCGTCGTCGTTTGACTGAG"},
  {name:"Spot-tag", seq:"CCAGACCGTGTTCGTGCCGTTTCTCACTGGTCTTCT"},
  {name:"E-tag", seq:"GGAGCCCCAGTTCCATACCCAGACCCATTGGAGCCACGT"},
  {name:"SpyTag003", seq:"AGAGGTGTTCCACACATCGTTATGGTTGACGCCTACAAGCGTTACAAG"},
  {name:"SnoopTag", seq:"AAGTTGGGTGACATCGAGTTCATCAAGGTTAATAAG"},
  {name:"Sortase motif", seq:"TTACCAGAGACCGGT"},
  {name:"ybbR", seq:"GACTCTTTGGAGTTCATCGCCTCTAAGTTGGCT"},
  {name:"TEV", seq:"GAGAACTTGTACTTCCAGGGT"},
  {name:"PreScission/3C", seq:"TTAGAGGTTTTGTTCCAGGGTCCT"},
  {name:"Thrombin", seq:"TTAGTTCCACGTGGTTCT"},
  {name:"Factor Xa", seq:"ATCGAGGGTCGT"},
];
// Oligo viability: with 8-nt overhangs anchored at the terminal A (5') and T
// (3'), the double-stranded core between them must be >= 6 bp for the two
// synthetic strands to anneal into an orderable oligo. Returns false otherwise
// (then the "Clone as oligo" option is not shown at all).
var OLIGO_OVH = 8, OLIGO_MIN_CORE = 6;
function oligoViable(seq){
  var s = normSeq(seq); var a = s.indexOf('A'), t = s.lastIndexOf('T');
  if(a < 0 || t < 0) return false;
  return (t - a - (2*OLIGO_OVH - 1)) >= OLIGO_MIN_CORE;   // core bp between the overhangs
}
function tagSeqByName(n){ for(var i=0;i<TAGS.length;i++) if(TAGS[i].name===n) return TAGS[i].seq; return ''; }
// Resolve a selector value + custom boxes to {name, seq} (or null if none).
function resolveTag(sel, custom, customname){
  if(!sel || sel==='none') return null;
  if(sel==='__custom__'){
    var s=(custom||'').replace(/\s/g,'');
    if(!s) return null;
    var nm=(customname||'').trim() || 'custom tag';
    return { name:nm, seq:s };
  }
  return { name: sel, seq: tagSeqByName(sel) };
}

function normSeq(s){ return String(s == null ? '' : s).toUpperCase().replace(/[^A-Z]/g, ''); }
function startsATG(s){ return normSeq(s).slice(0,3) === 'ATG'; }
function endsStop(s){ var n = normSeq(s); return n.length >= 3 && !!STOPS[n.slice(-3)]; }
function getOpt(name){ return fragOpts[name] || (fragOpts[name] = { trimStart:false, trimStop:false, autofrag:false, tag5:'none', tag5custom:'', tag5customname:'', tag5atg:false, tag5oligo:false, tag5linker:'', tag5linkoligo:false, tag3:'none', tag3custom:'', tag3customname:'', tag3stop:false, tag3oligo:false, tag3linker:'', tag3linkoligo:false }); }
// A linker is DNA only.
function cleanLinker(s){ return normSeq(s).replace(/[^ACGT]/g, ''); }
// Mirrors the backend linker_in_oligo(): true when the linker can sit INSIDE the
// oligo, i.e. the oligo's gene-facing overhang is derived from the linker. For a
// 5' tag that is the 3' end (terminal T must be in the linker); for a 3' tag the
// 5' end (leading A must be in the linker). Also needs the combined tag+linker
// to be oligo-viable. When false the linker can only ride on the gene primer.
function linkerInOligo(tagseq, linker, end){
  var tg = normSeq(tagseq), lk = cleanLinker(linker);
  if(!lk) return false;
  var five = (end === 5 || end === '5');
  var combined = five ? (tg + lk) : (lk + tg);
  if(!oligoViable(combined)) return false;
  if(five){ return combined.lastIndexOf('T') >= tg.length; }
  var a = combined.indexOf('A'); return (a >= 0 && a < lk.length);
}
// Single source of truth for the oligo decision, used by the UI, badges and the
// payload so they never disagree. A tag can be cloned as an oligo when it is
// viable on its own OR when an included linker makes the combination viable and
// provides the gene-facing overhang (anchors). When the tag is viable only
// because of the linker, "include linker in oligo" is forced (the oligo can't
// exist without it). Returns: possible (offer "Clone as oligo"), oligo (it will
// be an oligo given the current flags), anchors (linker can sit in the oligo),
// showLink (offer the include-linker box), forced (include-linker is mandatory),
// linkoligo (effective include-linker value for the payload).
function oligoDecision(tagseq, linker, end, oligoOn, linkoligoOn){
  var tagViable = oligoViable(tagseq);
  var lk = cleanLinker(linker);
  var anchors = !!(lk && linkerInOligo(tagseq, lk, end));
  var possible = tagViable || anchors;
  var oligo = (oligoOn && possible) ? 1 : 0;
  var forced = !!(oligo && anchors && !tagViable);
  var showLink = !!(oligo && anchors);
  var linkoligo = (showLink && (forced || linkoligoOn)) ? 1 : 0;
  return { tagViable:tagViable, anchors:anchors, possible:possible,
           oligo:oligo, forced:forced, showLink:showLink, linkoligo:linkoligo };
}
// Add a start codon (ATG) before a 5' tag, or a stop (TAA) after a 3' tag, when opted in.
function applyTagCodon(t, end, on){
  if(!t || !on) return t;
  return (end === 5 || end === '5')
    ? { name:t.name, seq:'ATG' + t.seq }
    : { name:t.name, seq:t.seq + 'TAA' };
}

// Extract a GenBank record's ORIGIN sequence (uppercase, letters only).
function gbOriginSeq(block){
  var m = String(block).replace(/\r/g,'').match(/^ORIGIN[^\n]*\n([\s\S]*?)^\/\//m);
  return m ? normSeq(m[1]) : '';
}
// Ordered [{name, seq}] for the current fragments (both modes).
function fragmentSeqs(){
  if(inGbMode()) return gbRecords.map(function(r){ return { name:r.name, seq:gbOriginSeq(r.block) }; });
  return getFastaRecords(document.getElementById('fasta').value).map(function(r){ return { name:r.name, seq:normSeq(r.body) }; });
}

function renderFragOpts(){
  var frags = fragmentSeqs();
  currentFragOrder = frags.map(function(f){ return f.name; });
  var box = document.getElementById('fragopts');
  if(frags.length === 0){ box.style.display='none'; box.innerHTML=''; return; }
  var html = '<div class="cap">Per-sequence options &mdash; click a sequence to expand.</div>';
  for(var i=0;i<frags.length;i++){
    var name = frags[i].name, seq = frags[i].seq;
    var o = getOpt(name);
    var canATG = startsATG(seq), canStop = endsStop(seq);
    var canFrag = seq.length > 2000;      // offer auto-fragmentation only for large sequences
    if(!canATG) o.trimStart = false;      // keep state consistent with what's possible
    if(!canStop) o.trimStop = false;
    if(!canFrag) o.autofrag = false;
    var endL = seq.slice(0,3), endR = seq.slice(-3);
    var badges = [];
    if(o.trimStart) badges.push('&minus;ATG');
    if(o.trimStop)  badges.push('&minus;stop');
    if(o.autofrag)  badges.push('split &sim;1&nbsp;kb');
    var t5 = applyTagCodon(resolveTag(o.tag5, o.tag5custom, o.tag5customname), 5, o.tag5atg);
    var t3 = applyTagCodon(resolveTag(o.tag3, o.tag3custom, o.tag3customname), 3, o.tag3stop);
    if(t5){ var d5=oligoDecision(t5.seq,o.tag5linker,5,o.tag5oligo,o.tag5linkoligo); if(!d5.possible) o.tag5oligo=false; badges.push('5&prime;:'+escapeHtml(t5.name)+(o.tag5atg?'+ATG':'')+(cleanLinker(o.tag5linker)?'+linker':'')+(d5.oligo?' (oligo)':'')); }
    if(t3){ var d3=oligoDecision(t3.seq,o.tag3linker,3,o.tag3oligo,o.tag3linkoligo); if(!d3.possible) o.tag3oligo=false; badges.push('3&prime;:'+escapeHtml(t3.name)+(o.tag3stop?'+stop':'')+(cleanLinker(o.tag3linker)?'+linker':'')+(d3.oligo?' (oligo)':'')); }
    html += '<details class="fragopt" data-i="'+i+'"'+(openFragCards[name] ? ' open' : '')+'>'
      + '<summary><span class="idx">'+(i+1)+'</span><span class="nm" title="'+escapeHtml(name)+'">'+escapeHtml(name)+'</span>'
      + '<span class="badges">'+badges.join(' &middot; ')+'</span>'
      + '<span class="ends">'+escapeHtml(endL)+(seq.length>6?'&hellip;':'&middot;')+escapeHtml(endR)+'</span></summary>'
      + '<div class="optbody">'
      // only show a trim control when its codon is actually present
      + (canATG ? optChk(i,'trimStart',o.trimStart,'Trim start codon','ATG') : '')
      + (canStop ? optChk(i,'trimStop', o.trimStop, 'Trim stop codon', escapeHtml(endR)) : '')
      // auto-fragmentation: offered only for >2 kb sequences
      + (canFrag ? '<label class="chk tagoligo" title="Split this large sequence into ~1 kb sub-fragments that are USER-fused back into the full-length sequence. More reliable PCR and cloning for >2 kb fragments; the final assembly regenerates the original sequence.">'
          + '<input type="checkbox" data-i="'+i+'" data-key="autofrag"'+(o.autofrag?' checked':'')+'> Fragment for USER assembly '
          + '<span class="codon">'+seq.length+'&nbsp;bp &rarr; &sim;1&nbsp;kb pieces</span></label>' : '')
      + '</div>'
      + '<div class="optbody tagbody">'
      + tagSel(i,'5',o.tag5,o.tag5custom,o.tag5customname,o.tag5atg,o.tag5oligo,o.tag5linker,o.tag5linkoligo)
      + tagSel(i,'3',o.tag3,o.tag3custom,o.tag3customname,o.tag3stop,o.tag3oligo,o.tag3linker,o.tag3linkoligo)
      + '</div></details>';
  }
  box.innerHTML = html; box.style.display='block';
  // wire handlers via data-index (never embed arbitrary fragment names in markup/JS)
  box.querySelectorAll('input[type=checkbox][data-i]').forEach(function(cb){
    cb.addEventListener('change', function(){
      var nm = currentFragOrder[+cb.dataset.i]; if(nm == null) return;
      getOpt(nm)[cb.dataset.key] = cb.checked;
      renderFragOpts();
    });
  });
  box.querySelectorAll('details.fragopt[data-i]').forEach(function(d){
    d.addEventListener('toggle', function(){ var nm = currentFragOrder[+d.dataset.i]; if(nm != null) openFragCards[nm] = d.open; });
  });
  box.querySelectorAll('select[data-tagend]').forEach(function(sel){
    sel.addEventListener('change', function(){
      var nm = currentFragOrder[+sel.dataset.i]; if(nm == null) return;
      getOpt(nm)['tag'+sel.dataset.tagend] = sel.value;
      renderFragOpts();   // re-render to toggle the custom box and refresh badges
    });
  });
  box.querySelectorAll('input.tagcustom[data-tagend]').forEach(function(inp){
    inp.addEventListener('input', function(){
      var nm = currentFragOrder[+inp.dataset.i]; if(nm == null) return;
      getOpt(nm)['tag'+inp.dataset.tagend+'custom'] = inp.value;   // no re-render: keep focus
    });
    // On blur (change fires when the value changed and focus leaves), re-render so
    // the "Clone as oligo" checkbox appears/disappears as the custom sequence crosses
    // the viability threshold (oligoViable) - focus has already left, so no disruption.
    inp.addEventListener('change', function(){
      var nm = currentFragOrder[+inp.dataset.i]; if(nm == null) return;
      getOpt(nm)['tag'+inp.dataset.tagend+'custom'] = inp.value;
      renderFragOpts();
    });
  });
  box.querySelectorAll('input.tagcustomname[data-tagend]').forEach(function(inp){
    inp.addEventListener('input', function(){
      var nm = currentFragOrder[+inp.dataset.i]; if(nm == null) return;
      getOpt(nm)['tag'+inp.dataset.tagend+'customname'] = inp.value;   // no re-render: keep focus
    });
  });
  box.querySelectorAll('input.taglink[data-tagend]').forEach(function(inp){
    inp.addEventListener('input', function(){
      var nm = currentFragOrder[+inp.dataset.i]; if(nm == null) return;
      getOpt(nm)['tag'+inp.dataset.tagend+'linker'] = inp.value;   // no re-render: keep focus
    });
    // On blur, re-render so the "include linker in oligo" checkbox and the badge
    // appear/disappear as the linker crosses the anchor/viability threshold.
    inp.addEventListener('change', function(){
      var nm = currentFragOrder[+inp.dataset.i]; if(nm == null) return;
      getOpt(nm)['tag'+inp.dataset.tagend+'linker'] = inp.value;
      renderFragOpts();
    });
  });
}
function optChk(i,key,checked,label,codon){
  return '<label class="chk">'
    + '<input type="checkbox" data-i="'+i+'" data-key="'+key+'"'+(checked?' checked':'')+'> '
    + label+' <span class="codon">'+codon+'</span></label>';
}
function tagSel(i, end, sel, custom, customname, codonOn, oligoOn, linkerVal, linkoligoOn){
  var opts = '<option value="none"'+((!sel||sel==='none')?' selected':'')+'>None</option>';
  for(var k=0;k<TAGS.length;k++){
    var nm = TAGS[k].name;
    opts += '<option value="'+escapeHtml(nm)+'"'+(sel===nm?' selected':'')+'>'+escapeHtml(nm)+'</option>';
  }
  opts += '<option value="__custom__"'+(sel==='__custom__'?' selected':'')+'>Custom&hellip;</option>';
  var box = (sel==='__custom__')
    ? '<input class="tagcustomname" data-i="'+i+'" data-tagend="'+end+'" placeholder="name (for the .gb)" value="'+escapeHtml(customname||'')+'">'
      + '<input class="tagcustom" data-i="'+i+'" data-tagend="'+end+'" placeholder="sequence (sense)" value="'+escapeHtml(custom||'')+'">'
    : '';
  // Codon checkbox, only when a tag is actually selected: ATG before a 5' tag, stop (TAA) after a 3' tag.
  var codon = '', oligo = '', linkfield = '', linkoligo = '';
  if(sel && sel !== 'none'){
    var ckey = (end === '5') ? 'tag5atg' : 'tag3stop';
    var clabel = (end === '5') ? 'Add start codon <span class="codon">ATG</span>' : 'Add stop codon <span class="codon">TAA</span>';
    codon = '<label class="chk tagcodon"><input type="checkbox" data-i="'+i+'" data-key="'+ckey+'"'+(codonOn?' checked':'')+'> '+clabel+'</label>';
    // Optional linker between the tag and the gene (gene-facing side of the tag).
    var lkTitle = (end === '5') ? 'Linker inserted after the 5\u2032 tag, before the gene (DNA, sense)'
                                : 'Linker inserted before the 3\u2032 tag, after the gene (DNA, sense)';
    linkfield = '<input class="taglink" data-i="'+i+'" data-tagend="'+end+'" placeholder="linker (optional, DNA)" title="'+lkTitle+'" value="'+escapeHtml(linkerVal||'')+'">';
    // "Clone as oligo": shown when the tag is oligo-viable on its own OR when an
    // included linker makes the combination viable (oligoDecision).
    var tg = applyTagCodon(resolveTag(sel, custom, customname), end, codonOn);
    var d = tg ? oligoDecision(tg.seq, linkerVal, end, oligoOn, linkoligoOn) : null;
    if(d && d.possible){
      var okey = (end === '5') ? 'tag5oligo' : 'tag3oligo';
      oligo = '<label class="chk tagcodon tagoligo" title="Order this tag as a separate annealed FW/RV oligo instead of adding it to the gene primer">'
            + '<input type="checkbox" data-i="'+i+'" data-key="'+okey+'"'+(oligoOn?' checked':'')+'> Clone as oligo</label>';
      // "include linker in oligo": shown once the tag is cloned as an oligo and the
      // linker can form the gene-facing overhang. When the oligo is viable ONLY
      // because of the linker it is forced on and locked (disabled) so it can't be
      // unchecked into a non-viable oligo.
      if(d.showLink){
        var lkey = (end === '5') ? 'tag5linkoligo' : 'tag3linkoligo';
        var lchecked  = (d.forced || linkoligoOn) ? ' checked' : '';
        var ldisabled = d.forced ? ' disabled' : '';
        linkoligo = '<label class="chk tagcodon tagoligo'+(d.forced?' locked':'')+'" title="Include the linker in the ordered oligo (the oligo\u2019s gene-facing overhang is then derived from the linker). Unchecked: the linker goes on the gene primer and the oligo stays reusable across genes.">'
            + '<input type="checkbox" data-i="'+i+'" data-key="'+lkey+'"'+lchecked+ldisabled+'> include linker in oligo</label>';
      }
    }
  }
  return '<span class="tagwrap"><span class="taglbl">'+end+'&prime; tag</span>'
    + '<select data-i="'+i+'" data-tagend="'+end+'">'+opts+'</select>'+box+linkfield+codon+oligo+linkoligo+'</span>';
}

// Apply the active trims to a plain sequence (sense strand).
function applySeqTrims(seq, o){
  var s = normSeq(seq);
  if(o && o.trimStart && s.slice(0,3) === 'ATG') s = s.slice(3);
  if(o && o.trimStop && s.length >= 3 && STOPS[s.slice(-3)]) s = s.slice(0, -3);
  return s;
}
// Build the -i FASTA with per-fragment trims applied (headers preserved).
function buildTrimmedFasta(){
  var recs = getFastaRecords(document.getElementById('fasta').value);
  return recs.map(function(r){
    var o = fragOpts[r.name];
    var body = (o && (o.trimStart || o.trimStop)) ? applySeqTrims(r.body, o) : r.body;
    var h = (r.header != null && r.header !== '') ? r.header : r.name;
    return '>'+h+'\n'+body;
  }).join('\n') + '\n';
}
// Build the -igb GenBank with per-fragment trims applied (seq + coordinates).
// Collect per-fragment tags (in assembly order) for the payload.
function collectTags(){
  var names = orderedNames(), out = [];
  for(var i=0;i<names.length;i++){
    var o = fragOpts[names[i]]; if(!o) continue;
    var t5 = applyTagCodon(resolveTag(o.tag5, o.tag5custom, o.tag5customname), 5, o.tag5atg);
    if(t5){ var d5=oligoDecision(t5.seq,o.tag5linker,5,o.tag5oligo,o.tag5linkoligo);
      out.push({ frag:i+1, end:5, name:t5.name, seq:t5.seq, oligo:d5.oligo, linker:cleanLinker(o.tag5linker), linkoligo:d5.linkoligo }); }
    var t3 = applyTagCodon(resolveTag(o.tag3, o.tag3custom, o.tag3customname), 3, o.tag3stop);
    if(t3){ var d3=oligoDecision(t3.seq,o.tag3linker,3,o.tag3oligo,o.tag3linkoligo);
      out.push({ frag:i+1, end:3, name:t3.name, seq:t3.seq, oligo:d3.oligo, linker:cleanLinker(o.tag3linker), linkoligo:d3.linkoligo }); }
  }
  return out;
}
// 1-based indices (assembly order) of fragments marked for auto-fragmentation,
// gated by the same >2 kb rule the checkbox uses so a stale flag can't leak through.
function collectAutofrag(){
  var frags = fragmentSeqs(), out = [];
  for(var i=0;i<frags.length;i++){
    var o = fragOpts[frags[i].name];
    if(o && o.autofrag && frags[i].seq.length > 2000) out.push(i+1);
  }
  return out;
}
function buildTrimmedIgb(){
  if(!gbRecords.length) return '';
  return gbRecords.map(function(r){
    var o = fragOpts[r.name];
    var b = (o && (o.trimStart || o.trimStop)) ? gbApplyTrims(r.block, o) : r.block;
    return b.replace(/\s+$/,'');
  }).join('\n') + '\n';
}

// Trim a GenBank record: drop k5 (5') and/or k3 (3') nt, shift/clamp features,
// update the LOCUS length, and regenerate the ORIGIN block.
function gbApplyTrims(block, o){
  block = String(block).replace(/\r/g,'');
  var seq = gbOriginSeq(block);
  if(!seq) return block;
  var L = seq.length;
  var k5 = (o.trimStart && seq.slice(0,3) === 'ATG') ? 3 : 0;
  var k3 = (o.trimStop && L >= 3 && STOPS[seq.slice(-3)]) ? 3 : 0;
  if(!k5 && !k3) return block;
  var newSeq = seq.slice(k5, L - k3);
  var newLen = newSeq.length;
  block = gbEditFeatures(block, k5, k3, L, newLen);
  block = block.replace(/^(LOCUS\s+\S+\s+)(\d+)(\s+bp)/m, function(m,a,n,c){ return a + newLen + c; });
  block = block.replace(/^ORIGIN[^\n]*\n[\s\S]*?^\/\/[ \t]*$/m, gbOriginText(newSeq));
  return block;
}
function gbOriginText(seq){
  var out = 'ORIGIN\n', s = seq.toLowerCase(), i, j;
  for(i=0;i<s.length;i+=60){
    var chunk = s.slice(i,i+60), groups = [];
    for(j=0;j<chunk.length;j+=10) groups.push(chunk.slice(j,j+10));
    var num = String(i+1); while(num.length < 9) num = ' ' + num;
    out += num + ' ' + groups.join(' ') + '\n';
  }
  return out + '//';
}
// Rewrite each feature's location line within FEATURES..ORIGIN; drop features
// that fall entirely inside a trimmed region. (Single-line locations, matching
// what AMUSER's own GenBank parser accepts.)
function gbEditFeatures(block, k5, k3, L, newLen){
  var lines = block.split('\n'), out = [], i = 0;
  var featRe = /^(\s{1,10})(\S+)(\s+)((?:complement\(|join\(|order\(|<|>|\d).*)$/;
  for(; i<lines.length; i++){ out.push(lines[i]); if(/^FEATURES/.test(lines[i])){ i++; break; } }
  while(i < lines.length && !/^ORIGIN/.test(lines[i])){
    var fm = lines[i].match(featRe);
    if(fm){
      var res = transformLoc(fm[4], k5, k3, L, newLen);
      if(res === null){           // drop this feature and its qualifier lines
        i++;
        while(i < lines.length && !/^ORIGIN/.test(lines[i]) && !featRe.test(lines[i])) i++;
        continue;
      }
      out.push(fm[1] + fm[2] + fm[3] + res); i++;
    } else { out.push(lines[i]); i++; }
  }
  for(; i<lines.length; i++) out.push(lines[i]);
  return out.join('\n');
}
function transformLoc(loc, k5, k3, L, newLen){
  var nums = (loc.match(/\d+/g) || []).map(Number);
  if(!nums.length) return loc;
  var mn = Math.min.apply(null, nums), mx = Math.max.apply(null, nums);
  if(mx <= k5) return null;         // entirely within the 5' removed region
  if(mn > L - k3) return null;      // entirely within the 3' removed region
  return loc.replace(/\d+/g, function(d){
    var x = Number(d) - k5;
    if(x < 1) x = 1;
    if(x > newLen) x = newLen;
    return String(x);
  });
}

// ---- the ordered-fragment diagram (drag + buttons to reorder) -----
function renderDiagram(){
  var names = orderedNames();
  var box = document.getElementById('orderdiagram');
  if(names.length < 2){ box.style.display='none'; box.innerHTML=''; return; }
  var html = '<div class="chainlbl">Assembly order (5&prime; &rarr; 3&prime;) &mdash; drag a block, or use &#9664;/&#9654;, to reorder. This is the order the sequences are cloned in.</div><div class="chain">';
  for(var i=0;i<names.length;i++){
    if(i>0) html += '<span class="conn">&mdash;</span>';
    html += '<span class="chip" draggable="true" data-i="'+i+'"'
      + ' ondragstart="onDragStart(event,'+i+')" ondragover="onDragOver(event,'+i+')"'
      + ' ondragleave="onDragLeave(event)" ondrop="onDrop(event,'+i+')" ondragend="onDragEnd(event)">'
      + '<button type="button" class="mv" title="Move left"'+(i===0?' disabled':'')+' onclick="moveFragment('+i+','+(i-1)+')">&#9664;</button>'
      + '<span class="nm" title="'+escapeHtml(names[i])+'">'+escapeHtml(names[i])+'</span>'
      + '<button type="button" class="mv" title="Move right"'+(i===names.length-1?' disabled':'')+' onclick="moveFragment('+i+','+(i+1)+')">&#9654;</button>'
      + '<button type="button" class="rm" title="Remove this sequence" onclick="removeFragment('+i+')">&times;</button>'
      + '</span>';
  }
  html += '</div>';
  box.innerHTML = html; box.style.display='block';
}

var dragFrom = null;
function onDragStart(e, i){
  dragFrom = i;
  e.dataTransfer.effectAllowed = 'move';
  try { e.dataTransfer.setData('text/plain', String(i)); } catch(_){}
  e.currentTarget.classList.add('dragging');
}
function onDragOver(e, i){
  e.preventDefault(); e.dataTransfer.dropEffect = 'move';
  if(dragFrom !== null && i !== dragFrom) e.currentTarget.classList.add('over');
}
function onDragLeave(e){ e.currentTarget.classList.remove('over'); }
function onDrop(e, i){
  e.preventDefault();
  e.currentTarget.classList.remove('over');
  var from = dragFrom;
  if(from === null){ var d = parseInt(e.dataTransfer.getData('text/plain'),10); if(!isNaN(d)) from = d; }
  if(from === null) return;
  moveFragment(from, i);
  dragFrom = null;
}
function onDragEnd(e){
  e.currentTarget.classList.remove('dragging');
  var c = document.querySelectorAll('.chip.over'); for(var k=0;k<c.length;k++) c[k].classList.remove('over');
  dragFrom = null;
}

// Move the fragment at index `from` to index `to`, in whichever list is active.
function moveFragment(from, to){
  var names = orderedNames();
  if(from === to || from < 0 || to < 0 || from >= names.length || to >= names.length) return;
  if(inGbMode()){
    var it = gbRecords.splice(from,1)[0]; gbRecords.splice(to,0,it);
    updateGbStatus();
  } else {
    var recs = getFastaRecords(document.getElementById('fasta').value);
    var r = recs.splice(from,1)[0]; recs.splice(to,0,r);
    document.getElementById('fasta').value = fastaText(recs);
  }
  renderAll();
}
function removeFragment(i){
  if(inGbMode()){
    gbRecords.splice(i,1);
    if(!gbRecords.length){
      var fa = document.getElementById('fasta'); fa.disabled=false; fa.style.opacity=1;
      document.getElementById('gbfragstatus').innerHTML='';
    } else { updateGbStatus(); }
  } else {
    var recs = getFastaRecords(document.getElementById('fasta').value);
    recs.splice(i,1);
    document.getElementById('fasta').value = recs.length ? fastaText(recs) : '';
  }
  renderAll();
}

function rebuildJunctions(){
  var names = orderedNames();
  var nj = Math.max(names.length - 1, 0);
  var box = document.getElementById('fo');
  if(nj === 0){
    box.innerHTML = '<div class="foempty">Add two or more sequences to choose forced overhangs.</div>';
    return;
  }
  var html = '';
  for(var i=0;i<nj;i++){
    var sel = foState[i] || 'none';
    var up = names[i], dn = names[i+1];
    html += '<div class="juncrow">'
      + '<span class="jl">Junction '+(i+1)+' &nbsp;('+escapeHtml(up)+' &rarr; '+escapeHtml(dn)+')</span>'
      + '<select onchange="foState['+i+']=this.value">'
      + '<option value="none"'+(sel==='none'?' selected':'')+'>not forced</option>'
      + '<option value="A"'+(sel==='A'?' selected':'')+'>upstream seq. ('+escapeHtml(up)+')</option>'
      + '<option value="B"'+(sel==='B'?' selected':'')+'>downstream seq. ('+escapeHtml(dn)+')</option>'
      + '</select></div>';
  }
  box.innerHTML = html;
}

function polyChanged(){
  // Follow the polymerase's recommended primer concentration unless the user
  // has typed a value of their own.
  var defs = {phusion:'0.5', taq:'0.2'};
  var pc = document.getElementById('pc');
  if (pc.dataset.touched !== '1') pc.value = defs[document.getElementById('poly').value];
}
function cassetteChanged(){
  var custom = document.getElementById('cassette').value === '__CUSTOM__';
  document.getElementById('customcard').style.display = custom ? 'block' : 'none';
}

function tmModeChanged(){
  var on = document.querySelector('input[name=tmmode]:checked').value === 'value';
  document.getElementById('tmvalue').disabled = !on;
  if(on) document.getElementById('tmvalue').focus();
}
function selectTmValue(){
  document.querySelector('input[name=tmmode][value=value]').checked = true;
  document.getElementById('tmvalue').disabled = false;
}

// Read an uploaded file into a textarea so the user sees the sequence,
// then (for -i) refresh the junction list.
function loadFile(input, targetId){
  var f = input.files && input.files[0];
  if(!f) return;
  var rd = new FileReader();
  rd.onload = function(){
    document.getElementById(targetId).value = rd.result;
    if(targetId === 'fasta') rebuildJunctions();
  };
  rd.readAsText(f);
  input.value = '';   // allow re-selecting the same file later
}

function showMsg(cls, text){
  document.getElementById('msg').innerHTML =
    text ? '<div class="msg '+cls+'">'+escapeHtml(text)+'</div>' : '';
}
function escapeHtml(s){return String(s).replace(/[&<>]/g,function(c){
  return {'&':'&amp;','<':'&lt;','>':'&gt;'}[c];});}

function run(){
  var n = orderedNames().length;
  if(n === 0){ showMsg('err','Please enter at least one sequence.'); return; }
  showMsg('', '');

  // Temperature optimization -> -tm value
  var mode = document.querySelector('input[name=tmmode]:checked').value;
  var tm = '';
  if(mode === 'eq'){ tm = 'eq'; }
  else if(mode === 'value'){
    tm = document.getElementById('tmvalue').value.trim();
    if(tm === ''){ showMsg('err','Enter a temperature (40-72 \u00b0C), or choose another Tm option.'); return; }
    var v = parseFloat(tm);
    if(isNaN(v) || v < 40 || v > 72){ showMsg('err','Temperature must be between 40 and 72 \u00b0C.'); return; }
  }

  var nj = Math.max(n-1,0);
  var fo = [];
  for(var i=0;i<nj;i++){
    var sv = foState[i] || 'none';
    if(sv === 'A' || sv === 'B') fo.push({junction:i+1, side:sv});
  }
  var cassette = document.getElementById('cassette').value;
  if(cassette === '__CUSTOM__'){
    var cf = document.getElementById('ccustfull').value.trim();
    var c5 = document.getElementById('ccust5').value.trim();
    var c3 = document.getElementById('ccust3').value.trim();
    if(cf === '' || c5 === '' || c3 === ''){
      showMsg('err','Custom cassette: fill in the full cassette and both overhangs.'); return;
    }
  }

  var payload = {
    fasta: inGbMode() ? '' : buildTrimmedFasta(),
    igb: inGbMode() ? buildTrimmedIgb() : '',
    tags: collectTags(),
    autofrag: collectAutofrag(),
    cassette: cassette,
    ccustfull: document.getElementById('ccustfull').value,
    ccust5: document.getElementById('ccust5').value,
    ccust3: document.getElementById('ccust3').value,
    cname: document.getElementById('cname').value,
    tm: tm,
    poly: document.getElementById('poly').value,
    pc: document.getElementById('pc').value,
    fo: fo,
    dp: document.getElementById('dp').value,
    dpgb: gbPlasmid || '',
    oname: document.getElementById('oname').value
  };
  // Open the report tab NOW, inside the click gesture, so the browser does
  // not block it as a popup. It is pointed at the report once ready. The
  // report opens full-width in its own tab, avoiding the narrow-panel wrap.
  var rpt = window.open('about:blank', '_blank');

  var btn = document.getElementById('go');
  btn.disabled = true; btn.textContent = 'Designing...';
  fetch('/run',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)})
    .then(function(r){return r.json();})
    .then(function(d){
      btn.disabled=false; btn.textContent='Design primers';
      if(!d.ok){ if(rpt) rpt.close(); showMsg('err', d.error || 'Something went wrong.'); return; }
      var url = '/report/' + encodeURIComponent(d.name);
      if(rpt){ rpt.location.href = url; }
      var warn = (d.warn && d.warn.trim())
        ? '<div class="msg warn">'+escapeHtml(d.warn.trim())+'</div>' : '';
      var link = '<a href="'+url+'" target="_blank" rel="noopener">open the report</a>';
      document.getElementById('msg').innerHTML = warn +
        '<div class="msg ok">Report opened in a new tab' +
        (rpt ? '' : ' &mdash; your browser blocked the pop-up, click to '+link) +
        '.<div class="saved">Saved: '+escapeHtml(d.saved||'')+'</div></div>';
    })
    .catch(function(e){
      btn.disabled=false; btn.textContent='Design primers';
      if(rpt) rpt.close();
      showMsg('err','Could not reach the local AMUSER server. Is the launcher window still open?');
    });
}

renderAll();   // reflect any pre-filled content on load
</script>
</body>
</html>
HTML
    $html =~ s/__CASSETTE_OPTIONS__/$opts/;
    return $html;
}
