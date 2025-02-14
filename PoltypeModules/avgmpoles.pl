#!/usr/bin/perl 

if ($#ARGV != 5)
{
   die
     "usage $0 input_key_file input_xyz_file input_groups_file first_atom last_atom start_index_in_prm output_key_file output_xyz_file\n";
}

#This program takes the output of poledit.x (option 2) and averages multipoles.

$file      = $ARGV[0];
$xyzfile   = $ARGV[1];
$groupfile = $ARGV[2];

#Output file
$file_out     = $ARGV[3];
$xyz_file_out = $ARGV[4];

#Starting index of the first atom in parameter file
$index = $ARGV[5];

open(GFIL, "<$groupfile");

@groups = ();
$ii     = 0;
$nextid = 0;
my @grp_heads = ();
while (<GFIL>)
{
   $line = $_;
   @tmp = split(' ', $line);
#   push(@groups, [sort { $a <=> $a } @tmp]);
   my $groupid = shift @tmp;
   push @grp_heads, $tmp[0];
   $nextid += $#tmp + 1;
   foreach $atmidx (@tmp) {
      $typeid[$atmidx] = $groupid;
   }
}
close(GFIL);


#@groups = sort { ${$a}[0] <=> ${$b}[0] } @groups;
#
## @grp_heads maps group indices to the first member of each group, including singleton groups.
#@subtract = ();
#foreach $group (@groups)
#{
#   @subtract = (@subtract, @$group[1 .. $#$group]);
#}
#@dup{@subtract} = ();
#@grp_heads = grep { not exists $dup{$_} } ($first .. $last);
#
#$nextid = 1;
#$grpidx = 0;
#@typeid = ();
#
## @typeid maps atom number to new type (prior to offset with $index)
#foreach $atmnum ($first .. $last)
#{
#   if (not $typeid[$atmnum])
#   {
#      while ($atmnum > $groups[$grpidx][0] and $grpidx < $#groups)
#      {
#         foreach my $grpmem (@{$groups[$grpidx]})
#         {
#            if ($typeid[$grpmem])
#            {
#               print "error: $grpmem is defined in more than one group\n";
#            }
#            else
#            {
#               $typeid[$grpmem] = $nextid;
#            }
#         }
#         $grpidx++;
#         $nextid++;
#      }
#      if ($atmnum == $groups[$grpidx][0])
#      {
#         foreach my $grpmem (@{$groups[$grpidx]})
#         {
#            if ($typeid[$grpmem])
#            {
#               print "error: $grpmem is defined in more than one group\n";
#            }
#            else
#            {
#               $typeid[$grpmem] = $nextid;
#            }
#         }
#         $grpidx++;
#      }
#      else
#      {
#         $typeid[$atmnum] = $nextid;
#      }
#      $nextid++;
#   }

   #	else {print "not empty: $atmnum $typeid[$atmnum]\n";}
#}

#print "\n";
#foreach $atmnum($first..$last) {
#	print "$atmnum: $typeid[$atmnum]\n";
#}
#print "\n";
#
#for (my $x = 0; $x <= $#groups; $x++) {
#   for (my $y = 0; $y <= $#{$groups[$x]}; $y++) {
#      print "$groups[$x][$y]:"; } print "\n";}

#$start=1;
#$end=$last-$first+$start;
#$ntypes=(($end-$start+1)*5)+3;

$start  = 1;
$end    = $nextid;
$ntypes = (($end - $start + 1) * 5) + 3;

#open(FILE,$file);
#@lines=(<FILE>);
#close FILE;

open(XYZIN,  "<$xyzfile")      or die "Error: Cannot open $xyzfile";
open(XYZOUT, ">$xyz_file_out") or die "Error: Cannot open $xyz_file_out";
$line = <XYZIN>;
print XYZOUT "$line";
my $first = 1;
my $last = 0;
foreach $line (<XYZIN>)
{
   chomp $line;
   $line =~ m/(.*\.\d+\s+)(\d+)/;
   $atype = $typeid[$2];
   $line =~ s/(.*\.\d+\s+)(\d+) /$1 $atype/;
   print XYZOUT "$line\n";
   $last++;
}
close(XYZIN);
close(XYZOUT);

@lines = `grep -A4 '^multipole ' $file`;

open(OUT, ">$file_out");
my %atmnumtozthenbisectbool;

foreach $atmnum ($first .. $last)
{
   $atmnumtozthenbisectbool{$atmnum} = 0;
   $atmnumtozthenbisectbool{$typeid[$atmnum]} = 0;
   foreach $linenum (0 .. $#lines)
   {
      $anum = (split ' ', $lines[$linenum])[1];
      if ((split ' ', $lines[$linenum])[1] == $atmnum)
      {
#Bugfix 2013-8-2:" multipole ID c1 c2"
#when both c1 and c2 are negative (same as either is negative, meaning bisector frame), the $c1 also needs to flip the sign of tmpval
#previously the code only handles C2; $tmpval2 is nowadded
         #NOTE: If $c1 and $c2 exist, check $typeid values with existing values
         my @array1 = split(' ',$lines[$linenum]);
         my $size = scalar(@array1);
         if ($size==6) # case where there is bisector then z, otherwise it is normal
         {
            my $tmpval3 = (split ' ', $lines[$linenum])[4];

            $atmnumtozthenbisectbool{$atmnum} = 1;
         if ($tmpval3 < 0)
         {
            $c3[$typeid[$atmnum]] = -$typeid[-$tmpval3];
            #print "Content converted c3: $c3[$typeid[$atmnum]]\n";

         }
         else
         {
            $c3[$typeid[$atmnum]] = $typeid[$tmpval3];
            #print "Content converted c3: $c3[$typeid[$atmnum]]\n";
         }


            $e11[$typeid[$atmnum]] += (split ' ', $lines[$linenum])[5]; # this goes with if size==6
         }


         elsif ($size==4) # case where z only c3 is blank
         {
            $e11[$typeid[$atmnum]] += (split ' ', $lines[$linenum])[3]
         }

         else # 
         {
            $e11[$typeid[$atmnum]] += (split ' ', $lines[$linenum])[4]
         }

         $c1[$typeid[$atmnum]] = $typeid[(split ' ', $lines[$linenum])[2]];
         #print "Content converted c1: $c1[$typeid[$atmnum]]\n";
         my $tmpval = (split ' ', $lines[$linenum])[3]; # this corresponds to C2
         my $tmpval2 = (split ' ', $lines[$linenum])[2]; # corresponds to C1
         #print "line split ", (split ' ', $lines[$linenum]);
         #print "tmpval2 before ", $tmpval2;
         $c2[$typeid[$atmnum]] = $typeid[(split ' ', $lines[$linenum])[3]];
         if ($tmpval < 0)
         {
            $c2[$typeid[$atmnum]] = -$typeid[-$tmpval];
#            print "tmpval", $tmpval," ", $c1[$typeid[$atmnum]]," ",$c2[$typeid[$atmnum]],"\n";
         }
         else
         {
            $c2[$typeid[$atmnum]] = $typeid[$tmpval];
         }
         if ($tmpval2 < 0)
         {
            $c1[$typeid[$atmnum]] = -$typeid[-$tmpval2];
         }

         #$axdef[$typeid[$atmnum]]=(split ' ',$lines[$linenum])[4];
         
         $e21[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 1])[0];
         $e22[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 1])[1];
         $e23[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 1])[2];
         $e31[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 2])[0];
         $e41[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 3])[0];
         $e42[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 3])[1];
         $e51[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 4])[0];
         $e52[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 4])[1];
         $e53[$typeid[$atmnum]] += (split ' ', $lines[$linenum + 4])[2];
         $N[$typeid[$atmnum]]   += 1;
         #print "$atmnum $typeid[$atmnum] $N[$typeid[$atmnum]]\n";
      }
   }
}

@lines = `grep '^atom ' $file`;
foreach $line (@lines)
{
   chomp $line;
   $line =~ m/^atom\s+(\d+)/;

   #        my $atmnum = $1;
   #	print "$1 $typeid[$1] $groups[$1-1][0]\n";
   if (grep $1 eq $_, @grp_heads)
   {
      $atype = $typeid[$1];
      $line =~ s/^(atom\s+)(\d+)(\s+)(\d+)/$1$atype$3$atype/;
      print OUT "$line\n";
   }
}

print OUT "\n";
#foreach $typenum ($start .. $end)
foreach $typenum (@grp_heads)
{
   if ($N[$typeid[$typenum]])
   {
      $c0 = $typeid[$typenum]; # this is the same as anum in the last for loop

      $e11[$typeid[$typenum]] = $e11[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e21[$typeid[$typenum]] = $e21[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e22[$typeid[$typenum]] = $e22[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e23[$typeid[$typenum]] = $e23[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e31[$typeid[$typenum]] = $e31[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e41[$typeid[$typenum]] = $e41[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e42[$typeid[$typenum]] = $e42[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e51[$typeid[$typenum]] = $e51[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e52[$typeid[$typenum]] = $e52[$typeid[$typenum]] / $N[$typeid[$typenum]];
      $e53[$typeid[$typenum]] = $e53[$typeid[$typenum]] / $N[$typeid[$typenum]];
      if ($c2[$typeid[$typenum]] < 0 and $atmnumtozthenbisectbool{$typenum} == 0)  # if the C2 is negative due to bisector this is okay (why dont we also make y zero??) not for z-then-bisector or trisector
      {
         $e21[$typeid[$typenum]] = 0.0;
         $e51[$typeid[$typenum]] = 0.0;
      }
      if ($atmnumtozthenbisectbool{$typenum} == 1)
      {
      print (
                   "multipole",
                  "$c0",                 "$c1[$typeid[$typenum]]",
                  "$c2[$typeid[$typenum]]","$c3[$typeid[$typenum]]",
                 );
      printf OUT (
                  "%-9s%7d%7d%7d%7d%17f", "multipole",
                  "$c0",                 "$c1[$typeid[$typenum]]",
                  "$c2[$typeid[$typenum]]","$c3[$typeid[$typenum]]",       "$e11[$typeid[$typenum]]"
                 );
      
      }
      else
      {
      printf OUT (
                  "%-9s%7d%7d%7d%21.5f", "multipole",
                  "$c0",                 "$c1[$typeid[$typenum]]",
                  "$c2[$typeid[$typenum]]",       "$e11[$typeid[$typenum]]"
                 );
      }
      print OUT "\n";

      #print "multipole   $c0  $c1[$typeid[$typenum]]  $c2[$typeid[$typenum]]\t\t";
      #print "$e11[$typeid[$typenum]]\n";
      printf OUT (
                  "%37s%8.5f%11.5f%11.5f", " ",
                  "$e21[$typeid[$typenum]]",        "0.00000",
                  "$e23[$typeid[$typenum]]"
                 ); # makes the y component of dipole 0 because the dipole vector only needs to be described in a plane (the x and z plane)
      print OUT "\n";
      printf OUT ("%37s%8.5f", " ", "$e31[$typeid[$typenum]]"); # xx component of quadrupole
      print OUT "\n";
      printf OUT ("%37s%8.5f%11.5f", " ", "0.000", "$e42[$typeid[$typenum]]"); # makes the xy component of quadrupole zero
      print OUT "\n";
      printf OUT (
                  "%37s%8.5f%11.5f%11.5f", " ",
                  "$e51[$typeid[$typenum]]",        "0.000",
                  "$e53[$typeid[$typenum]]"
                 );
      print OUT "\n"; # makes the yz component of the quadrupole zero

      #print "\t\t\t\t\t$e21[$typenum]  0.000  $e23[$typenum]\n";
      #print "\t\t\t\t\t$e31[$typenum]\n";
      #print "\t\t\t\t\t0.000  $e42[$typenum]\n";
      #print "\t\t\t\t\t$e51[$typenum]  0.000  $e53[$typenum]\n";
   }
}

print OUT "\n";

@lines = `grep '^polarize ' $file`;
foreach $line (@lines)
{
   chomp $line;
   $line =~ m/^polarize\s+(\d+).*\.\d+((\s+\d+)*)/;
   my $atmnum = $1;
   my @grpwith = split /\s+/, $2;
   if (grep $1 eq $_, @grp_heads)
   {
      $atype = $typeid[$1];

      #Remove duplicate from @grpwith
      #Convert to new index
      foreach my $ii (1 .. $#grpwith)
      {
         $grpwith[$ii] = $typeid[$grpwith[$ii]];
      }
      my %seen = ();
      @grpwith = grep { !$seen{$_}++ } @grpwith;

      #print "@grpwith\n";

      $grpwithstr = join ' ', @grpwith;
      $line =~ s/^(polarize\s+)(\d+)(.*\.\d+)(\s+\d+)*/$1$atype$3$grpwithstr/;
      print OUT "$line\n";
   }
}

#NOTE: Check averages against each value for errors and print warnings.
