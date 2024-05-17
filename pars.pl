open(IN, 'pars.txt');
open(OUT, '>globalpars_c.h');
open(OUT2, '>globalpars.py');
print OUT "/* Global parameter file */\n";
print OUT2 "# Global parameter file\n";
while ($line = <IN>) {
  if ($line!~m/^\#/) {
    @data = split ' ', $line;
    if ((scalar @data)==2) {
      print OUT "\#define GLOBALPAR_$data[0] $data[1]\n";
      print OUT2 "GLOBALPAR_$data[0] = $data[1]\n";
    }
  }
}
close IN;
close OUT;
close OUT2;