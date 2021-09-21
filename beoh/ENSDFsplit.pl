#!/usr/bin/perl

#  Perl Script to split a big single ENSDF file into each nucleus
#  This creates files those names are ENSDFzzzaaa.dat, where zzz and aaa
#  are the Z and A numbers. These files should be installed under 
#  /usr/local/share/coh/ENSDF
#  unless the default location was changed by the configure script.
#
#  usage:
#         ENSDFsplit.pl ENSDFfile
#

$

@element = (
  "NN",
  "H ","HE","LI","BE","B ","C ","N ","O ","F ","NE",
  "NA","MG","AL","SI","P ","S ","CL","AR","K ","CA",
  "SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN",
  "GA","GE","AS","SE","BR","KR","RB","SR","Y ","ZR",
  "NB","MO","TC","RU","RH","PD","AG","CD","IN","SN",
  "SB","TE","I ","XE","CS","BA","LA","CE","PR","ND",
  "PM","SM","EU","GD","TB","DY","HO","ER","TM","YB",
  "LU","HF","TA","W ","RE","OS","IR","PT","AU","HG",
  "TL","PB","BI","PO","AT","RN","FR","RA","AC","TH",
  "PA","U ","NP","PU","AM","CM","BK","CF","ES","FM",
  "MD","NO","LR","RF","DB","SG","BH","HS","MT","DS",
  "RG","12","Nh","14","15","16","17","18","19","20");

$i = 0;
while(<>){
    $line[$i++] = $_;
#   if(/^ \-/){
    if(/^ *$/){
        $mass = substr($line[0],0,3);
        $name = substr($line[0],3,2);
        for($j=0 ; $j<=$#element ; $j++){
            if($name eq $element[$j]){
                last;
            }
        }

        $fname = "ENSDF".sprintf("%03d%03d",$j,int($mass)).".dat";

        if($line[0] =~ /B\- DECAY/){
#           print $fname,$line[0],"\n";
            open(FILE,">$fname");
            for($j=0 ; $j<$i-1 ; $j++){
                $line[$j] =~ s/\s*$//;
                print FILE $line[$j],"\n";
            }
            close(FILE);
        }

        $i=0;
    }
}
