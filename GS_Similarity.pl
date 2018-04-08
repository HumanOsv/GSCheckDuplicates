#!/usr/bin/env perl

# Author Code: Osvaldo Yañez Osses


use strict;
use warnings;
# install: 
#          sudo cpan Parallel::ForkManager
use Parallel::ForkManager;
use Benchmark;




#
# Global variable
my $numb_atoms;
#
my %Atomic_number = ( '89'  => 'Ac', '13' => 'Al', '95' => 'Am', '51'  => 'Sb',	
	              '18'  => 'Ar', '33' => 'As', '85' => 'At', '16'  => 'S',  
                      '56'  => 'Ba', '4'  => 'Be', '97' => 'Bk', '83'  => 'Bi',
		      '107' => 'Bh', '5'  => 'B' , '35' => 'Br', '48'  => 'Cd',	
	              '20'  => 'Ca', '98' => 'Cf', '6'  => 'C',  '58'  => 'Ce',
		      '55'  => 'Cs', '17'  => 'Cl','27'  => 'Co', '29'  => 'Cu',
		      '24'  => 'Cr', '96'  => 'Cm', '110' => 'Ds', '66'  => 'Dy',
		      '105' => 'Db', '99'  => 'Es', '68'  => 'Er', '21'  => 'Sc',
		      '50'  => 'Sn', '38'  => 'Sr', '63'  => 'Eu', '100' => 'Fm',
		      '9'   => 'F',  '15'  => 'P',  '87'  => 'Fr', '64'  => 'Gd',
		      '31'  => 'Ga', '32'  => 'Ge', '72'  => 'Hf', '108' => 'Hs',	
                      '2'   => 'He', '1'   => 'H',  '26'  => 'Fe', '67'  => 'Ho',
		      '49'  => 'In', '53'  => 'I',  '77'  => 'Ir', '70'  => 'Yb',
		      '39'  => 'Y',  '36'  => 'Kr', '57'  => 'La', '103' => 'Lr',
		      '3'   => 'Li', '71'  => 'Lu', '12'  => 'Mg', '25'  => 'Mn',	
                      '109' => 'Mt', '101' => 'Md', '80'  => 'Hg', '42'  => 'Mo',
		      '60'  => 'Nd', '10'  => 'Ne', '93'  => 'Np', '41'  => 'Nb',
		      '28'  => 'Ni', '7'   => 'N',  '102' => 'No', '79'  => 'Au',
		      '76'  => 'Os', '8'   => 'O',  '46'  => 'Pd', '47'  => 'Ag',
		      '78'  => 'Pt', '82'  => 'Pb', '94'  => 'Pu', '84'  => 'Po',
		      '19'  => 'K',  '59'  => 'Pr', '61'  => 'Pm', '91'  => 'Pa',
		      '88'  => 'Ra', '86'  => 'Rn', '75'  => 'Re', '45'  => 'Rh',
		      '37'  => 'Rb', '44'  => 'Ru', '104' => 'Rf', '62'  => 'Sm',
		      '106' => 'Sg', '34'  => 'Se', '14'  => 'Si', '11'  => 'Na',
		      '81'  => 'Tl', '73'  => 'Ta', '43'  => 'Tc', '52'  => 'Te',
		      '65'  => 'Tb', '22'  => 'Ti', '90'  => 'Th', '69'  => 'Tm',
		      '112' => 'Uub','116' => 'Uuh','111' => 'Uuu','118' => 'Uuo',
		      '115' => 'Uup','114' => 'Uuq','117' => 'Uus','113' => 'Uut',
		      '92'  => 'U',  '23'  => 'V',  '74'  => 'W',  '54'  => 'Xe',
                      '30'  => 'Zn', '40'  => 'Zr' );
					  
###################################
# read files
sub read_file {
	# filename
	my ($input_file) = @_;
	my @array = ();
	# open file
	open(FILE, "<", $input_file ) || die "Can't open $input_file: $!";
	while (my $row = <FILE>) {
		chomp($row);
		push (@array,$row);
	}
	close (FILE);
	# return array	
	return @array;
}
###################################
# Euclidean distance between points
sub Euclidean_distance {
	# array coords basin 1 and basin 2
	my ($p1,$p2,$p3, $axis_x, $axis_y, $axis_z) = @_;
	# variables
	my $x1 = $axis_x;
	my $y1 = $axis_y;
	my $z1 = $axis_z;
	# measure distance between two point
	my $dist = sqrt(
					($x1-$p1)**2 +
					($y1-$p2)**2 +
					($z1-$p3)**2
					); 
	return $dist;
}
###################################
# Between points
sub Mult_coords {
	# array coords basin 1 and basin 2
	my ($p1,$p2,$p3, $axis_x, $axis_y, $axis_z) = @_;
	# variables
	my $x1 = $axis_x;
	my $y1 = $axis_y;
	my $z1 = $axis_z;
	# measure distance between two point
	my $dist = 	$x1*$p1 +
				$y1*$p2 +
				$z1*$p3 ; 
	return $dist;
}
###################################
# Promedio
sub promedio {
	my ($num,$data) = @_;
	# write file
	my $sum = 0;
	for ( my $i = 0 ; $i < $num ; $i = $i + 1 ){
		$sum+= @$data[$i];
	}
	my $div = $sum / $num;
	return $div; 
}
###################################
# Grigoryan Springborg similitud
sub Grigoryan_Springborg {
	my ($numb_atoms,$array_coord_x_1,$array_coord_y_1, $array_coord_z_1,
	                $array_coord_x_2,$array_coord_y_2, $array_coord_z_2) = @_;
	#
	my @distance_alpha = ();
	my @distance_beta  = ();
	#
	my $sum_1 = 0;
	for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
		for ( my $j = 0 ; $j < $numb_atoms ; $j = $j + 1 ){
			if ( $i < $j ){
				my $distance = Euclidean_distance (@$array_coord_x_1[$i],@$array_coord_y_1[$i],@$array_coord_z_1[$i],
												@$array_coord_x_1[$j],@$array_coord_y_1[$j],@$array_coord_z_1[$j]);
				push (@distance_alpha,$distance);
			}
		}
	}
	#
	my $sum_2 = 0;
	for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
		for ( my $j = 0 ; $j < $numb_atoms ; $j = $j + 1 ){
			if ( $i < $j ){
				my $distance = Euclidean_distance (@$array_coord_x_2[$i],@$array_coord_y_2[$i],@$array_coord_z_2[$i],
												@$array_coord_x_2[$j],@$array_coord_y_2[$j],@$array_coord_z_2[$j]);
				push (@distance_beta,$distance);
			}
		}
	}
	#
	my $InterDist_1 = (2/($numb_atoms*($numb_atoms-1)));
	my $InterDist_2 = (($numb_atoms*($numb_atoms-1))/2);  
	#
	my @mol_alpha = ();
	my @mol_beta  = ();
	my @idx_1 = sort { $distance_alpha[$a] <=> $distance_alpha[$b] } 0 .. $#distance_alpha;
	my @idx_2 = sort { $distance_beta[$a]  <=> $distance_beta[$b]  } 0 .. $#distance_beta;
	@mol_alpha = @distance_alpha[@idx_1];
	@mol_beta  = @distance_beta[@idx_2];
	#
	my $num_1 = scalar (@mol_alpha);
	my $num_2 = scalar (@mol_beta);
	my $dim_alpha =  promedio ($num_1,\@mol_alpha);
	my $dim_beta  =  promedio ($num_2,\@mol_beta);
	#
	my $sumX;
	my $sumY;
	# Sin normalizar
	for ( my $i = 0 ; $i < $InterDist_2 ; $i = $i + 1 ){
		my $mult = ( $mol_alpha[$i] - $mol_beta[$i] )**2;
		$sumX+=$mult;
	}
	my $Springborg_1 = sqrt( $InterDist_1 * $sumX );
	# Normalizado
	for ( my $i = 0 ; $i < $InterDist_2 ; $i = $i + 1 ){
		my $mult = ( ($mol_alpha[$i]/$dim_alpha) - ($mol_beta[$i]/$dim_beta) )**2;
		$sumY+=$mult;
	}
	my $Springborg_2 = sqrt( $InterDist_1 * $sumY );
	#
	return $Springborg_2; 
}
###################################
# delete repeat data
sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}
###################################
# index duplicate data
sub index_elements {
	my ($duplicate_name,$files_name) = @_;
	# reference arrays	
	my @array_1     = @{$duplicate_name}; 
	my @array_2     = @{$files_name};
	my @array_index = ();
	#
	my @filtered = uniq(@array_1);
	foreach my $u (@filtered){
		my @del_indexes = reverse( grep { $array_2[$_] eq "$u" } 0..$#array_2);
		foreach my $k (@del_indexes) {
			push (@array_index,$k);
		}
	}
	return @array_index;
}
#
my ($threshold_duplicate) = @ARGV;
if (not defined $threshold_duplicate) {
	die "\nGrigoryan Springborg Similarity must be run with:\n\nUsage:\n\tGS_Similarity.pl [threshold duplicate]\n";
	exit(1);  
}
###################################
# MAIN
#
my $tiempo_inicial  = new Benchmark; #funcion para el tiempo de ejecucion del programa
#
my %Info_Coords = ();
my @array_keys  = ();
#
my @files_xyz   = glob "./XYZ_data/*.xyz";
#
for (my $i = 0 ; $i < scalar (@files_xyz) ; $i++) {
	my @data    = read_file ($files_xyz[$i]);
	$numb_atoms = $data[0];
	shift(@data);
	shift(@data);
	my $id      = sprintf("%06d",$i);
	$Info_Coords{$id} = \@data;
	push(@array_keys,$id);	
}
#
my $ncpus     = 100;
my $pm        = new Parallel::ForkManager($ncpus);
my $iteration = 0;
#
my $file_tmp = "Dupli.tmp";
open (FILE, ">$file_tmp") or die "Unable to open XYZ file: $file_tmp";
my $file_log = "Info_Duplicates.txt";
open (LOGDUPLI, ">$file_log") or die "Unable to open XYZ file: $file_log"; 
print LOGDUPLI "\n# # # SUMMARY SIMILAR STRUCTURES # # #\n\n";
for ( my $x = 0 ; $x < scalar (@array_keys); $x = $x + 1 ) {
	$pm->start($iteration) and next;
	# All children process havee their own random.			
	srand();
	for ( my $y = 0 ; $y < scalar (@array_keys); $y = $y + 1 ) {
		if ( $x < $y ){
			#
			my @matrix_1 = @{$Info_Coords{$array_keys[$x]}};
			my @matrix_2 = @{$Info_Coords{$array_keys[$y]}};
			# # # # # # # # # # # # # # # # #
			#
			my @array_name_atoms_1 = ();
			my @array_coord_x_1    = ();
			my @array_coord_y_1    = ();
			my @array_coord_z_1    = ();
			#
			my @array_name_atoms_2 = ();
			my @array_coord_x_2    = ();
			my @array_coord_y_2    = ();
			my @array_coord_z_2    = ();
			#
			for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
				my @array_tabs_1  = split (/\s+/,$matrix_1[$i]);
				#
				my $radii_val;
				my $other_element = 0;
				if ( exists $Atomic_number{$array_tabs_1[0]} ) {
					# exists
					$radii_val = $Atomic_number{$array_tabs_1[0]};
					$array_name_atoms_1[++$#array_name_atoms_1] = $radii_val;
				} else {
					# not exists
					$radii_val = $array_tabs_1[0] ;
					$array_name_atoms_1[++$#array_name_atoms_1] = $radii_val;
				}
				$array_coord_x_1[++$#array_coord_x_1]   = $array_tabs_1[1];
				$array_coord_y_1[++$#array_coord_y_1]   = $array_tabs_1[2];
				$array_coord_z_1[++$#array_coord_z_1]   = $array_tabs_1[3];
			}
			#
			for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
				my @array_tabs_2 = split (/\s+/,$matrix_2[$i]);
				#
				my $radii_val;
				my $other_element = 0;
				if ( exists $Atomic_number{$array_tabs_2[0]} ) {
					# exists
					$radii_val = $Atomic_number{$array_tabs_2[0]};
					$array_name_atoms_2[++$#array_name_atoms_2] = $radii_val;
				} else {
					# not exists
					$radii_val = $array_tabs_2[0] ;
					$array_name_atoms_2[++$#array_name_atoms_2] = $radii_val;
				}
				$array_coord_x_2[++$#array_coord_x_2]   = $array_tabs_2[1];
				$array_coord_y_2[++$#array_coord_y_2]   = $array_tabs_2[2];
				$array_coord_z_2[++$#array_coord_z_2]   = $array_tabs_2[3];
			}
			my $Springborg = Grigoryan_Springborg ($numb_atoms,\@array_coord_x_1 ,\@array_coord_y_1 ,\@array_coord_z_1 
			                                                  ,\@array_coord_x_2 ,\@array_coord_y_2 ,\@array_coord_z_2 );												  
			#
			if ( $Springborg < $threshold_duplicate ) {
				my $number      = sprintf '%.6f', $Springborg;
				print FILE "$array_keys[$y]\n";
				print FILE "Value = $number\n";
				print LOGDUPLI "# $files_xyz[$x] ~= $files_xyz[$y]\n";
				print LOGDUPLI "# Value = $number\n";
				print LOGDUPLI "------------------------\n";
			}
			#
		}
		$iteration++;
	}
	$pm->finish;	
}
close (FILE);
# Paralel
$pm->wait_all_children;
# # #
my @data_tmp = read_file ($file_tmp);
my @duplicates_name = ();
my @Value_simi      = ();
foreach my $info (@data_tmp) {
	if ( ($info =~ m/Value/) ) {
		my @array_tabs = ();
		@array_tabs    = split ('\s+',$info);
		push (@Value_simi,$array_tabs[2]);
	} else {
		push (@duplicates_name,$info);
	}
}
# Delete similar structures
my @index_files = index_elements (\@duplicates_name,\@array_keys);
my $file_xyz = "02Duplicates_coords.xyz";
open (DUPLIXYZ, ">$file_xyz") or die "Unable to open XYZ file: $file_log"; 
my $count_sim_struc = 0;
foreach my $id_index (@index_files) {
	my @coords_dup = @{$Info_Coords{$array_keys[$id_index]}};
	print DUPLIXYZ scalar (@coords_dup);
	print DUPLIXYZ "\n";
	print DUPLIXYZ "Duplicate Structure $files_xyz[$id_index]\n";
	foreach (@coords_dup) {
		print DUPLIXYZ "$_\n";
	}
	$count_sim_struc++;
}
print LOGDUPLI "\nNumber of Similar Structures = $count_sim_struc\n";
close (LOGDUPLI);  
close (DUPLIXYZ);
#
# Delete similar structures
for my $k (@index_files) {
	delete $Info_Coords{$array_keys[$k]};
}
my $file_x = "01Clean_Duplicates_coords.xyz";
open (CLEANDUPLIXYZ, ">$file_x") or die "Unable to open XYZ file: $file_x"; 
foreach my $key (sort(keys %Info_Coords)) {
	#
	my @new_matrix = @{$Info_Coords{$key}};
	print CLEANDUPLIXYZ scalar (@new_matrix);
	print CLEANDUPLIXYZ "\n";
	print CLEANDUPLIXYZ "Unique Structure $files_xyz[$key]\n";
	foreach my $u (@new_matrix) {
		print CLEANDUPLIXYZ "$u\n";
	}
}
close (CLEANDUPLIXYZ);
#
my $tiempo_final  = new Benchmark;
my $tiempo_total  = timediff($tiempo_final, $tiempo_inicial);
print "\n\tExecution Time: ",timestr($tiempo_total),"\n";
print "\n";
#
unlink ($file_tmp);
exit 0;
