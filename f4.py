#!/usr/local/bin/python3

# Michael Matschiner, 2015-05-17
# michaelmatschiner@mac.com

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse, textwrap, os, scipy, numpy, tempfile, random
from scipy import stats
from subprocess import call

# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    ----------------------------------------------------------------------
      Calculate the f4 statistic from SNP data in Treemix format and assess
      support for introgression by jackknifing and coalescent simulations.
      For format specification, see the Treemix manual at
      https://bitbucket.org/nygcresearch/treemix/downloads.
      The f4 statistic was originally described by Reich et al. (2009): 
      http://www.nature.com/nature/journal/v461/n7263/abs/nature08365.html.
      Simulations are performed using fastsimcoal2, which must be installed
      and executable with 'fsc252': http://cmpg.unibe.ch/software/fastsimcoal2.
    '''))
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.9'
    )
parser.add_argument(
    '-k',
    nargs=1,
    type=int,
    default=[-1],
    dest='jackknife_k',
    help="Number of SNPs per block for estimation of standard errors by jackknifing (default: off)."
    )
parser.add_argument(
    '-s',
    nargs=1,
    type=int,
    default=[-1],
    dest='number_of_simulations',
    help="Number of simulations (default: off)."
    )
parser.add_argument(
    '-o',
    nargs=1,
    type=str,
    default=["-1"],
    dest='output_dir',
    help="Name of a directory to which simulated data sets should be written (in treemix format) (default: off)."
    )
parser.add_argument(
    'infile',
    nargs='?',
    type=argparse.FileType('r'),
    default='-',
    help='The input file name.'
    )
parser.add_argument(
    'outfile',
    nargs='?',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='The output file name.'
    )

# Get the command line arguments.
args = parser.parse_args()
jackknife_k = args.jackknife_k[0]
number_of_simulations = args.number_of_simulations[0]
output_dir = args.output_dir[0]
infile = args.infile
outfile = args.outfile
if infile.isatty():
    print('No input file specified, and no input piped through stdin!')
    sys.exit(0)
if output_dir != "-1" and number_of_simulations < 1:
    print("")
    print("WARNING: An output directory for simulated data sets has been specified, but no simulations are performed!")

# Prepare the output string.
outfile.write("\n")
outfile.write("  f4.py\n")
outfile.write("----------------------------------------------------------------------\n")
outfile.write("  Michael Matschiner | 2015-05-17 | evoinformatics.eu\n")
outfile.write("\n")
outfile.write("Input\n")
if infile.name == '<stdin>':
    outfile.write("  File name: STDIN\n")
else:
    outfile.write("  File name: " + str(infile.name) + "\n")

# Read the input and calculate the f4 statistic.
instring = infile.read()
inlines = instring.split('\n')
header_line = inlines[0]
body_lines = inlines[1:]
# Remove empty lines and SNPs with completely missing data in one population, but count them before.
valid_body_lines = []
number_of_snps = 0
for body_line in body_lines:
    if body_line.strip() != "":
        number_of_snps += 1
        if "0,0" not in body_line:
            valid_body_lines.append(body_line)
body_lines = valid_body_lines
pops = header_line.strip().split()
if len(pops) != 4:
    print("ERROR: The input file should contain data for exactly four populations. The header line has only " + str(len(pops)) + " columns!")
    sys.exit(1)
outfile.write("  Assumed population topology: (" +  pops[0] + "," + pops[1] + "),(" + pops[2] + "," + pops[3] + ")\n")

f4s = []
number_of_valid_snps = 0
number_of_alleles_per_snp_a = []
number_of_alleles_per_snp_b = []
number_of_alleles_per_snp_c = []
number_of_alleles_per_snp_d = []
number_of_snps_variable_in_more_than_one_population = 0
number_of_snps_variable_on_both_sides_of_the_root = 0
for x in range(len(body_lines)):
    line_ary = body_lines[x].strip().split()
    if len(line_ary) != 4:
        print("ERROR: The input file should contain data for exactly four populations. Line " + str(x+2) + " has only " + str(len(line_ary)) + " columns!")
        sys.exit(1)
    a_ary = line_ary[0].split(",")
    b_ary = line_ary[1].split(",")
    c_ary = line_ary[2].split(",")
    d_ary = line_ary[3].split(",")
    number_of_alleles_a = int(a_ary[0])+int(a_ary[1])
    number_of_alleles_b = int(b_ary[0])+int(b_ary[1])
    number_of_alleles_c = int(c_ary[0])+int(c_ary[1])
    number_of_alleles_d = int(d_ary[0])+int(d_ary[1])
    if min(number_of_alleles_a,number_of_alleles_b,number_of_alleles_c,number_of_alleles_d) > 0:
        number_of_alleles_per_snp_a.append(number_of_alleles_a)
        number_of_alleles_per_snp_b.append(number_of_alleles_b)
        number_of_alleles_per_snp_c.append(number_of_alleles_c)
        number_of_alleles_per_snp_d.append(number_of_alleles_d)
        number_of_populations_with_more_than_one_allele_for_this_snp = 0
        variable_in_a_or_b = False
        a = int(a_ary[0])/number_of_alleles_a
        if a > 0 and a < 1:
            number_of_populations_with_more_than_one_allele_for_this_snp += 1
            variable_in_a_or_b = True
        b = int(b_ary[0])/number_of_alleles_b
        if b > 0 and b < 1:
            number_of_populations_with_more_than_one_allele_for_this_snp += 1
            variable_in_a_or_b = True
        variable_in_c_or_d = False
        c = int(c_ary[0])/number_of_alleles_c
        if c > 0 and c < 1:
            number_of_populations_with_more_than_one_allele_for_this_snp += 1
            variable_in_c_or_d = True
        d = int(d_ary[0])/number_of_alleles_d
        if d > 0 and d < 1:
            number_of_populations_with_more_than_one_allele_for_this_snp += 1
            variable_in_c_or_d = True
        f4s.append((a-b)*(c-d))
        number_of_valid_snps += 1
        if number_of_populations_with_more_than_one_allele_for_this_snp > 1:
            number_of_snps_variable_in_more_than_one_population += 1
        if variable_in_a_or_b and variable_in_c_or_d:
            number_of_snps_variable_on_both_sides_of_the_root += 1
    else:
        print("ERROR: One of the four populations has no allele in this line: ")
        print(lines[x].strip())
        sys.exit(1)

# Report data set statistics.
max_number_of_alleles_a = max(number_of_alleles_per_snp_a)
max_number_of_alleles_b = max(number_of_alleles_per_snp_b)
max_number_of_alleles_c = max(number_of_alleles_per_snp_c)
max_number_of_alleles_d = max(number_of_alleles_per_snp_d)
outfile.write("  Number of sequences per population: ")
outfile.write(str(max_number_of_alleles_a) + ", ")
outfile.write(str(max_number_of_alleles_b) + ", ")
outfile.write(str(max_number_of_alleles_c) + ", ")
outfile.write(str(max_number_of_alleles_d) + "\n")
outfile.write("  Number of SNPs: " + str(number_of_snps) + "\n")
outfile.write("  Number of SNPs excluded due to missing data: " + str(number_of_snps-number_of_valid_snps) + "\n")
outfile.write("  Number of SNPs used: " + str(number_of_valid_snps) + "\n")
proportion_missing_per_snp_a = []
for n in number_of_alleles_per_snp_a:
    proportion_missing_per_snp_a.append(1-(n/max_number_of_alleles_a))
proportion_missing_a = sum(proportion_missing_per_snp_a)/max(len(proportion_missing_per_snp_a),1)
proportion_missing_per_snp_b = []
for n in number_of_alleles_per_snp_b:
    proportion_missing_per_snp_b.append(1-(n/max_number_of_alleles_b))
proportion_missing_b = sum(proportion_missing_per_snp_b)/max(len(proportion_missing_per_snp_b),1)
proportion_missing_per_snp_c = []
for n in number_of_alleles_per_snp_c:
    proportion_missing_per_snp_c.append(1-(n/max_number_of_alleles_c))
proportion_missing_c = sum(proportion_missing_per_snp_c)/max(len(proportion_missing_per_snp_c),1)
proportion_missing_per_snp_d = []
for n in number_of_alleles_per_snp_d:
    proportion_missing_per_snp_d.append(1-(n/max_number_of_alleles_d))
proportion_missing_d = sum(proportion_missing_per_snp_d)/max(len(proportion_missing_per_snp_d),1)
proportion_missing_max = max(proportion_missing_a,proportion_missing_b,proportion_missing_c,proportion_missing_d)
outfile.write("  Proportion of missing data per population: ")
outfile.write("{0:.2f}".format(proportion_missing_a) + ", ")
outfile.write("{0:.2f}".format(proportion_missing_b) + ", ")
outfile.write("{0:.2f}".format(proportion_missing_c) + ", ")
outfile.write("{0:.2f}".format(proportion_missing_d) + "\n")
proportion_of_snps_variable_in_more_than_one_population = number_of_snps_variable_in_more_than_one_population/number_of_valid_snps
outfile.write("  Proportion of SNPs variable in more than one population: " + "{0:.2f}".format(proportion_of_snps_variable_in_more_than_one_population) + "\n")
proportion_of_snps_variable_on_both_sides_of_the_root = number_of_snps_variable_on_both_sides_of_the_root/number_of_valid_snps
outfile.write("  Proportion of SNPs variable on both sides of the root: " + "{0:.2f}".format(proportion_of_snps_variable_on_both_sides_of_the_root) + "\n")
observed_f4 = numpy.mean(f4s)
outfile.write("  Observed f4: " + "{0:.5f}".format(observed_f4) + "\n")
outfile.write("\n")

# If a k value has been specified, use it for block jack-knifing.
if jackknife_k != -1:
    outfile.write("Block jackknifing\n")
    if jackknife_k > len(body_lines):
        print("ERROR: The specified jackknifing block size (option 'k') is greater than the number of SNPs!")
        sys.exit(1)
    elif jackknife_k < 1:
        print("ERROR: The specified jackknifing block size (option '-k') is smaller than 1!")
        sys.exit(1)
    body_lines_per_block = [[]]
    for x in range(len(body_lines)):
        if len(body_lines_per_block[-1]) < jackknife_k:
            body_lines_per_block[-1].append(body_lines[x])
        else:
            body_lines_per_block.append([body_lines[x]])
    body_lines_per_valid_block = []
    for item in body_lines_per_block:
        if len(item) == jackknife_k:
            body_lines_per_valid_block.append(item)
    if len(body_lines_per_valid_block) < 2:
        print("ERROR: he specified jack-knifing block size (k) allows only a single block of SNPs!")
        sys.exit(1)
    f4_per_block = []
    for block_body_lines in body_lines_per_valid_block:
        block_f4s = []
        for x in range(len(block_body_lines)):
            if block_body_lines[x].strip() != "":
                line_ary = block_body_lines[x].strip().split()
                a_ary = line_ary[0].split(",")
                b_ary = line_ary[1].split(",")
                c_ary = line_ary[2].split(",")
                d_ary = line_ary[3].split(",")
                number_of_alleles_a = int(a_ary[0])+int(a_ary[1])
                number_of_alleles_b = int(b_ary[0])+int(b_ary[1])
                number_of_alleles_c = int(c_ary[0])+int(c_ary[1])
                number_of_alleles_d = int(d_ary[0])+int(d_ary[1])
                if min(number_of_alleles_a,number_of_alleles_b,number_of_alleles_c,number_of_alleles_d) > 0:
                    a = int(a_ary[0])/number_of_alleles_a
                    b = int(b_ary[0])/number_of_alleles_b
                    c = int(c_ary[0])/number_of_alleles_c
                    d = int(d_ary[0])/number_of_alleles_d
                    block_f4s.append((a-b)*(c-d))
                else:
                    print("ERROR: One of the four populations has no allele in this line: ")
                    print(block_body_lines[x].strip())
                    sys.exit(1)
        f4_per_block.append(numpy.mean(block_f4s))
    f4_per_block_mean = numpy.mean(f4_per_block)

    # Report Jack-knife results.
    outfile.write("  Length of jackknife blocks: " + str(jackknife_k) + "\n")
    outfile.write("  Number of jackknife blocks: " + str(len(body_lines_per_valid_block)) + "\n")
    outfile.write("  Jackknife blocks f4 mean: " + "{0:.5f}".format(f4_per_block_mean) + "\n")
    observed_jackknife_f4_standard_error = stats.sem(f4_per_block)
    outfile.write("  Jackknife blocks f4 standard error: " + "{0:.5f}".format(observed_jackknife_f4_standard_error) + "\n")
    observed_jackknife_f4_z_zcore = observed_f4/observed_jackknife_f4_standard_error
    outfile.write("  Jackknife blocks f4 z-value: " + "{0:.4f}".format(observed_jackknife_f4_z_zcore) + "\n")
    observed_jackknife_f4_p_value = scipy.stats.norm.cdf(observed_jackknife_f4_z_zcore)
    outfile.write("  Jackknife blocks f4 test for z-value = 0, p (assumes normality): " + "{0:.4f}".format(observed_jackknife_f4_p_value) + "\n")
    if len(f4_per_block) > 2:
        sw_normality_test_results = scipy.stats.shapiro(f4_per_block)
        outfile.write("  Shapiro-Wilk test for normality of jackknife block f4 values, W: " + "{0:.4f}".format(sw_normality_test_results[0]) + "\n")
        outfile.write("  Shapiro-Wilk test for normality of jackknife block f4 values, p: " + "{0:.4f}".format(sw_normality_test_results[1]) + "\n")
    else:
        outfile.write("  Shapiro-Wilk test for normality of Jackknife block f4 values, W: NA\n")
        outfile.write("  Shapiro-Wilk test for normality of Jackknife block f4 values, p: NA\n")
    outfile.write("\n")

# Run simulations to assess support for f4 < 0.
if number_of_simulations != -1:
    if number_of_simulations < 1:
        print("ERROR: The specified number of simulations (option '-s') is smaller than 1!")
        sys.exit(1)
    converged = False
    simulated_f4s = []
    if jackknife_k != -1:
        simulated_jackknife_f4_z_zcores = []
    simulation_proportion_of_snps_variable_in_more_than_one_population = []
    simulation_proportion_of_snps_variable_on_both_sides_of_the_root = []
    effective_population_size = 100
    time_of_second_divergence_x = 100
    effective_population_size_has_been_too_large = False
    effective_population_size_has_been_too_small = False
    time_of_second_divergence_has_been_too_large = False
    time_of_second_divergence_has_been_too_small = False
    while len(simulated_f4s) < number_of_simulations:
        if len(simulated_f4s) == 0:
            outfile.write("\rRunning simulations (in burnin) ...")
        else:
            outfile.write("\rRunning simulations (" + str(len(simulated_f4s)) + "/" + str(number_of_simulations) + ") ...                   ")
        time_of_second_divergence = 1000*((time_of_second_divergence_x/100)/(1+(time_of_second_divergence_x/100)))
        # Write a temporary file for fastsimcoal2.
        fsc_input_string = ""
        fsc_input_string += "//Number of population samples (demes)\n"
        fsc_input_string += "4\n"
        fsc_input_string += "//Population effective sizes (number of genes)\n"
        fsc_input_string += str(effective_population_size) + "\n"
        fsc_input_string += str(effective_population_size) + "\n"
        fsc_input_string += str(effective_population_size) + "\n"
        fsc_input_string += str(effective_population_size) + "\n"
        fsc_input_string += "//Sample sizes\n"
        fsc_input_string += str(max_number_of_alleles_a) + "\n"
        fsc_input_string += str(max_number_of_alleles_b) + "\n"
        fsc_input_string += str(max_number_of_alleles_c) + "\n"
        fsc_input_string += str(max_number_of_alleles_d) + "\n"
        fsc_input_string += "//Growth rates  : negative growth implies population expansion\n"
        fsc_input_string += "0\n"
        fsc_input_string += "0\n"
        fsc_input_string += "0\n"
        fsc_input_string += "0\n"
        fsc_input_string += "//Number of migration matrices : 0 implies no migration between demes\n"
        fsc_input_string += "0 migration matrices\n"
        fsc_input_string += "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix\n"
        fsc_input_string += "3 historical events\n"
        fsc_input_string += str(int(time_of_second_divergence)) + " 0 1 1 1 0 0\n"
        fsc_input_string += str(int(time_of_second_divergence)) + " 2 3 1 1 0 0\n"
        fsc_input_string += "1000 1 3 1 1 0 0\n"
        # fsc_input_string += "2000 3 3 1 " + "{0:.4f}".format(4/effective_population_size) + " 0 0\n"
        fsc_input_string += "//Number of independent loci [chromosome]\n"
        fsc_input_string += str(int(number_of_valid_snps/proportion_missing_max)) + " 0\n"
        fsc_input_string += "//Per chromosome: Number of linkage blocks\n"
        fsc_input_string += "1\n"
        fsc_input_string += "//per Block: data type, num loci, rec. rate and mut rate + optional parameters\n"
        fsc_input_string += "SNP 1 0 0\n"
        tmp_in = tempfile.NamedTemporaryFile(delete=False)
        tmp_in.write(fsc_input_string.encode('utf-8'))
        tmp_in.close()
        FNULL = open(os.devnull, 'w')
        call_list = ["fsc252", "-i", tmp_in.name, "-n", "1"]
        call(call_list, stdout=FNULL, stderr=FNULL)

        # Read the fastsimcoal2 output file.
        if "/" in tmp_in.name:
            tmp_file_name = tmp_in.name.split("/")[-1]
        else:
            tmp_file_name = tmp_in.name
        if os.path.isfile(tmp_file_name + "/" + tmp_file_name + "_1_1.arp"):
            tmp_out = open(tmp_file_name + "/" + tmp_file_name + "_1_1.arp")
            fsc_output_string = tmp_out.read()
            tmp_out.close()
            os.remove(tmp_file_name + "/" + tmp_file_name + "_1_1.arp")
            if os.path.isfile(tmp_file_name + "/" + tmp_file_name + "_1.arb"):
                os.remove(tmp_file_name + "/" + tmp_file_name + "_1.arb")
            if os.path.isfile(tmp_file_name + "/" + tmp_file_name + "_1.simparam"):
                os.remove(tmp_file_name + "/" + tmp_file_name + "_1.simparam")
            os.rmdir(tmp_file_name)
            if os.path.isfile("seed.txt"):
                os.remove("seed.txt")
        else:
            print("ERROR: Fastsimcoal2 output file " + tmp_file_name + "/" + tmp_file_name + "_1_1.arp" + " can not be found!")
            sys.exit(1)

        # Parse the fastsimcoal2 output file and create treemix format lines (without actually writing them to file).
        body_lines = []
        fsc_output_lines = fsc_output_string.split("\n")
        in_record = False
        in_samples = False
        simulated_ids = []
        simulated_snps = []
        for line in fsc_output_lines:
            if line.strip("") != "":
                if line[0] != "#":
                    if "{" in line:
                        in_record = True
                    elif "}" in line:
                        in_record = False
                    elif "[[Samples]]" in line:
                        in_samples = True
                    elif "[[Structure]]" in line:
                        in_samples = False
                    elif in_record and in_samples:
                        line_ary = line.split()
                        simulated_ids.append(line_ary[0])
                        simulated_snps.append(line_ary[2])
        continue_searching = True
        enough_snps_found = False
        pos = 0
        while continue_searching:
            if len(body_lines) == number_of_valid_snps:
                enough_snps_found = True
                continue_searching = False
            elif pos == len(simulated_snps[0]):
                continue_searching = False
            else:
                pop_a_alleles = []
                pop_b_alleles = []
                pop_c_alleles = []
                pop_d_alleles = []
                for x in range(len(simulated_ids)):
                    if "1_" in simulated_ids[x]:
                        pop_a_alleles.append(simulated_snps[x][pos])
                    elif "2_" in simulated_ids[x]:
                        pop_b_alleles.append(simulated_snps[x][pos])
                    elif "3_" in simulated_ids[x]:
                        pop_c_alleles.append(simulated_snps[x][pos])
                    elif "4_" in simulated_ids[x]:
                        pop_d_alleles.append(simulated_snps[x][pos])
                    else:
                        print("ERROR: Id of simulated individual " + simulated_ids[x] + " could not be assigned to a population!")
                        sys.exit(1)
                # Mask alleles according to the missing data in the observed SNPs.
                number_of_alleles_to_mask_a = len(pop_a_alleles) - number_of_alleles_per_snp_a[len(body_lines)]
                number_of_alleles_to_mask_b = len(pop_b_alleles) - number_of_alleles_per_snp_b[len(body_lines)]
                number_of_alleles_to_mask_c = len(pop_c_alleles) - number_of_alleles_per_snp_c[len(body_lines)]
                number_of_alleles_to_mask_d = len(pop_d_alleles) - number_of_alleles_per_snp_d[len(body_lines)]
                if min(number_of_alleles_to_mask_a,number_of_alleles_to_mask_b,number_of_alleles_to_mask_c,number_of_alleles_to_mask_d) < 0:
                    print("ERROR: Something went wrong with masking of simulated data!")
                    sys.exit(1)
                indices_of_alleles_to_mask_a = random.sample(range(len(pop_a_alleles)),number_of_alleles_to_mask_a)
                for x in indices_of_alleles_to_mask_a:
                    pop_a_alleles[x] = "N"
                indices_of_alleles_to_mask_b = random.sample(range(len(pop_b_alleles)),number_of_alleles_to_mask_b)
                for x in indices_of_alleles_to_mask_b:
                    pop_b_alleles[x] = "N"
                indices_of_alleles_to_mask_c = random.sample(range(len(pop_c_alleles)),number_of_alleles_to_mask_c)
                for x in indices_of_alleles_to_mask_c:
                    pop_c_alleles[x] = "N"
                indices_of_alleles_to_mask_d = random.sample(range(len(pop_d_alleles)),number_of_alleles_to_mask_d)
                for x in indices_of_alleles_to_mask_d:
                    pop_d_alleles[x] = "N"

                # Make sure that this is still a bi-allelic SNP after masking.
                allele_1_present = False
                allele_2_present = False
                if "0" in pop_a_alleles:
                    allele_1_present = True
                if "1" in pop_a_alleles:
                    allele_2_present = True
                if "0" in pop_b_alleles:
                    allele_1_present = True
                if "1" in pop_b_alleles:
                    allele_2_present = True
                if "0" in pop_c_alleles:
                    allele_1_present = True
                if "1" in pop_c_alleles:
                    allele_2_present = True
                if "0" in pop_d_alleles:
                    allele_1_present = True
                if "1" in pop_d_alleles:
                    allele_2_present = True
                
                # If this is the case, count allele frequencies.
                if allele_1_present and allele_2_present:
                    body_line = ""
                    body_line += str(pop_a_alleles.count("0")) + "," + str(pop_a_alleles.count("1")) + " "
                    body_line += str(pop_b_alleles.count("0")) + "," + str(pop_b_alleles.count("1")) + " "
                    body_line += str(pop_c_alleles.count("0")) + "," + str(pop_c_alleles.count("1")) + " "
                    body_line += str(pop_d_alleles.count("0")) + "," + str(pop_d_alleles.count("1"))
                    body_lines.append(body_line)

                pos += 1

        # If enough bi-allelic SNPs were present after masking calculate the f4 statistic (and write data sets to files).
        if enough_snps_found:

            # If an output directory has been specified, write this simulated data set in treemix format.
            if converged:
                if output_dir != "-1":
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                    treemix_string = header_line + "\n"
                    for line in body_lines:
                        treemix_string += line + "\n"
                    simulation_id_number_of_digits = len(str(number_of_simulations))
                    simulation_id = str(len(simulated_f4s)+1).rjust(simulation_id_number_of_digits).replace(" ","0")
                    treemix_file_name = "s" + simulation_id + ".txt"
                    treemix_file_name_w_path = output_dir + "/" + treemix_file_name
                    treemix_file = open(treemix_file_name_w_path, 'w')
                    treemix_file.write(treemix_string)

            # Calculate the f4 statistic.
            f4s = []
            number_of_snps_variable_in_more_than_one_population_this_simulation = 0
            number_of_snps_variable_on_both_sides_of_the_root_this_simulation = 0
            for x in range(len(body_lines)):
                line_ary = body_lines[x].strip().split()
                a_ary = line_ary[0].split(",")
                b_ary = line_ary[1].split(",")
                c_ary = line_ary[2].split(",")
                d_ary = line_ary[3].split(",")
                number_of_alleles_a = int(a_ary[0])+int(a_ary[1])
                number_of_alleles_b = int(b_ary[0])+int(b_ary[1])
                number_of_alleles_c = int(c_ary[0])+int(c_ary[1])
                number_of_alleles_d = int(d_ary[0])+int(d_ary[1])
                number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation = 0
                variable_in_a_or_b = False
                a = int(a_ary[0])/number_of_alleles_a
                if a > 0 and a < 1:
                    number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation += 1
                    variable_in_a_or_b = True
                b = int(b_ary[0])/number_of_alleles_b
                if b > 0 and b < 1:
                    number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation += 1
                    variable_in_a_or_b = True
                variable_in_c_or_d = False
                c = int(c_ary[0])/number_of_alleles_c
                if c > 0 and c < 1:
                    number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation += 1
                    variable_in_c_or_d = True
                d = int(d_ary[0])/number_of_alleles_d
                if d > 0 and d < 1:
                    number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation += 1
                    variable_in_c_or_d = True
                f4s.append((a-b)*(c-d))
                if number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation > 1:
                    number_of_snps_variable_in_more_than_one_population_this_simulation += 1
                if variable_in_a_or_b and variable_in_c_or_d:
                    number_of_snps_variable_on_both_sides_of_the_root_this_simulation += 1
            simulated_f4 = numpy.mean(f4s)

            # Jackknife the simulated data set.
            if converged:
                if jackknife_k != -1:
                    body_lines_per_block = [[]]
                    for x in range(len(body_lines)):
                        if len(body_lines_per_block[-1]) < jackknife_k:
                            body_lines_per_block[-1].append(body_lines[x])
                        else:
                            body_lines_per_block.append([body_lines[x]])
                    body_lines_per_valid_block = []
                    for item in body_lines_per_block:
                        if len(item) == jackknife_k:
                            body_lines_per_valid_block.append(item)
                    f4_per_block = []
                    for block_body_lines in body_lines_per_valid_block:
                        block_f4s = []
                        for x in range(len(block_body_lines)):
                            if block_body_lines[x].strip() != "":
                                line_ary = block_body_lines[x].strip().split()
                                a_ary = line_ary[0].split(",")
                                b_ary = line_ary[1].split(",")
                                c_ary = line_ary[2].split(",")
                                d_ary = line_ary[3].split(",")
                                number_of_alleles_a = int(a_ary[0])+int(a_ary[1])
                                number_of_alleles_b = int(b_ary[0])+int(b_ary[1])
                                number_of_alleles_c = int(c_ary[0])+int(c_ary[1])
                                number_of_alleles_d = int(d_ary[0])+int(d_ary[1])
                                if min(number_of_alleles_a,number_of_alleles_b,number_of_alleles_c,number_of_alleles_d) > 0:
                                    a = int(a_ary[0])/number_of_alleles_a
                                    b = int(b_ary[0])/number_of_alleles_b
                                    c = int(c_ary[0])/number_of_alleles_c
                                    d = int(d_ary[0])/number_of_alleles_d
                                    block_f4s.append((a-b)*(c-d))
                                else:
                                    print("ERROR: One of the four populations has no allele in this line: ")
                                    print(block_body_lines[x].strip())
                                    sys.exit(1)
                        f4_per_block.append(numpy.mean(block_f4s))
                    simulated_jackknife_f4_standard_error = stats.sem(f4_per_block)
                    simulated_jackknife_f4_z_zcore = simulated_f4/simulated_jackknife_f4_standard_error

            # Check the fit of the parameters of effective population size and time of second divergence.
            # Let random chance decide which of the two variable is adjusted:
            effective_population_size_scaler = 1 + random.random()/10
            if random.randint(0,1) == 0:
                if number_of_snps_variable_in_more_than_one_population_this_simulation > number_of_snps_variable_in_more_than_one_population:
                    effective_population_size_has_been_too_large = True
                    effective_population_size = max(4, int(effective_population_size/effective_population_size_scaler))
                elif number_of_snps_variable_in_more_than_one_population_this_simulation < number_of_snps_variable_in_more_than_one_population:
                    effective_population_size_has_been_too_small = True
                    effective_population_size = int(effective_population_size*effective_population_size_scaler)
            else:
                if number_of_snps_variable_on_both_sides_of_the_root_this_simulation > number_of_snps_variable_on_both_sides_of_the_root:
                    time_of_second_divergence_has_been_too_large = True
                    effective_population_size = max(4, int(effective_population_size/effective_population_size_scaler))
                    time_of_second_divergence_x -= random.randint(0,2)
                elif number_of_snps_variable_on_both_sides_of_the_root_this_simulation < number_of_snps_variable_on_both_sides_of_the_root:
                    time_of_second_divergence_has_been_too_small = True
                    effective_population_size = int(effective_population_size*effective_population_size_scaler)
                    time_of_second_divergence_x += random.randint(0,2)

            if converged:
                simulated_f4s.append(simulated_f4)
                if jackknife_k != -1:
                    simulated_jackknife_f4_z_zcores.append(simulated_jackknife_f4_z_zcore)
                simulation_proportion_of_snps_variable_in_more_than_one_population.append(number_of_snps_variable_in_more_than_one_population_this_simulation/number_of_valid_snps)
                simulation_proportion_of_snps_variable_on_both_sides_of_the_root.append(number_of_snps_variable_on_both_sides_of_the_root_this_simulation/number_of_valid_snps)

            # Save the results for this simulation: the f4 statistic, as well as the
            # proportion of SNPs variable in more than one population and
            # the proportion of SNPs variable on both sides of the root.
            if effective_population_size_has_been_too_large and effective_population_size_has_been_too_small:
                if time_of_second_divergence_has_been_too_large and time_of_second_divergence_has_been_too_small:
                    converged = True

    outfile.write("\rSimulations                                    \n")
    outfile.write("  Number of simulations: " + str(number_of_simulations) + "\n")
    outfile.write("  Proportion of SNPs variable in more than one population: ")
    outfile.write("{0:.2f}".format(numpy.mean(simulation_proportion_of_snps_variable_in_more_than_one_population)))
    outfile.write(" (+/-")
    outfile.write("{0:.2f}".format(numpy.std(simulation_proportion_of_snps_variable_in_more_than_one_population)))
    outfile.write(")\n")
    outfile.write("  Proportion of SNPs variable on both sides of the root: ")
    outfile.write("{0:.2f}".format(numpy.mean(simulation_proportion_of_snps_variable_on_both_sides_of_the_root)))
    outfile.write(" (+/-")
    outfile.write("{0:.2f}".format(numpy.std(simulation_proportion_of_snps_variable_on_both_sides_of_the_root)))
    outfile.write(")\n")
    outfile.write("  Simulated f4: ")
    outfile.write("{0:.5f}".format(numpy.mean(simulated_f4s)))
    outfile.write(" (+/-")
    outfile.write("{0:.5f}".format(numpy.std(simulated_f4s)))
    outfile.write(")\n")
    number_of_simulations_with_f4_more_extreme_than_observed = 0
    for simulated_f4 in simulated_f4s:
        if observed_f4 < 0 and simulated_f4 < observed_f4:
            number_of_simulations_with_f4_more_extreme_than_observed += 1
        elif observed_f4 > 0 and simulated_f4 > observed_f4:
            number_of_simulations_with_f4_more_extreme_than_observed += 1
        elif observed_f4 == 0 and simulated_f4 != observed_f4:
            number_of_simulations_with_f4_more_extreme_than_observed += 1
    if observed_f4 < 0:
        outfile.write("  Proportion of simulated f4 values smaller than the observed: ")
    elif observed_f4 > 0:
        outfile.write("  Proportion of simulated f4 values larger than the observed: ")
    elif observed_f4 == 0:
        outfile.write("  Proportion of simulated f4 values different from the observed: ")
    outfile.write("{0:.4f}".format(number_of_simulations_with_f4_more_extreme_than_observed/number_of_simulations))
    outfile.write("\n")
    if jackknife_k != -1:
        number_of_simulations_with_jackknife_f4_z_score_more_extreme_than_observed = 0
        for simulated_jackknife_f4_z_zcore in simulated_jackknife_f4_z_zcores:
            if observed_jackknife_f4_z_zcore < 0 and simulated_jackknife_f4_z_zcore < observed_jackknife_f4_z_zcore:
                number_of_simulations_with_jackknife_f4_z_score_more_extreme_than_observed += 1
            elif observed_jackknife_f4_z_zcore > 0 and simulated_jackknife_f4_z_zcore > observed_jackknife_f4_z_zcore:
                number_of_simulations_with_jackknife_f4_z_score_more_extreme_than_observed += 1
            elif observed_jackknife_f4_z_zcore == 0 and simulated_jackknife_f4_z_zcore != observed_jackknife_f4_z_zcore:
                number_of_simulations_with_jackknife_f4_z_score_more_extreme_than_observed += 1
        if observed_jackknife_f4_z_zcore < 0:
            outfile.write("  Proportion of simulated jackknife blocks f4 z-values smaller than the observed: ")
        elif observed_jackknife_f4_z_zcore > 0:
            outfile.write("  Proportion of simulated jackknife blocks f4 z-values larger than the observed: ")
        elif observed_jackknife_f4_z_zcore == 0:
            outfile.write("  Proportion of simulated jackknife blocks f4 z-values different from the observed: ")
        outfile.write("{0:.4f}".format(number_of_simulations_with_jackknife_f4_z_score_more_extreme_than_observed/number_of_simulations))
        outfile.write("\n")
    outfile.write("\n")
