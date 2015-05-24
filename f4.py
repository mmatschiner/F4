#!/usr/local/bin/python3

# Michael Matschiner, 2015-05-21
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
    version='%(prog)s 0.91'
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
    '-l',
    nargs=1,
    type=str,
    default=["-1"],
    dest='log_file_name',
    help="Name of a file to which simulation parameters should be written, to assess convergence (mostly for debugging) (default: off)."
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
log_file_name = args.log_file_name[0]
infile = args.infile
outfile = args.outfile
if infile.isatty():
    print('No input file specified, and no input piped through stdin!')
    sys.exit(0)
if output_dir != "-1" and number_of_simulations < 1:
    print("")
    print("WARNING: An output directory for simulated data sets has been specified, but no simulations are performed!")

# If a log file should be written, prepare it.
if log_file_name != "-1":
    if number_of_simulations < 1:
        print("")
        print("WARNING: An output directory for simulated data sets has been specified, but no simulations are performed!")
    else:
        log_file = open(log_file_name, "w")
        log_header_line = "State"
        log_header_line += "\t"
        log_header_line += "Effective_population_size"
        log_header_line += "\t"
        log_header_line += "Time_of_second_divergence"
        log_header_line += "\t"
        log_header_line += "Proportion_of_SNPs_variable_in_more_than_one_population"
        log_header_line += "\t"
        log_header_line += "Proportion_of_SNPs_variable_on_both_sides_of_the_root"
        log_header_line += "\t"
        log_header_line += "f4"
        log_header_line += "\t"
        log_header_line += "f4_standard_error"
        log_header_line += "\n"
        log_file.write(log_header_line)
        log_file.flush()
        number_of_log_body_lines_written = 0

# Prepare the output string.
outfile.write("\n")
outfile.write("  f4.py\n")
outfile.write("----------------------------------------------------------------------\n")
outfile.write("  Michael Matschiner | 2015-05-21 | evoinformatics.eu\n")
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
original_body_lines = inlines[1:]

# Remove empty lines and SNPs with completely missing data in one population, but count them before.
valid_body_lines = []
number_of_snps = 0
linkage_break_points = []
for body_line in original_body_lines:
    body_line = body_line.strip()
    if body_line == "":
        linkage_break_points.append(len(valid_body_lines))
    else:
        number_of_snps += 1
        if " 0,0" not in body_line and body_line[:3] != "0,0":
            valid_body_lines.append(body_line)
linkage_break_points = sorted(list(set(linkage_break_points)))
if len(valid_body_lines) in linkage_break_points:
    linkage_break_points.remove(len(valid_body_lines))
linkage_block_sizes = []
for linkage_break_point in linkage_break_points:
    linkage_block_sizes.append(linkage_break_point-sum(linkage_block_sizes))
linkage_block_sizes.append(len(valid_body_lines)-sum(linkage_block_sizes))
body_lines = valid_body_lines
pops = header_line.strip().split()
if len(pops) != 4:
    print("ERROR: The input file should contain data for exactly four populations. The header line has only " + str(len(pops)) + " columns!")
    sys.exit(1)
outfile.write("  Assumed population topology: (" +  pops[0] + "," + pops[1] + "),(" + pops[2] + "," + pops[3] + ")\n")

observed_f4s = []
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
        a = int(a_ary[0])/number_of_alleles_a
        if a > 0 and a < 1:
            number_of_populations_with_more_than_one_allele_for_this_snp += 1
        b = int(b_ary[0])/number_of_alleles_b
        if b > 0 and b < 1:
            number_of_populations_with_more_than_one_allele_for_this_snp += 1

        variable_among_a_and_b = False
        if round(a,2) != round(b,2):
            variable_among_a_and_b = True

        c = int(c_ary[0])/number_of_alleles_c
        if c > 0 and c < 1:
            number_of_populations_with_more_than_one_allele_for_this_snp += 1
        d = int(d_ary[0])/number_of_alleles_d
        if d > 0 and d < 1:
            number_of_populations_with_more_than_one_allele_for_this_snp += 1

        variable_among_c_and_d = False
        if round(c,2) != round(d,2):
            variable_among_c_and_d = True

        observed_f4s.append((a-b)*(c-d))
        number_of_valid_snps += 1
        if number_of_populations_with_more_than_one_allele_for_this_snp > 1:
            number_of_snps_variable_in_more_than_one_population += 1
        if variable_among_a_and_b and variable_among_c_and_d:
            number_of_snps_variable_on_both_sides_of_the_root += 1
    else:
        print("ERROR: One of the four populations has no allele in this line: ")
        print(body_lines[x].strip())
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
if len(linkage_block_sizes) > 1:
    outfile.write("  Number of linkage blocks detected (separated by empty lines): " + str(len(linkage_block_sizes)) + "\n")
    if len(linkage_block_sizes) < 7:
        outfile.write("  Numbers of SNPs used per linkage block: ")
        for x in range(len(linkage_block_sizes)-1):
            outfile.write(str(linkage_block_sizes[x]) + ", ")
        outfile.write(str(linkage_block_sizes[-1]) + "\n")
    else:
        outfile.write("  Numbers of SNPs used per linkage block: between ")
        outfile.write(str(min(linkage_block_sizes)))
        outfile.write(" and ")
        outfile.write(str(max(linkage_block_sizes)))
        outfile.write(", with median ")
        outfile.write(str(numpy.median(linkage_block_sizes)))
        outfile.write("\n")
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
outfile.write("  Proportions of missing data per population: ")
outfile.write("{0:.2f}".format(proportion_missing_a) + ", ")
outfile.write("{0:.2f}".format(proportion_missing_b) + ", ")
outfile.write("{0:.2f}".format(proportion_missing_c) + ", ")
outfile.write("{0:.2f}".format(proportion_missing_d) + "\n")
proportion_of_snps_variable_in_more_than_one_population = number_of_snps_variable_in_more_than_one_population/number_of_valid_snps
outfile.write("  Proportion of SNPs variable in more than one population: " + "{0:.3f}".format(proportion_of_snps_variable_in_more_than_one_population) + "\n")
proportion_of_snps_variable_on_both_sides_of_the_root = number_of_snps_variable_on_both_sides_of_the_root/number_of_valid_snps
outfile.write("  Proportion of SNPs variable on both sides of the root: " + "{0:.3f}".format(proportion_of_snps_variable_on_both_sides_of_the_root) + "\n")
observed_f4 = numpy.mean(observed_f4s)
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
        print("ERROR: The specified jack-knifing block size (k) allows only a single block of SNPs!")
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
    FNULL = open(os.devnull, 'w')
    converged = False
    simulated_f4s = []
    simulated_f4s_including_burnin = []
    if jackknife_k != -1:
        simulated_jackknife_f4_z_zcores = []
    simulation_proportion_of_snps_variable_in_more_than_one_population = []
    simulation_proportion_of_snps_variable_on_both_sides_of_the_root = []
    simulation_effective_population_sizes = []
    simulation_times_of_second_divergence = []
    effective_population_size = 100
    time_of_second_divergence = 500
    effective_population_size_has_been_too_large = False
    effective_population_size_has_been_too_small = False
    time_of_second_divergence_has_been_too_large = False
    time_of_second_divergence_has_been_too_small = False
    number_of_burnin_simulations = 0
    while len(simulated_f4s) < number_of_simulations:
        if len(simulated_f4s) == 0:
            outfile.write("\rRunning simulations (burnin " + str(number_of_burnin_simulations) + ")...")
            number_of_burnin_simulations += 1
        elif len(simulated_f4s) == 1:
            outfile.write("\r                                                                                ")
            outfile.write("\rRunning simulations (1/" + str(number_of_simulations) + ")...")
        else:
            outfile.write("\rRunning simulations (" + str(len(simulated_f4s)) + "/" + str(number_of_simulations) + ")...")
        # time_of_second_divergence = 1000*((time_of_second_divergence_x/100)/(1+(time_of_second_divergence_x/100)))

        # If just a single block of input lines was found, assume that all SNPs are unlinked. If however,
        # multiple blocks separated by empty lines were found, assume these to represent linkage blocks.
        if len(linkage_block_sizes) == 1:
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
            fsc_input_string += "//Number of independent loci [chromosome]\n"
            fsc_input_string += str(int(number_of_valid_snps/(1-proportion_missing_max))+1) + " 0\n"
            fsc_input_string += "//Per chromosome: Number of linkage blocks\n"
            fsc_input_string += "1\n"
            fsc_input_string += "//per Block: data type, num loci, rec. rate and mut rate + optional parameters\n"
            fsc_input_string += "SNP 1 0 0\n"
            tmp_in = tempfile.NamedTemporaryFile(delete=False)
            tmp_in.write(fsc_input_string.encode('utf-8'))
            tmp_in.close()
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
            not_enough_snps_found = False
            pos = 0
            while continue_searching:
                if len(body_lines) == number_of_valid_snps:
                    continue_searching = False
                elif pos == len(simulated_snps[0]):
                    continue_searching = False
                    not_enough_snps_found = True
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

        else:

            body_lines = []
            not_enough_snps_found = False
            for linkage_block_size in linkage_block_sizes:
                # Write a temporary file for fastsimcoal2, simulate a block of linked SNPs.
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
                fsc_input_string += "//Number of independent loci [chromosome]\n"
                fsc_input_string += "1 0\n"
                fsc_input_string += "//Per chromosome: Number of linkage blocks\n"
                fsc_input_string += "1\n"
                fsc_input_string += "//per Block: data type, num loci, rec. rate and mut rate + optional parameters\n"
                fsc_input_string += "SNP " + str(int(linkage_block_size/(1-proportion_missing_max))+1) + " 0 0\n"
                tmp_in = tempfile.NamedTemporaryFile(delete=False)
                tmp_in.write(fsc_input_string.encode('utf-8'))
                tmp_in.close()
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
                number_of_body_lines_this_linkage_block = 0
                continue_searching = True
                pos = 0
                while continue_searching:
                    if number_of_body_lines_this_linkage_block == linkage_block_size:
                        continue_searching = False
                    elif pos == len(simulated_snps[0]):
                        continue_searching = False
                        not_enough_snps_found = True
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
                            number_of_body_lines_this_linkage_block += 1

                        pos += 1

        # If enough bi-allelic SNPs were present after masking calculate the f4 statistic (and write data sets to files).
        if not_enough_snps_found == False:

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

                a = int(a_ary[0])/number_of_alleles_a
                if a > 0 and a < 1:
                    number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation += 1
                b = int(b_ary[0])/number_of_alleles_b
                if b > 0 and b < 1:
                    number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation += 1

                variable_among_a_and_b = False
                if round(a,2) != round(b,2):
                    variable_among_a_and_b = True

                c = int(c_ary[0])/number_of_alleles_c
                if c > 0 and c < 1:
                    number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation += 1
                d = int(d_ary[0])/number_of_alleles_d
                if d > 0 and d < 1:
                    number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation += 1

                variable_among_c_and_d = False
                if round(c,2) != round(d,2):
                    variable_among_c_and_d = True

                f4s.append((a-b)*(c-d))
                if number_of_populations_with_more_than_one_allele_for_this_snp_this_simulation > 1:
                    number_of_snps_variable_in_more_than_one_population_this_simulation += 1
                if variable_among_a_and_b and variable_among_c_and_d:
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
            # Decide which of the two variable is adjusted, based on how far each of the two parameters is from its optimum.
            effective_population_size_scaler = 1 + random.random()/20
            absolute_diff1 = numpy.absolute(number_of_snps_variable_in_more_than_one_population_this_simulation-number_of_snps_variable_in_more_than_one_population)
            proportional_diff1 = absolute_diff1/number_of_snps_variable_in_more_than_one_population
            absolute_diff2 = numpy.absolute(number_of_snps_variable_on_both_sides_of_the_root_this_simulation-number_of_snps_variable_on_both_sides_of_the_root)
            proportional_diff2 = absolute_diff2/number_of_snps_variable_on_both_sides_of_the_root
            change_effective_population_size = False
            change_second_divergence_time = False
            if proportional_diff1 > proportional_diff2:
                if random.random() < proportional_diff1/(proportional_diff1+proportional_diff2):
                    change_effective_population_size = True
                else:
                    change_second_divergence_time = True
            elif proportional_diff1 < proportional_diff2:
                if random.random() < proportional_diff2/(proportional_diff1+proportional_diff2):
                    change_second_divergence_time = True
                else:
                    change_effective_population_size = True
            else:
                if random.randint(0,1) == 0:
                    change_effective_population_size = True
                else:
                    change_second_divergence_time = True
            print("change_effective_population_size: " + str(change_effective_population_size) + " change_second_divergence_time:" + str(change_second_divergence_time))
            if change_effective_population_size:
                if number_of_snps_variable_in_more_than_one_population_this_simulation > number_of_snps_variable_in_more_than_one_population:
                    effective_population_size_has_been_too_large = True
                    effective_population_size = max(4, int(effective_population_size/effective_population_size_scaler))
                elif number_of_snps_variable_in_more_than_one_population_this_simulation < number_of_snps_variable_in_more_than_one_population:
                    effective_population_size_has_been_too_small = True
                    effective_population_size = int(effective_population_size*effective_population_size_scaler)
            elif change_second_divergence_time:
                if number_of_snps_variable_on_both_sides_of_the_root_this_simulation > number_of_snps_variable_on_both_sides_of_the_root:
                    time_of_second_divergence_has_been_too_large = True
                    # effective_population_size = max(4, int(effective_population_size/effective_population_size_scaler))
                    # time_of_second_divergence_x -= random.randint(0,2)
                    time_of_second_divergence -= 10*random.random()
                    if time_of_second_divergence < 1:
                        time_of_second_divergence = 1
                elif number_of_snps_variable_on_both_sides_of_the_root_this_simulation < number_of_snps_variable_on_both_sides_of_the_root:
                    time_of_second_divergence_has_been_too_small = True
                    # effective_population_size = int(effective_population_size*effective_population_size_scaler)
                    # time_of_second_divergence_x += random.randint(0,2)
                    time_of_second_divergence += 10*random.random()
                    if time_of_second_divergence > 999:
                        time_of_second_divergence = 999
            else:
                print("ERROR: None of the two parameters effective population size and second divergence time could be optimized in the last step.")
                sys.exit(1)

            if converged:
                simulated_f4s.append(simulated_f4)
                if jackknife_k != -1:
                    simulated_jackknife_f4_z_zcores.append(simulated_jackknife_f4_z_zcore)
                simulation_proportion_of_snps_variable_in_more_than_one_population.append(number_of_snps_variable_in_more_than_one_population_this_simulation/number_of_valid_snps)
                simulation_proportion_of_snps_variable_on_both_sides_of_the_root.append(number_of_snps_variable_on_both_sides_of_the_root_this_simulation/number_of_valid_snps)
                simulation_effective_population_sizes.append(effective_population_size)
                # simulation_times_of_second_divergence.append(1000*((time_of_second_divergence_x/100)/(1+(time_of_second_divergence_x/100))))
                simulation_times_of_second_divergence.append(time_of_second_divergence)
            simulated_f4s_including_burnin.append(simulated_f4)

            if log_file_name != "-1":
                log_body_line = str(number_of_log_body_lines_written)
                log_body_line += "\t"
                log_body_line += str(effective_population_size)
                log_body_line += "\t"
                # log_body_line += str(1000*((time_of_second_divergence_x/100)/(1+(time_of_second_divergence_x/100))))
                log_body_line += str(time_of_second_divergence)
                log_body_line += "\t"
                log_body_line += str(number_of_snps_variable_in_more_than_one_population_this_simulation/number_of_valid_snps)
                log_body_line += "\t"
                log_body_line += str(number_of_snps_variable_on_both_sides_of_the_root_this_simulation/number_of_valid_snps)
                log_body_line += "\t"
                log_body_line += str(simulated_f4)
                if len(simulated_f4s_including_burnin) == 1:
                    simulated_f4_standard_error_last_10 = 0
                elif len(simulated_f4s_including_burnin) < 10:
                    simulated_f4_standard_error_last_10 = stats.sem(simulated_f4s_including_burnin)
                else:
                    simulated_f4_standard_error_last_10 = stats.sem(simulated_f4s_including_burnin[-10:])
                log_body_line += "\t"
                log_body_line += str(simulated_f4_standard_error_last_10)
                log_body_line += "\n"
                number_of_log_body_lines_written += 1
                log_file.write(log_body_line)
                log_file.flush()

            # Save the results for this simulation: the f4 statistic, as well as the
            # proportion of SNPs variable in more than one population and
            # the proportion of SNPs variable on both sides of the root.
            if effective_population_size_has_been_too_large and effective_population_size_has_been_too_small:
                if time_of_second_divergence_has_been_too_large and time_of_second_divergence_has_been_too_small:
                    if proportional_diff1 < 0.05 and proportional_diff2 < 0.05:
                        if converged == False:
                            converged = True

    # Check whether there may have been convergence issues. If both the proportion of SNPs variable in more than one
    # population and the proportion of SNPs variable on both sides of the root lie within the standard deviation of
    # the values resulting from simulations, there seem to be no issues.
    simulation_proportion_of_snps_variable_in_more_than_one_population_mean = numpy.mean(simulation_proportion_of_snps_variable_in_more_than_one_population)
    simulation_proportion_of_snps_variable_in_more_than_one_population_std = numpy.std(simulation_proportion_of_snps_variable_in_more_than_one_population)
    lower1 = simulation_proportion_of_snps_variable_in_more_than_one_population_mean-simulation_proportion_of_snps_variable_in_more_than_one_population_std
    upper1 = simulation_proportion_of_snps_variable_in_more_than_one_population_mean+simulation_proportion_of_snps_variable_in_more_than_one_population_std
    simulation_proportion_of_snps_variable_on_both_sides_of_the_root_mean = numpy.mean(simulation_proportion_of_snps_variable_on_both_sides_of_the_root)
    simulation_proportion_of_snps_variable_on_both_sides_of_the_root_std = numpy.std(simulation_proportion_of_snps_variable_on_both_sides_of_the_root)
    lower2 = simulation_proportion_of_snps_variable_on_both_sides_of_the_root_mean-simulation_proportion_of_snps_variable_on_both_sides_of_the_root_std
    upper2 = simulation_proportion_of_snps_variable_on_both_sides_of_the_root_mean+simulation_proportion_of_snps_variable_on_both_sides_of_the_root_std
    convergence_issues = False
    if proportion_of_snps_variable_in_more_than_one_population < lower1:
        convergence_issues = True
    elif proportion_of_snps_variable_in_more_than_one_population > upper1:
        convergence_issues = True
    elif proportion_of_snps_variable_on_both_sides_of_the_root < lower2:
        convergence_issues = True
    elif proportion_of_snps_variable_on_both_sides_of_the_root > upper2:
        convergence_issues = True

    outfile.write("\rSimulations                                    \n")
    outfile.write("  Number of burnin simulations: " + str(number_of_burnin_simulations) + "\n")
    outfile.write("  Number of post-burnin simulations: " + str(number_of_simulations) + "\n")
    outfile.write("  Effective population size: " + "{0:.2f}".format(numpy.mean(simulation_effective_population_sizes)))
    outfile.write(" (+/-")
    outfile.write("{0:.2f}".format(numpy.std(simulation_effective_population_sizes)))
    outfile.write(")\n")
    outfile.write("  Divergence time ratio: " + "{0:.4f}".format(numpy.mean(simulation_times_of_second_divergence)/1000))
    outfile.write(" (+/-")
    outfile.write("{0:.4f}".format(numpy.std(simulation_times_of_second_divergence)/1000))
    outfile.write(")\n")
    outfile.write("  Proportion of SNPs variable in more than one population: ")
    outfile.write("{0:.3f}".format(simulation_proportion_of_snps_variable_in_more_than_one_population_mean))
    outfile.write(" (+/-")
    outfile.write("{0:.3f}".format(simulation_proportion_of_snps_variable_in_more_than_one_population_std))
    outfile.write(")\n")
    outfile.write("  Proportion of SNPs variable on both sides of the root: ")
    outfile.write("{0:.3f}".format(simulation_proportion_of_snps_variable_on_both_sides_of_the_root_mean))
    outfile.write(" (+/-")
    outfile.write("{0:.3f}".format(simulation_proportion_of_snps_variable_on_both_sides_of_the_root_std))
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
    proportion_of_simulations_with_f4_more_extreme_than_observed = number_of_simulations_with_f4_more_extreme_than_observed/number_of_simulations
    outfile.write("{0:.4f}".format(proportion_of_simulations_with_f4_more_extreme_than_observed))
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

outfile.write("Interpretation\n")
if observed_f4 == 0:
    outfile.write("  The observed f4 value is exactly zero, there is no sign of introgression.\n")
else:
    if number_of_simulations == -1:
        outfile.write("  The observed f4 value differs from zero, but whether the observed difference\n")
        outfile.write("  could result from incomplete lineage sorting alone should be tested with \n")
        outfile.write("  simulations (option -s).\n")
    else:
        if proportion_of_simulations_with_f4_more_extreme_than_observed < 0.01:
            outfile.write("  The observed f4 value differs from zero, and there is a very low\n")
            outfile.write("  probability that the observed difference can result from incomplete\n")
            outfile.write("  lineage sorting alone.\n")
        elif proportion_of_simulations_with_f4_more_extreme_than_observed < 0.05:
            outfile.write("  The observed f4 value differs from zero, and there is a low probability\n")
            outfile.write("  that the observed difference can result from incomplete lineage\n")
            outfile.write("  sorting alone.\n")
        else:
            outfile.write("  The observed f4 value differs from zero, but the observed difference\n")
            outfile.write("  could result from incomplete lineage sorting alone.\n")

        if convergence_issues:
            outfile.write("  Note that the simulated proportion of SNPs that are variable in one or\n")
            outfile.write("  more populations and/or the simulated proportion of SNPs that are variable\n")
            outfile.write("  on both sides of the root differ from the observed, which indicates\n")
            outfile.write("  convergence issues. This might be solved by increasing the number of\n")
            if log_file_name == "-1":
                outfile.write("  simulations. Also consider using option '-l' to write a log file that\n")
                outfile.write("  could be more informative regarding convergence.\n")
            else:
                outfile.write("  simulations. To assess convergence in more detail, use Tracer to open the\n")
                outfile.write("  log file '" + str(log_file_name) + "'.\n")

    # Check which and how many SNPs are driving this pattern.
    outlier_valid_body_lines = []
    outlier_f4s = []
    reduced_f4 = observed_f4
    reduced_f4s = observed_f4s
    if observed_f4 < 0:
        while reduced_f4 < 0:
            # Identify the most negative SNP.
            min_f4 = 0
            min_f4_index = 0
            for x in range(len(reduced_f4s)):
                if reduced_f4s[x] < min_f4:
                    min_f4 = reduced_f4s[x]
                    min_f4_index = x
            # Replace it with 0.
            reduced_f4s[min_f4_index] = 0
            reduced_f4 = numpy.mean(reduced_f4s)
            outlier_valid_body_lines.append(valid_body_lines[min_f4_index])
            outlier_f4s.append(min_f4)
    if observed_f4 > 0:
        while reduced_f4 > 0:
            # Identify the most positive SNP.
            max_f4 = 0
            max_f4_index = 0
            for x in range(len(reduced_f4s)):
                if reduced_f4s[x] > max_f4:
                    max_f4 = reduced_f4s[x]
                    max_f4_index = x
            # Replace it with 0.
            reduced_f4s[max_f4_index] = 0
            reduced_f4 = numpy.mean(reduced_f4s)
            outlier_valid_body_lines.append(valid_body_lines[max_f4_index])
            outlier_f4s.append(max_f4)

    printed_outlier_snps_body_line_indices = []
    index_string_length = max(len("Line"),len(str(len(original_body_lines))))
    outlier_lines = []
    for x in range(len(outlier_valid_body_lines)):
        outlier_snps_body_line_indices = []
        for y in range(len(original_body_lines)):
            if original_body_lines[y] == outlier_valid_body_lines[x]:
                outlier_snps_body_line_indices.append(y)
        # Select the first line index that has not been printed yet.
        selected_outlier_snps_body_line_index = 0
        for z in range(len(outlier_snps_body_line_indices)):
            if outlier_snps_body_line_indices[z] not in printed_outlier_snps_body_line_indices:
                selected_outlier_snps_body_line_index = outlier_snps_body_line_indices[z]
                printed_outlier_snps_body_line_indices.append(selected_outlier_snps_body_line_index)
                break
        selected_outlier_snps_body_line_index_string = str(selected_outlier_snps_body_line_index+2)
        outlier_lines.append("  " + selected_outlier_snps_body_line_index_string.rjust(index_string_length) + " | " + outlier_valid_body_lines[x] + " | " + "{0:.4f}".format(outlier_f4s[x]) + "\n")

    # Check wether the indices of the outlier SNPs are more clustered than expected by chance.
    linked = False
    very_linked = False

    if linkage_break_points == [] and len(outlier_valid_body_lines) > 1:
        observed_index_distances = []
        observed_cluster_measure = 0
        number_of_reshuffled_indices_that_are_more_clustered_than_the_observed = 0
        sorted_printed_outlier_snps_body_line_indices = sorted(printed_outlier_snps_body_line_indices)
        for x in range(len(sorted_printed_outlier_snps_body_line_indices)-1):
            observed_index_distances.append(sorted_printed_outlier_snps_body_line_indices[x+1]-sorted_printed_outlier_snps_body_line_indices[x])
        for observed_index_distance in observed_index_distances:
            observed_cluster_measure += 1/observed_index_distance
        shuffled_cluster_measures = []
        number_of_shuffle_repetitions = 1000
        for x in range(number_of_shuffle_repetitions):
            shuffled_printed_outlier_snps_body_line_indices = random.sample(range(len(original_body_lines)),len(printed_outlier_snps_body_line_indices))
            shuffled_index_distances = []
            shuffled_cluster_measure = 0
            sorted_shuffled_printed_outlier_snps_body_line_indices = sorted(shuffled_printed_outlier_snps_body_line_indices)
            for x in range(len(sorted_shuffled_printed_outlier_snps_body_line_indices)-1):
                shuffled_index_distances.append(sorted_shuffled_printed_outlier_snps_body_line_indices[x+1]-sorted_shuffled_printed_outlier_snps_body_line_indices[x])
            for shuffled_index_distance in shuffled_index_distances:
                shuffled_cluster_measure += 1/shuffled_index_distance
            shuffled_cluster_measures.append(shuffled_cluster_measure)
        for shuffled_cluster_measure in shuffled_cluster_measures:
            if shuffled_cluster_measure > observed_cluster_measure:
                number_of_reshuffled_indices_that_are_more_clustered_than_the_observed += 1
        if number_of_reshuffled_indices_that_are_more_clustered_than_the_observed < 0.05 * number_of_shuffle_repetitions:
            linked = True
        if number_of_reshuffled_indices_that_are_more_clustered_than_the_observed < 0.01 * number_of_shuffle_repetitions:
            very_linked = True

        if linked:
            # Adjust weights of all positions with clustered outliers.
            index_weights = []
            for x in range(len(original_body_lines)):
                if x not in sorted_printed_outlier_snps_body_line_indices:
                    index_weights.append("N")
                else:
                    if len(sorted_printed_outlier_snps_body_line_indices) == 1:
                        dist_left = x+1
                        dist_right = len(original_body_lines)-x
                    else:                        
                        sorted_printed_outlier_snps_body_line_indices_index = sorted_printed_outlier_snps_body_line_indices.index(x)
                        if sorted_printed_outlier_snps_body_line_indices_index == 0:
                            dist_left = x+1
                            dist_right = observed_index_distances[sorted_printed_outlier_snps_body_line_indices_index]
                        elif sorted_printed_outlier_snps_body_line_indices_index == len(sorted_printed_outlier_snps_body_line_indices)-1:
                            dist_left = observed_index_distances[sorted_printed_outlier_snps_body_line_indices_index-1]
                            dist_right = len(original_body_lines)-x
                        else:
                            dist_left = observed_index_distances[sorted_printed_outlier_snps_body_line_indices_index-1]
                            dist_right = observed_index_distances[sorted_printed_outlier_snps_body_line_indices_index]
                    weight_component_left = 1-(1/dist_left)
                    weight_component_right = 1-(1/dist_right)
                    index_weight = 0.5*(weight_component_left+weight_component_right)
                    index_weights.append(index_weight)

            # Add the weights to the outlier_lines.
            for x in range(len(outlier_lines)):
                line_number = int(outlier_lines[x].split("|")[0].strip())-2
                line_weight = index_weights[line_number]
                outlier_lines[x] = outlier_lines[x].replace("\n","") + " | " + "{0:.2f}".format(line_weight) + "\n"

            # Since the positions referred to the original data set (before removing lines with "0,0" missing data),
            # the original lines have to be read once more.
            valid_body_lines2 = []
            valid_index_weights = []
            weighted_f4s = []
            for x in range(len(original_body_lines)):
                if original_body_lines[x].strip() != "":
                    if "0,0" not in original_body_lines[x]:
                        valid_body_lines2.append(original_body_lines[x])
                        valid_index_weights.append(index_weights[x])
            # Replace "N"s in valid_index_weights so that the mean weight is 1.
            number_of_ns = 0
            weight_sum = 0
            for x in range(len(valid_index_weights)):
                if valid_index_weights[x] == "N":
                    number_of_ns += 1
                else:
                    weight_sum += valid_index_weights[x]
            for x in range(len(valid_index_weights)):
                if valid_index_weights[x] == "N":
                    valid_index_weights[x] = (len(valid_index_weights)-weight_sum)/number_of_ns
            for x in range(len(valid_body_lines2)):
                line_ary = valid_body_lines2[x].strip().split()
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
                    weighted_f4s.append((a-b)*(c-d)*valid_index_weights[x])
            weighted_f4 = numpy.mean(weighted_f4s)

            if number_of_simulations != -1:
                number_of_simulations_with_f4_more_extreme_than_observed_weighted = 0
                for simulated_f4 in simulated_f4s:
                    if weighted_f4 < 0 and simulated_f4 < weighted_f4:
                        number_of_simulations_with_f4_more_extreme_than_observed_weighted += 1
                    elif weighted_f4 > 0 and simulated_f4 > weighted_f4:
                        number_of_simulations_with_f4_more_extreme_than_observed_weighted += 1
                    elif weighted_f4 == 0 and simulated_f4 != weighted_f4:
                        number_of_simulations_with_f4_more_extreme_than_observed_weighted += 1
                proportion_of_simulations_with_f4_more_extreme_than_observed_weighted = number_of_simulations_with_f4_more_extreme_than_observed_weighted/number_of_simulations


    if observed_f4 < 0:
        if len(outlier_valid_body_lines) == 1:
            outfile.write("  The negative f4 value is primarily driven by a single SNP.\n")
            outfile.write("  Removing this SNP from the data set brings the f4 value back to zero.\n")
            outfile.write("  More data will be needed to corroborate support for introgression.\n")
            outfile.write("  More information for this SNP is listed below.\n")
        elif len(outlier_valid_body_lines) > 1:
            outfile.write("  The negative f4 value is primarily driven by " + str(len(outlier_valid_body_lines)) + " SNPs.\n")
            outfile.write("  Removing these SNPs from the data set brings the f4 value back to zero.\n")
            if very_linked:
                outfile.write("  Note that the positions of these " + str(len(outlier_valid_body_lines)) + " SNPs are much more clustered than\n")
                outfile.write("  expected by chance, indicating linkage and thus non-independence of their\n")
                outfile.write("  f4 values. After downweighting clustered f4 outliers, the observed f4 \n")
                outfile.write("  is ")
                outfile.write("{0:.5f}".format(weighted_f4))
                if number_of_simulations != -1:
                    outfile.write(", and the proportion of simulated f4 values smaller than\n")
                    outfile.write("  this is ")
                    outfile.write("{0:.4f}".format(proportion_of_simulations_with_f4_more_extreme_than_observed_weighted))
                outfile.write(".\n")
            elif linked:
                outfile.write("  Note that the positions of these " + str(len(outlier_valid_body_lines)) + " SNPs are more clustered than\n")
                outfile.write("  expected by chance, indicating linkage and thus non-independence of their\n")
                outfile.write("  f4 values. After downweighting clustered f4 outliers, the observed f4 \n")
                outfile.write("  is ")
                outfile.write("{0:.5f}".format(weighted_f4))
                outfile.write(".\n")
            if number_of_simulations != -1:
                if proportion_of_simulations_with_f4_more_extreme_than_observed < 0.05:
                    if len(outlier_valid_body_lines) < 100 and len(outlier_valid_body_lines)/len(valid_body_lines) < 0.1:
                        outfile.write("  It might be worth checking whether these " + str(len(outlier_valid_body_lines)) + " SNPs are more often in or \n")
                        outfile.write("  near genes than expected. If so, this would suggest selection as an \n")
                        outfile.write("  alternative explanation for the negative f4 value.\n")
            outfile.write("  More information for these " + str(len(outlier_valid_body_lines)) + " SNP is listed below.\n")
        else:
            print("ERROR: Unexpected result during search for outlier f4 values")

    if observed_f4 > 0:
        if len(outlier_valid_body_lines) == 1:
            outfile.write("  The positive f4 value is primarily driven by a single SNP.\n")
            outfile.write("  Removing this SNP from the data set brings the f4 value back to zero.\n")
            outfile.write("  More data will be needed to corroborate support for introgression.\n")
            outfile.write("  More information for this SNP is listed below.\n")
        elif len(outlier_valid_body_lines) > 1:
            outfile.write("  The positive f4 value is primarily driven by " + str(len(outlier_valid_body_lines)) + " SNPs.\n")
            outfile.write("  Removing these SNPs from the data set brings the f4 value back to zero.\n")
            if very_linked:
                outfile.write("  Note that the positions of these " + str(len(outlier_valid_body_lines)) + " SNPs are much more clustered than\n")
                outfile.write("  expected by chance, indicating linkage and thus non-independence of their\n")
                outfile.write("  f4 values. After downweighting clustered f4 outliers, the observed f4 \n")
                outfile.write("  is ")
                outfile.write("{0:.5f}".format(weighted_f4))
                if number_of_simulations != -1:
                    outfile.write(", and the proportion of simulated f4 values larger than\n")
                    outfile.write("  this is ")
                    outfile.write("{0:.4f}".format(proportion_of_simulations_with_f4_more_extreme_than_observed_weighted))
                outfile.write(".\n")
            elif linked:
                outfile.write("  Note that the positions of these " + str(len(outlier_valid_body_lines)) + " SNPs are more clustered than\n")
                outfile.write("  expected by chance, indicating linkage and thus non-independence of their\n")
                outfile.write("  f4 values. After downweighting clustered f4 outliers, the observed f4 \n")
                outfile.write("  is ")
                outfile.write("{0:.5f}".format(weighted_f4))
                outfile.write(".\n")
            if number_of_simulations != -1:
                if proportion_of_simulations_with_f4_more_extreme_than_observed < 0.05:
                    if len(outlier_valid_body_lines) < 100 and len(outlier_valid_body_lines)/len(valid_body_lines) < 0.1:
                        outfile.write("  It might be worth checking whether these " + str(len(outlier_valid_body_lines)) + " SNPs are more often in or \n")
                        outfile.write("  near genes than expected. If so, this would suggest selection as an \n")
                        outfile.write("  alternative explanation for the positive f4 value.\n")
            outfile.write("  More information for these " + str(len(outlier_valid_body_lines)) + " SNP is listed below.\n")
        else:
            print("ERROR: Unexpected result during search for outlier f4 values")
    outfile.write("\n")
    outfile.write("Outliers\n")
    allele_freq_string_lengths = []
    for outlier_line in outlier_lines:
        allele_freq_string_lengths.append(len(outlier_line.split("|")[1].strip()))
    allele_freq_col_width = max(len("Allele frequencies"),max(allele_freq_string_lengths))
    if outlier_lines[0].count("|") == 2:
        outfile.write("  " + "Line".rjust(index_string_length) + " | " + "Allele frequencies".ljust(allele_freq_col_width) + " | f4\n")
    elif outlier_lines[0].count("|") == 3:
        f4_string_length = len(outlier_lines[0].split("|")[2].strip())
        outfile.write("  " + "Line".rjust(index_string_length) + " | " + "Allele frequencies".ljust(allele_freq_col_width) + " | " + "f4".ljust(f4_string_length) + " | Adjusted weight\n")
    if outlier_lines[0].count("|") == 2:
        for x in range(len(outlier_lines)):
            outlier_ary = outlier_lines[x].split(" | ")
            outlier_line = outlier_ary[0] + " | " + outlier_ary[1].ljust(allele_freq_col_width) + " | " + outlier_ary[2]
            outfile.write(outlier_line)
    elif outlier_lines[0].count("|") == 3:
        for x in range(len(outlier_lines)):
            outlier_ary = outlier_lines[x].split(" | ")
            outlier_line = outlier_ary[0] + " | " + outlier_ary[1].ljust(allele_freq_col_width) + " | " + outlier_ary[2] + " | " + outlier_ary[3]
            outfile.write(outlier_line)

outfile.write("\n")

if number_of_simulations != -1:
    if log_file_name != "-1" or output_dir != "-1":
        outfile.write("Output\n")
    if log_file_name != "-1":
        outfile.write("  Log file (readable with Tracer): " + str(log_file_name) + "\n")
    if output_dir != "-1":
        outfile.write("  Post-burnin simulated data sets (readable with Treemix): " + str(number_of_simulations) + " files in " + str(output_dir) + "\n")
    if log_file_name != "-1" or output_dir != "-1":
        outfile.write("\n")

