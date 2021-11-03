import argparse
import subprocess

parser = argparse.ArgumentParser(description="SHM_CSR")
parser.add_argument("input", help='Path of the input file')
parser.add_argument('-m', '--method', help='Specify the method for the run: default= custom', default='custom')
parser.add_argument('result_html', help='Path of the output file')
parser.add_argument('result_dir', help='Path where the results will be stored')
parser.add_argument('title', help='Specify the title for the run')
parser.add_argument('-if', '--inc-fr', help='Specify ', default="-")
parser.add_argument('func_filter', help='functionality filter input. Options: productive, unproductive, '
                                        'remove_unknown')
parser.add_argument('unique', help='Options: "VGene,CDR3.IMGT.AA,best_match_class", "VGene,CDR3.IMGT.AA", '
                                   '"productive CDR3.IMGT.AA,best_match_class", "CDR3.IMGT.AA", "VGene,'
                                   'CDR3.IMGT.seq,best_match_class", "VGene,CDR3.IMGT.seq",'
                                   ' "CDR3.IMGT.seq,best_match_class", "CDR3.IMGT.seq", Sequence.ID ')
parser.add_argument('-no', '--naive-output', help='New IMGT archives per class in history. Option: no, yes. Default=no',
                    default='no')
parser.add_argument('-noa', '--naive-output-ca', help='input naive ca: options path to output history ', default='XXX')
parser.add_argument('-nog', '--naive-output-cg', help=' input naive cg: options path to output history ', default='XXX')
parser.add_argument('-nom', '--naive-output-cm', help=' input naive cm: options path to output history ', default='XXX')
parser.add_argument('-noe', '--naive-output-ce', help=' input naive ce: options path to output history ', default='XXX')
parser.add_argument('-nol', '--naive-output-all', help=' input naive all: options path to output history ',
                    default='XXX')
parser.add_argument('filter_unique', help='input for unique filter. Options: remove, remove_vjaa, keep, no')
parser.add_argument('-uc', '--unique-count', help='input for the unique count. Default=2', default='2')
parser.add_argument('class_filter', help='input for class filter. Options: 70_70, 60_55, 70_0, 60_0, 19_0, '
                                         '101_101 ')
parser.add_argument('region_filter', help='input for region filter. Options: leader, CDR1, FR1, FR2')
parser.add_argument('-fa', '--fast', help='input fast process. input: yes, no. Default=no', default='no')

args = parser.parse_args()


subprocess.call(["bash", "wrapper.sh", args.input, args.method, args.result_html, args.result_dir,
                args.title, args.inc_fr, args.func_filter, args.unique, args.naive_output, args.naive_output_ca,
                args.naive_output_cg, args.naive_output_cm, args.naive_output_ce, args.naive_output_all,
                 args.filter_unique, args.unique_count, args.class_filter, args.region_filter, args.fast])

