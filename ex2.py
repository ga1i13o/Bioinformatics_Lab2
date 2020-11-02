from genominator import aligners
import argparse
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(description="performs global alignment on 2 sequences")
    parser.add_argument("seq1", help="first sequence",
                        type=str)
    parser.add_argument("seq2", help="second sequence",
                        type=str)
    parser.add_argument('scores', nargs='+', help='match, mismatch, and gap score',
                        type=int)

    args = parser.parse_args()
    scores = args.scores
    if len(scores) != 3:
        print("you have to enter 3 scores for match, mismatch and gap")
        sys.exit(-1)

    settings = {"seq1": args.seq1, "seq2": args.seq2, "scores": scores}

    return settings


def main(settings):
    loc_align = aligners.LocalAligner(settings['scores'][0],
                                        settings['scores'][1],
                                        settings['scores'][2])
    loc_align.align(settings['seq1'], settings['seq2'])


if __name__ == '__main__':
    settings = parse_arguments()
    main(settings)

