#!/usr/bin/env python3

import os
import logging
import argparse
import matrixcomponent

import matrixcomponent.JSONparser as JSONparser

LOGGER = logging.getLogger(__name__)
"""logging.Logger: The logger for this module"""


def setup_logging(output_dir):
    """Setup the logging, add a log file"""
    log_name = os.path.join(output_dir, 'log')
    handler = logging.FileHandler(log_name)
    handler.setLevel(args.log_level)
    handler.setFormatter(logging.Formatter(matrixcomponent.LOGGING_FORMAT_STR,
                                           datefmt=matrixcomponent.LOGGING_DATE_FORMAT))
    logging.getLogger().addHandler(handler)


# Helper class to allow multi-line help messages for argparse user parameters:
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def get_arguments():
    """Create the command line interface and return the command line arguments

    Returns
    -------
    Namespace
        The command line arguments

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="")

    parser.add_argument('-j', '--json-file',
                            dest='json_file',
                            required=True,
                            help='input JSON file')

    parser.add_argument('-o', '--out-folder',
                        dest='output_folder',
                        required=True,
                        help='output folder')

    parser.add_argument('-l', '--log-level',
                        default='DEBUG',
                        choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
                        help='level of logging verbosity. DEBUG is most verbose')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = get_arguments()

    setup_logging(args.output_folder)

    LOGGER.info("starting...\n")

    Paths = JSONparser.parse(args.json_file)


#Tasks
# Get example output for tests - Ted
# Parser for ODGI Bin file format - Ted
# Component Segmentation Detection - Josiah and Joerg
#   Python memory object model
# Output format


# For each Path
#     For each node
#         Is there a gap?
#             Is the gap range anywhere else in this individual?
#                 Divider should be somewhere in here
#                 Tolerable range?
#                 Stack up others using the same Link
