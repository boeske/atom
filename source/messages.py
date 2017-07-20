"""
ATOM: Terminal color messages
Version 2.0
"""

import sys


class tc:
    RED = '\033[91m'
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    SUBHEADER = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# ------------------------------------------------------------------------------
#                   Informative terminal messages with colors
# ------------------------------------------------------------------------------
def bold(string):
    return tc.BOLD + tc.RED + string + tc.ENDC


def boldgreen(string):
    return tc.BOLD + tc.OKGREEN + string + tc.ENDC


def message(string):
    print("")
    print("---------------------------------------")
    print(tc.SUBHEADER + "* " + string + tc.ENDC)
    print("---------------------------------------")


def message_info(string):
    print(tc.OKBLUE + "   - " + string + tc.ENDC)


def message_done():
    print(tc.OKGREEN + "     --> Done" + tc.ENDC)


def message_warning(string):
    print(tc.FAIL + "   * " + string + tc.ENDC)


def message_error(string):
    print("")
    print(tc.FAIL + "---------------------------------------")
    print(tc.BOLD + "### " + string + tc.ENDC)
    print(tc.FAIL + "---------------------------------------" + tc.ENDC)
    sys.exit()


def message_ok(string):
    print("")
    print("---------------------------------------")
    print(tc.OKGREEN + "* " + string + tc.ENDC)
    print("---------------------------------------")


def startmessage():
    print(",--------------------------------------,")
    print("| " + tc.BOLD + tc.OKBLUE + "              ATOM 2.0      " + tc.ENDC + "         |")
    print("|--------------------------------------|")
    print("| " + tc.OKGREEN + "  The " +
          tc.BOLD + tc.OKBLUE + "A" +
          tc.ENDC + tc.OKGREEN + "lmighty " +
          tc.BOLD + tc.OKBLUE + "TO" +
          tc.ENDC + tc.OKGREEN + "pography " +
          tc.BOLD + tc.OKBLUE + "M" +
          tc.ENDC + tc.OKGREEN + "odifier  " + tc.ENDC + " |")
    print("'--------------------------------------'")


startmessage()
