#!/usr/bin/env python
import logging
import argh
import sys
import os
import inspect
import datetime

from baghera_tool.logging import setup_console_logger
from baghera_tool import regression as regression
from baghera_tool import preprocess as pre
from baghera_tool import heritability as her


def main():
    setup_console_logger()

    """
    command dispatcher
    """
    argh.dispatch_commands([pre.create_files,
                            pre.generate_snp_file,
                            regression.regression,
                            her.heritability

                            ])

if __name__ == "__main__":
    """
    MAIN
    """
    main()
