import logging
import datetime

def setup_console_logger():
    logging.basicConfig(level=logging.INFO, format=" [%(levelname)s]: %(message)s")

def setup_logger(name, log_file, level=logging.INFO):
    """Function setup as many loggers as you want"""

    handler = logging.FileHandler(log_file)
    formatter = logging.Formatter("[%(levelname)s]: %(message)s")
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def initialise_log(logger_object: 'pass logger object',
                    analysis: 'string for the analysis type', 
                    input_filenames: 'pass list of input filenames',
                    output_filenames: 'pass list of output filenames',
                    sweeps, 
                    tune, 
                    chromosome='all', 
                    other_params_diz:'dictionary of other parameters' ={}  ):

    now = datetime.datetime.now()

    logger_object.info(
        "BAGHERA resutls\n "
        "--------------------------------\n"
        + "Current date & time "
        + now.strftime("%Y-%m-%d %H:%M \n")
    )

    logger_object.info("Analysis: %s \n" %str(analysis) )

    logger_object.info(" Input files: %s \n" %(','.join(input_filenames)))
    logger_object.info(" Output files: %s \n" %(','.join(output_filenames)))

    logger_object.info("- sweeps: %d, \n- burnin: %d, \n- chr: %s\nOther params:\n" %(sweeps, tune, str(chromosome)))
    for k,val in other_params_diz.items():
        logger_object.info("\t- %s: %s\n" %(str(k),str(val)))