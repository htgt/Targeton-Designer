import logging
import os


class CustomLogger:

    """
            This the CustomLogger class to use in our scripts where we need to log the output on 
            terminal and in log file. A log file is generated and placed at /logs/logs.log directory.

            Here, default logger level is DEBUG. 
            This means you can log INFO, WARNING, ERROR, and EXCEPTION.

            Example of use cases:

            from custom_logger.custom_logger import CustomLogger

            logger = CustomLogger(__name__)
            
            logger.info("Info")
            logger.warn("Warning!")
            logger.error("Error")
            logger.exception("Exception")
    """

    def __init__(self, module_name):
        self.module_name = module_name

        # Get directory for logs
        log_directory = os.path.join(os.getcwd(), "logs/")

        # Ensure the log directory exists
        if not os.path.exists(log_directory):
            os.makedirs(log_directory)
        log_file_path = os.path.join(log_directory, f"logs.log")


        # Set up logging to file and format
        logging.basicConfig(
            filename=log_file_path, 
            filemode='a', 
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            level=logging.DEBUG,
            force=True,
        )

        self.handler = logging.StreamHandler()
        self.handler.setLevel(logging.DEBUG)

        # Create logger instance

        self.logger = logging.getLogger(self.module_name)

        # Create custom formatters for different levels

        self.info_formatter = logging.Formatter("%(message)s")
        self.debug_formatter = logging.Formatter("DEBUG: %(message)s")
        self.warning_formatter = logging.Formatter("WARNING: %(message)s")
        self.error_formatter = logging.Formatter("%(asctime)s - ERROR - %(name)s - %(levelname)s - %(message)s")
        self.exception_formatter = logging.Formatter("%(asctime)s - EXCEPTION - %(name)s - %(levelname)s - %(message)s") 
        self.critical_formatter = logging.Formatter("%(asctime)s - CRITICAL - %(name)s - %(levelname)s - %(message)s")


        # Check if there no additional logger is being handled
        
        if not self.logger.handlers:
            self.handler = logging.StreamHandler()
            self.handler.setLevel(logging.DEBUG)
            self.logger.addHandler(self.handler)



    def setup_formatter(self, formatter) -> None:
        """
        
            Function setup_formatter() will format logs according to log level and given log format.

            This function takes one arguments formatter and returns None.  
            For example: self.setup_logger(self.debug_formatter) -> None

        """
        self.handler.setFormatter(formatter)

    
    # Customize and override log level methods

    def info(self, message: str) -> None:
        self.setup_formatter(self.info_formatter)
        self.logger.info(message)
    
    def debug(self, message: str) -> None:
        self.setup_formatter(self.debug_formatter)
        self.logger.debug(message)

    def warning(self, message: str) -> None:
        self.setup_formatter(self.warning_formatter)
        self.logger.warning(message)
    
    def error(self, message: str) -> None:
        self.setup_formatter(self.error_formatter)
        self.logger.error(message)

    def exception(self, message: str) -> None:
        self.setup_formatter(self.exception_formatter)
        self.logger.exception(message)

    def critical(self, message: str) -> None:
        self.setup_formatter(self.critical_formatter)
        self.logger.critical(message)


