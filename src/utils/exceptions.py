class FileFormatError(Exception):
    pass


class FileTypeError(Exception):
    pass


class FileValidationError(Exception):
    pass


class SlicerError(Exception):
    pass


class Primer3Error(Exception):
    pass


class FolderCreatorError(Exception):
    pass


class OutputError(Exception):
    pass


class InvalidConfigError(Exception):
    """Config file is not a valid JSON"""
    pass


class InputTypeError(Exception):
    pass
